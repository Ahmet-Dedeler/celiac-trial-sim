"""Variance model for celiac disease histology endpoints.

Two sources feed this model:

1. **Empirical between-patient SDs** for change-from-baseline, extracted from posted
   ClinicalTrials.gov results (see `ctsim.curate`). This is the quantity that actually
   determines trial power.

2. **Literature measurement-error estimates**, used to decompose that total variance
   into reader / sampling / biological components.

Every literature constant below carries its citation. None of them are tuned.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[2]
CURATED = ROOT / "data" / "curated" / "histology_endpoints.csv"


# ---------------------------------------------------------------------------
# Literature constants
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Constant:
    value: float
    units: str
    source: str
    note: str = ""


LIT = {
    # Taavela et al., PLoS One 2013, "Validation of Morphometric Analyses of
    # Small-Intestinal Biopsy Readouts in Celiac Disease" (PMID 24146832;
    # PLoS One 8(10):e76163). NB: 24098545 is a different paper entirely.
    #
    # These "error margins" are Bland-Altman repeatability coefficients: 2 x the SD of
    # the *paired difference* between two reads, not 2 x the SD of one read. The paper's
    # own interobserver limits of agreement pin this down — reported as -0.516 to +0.375,
    # a span of 0.891, and 2 * 1.96 * 0.227 = 0.890. Had 0.227 been a single-read SD the
    # span would have been 1.26. So the value stored here is SD(read_a - read_b).
    #
    # That distinction is the whole ballgame for the decomposition below, because a
    # change-from-baseline score is itself a difference of two independent reads. The
    # reader noise it carries is therefore this quantity directly, with no further
    # sqrt(2) inflation.
    "vhcd_interobserver_sd_of_difference": Constant(
        0.454 / 2, "VH:CD ratio units",
        "Taavela 2013 PLoS One", "error margin (2SD) = 0.454; two different readers",
    ),
    "vhcd_intraobserver_sd_of_difference": Constant(
        0.318 / 2, "VH:CD ratio units",
        "Taavela 2013 PLoS One", "error margin (2SD) = 0.318; one reader, re-read",
    ),
    # Read the provenance before quoting this one. It is NOT an anchor-based minimal
    # clinically important difference — no such study exists for VH:CD. Taavela derived
    # it from their own reader error: "twice the standard deviation was only 0.318 in
    # the intraobserver Bland-Altman analysis. A cautious new cut-off value of 0.4 could
    # thus be assigned to represent a clinically relevant difference between
    # measurements." It is a reading-reproducibility floor for one patient's paired
    # biopsies, rounded up, and the value it replaced (0.5) was itself a lab convention.
    #
    # Tampere consensus (Gut 2018;67:1410) adopted >0.4 citing only this paper, and
    # graded it D. FDA's 2022 draft guidance for celiac drug development names the
    # Marsh-Oberhuber scale and never mentions VH:CD or any effect size at all.
    #
    # Two consequences. It is a within-patient detection limit, so comparing it to a
    # between-patient SD is a ratio of two different things and the "noise exceeds
    # signal" phrasing oversells it. And it is not what trials target: ZED1227 powered
    # on 0.6, KAN-101 on 0.50 (see ctsim.published). Prefer each trial's own registered
    # target where one exists; this constant is the fallback for trials with none.
    "vhcd_clinically_significant": Constant(
        0.40, "VH:CD ratio units",
        "Taavela 2013 PLoS One (PMID 24146832); adopted by Tampere consensus, Gut 2018",
        "2SD intraobserver reading error, rounded up. Grade D. Not an anchor-based MCID",
    ),
    "iel_intraobserver_cv_cd3_paraffin": Constant(
        0.342, "fraction", "Taavela 2013 PLoS One", "34.2% error margin",
    ),
    "iel_intraobserver_cv_he": Constant(
        0.532, "fraction", "Taavela 2013 PLoS One", "53.2% error margin, H&E",
    ),
    # Taavela cite this rather than derive it: "These findings agree well with previous
    # studies in which a change of more than 30% in lymphocyte count has been considered
    # clinically relevant", pointing at Pollock MA et al., Ann Clin Biochem 1992;29:556
    # — a clinical-chemistry methods paper, not a celiac outcome study. Note also that
    # the H&E intraobserver error margin above (53.2%) is larger than this threshold, so
    # 30% is only reachable with CD3 immunostaining.
    "iel_clinically_significant": Constant(
        0.30, "fraction", "Pollock 1992 Ann Clin Biochem, via Taavela 2013 PLoS One",
        "30% change in T-cell IEL density; not celiac-specific in origin",
    ),
    # Taavela et al., Front Immunol 2021;12:713854 — same group, same SOP, 74 specimens
    # and 5 observers, reporting substantially worse reproducibility than the 2013 study.
    # Carried so the reader-error share can be quoted as a range rather than a point.
    "vhcd_intraobserver_sd_of_difference_2021_he": Constant(
        0.528 / 2, "VH:CD ratio units",
        "Taavela 2021 Front Immunol", "H&E error range (2SD) = 0.528",
    ),
    "vhcd_intraobserver_sd_of_difference_2021_apoa4": Constant(
        0.388 / 2, "VH:CD ratio units",
        "Taavela 2021 Front Immunol", "APOA4-stained error range (2SD) = 0.388",
    ),
    "vhcd_interobserver_sd_of_difference_2021_he": Constant(
        1.017 / 2, "VH:CD ratio units",
        "Taavela 2021 Front Immunol", "H&E average interobserver error range (2SD)",
    ),
    # Corazza et al., Clin Gastroenterol Hepatol 2007 (PMID 17544877).
    "marsh_interobserver_kappa": Constant(
        0.35, "kappa", "Corazza 2007 CGH", "Marsh-Oberhuber; 0.55 for a simplified scale",
    ),
    # Bonamico et al., Am J Gastroenterol 2010 (PMID 20372112).
    "patchiness_between_site_discordance": Constant(
        0.466, "fraction of patients",
        "Bonamico 2010 AJG", "different lesion grade at different duodenal sites",
    ),
    "patchiness_within_biopsy_discordance": Constant(
        0.169, "fraction of patients",
        "Bonamico 2010 AJG", "variable lesions within a single biopsy fragment",
    ),
}


# ---------------------------------------------------------------------------
# Empirical SDs from posted trial results
# ---------------------------------------------------------------------------

@dataclass
class EmpiricalSD:
    nct_id: str
    arm_label: str
    n_arm: int | None
    delta: float
    sd: float
    sd_source: str
    time_frame: str
    is_placebo: bool = False
    mechanism: str = ""
    design: str = "unknown"
    param_type: str = "MEAN"

    @property
    def sd_scale(self) -> str:
        """Which variance this SD actually is — they are not interchangeable.

        A raw change score has variance 2*sigma^2*(1-rho). An ANCOVA/MMRM least-squares
        mean adjusts for baseline, so its residual variance is sigma^2*(1-rho^2), which
        is smaller by a factor of (1+rho)/2. Pooling the two without saying so quietly
        averages two different quantities and lands between them.
        """
        return ("ancova_residual" if "SQUARES" in (self.param_type or "").upper()
                else "change_score")

    @property
    def is_unmedicated(self) -> bool:
        """Did this arm receive gluten but no investigational drug?

        These are the arms that measure what a challenge protocol does to the mucosa on
        its own: placebo arms, plus every arm of a challenge-methodology study, where
        nobody got a drug.
        """
        return self.is_placebo or self.mechanism == "gluten_challenge_methodology"


def load_empirical(endpoint: str = "VHCD", scale: str = "ratio",
                   changes_only: bool = True) -> list[EmpiricalSD]:
    """Load per-arm change-from-baseline SDs from the curated dataset."""
    out: list[EmpiricalSD] = []
    with CURATED.open() as fh:
        for r in csv.DictReader(fh):
            if r["endpoint"] != endpoint or r["scale"] != scale:
                continue
            if not r["sd"]:
                continue
            # Baseline rows carry absolute values, not changes. Separate them on the
            # outcome title, which tags baseline rows explicitly, rather than on the
            # value: a change score can legitimately be positive (mucosal healing), so
            # thresholding on the value would silently drop real recovery data.
            title = r["outcome_title"].lower()
            if changes_only and "[baseline]" in title:
                continue
            out.append(
                EmpiricalSD(
                    nct_id=r["nct_id"],
                    arm_label=r["arm_label"],
                    n_arm=int(r["n_arm"]) if r["n_arm"] else None,
                    delta=float(r["value"]),
                    sd=float(r["sd"]),
                    sd_source=r["sd_source"],
                    time_frame=r["time_frame"],
                    is_placebo=r["is_placebo"].strip().lower() == "true",
                    mechanism=r["mechanism"],
                    design=r["design"],
                    param_type=r.get("param_type", "MEAN"),
                )
            )
    return out


def load_baselines(endpoint: str = "VHCD", scale: str = "ratio") -> list[EmpiricalSD]:
    """Load the *baseline* (pre-treatment absolute value) rows, where posted.

    Only a few trials post these. They are what makes `baseline_correlation` possible.
    """
    out: list[EmpiricalSD] = []
    with CURATED.open() as fh:
        for r in csv.DictReader(fh):
            if r["endpoint"] != endpoint or r["scale"] != scale or not r["sd"]:
                continue
            if "[baseline]" not in r["outcome_title"].lower():
                continue
            out.append(
                EmpiricalSD(
                    nct_id=r["nct_id"], arm_label=r["arm_label"],
                    n_arm=int(r["n_arm"]) if r["n_arm"] else None,
                    delta=float(r["value"]), sd=float(r["sd"]),
                    sd_source=r["sd_source"], time_frame=r["time_frame"],
                    is_placebo=r["is_placebo"].strip().lower() == "true",
                    mechanism=r["mechanism"],
                    design=r["design"],
                    param_type=r.get("param_type", "MEAN"),
                )
            )
    return out


def pooled_sd(rows: list[EmpiricalSD]) -> float:
    """Sample-size-weighted pooled SD across arms.

    Uses the standard pooled-variance estimator, sum((n_i-1) s_i^2) / sum(n_i-1).
    Arms without an analysed N are given weight 1.
    """
    num = 0.0
    den = 0.0
    for r in rows:
        w = (r.n_arm - 1) if r.n_arm and r.n_arm > 1 else 1
        num += w * r.sd**2
        den += w
    return math.sqrt(num / den) if den else float("nan")


# ---------------------------------------------------------------------------
# Variance decomposition
# ---------------------------------------------------------------------------

@dataclass
class VarianceDecomposition:
    total_sd: float
    reader_sd: float
    residual_sd: float  # biology + sampling + orientation, not separable from these data
    reader_share: float
    residual_share: float

    def describe(self) -> str:
        return (
            f"total between-patient SD of ΔVH:CD = {self.total_sd:.3f}\n"
            f"  reader (intra-observer, single central reader): {self.reader_sd:.3f} "
            f"({self.reader_share:6.1%} of variance)\n"
            f"  residual (true biology + biopsy site + orientation): {self.residual_sd:.3f} "
            f"({self.residual_share:6.1%} of variance)"
        )


def decompose_vhcd(total_sd: float) -> VarianceDecomposition:
    """Split observed ΔVH:CD variance into reader vs everything-else.

    Two things have to be right here, and both cut the reader term down:

    1. Trials use a single blinded central reader, so the relevant term is
       *intra*-observer, not inter-observer. Using the interobserver figure would
       overstate reader noise — that is the error behind the folk claim that
       measurement error alone exceeds the effect size.

    2. The stored constant is already the SD of a difference between two reads (see the
       note on LIT). A change-from-baseline score is also a difference of two reads, so
       it carries exactly that SD. Multiplying by another sqrt(2) — tempting, because
       "a change score is a difference of two reads, so double the variance" is true of
       *single-read* variance — double counts, and inflates the reader share from
       ~5% to ~9%.

    Both errors point the same way: they make pathologist disagreement look like a
    bigger share of the problem than it is.
    """
    reader_sd = LIT["vhcd_intraobserver_sd_of_difference"].value
    reader_var = reader_sd**2
    total_var = total_sd**2
    residual_var = max(total_var - reader_var, 0.0)
    return VarianceDecomposition(
        total_sd=total_sd,
        reader_sd=reader_sd,
        residual_sd=math.sqrt(residual_var),
        reader_share=reader_var / total_var,
        residual_share=residual_var / total_var,
    )


# ---------------------------------------------------------------------------
# Baseline correlation and the analysis-method question
# ---------------------------------------------------------------------------

@dataclass
class BaselineCorrelation:
    nct_id: str
    arm_label: str
    sd_baseline: float
    sd_change: float
    rho: float


def baseline_correlation(nct_id: str | None = None) -> list[BaselineCorrelation]:
    """Recover the baseline-to-follow-up correlation of VH:CD from posted summaries.

    Nobody posts rho. But a trial that posts *both* the baseline SD and the
    change-from-baseline SD for the same arm has already told you what it is:

        Var(change) = sd_b^2 + sd_f^2 - 2*rho*sd_b*sd_f

    Assuming the follow-up SD equals the baseline SD (sd_f = sd_b = sigma) this
    collapses to Var(change) = 2*sigma^2*(1 - rho), so

        rho = 1 - Var(change) / (2 * Var(baseline))

    The equal-SD assumption is the weak link and it is not free: under a gluten
    challenge the mucosa flattens, which compresses the follow-up spread, so sd_f is
    probably a little *below* sd_b and these estimates are correspondingly rough.
    They are reported as a range, never as a single number.
    """
    changes = {(r.nct_id, r.arm_label): r for r in load_empirical()}
    out: list[BaselineCorrelation] = []
    for b in load_baselines():
        if nct_id and b.nct_id != nct_id:
            continue
        ch = changes.get((b.nct_id, b.arm_label))
        if ch is None or b.sd <= 0:
            continue
        out.append(
            BaselineCorrelation(
                nct_id=b.nct_id, arm_label=b.arm_label,
                sd_baseline=b.sd, sd_change=ch.sd,
                rho=1 - ch.sd**2 / (2 * b.sd**2),
            )
        )
    return out


def ancova_variance_ratio(rho: float) -> float:
    """Variance of an ANCOVA contrast relative to a change-score contrast.

    Analysing follow-up with baseline as a covariate has residual variance
    sigma^2*(1 - rho^2); analysing the change score has 2*sigma^2*(1 - rho). The ratio
    is (1 + rho) / 2, so ANCOVA is never worse and the gain grows as rho falls.

    This is the cheapest available fix — it costs nothing but a line in the statistical
    analysis plan — which is exactly why it is worth knowing how little it buys.
    """
    return (1.0 + rho) / 2.0


# ---------------------------------------------------------------------------
# Would a different endpoint help?
# ---------------------------------------------------------------------------

@dataclass
class EndpointComparison:
    endpoint: str
    scale: str
    meaningful: float
    meaningful_note: str
    pooled_sd: float
    n_arms: int
    standardized_effect: float  # meaningful / SD — this is what sets required N


def endpoint_comparison(population: set[str] | None = None) -> list[EndpointComparison]:
    """VH:CD versus IEL density on the only footing that matters for trial size.

    The obvious response to "VH:CD is too noisy" is "then use a different endpoint".
    Required N depends on the *standardized* effect (meaningful change / SD), not on
    either quantity alone, so a noisier endpoint with a proportionally larger meaningful
    change costs nothing. This computes both and lets them be compared.

    Restricted by default to gluten-challenge trials, because IEL density behaves
    differently in a challenged mucosa than in non-responsive celiac on a stable
    gluten-free diet, and pooling those two would compare populations, not endpoints.

    Caveat that limits how far this can be pushed: the IEL threshold is a *relative*
    30% change, so converting it to absolute cells/100 needs a baseline IEL density,
    and exactly one trial in the dataset posts one (NCT03409796, 26.7 cells/100).
    """
    population = population or {"NCT06001177", "NCT03409796"}

    vhcd = [r for r in load_empirical("VHCD", "ratio") if r.nct_id in population]
    iel = [r for r in load_empirical("IEL", "cells_per_100") if r.nct_id in population]

    iel_baselines = [b.delta for b in load_baselines("IEL", "cells_per_100")]
    if not iel_baselines or not vhcd or not iel:
        return []
    baseline_iel = mean(iel_baselines)
    iel_meaningful = LIT["iel_clinically_significant"].value * baseline_iel

    out = [
        EndpointComparison(
            endpoint="VH:CD", scale="ratio",
            meaningful=LIT["vhcd_clinically_significant"].value,
            meaningful_note="absolute threshold (Taavela 2013)",
            pooled_sd=pooled_sd(vhcd), n_arms=len(vhcd),
            standardized_effect=(LIT["vhcd_clinically_significant"].value
                                 / pooled_sd(vhcd)),
        ),
        EndpointComparison(
            endpoint="IEL density", scale="cells per 100 enterocytes",
            meaningful=iel_meaningful,
            meaningful_note=f"30% of a {baseline_iel:.1f} baseline (Taavela 2013)",
            pooled_sd=pooled_sd(iel), n_arms=len(iel),
            standardized_effect=iel_meaningful / pooled_sd(iel),
        ),
    ]
    return out


def reader_share_range(total_sd: float) -> list[tuple[str, float, float]]:
    """Reader share of variance under every published estimate of reader error.

    The 2013 and 2021 validation studies come from the same group using the same SOP and
    disagree by a factor of ~1.7 on intraobserver error. Quoting only the smaller one
    would be picking the estimate that flatters the conclusion. Reported as a range so
    the claim can be attacked at its weakest point rather than its strongest.

    Only the intraobserver rows are the honest comparator for a trial with one blinded
    central reader; the interobserver row is included to show what the number becomes if
    that assumption fails.
    """
    keys = [
        ("Taavela 2013, intraobserver", "vhcd_intraobserver_sd_of_difference"),
        ("Taavela 2021, intraobserver APOA4",
         "vhcd_intraobserver_sd_of_difference_2021_apoa4"),
        ("Taavela 2021, intraobserver H&E",
         "vhcd_intraobserver_sd_of_difference_2021_he"),
        ("Taavela 2021, INTERobserver H&E (single reader assumption fails)",
         "vhcd_interobserver_sd_of_difference_2021_he"),
    ]
    return [(label, LIT[k].value, LIT[k].value**2 / total_sd**2) for label, k in keys]


def summary() -> str:
    rows = load_empirical()
    sd = pooled_sd(rows)
    dec = decompose_vhcd(sd)
    lines = [
        "Empirical ΔVH:CD between-patient SDs from posted trial results",
        "-" * 62,
    ]
    for r in rows:
        lines.append(
            f"  {r.nct_id}  {r.arm_label[:34]:<36} n={r.n_arm!s:<4} "
            f"Δ={r.delta:+.3f}  SD={r.sd:.3f}  ({r.sd_source})"
        )
    lines += [
        "",
        (f"pooled SD = {sd:.3f}  (across {len(rows)} arms, "
         f"{len({r.nct_id for r in rows})} trials)"),
        f"clinically meaningful change = {LIT['vhcd_clinically_significant'].value}",
        f"  ratio SD / meaningful-change = {sd / LIT['vhcd_clinically_significant'].value:.2f}",
        "",
        dec.describe(),
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
