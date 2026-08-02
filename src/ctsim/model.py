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
    # Small-Intestinal Biopsy Readouts in Celiac Disease" (PMID 24098545).
    # Error margins are reported as 2SD; we store the implied 1SD.
    "vhcd_interobserver_sd": Constant(
        0.454 / 2, "VH:CD ratio units",
        "Taavela 2013 PLoS One", "reported error margin (2SD) = 0.454 on a single read",
    ),
    "vhcd_intraobserver_sd": Constant(
        0.318 / 2, "VH:CD ratio units",
        "Taavela 2013 PLoS One", "reported error margin (2SD) = 0.318",
    ),
    "vhcd_clinically_significant": Constant(
        0.40, "VH:CD ratio units",
        "Taavela 2013 PLoS One", "threshold for a clinically meaningful change",
    ),
    "iel_intraobserver_cv_cd3_paraffin": Constant(
        0.342, "fraction", "Taavela 2013 PLoS One", "34.2% error margin",
    ),
    "iel_intraobserver_cv_he": Constant(
        0.532, "fraction", "Taavela 2013 PLoS One", "53.2% error margin, H&E",
    ),
    "iel_clinically_significant": Constant(
        0.30, "fraction", "Taavela 2013 PLoS One", "30% change in T-cell IEL density",
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
            # Baseline rows carry absolute values, not changes. On the ratio scale a
            # baseline VH:CD is ~2-3 and a change is <=0; use that to separate them.
            val = float(r["value"])
            if changes_only and val > 1.0:
                continue
            out.append(
                EmpiricalSD(
                    nct_id=r["nct_id"],
                    arm_label=r["arm_label"],
                    n_arm=int(r["n_arm"]) if r["n_arm"] else None,
                    delta=val,
                    sd=float(r["sd"]),
                    sd_source=r["sd_source"],
                    time_frame=r["time_frame"],
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

    A change score is the difference of two independent reads, so reader variance on
    the change is 2x the single-read variance: sd_change = sqrt(2) * sd_single.

    Trials use a single blinded central reader, so the relevant single-read term is
    *intra*-observer, not inter-observer. Using the inter-observer figure here would
    overstate the reader contribution — which is exactly the error that makes the
    naive "measurement error exceeds the effect size" claim too strong.
    """
    reader_sd = math.sqrt(2) * LIT["vhcd_intraobserver_sd"].value
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
            f"  {r.nct_id}  {r.arm_label[:34]:<36} n={str(r.n_arm):<4} "
            f"Δ={r.delta:+.3f}  SD={r.sd:.3f}  ({r.sd_source})"
        )
    lines += [
        "",
        f"pooled SD = {sd:.3f}  (across {len(rows)} arms, "
        f"{len({r.nct_id for r in rows})} trials)",
        f"clinically meaningful change = {LIT['vhcd_clinically_significant'].value}",
        f"  ratio SD / meaningful-change = {sd / LIT['vhcd_clinically_significant'].value:.2f}",
        "",
        dec.describe(),
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
