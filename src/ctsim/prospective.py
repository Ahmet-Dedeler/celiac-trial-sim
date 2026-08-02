"""Predictions about celiac trials that have not reported yet.

Everything else in this repo is retrospective, which is the easy direction. A model that
explains finished trials has had every chance to be fitted to them. So this module states,
before the fact and with a date on it, what three currently-running trials can and cannot
resolve on the villous height : crypt depth endpoint.

They read out between September 2026 and mid-2027. When they do, these predictions are
either right or they are not, and there is no version of this file where that is
ambiguous.

Design parameters here are **not** from ClinicalTrials.gov, which carries none of the
load-bearing detail for any of the three. They come from the sponsors' own protocols
published on the EU Clinical Trials Information System. Two were partly redacted and
recovered another way: Sanofi's gluten dose is blacked out in the protocol but printed in
the patient information sheet, and Dr Falk's is blacked out in protocol v9.0 but readable
in the tracked-changes v8.0 in the same package.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ctsim.published import Source
from ctsim.simulate import mde
from ctsim.variance import fit_injury_variance

# CTIS does not serve its PDFs from ordinary links, which is most of the reason this
# corpus goes unused. A document is fetched in two steps, and the URL in a `Source`
# below is only the human-facing page. To pull the actual protocol:
#
#   curl -s https://euclinicaltrials.eu/ctis-public-api/documents/<ct>/<uuid>/download
#   # -> {"url": "<presigned S3 link, valid 5 minutes>"}, then curl that
#
# The uuids are pinned here so every protocol claim in this module can be taken back
# to the page it came from, by anyone, without clicking through the portal.
CTIS_DOCUMENTS = {
    "TEV-53408 protocol v5.0":
        ("2024-517081-42-00", "d14dcc6c-2eb2-4363-97a5-7b85af771a2a"),
    "ASPIRION protocol v5 (EN, redacted)":
        ("2024-511213-38-00", "69756e33-7854-41a1-95c8-d33e17b8c9f3"),
    "ASPIRION patient information sheet v2.0 (EN) — carries the gluten dose":
        ("2024-511213-38-00", "7cf31c12-df21-4aa6-ad18-4ccbf8fef268"),
    "CEC-013/CEL protocol, tracked changes v8.0":
        ("2023-506150-21-00", "c6bf4218-0796-4131-b4a1-2abe47a6f3d6"),
}


def ctis_download_endpoint(doc_key: str) -> str:
    """The first hop of the CTIS two-step download, for a key of CTIS_DOCUMENTS."""
    ct, uuid = CTIS_DOCUMENTS[doc_key]
    return f"https://euclinicaltrials.eu/ctis-public-api/documents/{ct}/{uuid}/download"


@dataclass(frozen=True)
class LiveTrial:
    name: str
    nct_id: str
    ct_number: str          # EU CTIS, where the real protocol lives
    sponsor: str
    n_per_arm: int
    n_control: int
    design: str             # prevention | restoration
    gluten_mg_per_day: float
    challenge_weeks: int
    entry_criterion: str
    vhcd_is_powered: bool
    assumed_sd: float | None
    readout: str
    sources: list[Source] = field(default_factory=list)

    @property
    def effective_n(self) -> float:
        """Harmonic-mean N per arm, which is what an unequal contrast really buys."""
        return 2.0 / (1.0 / self.n_per_arm + 1.0 / self.n_control)


# The gluten doses below are the single most important numbers in this file, and they
# span roughly 15-fold. Note which direction that runs relative to the entry criterion.
TEV_53408 = LiveTrial(
    name="TEV-53408 (Teva, anti-IL-15)",
    nct_id="NCT06807463",
    ct_number="2024-517081-42-00",
    sponsor="Teva",
    n_per_arm=20, n_control=20,          # 24/24 randomised, 40 evaluable targeted
    design="prevention",
    gluten_mg_per_day=3000.0,
    challenge_weeks=6,
    entry_criterion="Vh:Cd >= 2.0 (healed mucosa, then challenged)",
    vhcd_is_powered=False,
    assumed_sd=0.65,
    readout="2026-09",
    sources=[Source(
        citation="Teva. Clinical Trial Protocol with Amendment 05, VV-06538893 v8.0, "
                 "EU CTIS 2024-517081-42-00.",
        locator="Section 1.1 (challenge), 5.1 (entry), 9.16 (sample size)",
        url="https://euclinicaltrials.eu/ctis-public/view/2024-517081-42-00",
        quote="participants will begin a gluten challenge of 3 g gluten/day and will "
              "continue the gluten challenge for 6 weeks. ... The proposed sample size "
              "is typical for a cohort size in an early development trial... It is not "
              "driven by power considerations. Forty evaluable participants... will "
              "allow estimation of the true treatment effect within a margin of error "
              "of 0.40 using a 2-sided 95% CI, assuming... a common standard deviation "
              "of 0.65",
    )],
)

AMLITELIMAB = LiveTrial(
    name="Amlitelimab / ASPIRION (Sanofi, anti-OX40L)",
    nct_id="NCT06557772",
    ct_number="2024-511213-38-00",
    sponsor="Sanofi",
    n_per_arm=34, n_control=34,          # the SIGE contrast is a single 34-vs-34
    design="restoration",
    gluten_mg_per_day=250.0,
    # SIGE runs for the final 12 weeks of a 24-week controlled period, so this is the
    # exposure window, not the endpoint window: Vh:Cd is read from baseline to week 24.
    challenge_weeks=12,
    entry_criterion="Vh:Cd < 2.5 (active atrophy)",
    vhcd_is_powered=True,
    assumed_sd=None,                      # excised from the protocol; see note below
    readout="2026-08",
    sources=[Source(
        citation="Sanofi. ASPIRION protocol Amendment 05 (09-Dec-2025) and patient "
                 "information sheet, EU CTIS 2024-511213-38-00.",
        locator="Protocol section 9.4 (sample size); ICF (gluten dose, unredacted)",
        url="https://euclinicaltrials.eu/ctis-public/view/2024-511213-38-00",
        quote="The amount of gluten in the SIGE capsule is small, approximately 250 mg "
              "of vital wheat gluten. One capsule will be taken daily from Visit 6 to "
              "Visit 9 of the study for 12 weeks in total. ... A total sample size of "
              "204 participants (randomization ratio 1:1:1:1:1:1, ie, 34 per "
              "intervention group) was determined to demonstrate superiority of each of "
              "the amlitelimab groups versus the matching placebo group on change in "
              "Vh:Cd from baseline to Week 24 with at least 80% power in both the GF and "
              "SIGE populations at a 0.05 two-sided significance level based on the "
              "following assumptions: [blank] ... assuming a common standard [blank] the "
              "expected precision ... will be [blank].",
    )],
)

ZED1227_CEC013 = LiveTrial(
    name="ZED1227 CEC-013/CEL cohort 1 (Dr Falk)",
    nct_id="NCT07298343",
    ct_number="2023-506150-21-00",
    sponsor="Dr Falk Pharma",
    n_per_arm=72, n_control=48,
    design="restoration",
    gluten_mg_per_day=214.0,              # 0.5 g x 3/week, capped at 2.2 g/week
    challenge_weeks=15,
    entry_criterion="VH:CrD <= 2.5 (active atrophy)",
    vhcd_is_powered=False,
    assumed_sd=None,
    readout="2027-06",
    sources=[Source(
        citation="Dr Falk Pharma. CEC-013/CEL protocol, tracked-changes v8.0 "
                 "(07-Aug-2025), EU CTIS 2023-506150-21-00. The current v9.0 redacts "
                 "these; v8.0 in the same package does not.",
        locator="Synopsis (design, SIGE), section 9.4 (sample size)",
        url="https://euclinicaltrials.eu/ctis-public/view/2023-506150-21-00",
        quote="The SIGE bars contain approximately 0.5 g of gluten... approximately 500 "
              "mg three times a week [TIW], but not more than 2.2 g per week ... A "
              "sample size of 41 subjects per treatment group will provide 80% power to "
              "detect a standardised mean difference... of 0.6 ... in the change from "
              "baseline in CDSD GI Specific Symptom Score",
    )],
)

LIVE_TRIALS = [TEV_53408, AMLITELIMAB, ZED1227_CEC013]


# ---------------------------------------------------------------------------
# What a gluten dose does to a mucosa, from the dose-response literature
# ---------------------------------------------------------------------------

# The threshold literature is the only public basis for predicting the injury a given
# challenge will cause. It is unusually clean because the studies were built to find a
# safe dose, so they bracket the interesting range.
GLUTEN_DOSE_RESPONSE = [
    (10.0,   -0.01, "Catassi 2007: 10 mg/day for 90 days, ~1% VH:CD reduction"),
    (50.0,   -0.20, "Catassi 2007: 50 mg/day for 90 days, ~20% VH:CD reduction"),
    (3000.0, -0.61, "ZED1227 CEC-3 and KAN-101 SynCeD placebo arms: 3 g/day, -0.61"),
    (10000.0, -1.53, "Leonard 2021: 10 g/day for 14 days, -1.53"),
]

DOSE_RESPONSE_SOURCE = Source(
    citation="Catassi C et al. A prospective, double-blind, placebo-controlled trial to "
             "establish a safe gluten threshold for patients with celiac disease. "
             "Am J Clin Nutr 2007;85:160-6. PMID 17209192.",
    locator="Results: histologic change by dose group",
    url="https://pubmed.ncbi.nlm.nih.gov/17209192/",
    quote="10 mg/day produced ~1% reduction and 50 mg/day ~20% reduction in VH:CD over "
          "90 days.",
)


# Every VH:CD improvement a restoration-design celiac trial has produced, so the
# predictions below can be checked against something rather than asserted. The list is
# short because the design is: nobody has moved this endpoint upward by much.
#
# Two different quantities live here, and conflating them is the easy mistake. Entries
# marked `difference` are drug minus placebo, which is what an MDE has to clear.
# `single arm` is one arm's own change from baseline, which is a larger number for a
# reason that has nothing to do with a drug working.
RESTORATION_BENCHMARKS = [
    (+0.14, "difference: ZED1227 CEC-004/CEL, 50 mg QD vs placebo, 397 patients "
            "(0.26 vs 0.12, p<0.05) — the largest drug-minus-placebo VH:CD difference "
            "any restoration-design celiac trial has shown, and it came out of a trial "
            "that still missed its primary endpoint"),
    (+0.27, "single arm: CeliAction placebo arm, week 12 — improvement with no drug at "
            "all, which is what SIGE exists to suppress. Not an MDE comparator"),
    (-0.33, "difference: TAK-062 vs placebo, the only completed test of healing against "
            "ongoing SIGE — the drug arm went the wrong way"),
]

# The bar both restoration trials below have to clear. It is a between-arm difference,
# so only the `difference` benchmarks are eligible to set it.
BEST_RESTORATION_EFFECT = 0.14

CEC004_SOURCE = Source(
    citation="Dr Falk Pharma / Zedira. Clinical, histologic and safety outcomes of "
             "ZED1227 treatment in symptomatic celiac disease patients: results from "
             "the phase 2b CEC-004/CEL randomized clinical trial. UEG Week 2025.",
    locator="Results: mean change in VH:CrD by arm; primary endpoint",
    url="https://gutflix.eu/search/d/47fdfd28-a813-11f0-bca4-0242ac140006",
    quote="397 patients randomised. Mean change in VH:CrD from baseline was 0.21 (10 mg "
          "TID), 0.17 (25 mg TID) and 0.26 (50 mg QD) versus 0.12 for placebo; the "
          "between-group difference was significant only for 50 mg QD (p<0.05). The "
          "primary endpoint — the proportion improving both VH:CrD and symptoms — was "
          "met by 25.3% on 50 mg QD versus 19.3% on placebo, not significant.",
)


def best_known_protection() -> float:
    """Most of a 3 g/day gluten injury any drug has ever prevented.

    Derived from ZED1227 CEC-3 rather than written down, so it cannot drift: placebo
    fell 0.61, the best drug arm fell 0.12. The three doses span 72-80%, and it is the
    top of that range a new trial has to be able to see, not the middle.
    """
    from ctsim.published import ZED1227_ARMS

    placebo = next(a for a in ZED1227_ARMS if a[5])
    return max(1.0 - abs(a[2]) / abs(placebo[2]) for a in ZED1227_ARMS if not a[5])


@dataclass
class Prediction:
    trial: LiveTrial
    expected_injury: float | None
    sd_expected: float
    mde_80: float
    min_protection: float | None
    verdict: str


def predict(trial: LiveTrial) -> Prediction:
    """What this trial can resolve, given the empirical variance model."""
    model = fit_injury_variance()

    if trial.design == "prevention":
        # 3 g/day for 6 weeks is exactly the challenge that produced -0.61 in ZED1227's
        # placebo arm, and -0.61 again in KAN-101's. That is a direct read-across, not
        # an extrapolation.
        injury = 0.61
        sd = model.sd_at(injury)
        m = mde(int(trial.effective_n), sd)
        prot = m / injury
        best = best_known_protection()
        verdict = (
            f"needs {prot:.0%} protection for 80% power. The best ZED1227 arm "
            f"delivered {best:.0%}."
            + ("" if prot <= best else " Would miss a drug as good as the best one known.")
        )
    else:
        # Restoration: the control arm is not supposed to deteriorate, so there is no
        # injury to bound the effect. The ceiling is headroom to a normal mucosa, and
        # the entry criterion is what sets it.
        injury = None
        sd = model.floor          # a flat control arm sits at the model's floor
        m = mde(int(trial.effective_n), sd)
        prot = None
        ratio = m / BEST_RESTORATION_EFFECT
        verdict = (
            f"needs a healing difference of {m:.2f} to reach 80% power. The largest "
            f"drug-minus-placebo difference any restoration-design celiac trial has "
            f"produced is {BEST_RESTORATION_EFFECT:+.2f} (ZED1227 CEC-004, which still "
            f"missed its primary). This trial needs {ratio:.1f}x that."
        )
    return Prediction(trial, injury, sd, m, prot, verdict)


def report() -> str:
    model = fit_injury_variance()
    lines = [
        "PROSPECTIVE: what three running celiac trials can resolve",
        "=" * 78,
        "  Stated before readout. Design parameters from sponsor protocols on EU CTIS,",
        "  not ClinicalTrials.gov, which carries none of them.",
        "",
        f"  Variance model: SD = {model.floor:.3f} + {model.slope:.3f} x |injury|",
        "",
        f"  {'trial':<34}{'gluten/day':>12}{'weeks':>7}{'entry':>26}",
    ]
    for t in LIVE_TRIALS:
        lines.append(
            f"  {t.name[:32]:<34}{t.gluten_mg_per_day:>9.0f} mg{t.challenge_weeks:>7}"
            f"{t.entry_criterion.split('(')[0].strip():>26}"
        )

    lines += [
        "",
        "  The gluten dose spans 14-fold, and it runs OPPOSITE to the damage at entry:",
        "  the trial giving the largest challenge is the one enrolling healed patients,",
        "  and the two giving the smallest are enrolling patients already atrophic.",
        "",
        "-" * 78,
        "Per-trial predictions",
        "-" * 78,
    ]
    for t in LIVE_TRIALS:
        p = predict(t)
        lines.append(f"\n  {t.name}  [{t.nct_id}, reads out {t.readout}]")
        lines.append(f"    design      : {t.design}, n={t.n_per_arm} vs {t.n_control}")
        lines.append(f"    VH:CD powered: {'yes' if t.vhcd_is_powered else 'NO'}"
                     + ("" if t.assumed_sd is None else
                        f"   sponsor assumed SD {t.assumed_sd}"))
        lines.append(f"    expected SD : {p.sd_expected:.3f}   MDE(80%) = {p.mde_80:.3f}")
        lines.append(f"    -> {p.verdict}")

    lines += [
        "",
        "-" * 78,
        "Benchmark: what a restoration design has ever produced on this endpoint",
        "-" * 78,
    ]
    for value, note in RESTORATION_BENCHMARKS:
        lines.append(f"  {value:+.2f}   {note}")
    lines += [
        "",
        "  Only the `difference` rows are comparable to an MDE. Both restoration trials",
        "  need a larger drug-minus-placebo difference than any celiac trial has yet",
        "  delivered on VH:CD. That is the prediction, and it is checkable.",
        "",
        "-" * 78,
        "Where each one is exposed",
        "-" * 78,
        "  Teva is the only one of the three that manufactures a large, predictable",
        "  histologic insult (3 g/day is the dose behind every positive VH:CD result on",
        "  record), enrols a healed population that can show it, and assumes a",
        "  conservative SD (0.65 against an empirical 0.58). It is also the only one that",
        "  openly declines to power for the endpoint. Its exposure is arithmetic: at 20",
        "  evaluable per arm the 95% interval is about +/-0.40 on an effect whose ceiling",
        "  is 0.61, so the result will be compatible with anything from no effect to",
        "  complete protection.",
        "",
        "  Sanofi and Dr Falk run the opposite design: ~250 and ~214 mg/day into mucosas",
        "  that are already atrophic. At those doses the threshold literature expects",
        "  little to no further flattening (50 mg/day gives ~20% over 90 days; 10 mg/day",
        "  gives ~1%), which is the point — SIGE exists to stop the placebo arm from",
        "  *improving*, not to injure it. TAK-062 shows it works: its SIGE placebo arm",
        "  moved -0.006 over 24 weeks, dead flat.",
        "",
        "  But a flat control arm relocates the whole burden onto drug-driven healing,",
        "  and both trials admit patients up to a Vh:Cd of 2.5 — close enough to normal",
        "  that some enrollees have well under half a unit of headroom. The one completed",
        "  test of a drug healing against ongoing SIGE is TAK-062, where the drug arm",
        "  moved -0.338, worse than its flat placebo.",
        "",
        "  Falk's VH:CrD is a secondary endpoint riding on a sample size computed for a",
        "  symptom score; the protocol contains no VH:CrD effect-size or SD assumption",
        "  anywhere. Sanofi powered on VH:CD at 80% but its assumptions are physically",
        "  excised from the protocol, in every language version.",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
