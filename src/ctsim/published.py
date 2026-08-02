"""Trial data that exists only in papers and protocols, not in the registry API.

Two kinds of thing live here, both hand-entered, both carrying a verbatim quote:

1. **Results** from trials the ClinicalTrials.gov API cannot reach. The most important
   is ZED1227 — the only celiac drug trial to date that hit a histologic primary
   endpoint, registered in EudraCT with no NCT number at all.

2. **Design assumptions**: the effect size and SD each trial powered itself on, taken
   from its published statistical-analysis plan or protocol. These are what make it
   possible to ask whether the field's variance assumptions were right, rather than
   only whether its sample sizes were big enough.

Hand-entered data is a liability, so the rules here are stricter than for the scraped
dataset: every number carries the document it came from and the sentence it appeared
in, and anything requiring arithmetic shows the arithmetic in code rather than baking
in a computed constant.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ctsim.model import EmpiricalSD


@dataclass(frozen=True)
class Source:
    citation: str
    locator: str  # table / section, so a reader can find it
    url: str
    quote: str


# ---------------------------------------------------------------------------
# ZED1227 (CEC-3) — Schuppan et al., NEJM 2021
# ---------------------------------------------------------------------------

ZED1227_SOURCE = Source(
    citation="Schuppan D et al. A Randomized Trial of a Transglutaminase 2 Inhibitor "
             "for Celiac Disease. N Engl J Med 2021;385:35-45. PMID 34192430. "
             "EudraCT 2017-002241-30 (protocol CEC-3).",
    locator="Table 2, rows 'Change from baseline' (least-squares means, 95% CI)",
    url="https://www.nejm.org/doi/10.1056/NEJMoa2032441",
    quote="(estimated change, -0.61; 95% confidence interval [CI], -0.78 to -0.44)"
          " ... (estimated change, -0.12 [95% CI, -0.27 to 0.03] and -0.13 [95% CI,"
          " -0.28 to 0.03], respectively); efficacy in the 10-mg group was slightly"
          " less (estimated change, -0.17; 95% CI, -0.33 to -0.01).",
)

# Per-arm LS-mean change in VH:CD with its 95% CI, and the N contributing to the
# primary histology analysis. Table 2 of the paper.
#   (arm label, n analysed, LS-mean change, CI low, CI high, is_placebo)
ZED1227_ARMS = [
    ("Placebo",        30, -0.61, -0.78, -0.44, True),
    ("ZED1227 10 mg",  35, -0.17, -0.33, -0.01, False),
    ("ZED1227 50 mg",  39, -0.12, -0.27, +0.03, False),
    ("ZED1227 100 mg", 38, -0.13, -0.28, +0.03, False),
]

# Raw (unadjusted) baseline and week-6 levels, Table 2, "mean +/- SD". Kept because
# they are the only route to a change-score-scale SD for this trial, and because the
# baseline SD is what `baseline_correlation` needs.
#   (arm, n, baseline mean, baseline SD, week6 mean, week6 SD)
ZED1227_LEVELS = [
    ("Placebo",        30, 1.98, 0.33, 1.39, 0.61),
    ("ZED1227 10 mg",  35, 2.01, 0.30, 1.85, 0.53),
    ("ZED1227 50 mg",  39, 2.04, 0.32, 1.91, 0.44),
    ("ZED1227 100 mg", 38, 2.09, 0.35, 1.94, 0.48),
]

Z95 = 1.959963985


def zed1227_residual_sd() -> tuple[float, int]:
    """Residual SD behind ZED1227's LS-mean changes, with its degrees of freedom.

    Every arm's CI implies nearly the same SD (0.475-0.487) — which is the tell that
    these are **not four independent estimates**. They come from one generalized linear
    model with a common residual variance, so treating them as four arms would inflate
    the degrees of freedom roughly fourfold and make the pooled estimate look far more
    certain than it is.

    So this returns a single estimate at the model's own df: n_total - n_groups - 1
    (four group effects plus the baseline VH:CD covariate).
    """
    per_arm = []
    for _label, n, _delta, lo, hi in ((a[0], a[1], a[2], a[3], a[4]) for a in ZED1227_ARMS):
        se = (hi - lo) / (2 * Z95)
        per_arm.append(se * math.sqrt(n))
    n_total = sum(a[1] for a in ZED1227_ARMS)
    df = n_total - len(ZED1227_ARMS) - 1
    return sum(per_arm) / len(per_arm), df


def zed1227_control_injury() -> tuple[float, float, float]:
    """Placebo-arm mucosal injury under a 6-week 3 g/day gluten challenge, with CI."""
    label, _n, delta, lo, hi, _p = ZED1227_ARMS[0]
    assert label == "Placebo"
    return delta, lo, hi


def zed1227_rows() -> list[EmpiricalSD]:
    """ZED1227 as a single row on the ANCOVA-residual scale.

    Deliberately one row, not four — see `zed1227_residual_sd`. The n is set so that
    `pooled_sd`'s (n-1) weighting reproduces the model's real degrees of freedom.
    """
    sd, df = zed1227_residual_sd()
    delta, _lo, _hi = zed1227_control_injury()
    return [
        EmpiricalSD(
            nct_id="EudraCT2017-002241-30",
            arm_label="Placebo (ZED1227 CEC-3)",
            n_arm=df + 1,
            delta=delta,
            sd=sd,
            sd_source="from_lsmean_ci95",
            time_frame="Baseline to week 6",
            is_placebo=True,
            mechanism="tg2_inhibitor",
            design="prevention",
        )
    ]


# ---------------------------------------------------------------------------
# Change-from-baseline VH:CD arms published in papers but not the registry
# ---------------------------------------------------------------------------

# Genuine arithmetic means with genuine SDs of the change score — the same quantity the
# `posted_sd` registry rows carry, so these merge directly. LS-mean-only trials are
# deliberately absent: see ALV003-1021 under NOT_AVAILABLE.
#   (trial, id, arm, n, delta, sd, is_placebo, design)
PAPER_ARMS = [
    # Murray et al., Gastroenterology 2022;163:1510-21, Table 2 (mITT).
    ("IMGX003 CeliacShield", "NCT03585478", "Placebo",   22, -0.35, 0.616, True,  "prevention"),
    ("IMGX003 CeliacShield", "NCT03585478", "IMGX003",   21, -0.04, 0.466, False, "prevention"),
    # Murray et al., Gastroenterology 2017;152:787-98, Table 2 (MITT, week 12).
    # No gluten challenge: real-world exposure on a stable GFD, hence restoration.
    ("CeliAction", "NCT01917630", "Placebo",            125, 0.27, 0.401, True,  "restoration"),
    ("CeliAction", "NCT01917630", "Latiglutenase 100mg", 47, 0.12, 0.463, False, "restoration"),
    ("CeliAction", "NCT01917630", "Latiglutenase 300mg", 77, 0.15, 0.413, False, "restoration"),
    ("CeliAction", "NCT01917630", "Latiglutenase 450mg", 39, 0.05, 0.501, False, "restoration"),
    ("CeliAction", "NCT01917630", "Latiglutenase 600mg", 80, 0.14, 0.493, False, "restoration"),
    ("CeliAction", "NCT01917630", "Latiglutenase 900mg", 37, 0.11, 0.460, False, "restoration"),
]

PAPER_ARM_SOURCES = {
    "NCT03585478": Source(
        citation="Murray JA et al. Latiglutenase Protects the Mucosa and Attenuates "
                 "Symptom Severity in Patients With Celiac Disease Exposed to a Gluten "
                 "Challenge. Gastroenterology 2022;163:1510-1521.",
        locator="Table 2, 'Histologic Efficacy Analysis - mITT Population', "
                "row 'Change in Vh:Cd', mean (SD)",
        url="https://pmc.ncbi.nlm.nih.gov/articles/PMC9707643/",
        quote="Change in Vh:Cd: placebo -0.35 (0.616); IMGX003 -0.04 (0.466); "
              "between-group P = .0570.",
    ),
    "NCT01917630": Source(
        citation="Murray JA et al. No Difference Between Latiglutenase and Placebo in "
                 "Reducing Villous Atrophy or Improving Symptoms in Patients With "
                 "Symptomatic Celiac Disease. Gastroenterology 2017;152:787-798.",
        locator="Table 2, change in Vh:Cd from baseline to week 12, MITT, mean (SD)",
        url="https://www.gastrojournal.org/article/S0016-5085(16)35346-X/fulltext",
        quote="Change in Vh:Cd, mean (SD): placebo 0.27 (0.401); 100 mg 0.12 (0.463); "
              "300 mg 0.15 (0.413); 450 mg 0.05 (0.501); 600 mg 0.14 (0.493); "
              "900 mg 0.11 (0.460).",
    ),
}

# Trials whose histology is published but whose dispersion is not recoverable. Listed
# so the absence is a recorded fact rather than an oversight, and so nobody spends a
# second afternoon looking.
NOT_AVAILABLE = {
    "ALV003-1021 (Lahdeaho 2014, Gastroenterology)":
        "No SD, SE, CI or IQR for VH:CrD anywhere in the paper or its three "
        "supplements; Figure 4 plots mean +/- SE as graphics only. The paper also "
        "contradicts itself, calling the same 2.0 value a mean in the abstract and a "
        "median in the results. Not recoverable from the VCIEL re-analysis either: "
        "its effect size needs both an unpublished mean difference and an unpublished "
        "baseline SD, which is one equation in two unknowns.",
    "Nexvax2 RESET CeD (NCT03644069)":
        "Histology was an unregistered exploratory/safety endpoint, published only in "
        "the author manuscript. Table 4 gives SDs on the baseline and end-of-study "
        "*levels* and an ANCOVA between-arm mean difference of 0.22 (0.06, 0.39), but "
        "no within-arm SD of the change. The caption says 'Mean (standard deviation)' "
        "while the results text calls the same values medians.",
    "AMG 714 (NCT02637141, NCT02633020)":
        "Registry reports percent change only. Both Lancet Gastro Hepatol 2019 papers "
        "are closed-access with no repository copy, so whether they give raw ratio "
        "units is untested rather than known.",
    "TAK-062 (NCT05353985)":
        "No raw SD in any source; no per-arm baseline Vh:Cd; IELs measured but reported "
        "only qualitatively. The residual SD used here comes from the posted ANCOVA "
        "contrast, which is the only valid route.",
    "TAK-101 Phase 2b (NCT04530123)":
        "Has no histology endpoint at all — 'villous', 'crypt' and 'intraepithelial' "
        "return zero hits across the registry record. No VH:CD data will ever exist.",
}


def paper_rows() -> list[EmpiricalSD]:
    """Published arms as `EmpiricalSD`, mergeable with the scraped registry rows."""
    return [
        EmpiricalSD(
            nct_id=nct, arm_label=arm, n_arm=n, delta=delta, sd=sd,
            sd_source="paper_posted_sd", time_frame="", is_placebo=placebo,
            mechanism="glutenase", design=design, param_type="MEAN",
        )
        for _trial, nct, arm, n, delta, sd, placebo, design in PAPER_ARMS
    ]


# ---------------------------------------------------------------------------
# A measured variance decomposition
# ---------------------------------------------------------------------------

# Takeda reported a direct decomposition of Vh:Cd variability from the TAK-062 Phase 2
# — the experiment this repo previously said could not be done from public data.
# Retrieved from UEG Week 2025 abstract MP739 (Maxwell, Isola, Valimaki, Robert, Cheng,
# Zarei, Leffler), verbatim:
#
#   "Analysis of the source of variability in the Vh:Cd measurement showed that 52% of
#    the variability was at the patient level and 23% at the biopsy level. Only 1% of
#    the variability was due to reader effects."
#
# Caveat that has to travel with it: the same abstract states Vh:Cd "measurements were
# made by multiple readers ... and averaged", so 1% is the reader term *after* averaging
# readers, not the error of a single read. It is a floor, not a like-for-like comparison
# with Taavela's single-reader figures. The biopsy-level share is the useful number:
# it is a measured value for what averaging more fragments per timepoint can remove.
MEASURED_VARIANCE_SHARES = {
    "patient": 0.52,
    "biopsy": 0.23,
    "reader": 0.01,
}

MEASURED_VARIANCE_SOURCE = Source(
    citation="Maxwell JR, Isola J, Valimaki A, Robert M, Cheng J, Zarei M, Leffler DA. "
             "Assessment of histologic endpoints in a phase 2 trial of TAK-062 in "
             "celiac disease. UEG Week 2025, abstract MP739.",
    locator="Results, first sentence",
    url="https://ueg2025.abstract.documedias.systems/api/v1/manager/abstract/multi/"
        "html/id/4206/template/planner_preview",
    quote="Analysis of the source of variability in the Vh:Cd measurement showed that "
          "52% of the variability was at the patient level and 23% at the biopsy level. "
          "Only 1% of the variability was due to reader effects.",
)


# ---------------------------------------------------------------------------
# What each trial assumed when it sized itself
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DesignAssumption:
    trial: str
    nct_id: str
    target_effect: float   # between-arm difference the trial powered to detect
    assumed_sd: float
    n_per_arm: int
    claimed_power: float
    endpoint_role: str     # PRIMARY | SECONDARY — was VH:CD what it powered on?
    source: Source
    # "vhcd_ratio" is change in the raw VH:CD ratio (compatible with our observed
    # arms). "percent_change" is baseline-to-follow-up %-change — a different scale
    # that cannot be compared to our ratio-unit SDs without an unpublished baseline.
    units: str = "vhcd_ratio"


ASSUMPTIONS: list[DesignAssumption] = [
    DesignAssumption(
        trial="ZED1227 (CEC-3)",
        nct_id="EudraCT2017-002241-30",
        target_effect=0.6,
        assumed_sd=0.8,
        n_per_arm=34,
        claimed_power=0.80,
        endpoint_role="PRIMARY",
        source=Source(
            citation="Schuppan D et al. NEJM 2021;385:35-45.",
            locator="Statistical Analysis; identical text in Protocol section 9.4",
            url="https://www.nejm.org/doi/10.1056/NEJMoa2032441",
            quote="We estimated that a sample of 136 patients, or 34 patients per "
                  "group, would provide the trial with 80% power for the primary "
                  "analysis, assuming an alpha error of 0.05, an effect size of 0.6, "
                  "and standard deviation of 0.8.",
        ),
    ),
    DesignAssumption(
        trial="KAN-101 SynCeD",
        nct_id="NCT06001177",
        target_effect=0.50,
        assumed_sd=0.5,
        n_per_arm=26,
        claimed_power=0.94,
        endpoint_role="PRIMARY",
        source=Source(
            citation="Anokion/Pfizer. Protocol KAN-101-03 Statistical Analysis Plan, "
                     "2025-02-03, section 1.3 'Sample Size'.",
            locator="Section 1.3",
            url="https://cdn.clinicaltrials.gov/large-docs/77/NCT06001177/SAP_001.pdf",
            quote="the study will have approximately 94% power under a sample size of "
                  "26 participants per arm to detect a treatment difference in LSM of "
                  "0.50 in change in ratio of Vh:Cd between treatment arms, assuming "
                  "that the common standard deviation is 0.5 and the difference "
                  "between the means is 0.5.",
        ),
    ),
    DesignAssumption(
        trial="IMGX003 CeliacShield",
        nct_id="NCT03585478",
        target_effect=0.40,
        assumed_sd=0.45,
        n_per_arm=25,
        claimed_power=0.86,
        endpoint_role="PRIMARY",
        source=Source(
            citation="Murray JA et al. Gastroenterology 2022;163:1510-1521.",
            locator="Statistical analysis / sample size paragraph",
            url="https://pmc.ncbi.nlm.nih.gov/articles/PMC9707643/",
            quote="A sample size of 50 ITT patients (25 patients per treatment group) "
                  "provided 86% power and two-sided 5% Type 1 error to detect a 0.40 "
                  "between the treatment groups in change from baseline Vh:Cd at "
                  "Week 6 assuming a standard deviation of 0.45.",
        ),
    ),
    DesignAssumption(
        trial="AMG 714 CELIM-NRCD-001",
        nct_id="NCT02637141",
        target_effect=40.0,   # percentage points of %-change, not ratio units
        assumed_sd=36.0,
        n_per_arm=17,         # evaluable completers the power calc used
        claimed_power=0.888,
        endpoint_role="PRIMARY",
        units="percent_change",
        source=Source(
            citation="Celimmune. CELIM-NRCD-001 Statistical Analysis Plan v1.0, "
                     "19 Apr 2017, section 6 'Estimation of sample size'.",
            locator="Section 6",
            url="https://cdn.clinicaltrials.gov/large-docs/41/NCT02637141/SAP_001.pdf",
            quote="Common SD = 36 for the baseline to Week 12 %-change in VH:CD. "
                  "... close to 90% (88.8%) power to detect a 40-point difference "
                  "between the placebo arm and the 300mg high-dose arm ... "
                  "17 completed evaluable subjects/arm ... will be needed for the "
                  "primary analysis.",
        ),
    ),
]


# Trials that ran histology but sized themselves on a different primary. Recording
# the absence is the finding: their VH:CD result was never powered, so a miss (or a
# hit) on histology is not a statement about assay design the way the rows above are.
NOT_POWERED_ON_VHCD: dict[str, Source] = {
    "NCT03738475": Source(
        citation="COUR/Takeda. Protocol TGLIA-5.002 v2.0, section 8.1 Sample Size.",
        locator="Section 8.1",
        url="https://cdn.clinicaltrials.gov/large-docs/75/NCT03738475/Prot_001.pdf",
        quote="For the primary efficacy endpoint, using a 2-sided 0.05 significance "
              "level, a statistical power of ~70%, and assuming an increase in mean "
              "IFN-ɣ SFUs in the placebo group of 75 (standard deviation [SD] 100) "
              "and in the TIMP-GLIA group of 5 (SD 10), 15 subjects per group will "
              "allow detection of a difference of 70 in mean reductions in SFUs.",
    ),
    "NCT05353985": Source(
        citation="Takeda. TAK-062-2001 Statistical Analysis Plan, section 4.0 "
                 "Sample-Size Determination.",
        locator="Section 4.0",
        url="https://cdn.clinicaltrials.gov/large-docs/85/NCT05353985/SAP_001.pdf",
        quote="A sample size of 53 subjects per treatment group will provide 80% "
              "power to detect a standardized mean difference (mean difference/SD) "
              "of 0.55 between TAK-062 and placebo in the change from baseline "
              "CDSD weekly score (assuming common SD).",
    ),
    "NCT04424927": Source(
        citation="Provention/Sanofi. PRV-015-002b Statistical Analysis Plan v2.0, "
                 "section 7 Estimation of Sample Size.",
        locator="Section 7",
        url="https://cdn.clinicaltrials.gov/large-docs/27/NCT04424927/SAP_001.pdf",
        quote="A proposed sample size of approximately 50 evaluable subjects within "
              "each treatment group would provide approximately 80% power to detect "
              "a 0.40 difference from placebo and any given active treatment group "
              "... in the primary endpoint, change through Week 24 in the Abdominal "
              "Symptoms domain score in CeD PRO.",
    ),
}
