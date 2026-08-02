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
]
