"""The gluten challenge is a dose-*and-time* problem. The field designs it as a dose.

Every prevention-design celiac trial has to manufacture the injury it then tries to
prevent. How much injury a challenge delivers is therefore the first design decision,
and it fixes the ceiling on everything downstream: the effect size, the variance (see
`ctsim.variance`), and the sample size. Sponsors pick a dose. They state it in grams per
day, argue about whether it is tolerable, and then choose a duration almost as an
afterthought — 14 days here, 6 weeks there.

Assembled from primary protocols, the published record says the duration is doing more
of the work than the dose.

    3 g/day of gluten protein, three independent studies, one dose:

        14 days  ->  0.06 VH:CD lost   (Leonard 2021, n=7)
        42 days  ->  0.61              (ZED1227 CEC-3 placebo, n=30)
        78 days  ->  1.14              (Lahdeaho 2011, n=21)

That is close to linear in time at roughly 0.015 VH:CD units per day, with no sign of
flattening out to 78 days. Nobody has put those three numbers next to each other,
because they sit in three papers that do not cite one another for this purpose.

Meanwhile the dose lever is much weaker than it looks. Within Lahdeaho's 21 patients
the delivered dose varied more than fourfold (1.3-5.0 g/day) and correlated with injury
at r = -0.14: between-patient variation swamped it entirely.

Two consequences for design, and they point the same way:

1. A longer challenge buys injury roughly in proportion to time, while a harsher one
   saturates. Since required N falls with the square of the standardized effect, and
   `ctsim.variance` says the SD grows only sub-linearly with injury, **duration is the
   cheaper lever**.
2. It is also the kinder one. Lahdeaho lost 7 of 25 patients to symptoms at 3-5 g/day.
   A lower dose for longer delivers the same injury with less to withdraw from.

Units warning, and it is not pedantry. "Gluten" means gluten *protein* in some protocols
and vital wheat gluten *powder* in others, and the two differ by about 30%. TAK-101's
protocol is explicit that a packet holds "approximately 6 g gluten (approximately 8.5 g
of powder)"; Sanofi's ASPIRION states its SIGE dose as 250 mg of *vital wheat gluten*,
which is roughly 190 mg of protein. Every dose in this module is normalised to grams of
gluten **protein** per day, and where a conversion was applied it is recorded.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ctsim.published import Source
from ctsim.simulate import mde
from ctsim.variance import InjuryVarianceModel, fit_injury_variance

# Vital wheat gluten is not pure gluten protein. Two independent primary sources put the
# protein fraction in the same place, so a single constant is defensible:
#   - TAK-101 protocol: "approximately 6 g gluten (approximately 8.5 g of powder)" = 71%
#   - Leonard 2021 methods: "Flour protein fraction was 66% of the total flour content
#     and all protein was gluten (3 g dose = 4.5 g flour)"                       = 67%
# Used only to convert doses that a sponsor stated in powder rather than protein.
VITAL_WHEAT_GLUTEN_PROTEIN_FRACTION = 0.70


@dataclass(frozen=True)
class ChallengeArm:
    """One published gluten-challenge arm, on the raw VH:CD change scale."""

    label: str
    g_protein_per_day: float
    days: int
    injury: float               # |mean change| in VH:CD, positive = mucosa got worse
    sd: float | None            # SD of the change score, where published
    n: int
    baseline_vhcd: float | None
    dose_as_stated: str         # what the sponsor wrote, before normalisation
    source: Source

    @property
    def cumulative_g(self) -> float:
        return self.g_protein_per_day * self.days

    @property
    def injury_per_day(self) -> float:
        return self.injury / self.days


# ---------------------------------------------------------------------------
# The dataset
# ---------------------------------------------------------------------------
# Only arms where a *challenge* was administered and a VH:CD change was published.
# Restoration designs are excluded on purpose: with no challenge there is no dose, and
# their control arms are measuring something else entirely.

CHALLENGE_ARMS: list[ChallengeArm] = [
    ChallengeArm(
        label="Catassi 2007, 10 mg/day",
        g_protein_per_day=0.010, days=90, injury=0.02, sd=None, n=13,
        baseline_vhcd=2.20,
        dose_as_stated="10 mg gluten in a capsule, daily",
        source=Source(
            citation="Catassi C et al. A prospective, double-blind, placebo-controlled "
                     "trial to establish a safe gluten threshold for patients with "
                     "celiac disease. Am J Clin Nutr 2007;85:160-6. PMID 17209192.",
            locator="Results: Vh/Cd change by dose group",
            url="https://pubmed.ncbi.nlm.nih.gov/17209192/",
            quote="patients were assigned to ingest daily for 90 d a capsule containing "
                  "0, 10, or 50 mg gluten. ... baseline Vh/Cd 2.20; 10 mg: -1% change.",
        ),
    ),
    ChallengeArm(
        label="Catassi 2007, 50 mg/day",
        g_protein_per_day=0.050, days=90, injury=0.44, sd=None, n=13,
        baseline_vhcd=2.20,
        dose_as_stated="50 mg gluten in a capsule, daily",
        source=Source(
            citation="Catassi C et al. Am J Clin Nutr 2007;85:160-6. PMID 17209192.",
            locator="Results: Vh/Cd change by dose group",
            url="https://pubmed.ncbi.nlm.nih.gov/17209192/",
            quote="50 mg: -20% decline in Vh/Cd from a baseline of 2.20 over 90 days "
                  "(0.20 x 2.20 = 0.44 in ratio units).",
        ),
    ),
    ChallengeArm(
        label="Leonard 2021, 3 g/day",
        g_protein_per_day=3.0, days=14, injury=0.06, sd=0.516, n=7,
        baseline_vhcd=2.1,
        dose_as_stated="3 g gluten protein (4.5 g Bob's Red Mill flour)",
        source=Source(
            citation="Leonard MM et al. Evaluating Responses to Gluten Challenge: A "
                     "Randomized, Double-Blind, 2-Dose Gluten Challenge Trial. "
                     "Gastroenterology 2021;160:720-733. PMID 33130104. NCT03409796.",
            locator="Methods (gluten preparation); registry results, change from "
                    "baseline in VH:CD, Baseline and Day 15",
            url="https://pmc.ncbi.nlm.nih.gov/articles/PMC7878429/",
            quote="Flour protein fraction was 66% of the total flour content and all "
                  "protein was gluten (3 g dose = 4.5 g flour; 10 g dose = 15 g flour).",
        ),
    ),
    ChallengeArm(
        label="IMGX003 CeliacShield placebo, 2 g/day",
        g_protein_per_day=2.0, days=42, injury=0.35, sd=0.616, n=22,
        baseline_vhcd=2.95,
        dose_as_stated="2 g gluten per day for 6 weeks",
        source=Source(
            citation="Murray JA et al. Latiglutenase Protects the Mucosa and Attenuates "
                     "Symptom Severity in Patients With Celiac Disease Exposed to a "
                     "Gluten Challenge. Gastroenterology 2022;163:1510-1521. "
                     "NCT03585478.",
            locator="Methods (challenge); Table 2, change in Vh:Cd, mITT",
            url="https://pmc.ncbi.nlm.nih.gov/articles/PMC9707643/",
            quote="2 g of gluten per day for 6 weeks. ... Change in Vh:Cd: placebo "
                  "-0.35 (0.616).",
        ),
    ),
    ChallengeArm(
        label="ZED1227 CEC-3 placebo, 3 g/day",
        g_protein_per_day=3.0, days=42, injury=0.61, sd=0.481, n=30,
        baseline_vhcd=1.98,
        dose_as_stated="3 g gluten per day for 6 weeks",
        source=Source(
            citation="Schuppan D et al. A Randomized Trial of a Transglutaminase 2 "
                     "Inhibitor for Celiac Disease. N Engl J Med 2021;385:35-45. "
                     "PMID 34192430. EudraCT 2017-002241-30.",
            locator="Methods (3 g/day, 6 weeks); Table 2 and Figure 2, placebo arm",
            url="https://www.nejm.org/doi/full/10.1056/NEJMoa2032441",
            quote="daily gluten challenge (3 g/day) for 6 weeks; placebo-arm estimated "
                  "change in VH:CD -0.61 (95% CI -0.78 to -0.44).",
        ),
    ),
    ChallengeArm(
        label="TAK-101 Ph2a placebo, 12->6 g/day",
        # 12 g/day for days 15-17 then 6 g/day for days 18-28: 102 g over 14 days.
        g_protein_per_day=102.0 / 14.0, days=14, injury=0.63, sd=0.657, n=15,
        baseline_vhcd=3.01,
        dose_as_stated="12 g gluten x 3 days, then 6 g/day x 11 days",
        source=Source(
            citation="Kelly CP et al. TAK-101 Nanoparticles Induce Gluten-Specific "
                     "Tolerance in Celiac Disease. Gastroenterology 2021;161:66-80. "
                     "PMID 33722583. Protocol posted on NCT03738475, section 5.7.",
            locator="Protocol section 5.7 'Gluten Challenge'; registry results",
            url="https://cdn.clinicaltrials.gov/large-docs/75/NCT03738475/Prot_000.pdf",
            quote="Bob's Red Mill Vital Wheat Gluten (75-80% Protein) powder will be "
                  "used for the gluten challenge. It will be supplied by the sponsor in "
                  "individual foil packets, each containing approximately 6 g gluten "
                  "(approximately 8.5 g of powder) per packet. Subjects will consume "
                  "12 g of gluten for the first 3 days of the gluten challenge (Days "
                  "15-17); 2 packets per day. For the remainder of the gluten challenge "
                  "(Days 18-28), subjects will consume 6 g of gluten per day; 1 packet.",
        ),
    ),
    ChallengeArm(
        label="KAN-101 SynCeD placebo, 9 g/day",
        g_protein_per_day=9.0, days=14, injury=0.61, sd=0.614, n=25,
        baseline_vhcd=None,      # entry required >=2.3; per-arm baseline not posted
        dose_as_stated="9 g/day gluten protein as 12 g vital wheat gluten",
        source=Source(
            citation="Anokion/Kanyos Bio. Protocol KAN-101-03, Final Protocol Amendment "
                     "01, 20 March 2024, posted on NCT06001177.",
            locator="Synopsis, Overall Design (gluten challenge); Key Inclusion "
                    "Criterion 6",
            url="https://cdn.clinicaltrials.gov/large-docs/77/NCT06001177/Prot_000.pdf",
            quote="Participants will undergo a 2-week GC in which they will ingest 9 "
                  "g/day of gluten protein in the form of 12 g vital wheat gluten. ... "
                  "Screening intestinal biopsy demonstrating Vh:Cd ratio of 2.3 or "
                  "higher.",
        ),
    ),
    ChallengeArm(
        label="Lahdeaho 2011, 3.1 g/day mean",
        g_protein_per_day=3.11, days=78, injury=1.14, sd=1.053, n=21,
        baseline_vhcd=2.97,
        dose_as_stated="low (1-3 g) or moderate (3-5 g) gluten daily for 12 weeks",
        source=Source(
            citation="Lahdeaho ML et al. Small-bowel mucosal changes and antibody "
                     "responses after low- and moderate-dose gluten challenge in celiac "
                     "disease. BMC Gastroenterol 2011;11:129. PMID 22115041.",
            locator="Table 2, per-patient Vh/CrD before (I) and after (II)",
            url="https://pmc.ncbi.nlm.nih.gov/articles/PMC3240817/",
            quote="Twenty-five celiac disease adults were challenged with low (1-3 g) "
                  "or moderate (3-5g) doses of gluten daily for 12 weeks. ... the "
                  "gluten challenge lasted a median of 84 days (range 29-103 days) with "
                  "an average of 3.1 g daily gluten consumption (range 1.3-5.0 g/day).",
        ),
    ),
    ChallengeArm(
        label="Leonard 2021, 10 g/day",
        g_protein_per_day=10.0, days=14, injury=1.53, sd=0.941, n=7,
        baseline_vhcd=2.3,
        dose_as_stated="10 g gluten protein (15 g Bob's Red Mill flour)",
        source=Source(
            citation="Leonard MM et al. Gastroenterology 2021;160:720-733. "
                     "PMID 33130104. NCT03409796.",
            locator="Registry results, change from baseline in VH:CD, Baseline and "
                    "Day 15",
            url="https://pmc.ncbi.nlm.nih.gov/articles/PMC7878429/",
            quote="10 g dose = 15 g flour; change from baseline in VH:CD -1.53 (0.941), "
                  "n = 7.",
        ),
    ),
]


# ---------------------------------------------------------------------------
# Individual patient data: the only place dose and duration vary independently
# ---------------------------------------------------------------------------
# Across studies, dose and duration are confounded almost perfectly by design — the
# high-dose studies run 14 days and the trace-dose ones run 90. Nothing can be
# untangled from aggregate arms alone. Lahdeaho's Table 2 is the exception, because
# achieved dose (1.3-5.0 g/day) and achieved duration (29-103 days) both vary *within*
# one protocol, one centre and one reader.
#   (mean g/day, days, Vh/CrD before, Vh/CrD after)
LAHDEAHO_IPD: list[tuple[float, int, float, float]] = [
    (5.0, 29, 2.8, 0.8), (4.9, 84, 3.5, 3.4), (4.9, 38, 2.9, 1.3), (4.7, 45, 2.7, 0.2),
    (4.1, 61, 3.0, 3.5), (4.0, 91, 3.0, 0.6), (3.6, 91, 2.8, 1.4), (3.6, 84, 3.0, 2.6),
    (3.4, 88, 3.8, 2.3), (3.3, 86, 2.7, 0.6), (2.8, 89, 2.5, 3.1), (2.7, 84, 3.0, 0.6),
    (2.6, 81, 2.9, 3.1), (2.4, 85, 2.7, 1.9), (2.2, 103, 4.2, 1.3), (2.1, 93, 1.3, 0.1),
    (2.1, 85, 3.3, 2.4), (2.1, 84, 2.9, 1.7), (2.1, 83, 3.2, 3.4), (1.4, 77, 3.4, 3.0),
    (1.3, 78, 2.5, 0.8),
]

LAHDEAHO_IPD_SOURCE = Source(
    citation="Lahdeaho ML et al. BMC Gastroenterol 2011;11:129. PMID 22115041.",
    locator="Table 2, columns 'Mean daily gluten intake (g)', 'Duration of gluten "
            "challenge (days)', Vh/CrD I and II",
    url="https://pmc.ncbi.nlm.nih.gov/articles/PMC3240817/",
    quote="No | Mean daily gluten intake (g) | Duration of gluten challenge (days) | "
          "Vh/CrD I | Vh/CrD II  ...  1 | 5.0 | 29 | 2.8 | 0.8",
)


@dataclass(frozen=True)
class IPDFit:
    """How much of the injury a challenge causes is explained by how it was designed."""

    r_dose: float
    r_duration: float
    r_cumulative: float
    r2_joint: float
    injury_mean: float
    injury_sd: float
    n: int


def lahdeaho_ipd_fit() -> IPDFit:
    """Regress per-patient injury on the two things a protocol actually controls.

    The answer is the point: almost nothing. Within one protocol, a fourfold spread in
    dose and a threefold spread in duration jointly explain a few percent of who ends
    up with an atrophic mucosa. That is not an argument against a dose-response — the
    across-study numbers show one plainly — it is a statement about the range trials
    operate in, where between-patient variation dominates.
    """
    dose = np.array([r[0] for r in LAHDEAHO_IPD], dtype=float)
    days = np.array([float(r[1]) for r in LAHDEAHO_IPD])
    injury = np.array([r[2] - r[3] for r in LAHDEAHO_IPD])  # positive = deterioration
    cumulative = dose * days

    design = np.column_stack([np.ones_like(dose), dose, days])
    coef, *_ = np.linalg.lstsq(design, injury, rcond=None)
    resid = injury - design @ coef
    ss_tot = float(((injury - injury.mean()) ** 2).sum())
    r2 = 1.0 - float((resid**2).sum()) / ss_tot

    return IPDFit(
        r_dose=float(np.corrcoef(dose, injury)[0, 1]),
        r_duration=float(np.corrcoef(days, injury)[0, 1]),
        r_cumulative=float(np.corrcoef(cumulative, injury)[0, 1]),
        r2_joint=r2,
        injury_mean=float(injury.mean()),
        injury_sd=float(injury.std(ddof=1)),
        n=len(LAHDEAHO_IPD),
    )


# ---------------------------------------------------------------------------
# The one dose held constant across three studies
# ---------------------------------------------------------------------------

# Three studies, run years apart on three continents, all used 3 g/day of gluten
# protein and differed only in how long they ran it. That is the closest thing the
# published record has to a controlled experiment on challenge duration, and it is the
# backbone of everything below.
THREE_GRAM_SERIES = ("Leonard 2021, 3 g/day",
                     "ZED1227 CEC-3 placebo, 3 g/day",
                     "Lahdeaho 2011, 3.1 g/day mean")


@dataclass(frozen=True)
class DurationModel:
    """Injury accumulated per day of challenge at a fixed 3 g/day dose."""

    per_day: float
    per_day_se: float
    intercept: float
    r: float
    n_studies: int

    def injury_at(self, days: float) -> float:
        """Expected |change| in VH:CD after `days` at 3 g/day. Floored at zero."""
        return max(0.0, self.intercept + self.per_day * days)


def fit_duration_model() -> DurationModel:
    """Weighted straight line through the three 3 g/day studies.

    Weighted by sqrt(n), so ZED1227's 30 patients and Lahdeaho's 21 are not outvoted by
    Leonard's 7. Three points is three points — the standard error is reported so the
    fit is not mistaken for more than it is.
    """
    arms = [a for a in CHALLENGE_ARMS if a.label in THREE_GRAM_SERIES]
    assert len(arms) == 3, "the 3 g/day series must have exactly three studies"
    x = np.array([float(a.days) for a in arms])
    y = np.array([a.injury for a in arms])
    w = np.sqrt(np.array([float(a.n) for a in arms]))

    design = np.column_stack([np.ones_like(x), x]) * w[:, None]
    coef, *_ = np.linalg.lstsq(design, y * w, rcond=None)
    intercept, slope = float(coef[0]), float(coef[1])

    resid = (y - (intercept + slope * x)) * w
    dof = max(1, len(x) - 2)
    s2 = float((resid**2).sum()) / dof
    cov = s2 * np.linalg.inv(design.T @ design)
    slope_se = float(np.sqrt(cov[1, 1]))

    return DurationModel(
        per_day=slope, per_day_se=slope_se, intercept=intercept,
        r=float(np.corrcoef(x, y)[0, 1]), n_studies=len(arms),
    )


# ---------------------------------------------------------------------------
# Composing it: what a proposed challenge design costs in patients
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ChallengeDesign:
    days: int
    injury: float
    sd: float
    mde_80: float
    n_per_arm: int
    protection_needed: float


def required_n(
    days: int,
    protection: float = 0.50,
    model: InjuryVarianceModel | None = None,
    duration: DurationModel | None = None,
    max_n: int = 5000,
) -> ChallengeDesign:
    """Patients per arm to detect a drug preventing `protection` of a 3 g/day challenge.

    Chains the two empirical models this repo has: duration sets the injury, injury sets
    the SD, and the SD sets the sample size. Every link is fitted to published arms, and
    every link is the weakest part of the chain in a different way — see `report`.
    """
    model = fit_injury_variance() if model is None else model
    duration = fit_duration_model() if duration is None else duration

    injury = duration.injury_at(days)
    sd = model.sd_at(injury)
    effect = protection * injury

    n = 2
    while n < max_n and mde(n, sd) > effect:
        n += 1
    return ChallengeDesign(
        days=days, injury=injury, sd=sd, mde_80=mde(n, sd), n_per_arm=n,
        protection_needed=mde(n, sd) / injury if injury > 0 else float("inf"),
    )


# ---------------------------------------------------------------------------
# The SIGE contradiction
# ---------------------------------------------------------------------------
# Two live trials are betting that a few hundred mg of gluten a day for three months
# leaves the mucosa alone. Takeda's own protocol says the opposite is established, and
# cites the paper that established it. Both cannot be right, and both are checkable.

SIGE_DESIGN_INTENT = Source(
    citation="Takeda. Study TAK-062-2001 protocol incorporating Amendment No. 5, "
             "16 November 2023, posted on NCT05353985.",
    locator="Study summary / section 6.1 Study Design; section 4.1 Background",
    url="https://cdn.clinicaltrials.gov/large-docs/85/NCT05353985/Prot_001.pdf",
    quote="It is expected that SIGE will result in stable, rather than worsening, "
          "symptoms and enteropathy in the TAK-062 placebo groups. ... The threshold of "
          "daily gluten that will cause mucosal injury in both adults and children is "
          "thought to be 10 to 50 mg/day, or about one-hundredth of a slice of bread "
          "(Catassi et al. 2007; Gibert et al. 2013).",
)

# Doses normalised to gluten protein. Sanofi states vital wheat gluten, so it converts.
#   (label, g protein/day, days, entry criterion, what the sponsor wrote)
SIGE_ARMS = [
    ("TAK-062 placebo (observed)", None, 168, "Vh:Cd <2.5",
     "gluten content of the SIGE bar is redacted throughout the protocol"),
    ("Sanofi ASPIRION SIGE", 0.250 * VITAL_WHEAT_GLUTEN_PROTEIN_FRACTION, 84,
     "Vh:Cd <2.5", "approximately 250 mg of vital wheat gluten daily for 12 weeks"),
    ("Dr Falk CEC-013 SIGE", 0.214, 105, "Vh:CrD <=2.5",
     "approximately 500 mg three times a week, not more than 2.2 g per week"),
]

# The three long-challenge arms, all on a *healed* mucosa, that between them span
# 0.01 to 3.11 g/day at 78-90 days. This is the only dose-response the record supports
# at a duration comparable to a SIGE regimen.
LONG_CHALLENGE_POINTS = ("Catassi 2007, 10 mg/day",
                         "Catassi 2007, 50 mg/day",
                         "Lahdeaho 2011, 3.1 g/day mean")

TAK062_ENTRY = Source(
    citation="Takeda. Study TAK-062-2001 protocol incorporating Amendment No. 5, "
             "16 November 2023, posted on NCT05353985.",
    locator="Inclusion criterion 7; randomisation stratification",
    url="https://cdn.clinicaltrials.gov/large-docs/85/NCT05353985/Prot_001.pdf",
    quote="The subject has small intestinal villous atrophy on duodenal biopsy defined "
          "as Vh:Cd <2.5 at Week -4. ... stratified by ... mild to moderate (Vh:Cd 1.5 "
          "to <2.5) versus moderate to severe (Vh:Cd <1.5) histologic injury at "
          "baseline.",
)


def fit_long_challenge_dose_response() -> tuple[float, float, float]:
    """Hill curve through the three ~85-day arms: (Imax, D50 in g/day, hill slope).

    Three points and three parameters, so this interpolates rather than fits — there is
    no residual and no test of shape. It is used only *between* its anchors (0.01 to
    3.11 g/day), which is where the SIGE doses fall, and never beyond them.

    A plain saturating curve cannot pass through these points: 10 mg gives 0.02 and
    50 mg gives 0.44, so a fivefold dose rise multiplies injury twentyfold. That is a
    threshold, and it is why the Hill slope is free rather than fixed at 1.
    """
    from scipy.optimize import fsolve

    arms = [next(a for a in CHALLENGE_ARMS if a.label == lbl)
            for lbl in LONG_CHALLENGE_POINTS]
    d = np.array([a.g_protein_per_day for a in arms])
    y = np.array([a.injury for a in arms])

    def residuals(p):
        imax, log_d50, h = p
        d50 = math.exp(log_d50)
        return imax * d**h / (d50**h + d**h) - y

    imax, log_d50, h = fsolve(residuals, [1.3, math.log(0.05), 1.5], full_output=False)
    return float(imax), float(math.exp(log_d50)), float(h)


def sige_injury_on_healed_mucosa(g_per_day: float, days: int) -> float:
    """Injury a SIGE dose would cause in a patient who still had villi to lose.

    Two caveats that matter more than the number. The dose-response is anchored at
    78-90 days, so the duration scaling here is a straight proportion off ~85 days and
    is the weakest part. And every study behind the curve enrolled *healed* patients —
    which is exactly the population no SIGE trial enrols.
    """
    imax, d50, h = fit_long_challenge_dose_response()
    at_reference = imax * g_per_day**h / (d50**h + g_per_day**h)
    return at_reference * (days / 85.0)


def report() -> str:
    dur = fit_duration_model()
    ipd = lahdeaho_ipd_fit()
    var = fit_injury_variance()

    lines = [
        "CHALLENGE: how much injury a gluten challenge actually delivers",
        "=" * 78,
        "  Doses normalised to grams of gluten PROTEIN per day. Sponsors state some in",
        "  vital wheat gluten powder, which is about 30% heavier for the same protein.",
        "",
        f"  {'arm':34}{'g/day':>8}{'days':>6}{'cumul':>8}{'injury':>8}{'n':>5}",
    ]
    for a in sorted(CHALLENGE_ARMS, key=lambda a: a.injury):
        lines.append(
            f"  {a.label[:32]:34}{a.g_protein_per_day:>8.2f}{a.days:>6}"
            f"{a.cumulative_g:>8.1f}{a.injury:>8.2f}{a.n:>5}"
        )

    lines += [
        "",
        "-" * 78,
        "One dose, three studies, three durations",
        "-" * 78,
        "  Across studies, dose and duration are confounded by design: the high-dose",
        "  studies run 14 days and the trace-dose ones run 90. The exception is 3 g/day,",
        "  which three independent studies used at three different durations.",
        "",
        f"  {'study':34}{'days':>6}{'injury':>8}{'per day':>10}{'n':>5}",
    ]
    for label in THREE_GRAM_SERIES:
        a = next(x for x in CHALLENGE_ARMS if x.label == label)
        lines.append(
            f"  {a.label[:32]:34}{a.days:>6}{a.injury:>8.2f}"
            f"{a.injury_per_day:>10.4f}{a.n:>5}"
        )
    lines += [
        "",
        f"  injury = {dur.intercept:+.3f} + {dur.per_day:.5f} (+/-{dur.per_day_se:.5f})"
        f" x days     r = {dur.r:.3f}",
        "",
        "  Roughly 0.015 VH:CD units lost per day of challenge, and still climbing at",
        "  78 days. Three points cannot rule out curvature; they do rule out the",
        "  saturation a 14-day challenge would need for it to be a sensible default.",
        "",
        "-" * 78,
        "Meanwhile the dose lever, measured where it can be measured",
        "-" * 78,
        f"  Lahdeaho 2011, {ipd.n} patients, one protocol, one reader.",
        "  Achieved dose spans 1.3-5.0 g/day and duration 29-103 days.",
        "",
        f"    corr(injury, dose)        = {ipd.r_dose:+.3f}",
        f"    corr(injury, duration)    = {ipd.r_duration:+.3f}",
        f"    corr(injury, cumulative)  = {ipd.r_cumulative:+.3f}",
        f"    R^2, dose + duration      = {ipd.r2_joint:.3f}",
        "",
        f"  Mean injury {ipd.injury_mean:.2f}, SD {ipd.injury_sd:.2f}. Design explains"
        f" {ipd.r2_joint:.0%} of who ends up",
        "  atrophic; the patient explains the rest. A sponsor choosing between 3 and 5",
        "  g/day is tuning a knob that is not connected to much.",
        "",
        "-" * 78,
        "What that costs in patients, per arm, at 3 g/day, for a drug preventing 50%",
        "-" * 78,
        f"  Injury from the duration model; SD from the injury model"
        f" ({var.floor:.3f} + {var.slope:.3f} x injury).",
        "  Rows past 78 days are extrapolation — the longest challenge on record.",
        "",
        f"  {'challenge':>12}{'injury':>9}{'SD':>8}{'MDE(80%)':>11}{'N/arm':>8}",
    ]
    for weeks in (2, 4, 6, 8, 12, 16):
        d = required_n(weeks * 7)
        flag = "  (extrapolated)" if weeks * 7 > 78 else ""
        lines.append(
            f"  {str(weeks) + ' weeks':>12}{d.injury:>9.2f}{d.sd:>8.3f}"
            f"{d.mde_80:>11.3f}{d.n_per_arm:>8}{flag}"
        )
    six, twelve = required_n(42), required_n(84)
    lines += [
        "",
        f"  Doubling a 6-week challenge to 12 weeks takes N/arm from {six.n_per_arm} to"
        f" {twelve.n_per_arm}",
        f"  ({1 - twelve.n_per_arm / six.n_per_arm:.0%} fewer patients) at the same"
        " 3 g/day dose. Lahdeaho lost 7 of 25",
        "  patients to symptoms at 3-5 g/day, so raising the dose instead costs",
        "  dropouts that a longer, gentler challenge does not.",
        "",
        "-" * 78,
        "The SIGE arms are betting against a paper their own sponsor cites",
        "-" * 78,
        "  Takeda's protocol states the design intent and the threshold literature in",
        "  the same document:",
        "",
        f"    \"{SIGE_DESIGN_INTENT.quote[:72]}...\"",
        "",
        "  TAK-062's SIGE placebo arm then moved -0.006 over 24 weeks: dead flat, exactly",
        "  as designed. Yet the only dose-response measured over a comparable duration",
        "  says a few hundred mg a day is not a small exposure:",
        "",
        f"  {'arm':30}{'mg protein/d':>14}{'days':>6}{'entry':>14}"
        f"{'% of max':>10}{'implied':>9}",
    ]
    imax, d50, h = fit_long_challenge_dose_response()
    for label, g, days, entry, _stated in SIGE_ARMS:
        if g is None:
            lines.append(f"  {label[:28]:30}{'redacted':>14}{days:>6}{entry:>14}"
                         f"{'n/a':>10}{'n/a':>9}")
            continue
        frac = g**h / (d50**h + g**h)
        lines.append(
            f"  {label[:28]:30}{g * 1000:>14.0f}{days:>6}{entry:>14}"
            f"{frac:>10.0%}{sige_injury_on_healed_mucosa(g, days):>9.2f}"
        )

    lines += [
        "",
        "  Interpolated on a Hill curve through the three long-challenge arms:",
        f"  Imax {imax:.2f}, D50 {d50 * 1000:.0f} mg/day, slope {h:.2f}. Threshold-shaped,",
        "  because 10 mg gives 0.02 and 50 mg gives 0.44 — a fivefold dose for twenty",
        "  times the injury.",
        "",
        "  Read the '% of max' column twice. On this curve the half-maximal dose is",
        f"  {d50 * 1000:.0f} mg/day, so both SIGE regimens sit near the top of the"
        " dose-response, not",
        "  near the bottom. Whatever else SIGE is, on a mucosa with villi to lose it is",
        "  not a trace exposure. 'Simulated inadvertent' describes the intent behind the",
        "  dose, not its size.",
        "",
        "  But look at the entry column, because it is the whole argument. Every study",
        "  behind that curve challenged a *healed* mucosa: Catassi's patients started at",
        "  2.20, Lahdeaho's at 2.97. All three SIGE trials enrol patients who are already",
        "  atrophic, Vh:Cd under 2.5, and TAK-062 stratified its own enrolment at 1.5.",
        "",
        "  So TAK-062's flat control arm does not show that SIGE spares the mucosa. It is",
        "  equally consistent with a mucosa that had already fallen as far as it goes.",
        "  Those two readings are not distinguishable from anything published, and they",
        "  imply opposite things for the two trials now running. Nobody has separated",
        "  them, and nobody needs new data to try — TAK-062's own baseline strata would",
        "  do it, since the <1.5 group has less room to fall than the 1.5-2.5 group.",
        "",
        "  It matters for sample size either way. A restoration design assumes its",
        "  control arm sits at the variance model's floor. If SIGE injures at all, the",
        "  SD rises with it:",
        "",
        f"  {'assumed control injury':>24}{'SD':>8}{'MDE at n=34/arm':>18}",
    ]
    for inj in (0.0, 0.2, 0.4, 0.6):
        sd = var.sd_at(inj)
        lines.append(f"  {inj:>24.2f}{sd:>8.3f}{mde(34, sd):>18.3f}")
    lines += [
        "",
        "  Sanofi's ASPIRION reads a 34-vs-34 contrast. At a flat control arm it needs",
        "  a healing difference of 0.27; if its SIGE injures by 0.4 it needs 0.35. The",
        "  largest such difference any restoration-design celiac trial has produced is",
        "  0.14.",
        "",
        "-" * 78,
        "What would sink this",
        "-" * 78,
        "  1. The duration line is three points at one dose. r = 0.997 on three points",
        "     is not evidence of linearity, it is what three points do. A fourth study",
        "     at 3 g/day and any duration would be worth more than everything else here.",
        "  2. Those three studies differ in centre, decade, reader and population, so",
        "     the 'one dose held constant' framing controls for dose and nothing else.",
        "     Leonard's 14-day point is n = 7 and its CI spans 0.06 +/- 0.38, which is",
        "     wide enough to sit on a curve as easily as a line.",
        "  3. Injury is bounded below by zero villous height, so it cannot grow linearly",
        "     forever. On a baseline near 3.0 the ceiling is around 2.5. Extrapolating",
        "     past ~12 weeks runs into that whatever the fit says.",
        "  4. The dose-response is anchored on three arms and interpolated with three",
        "     free parameters, so it has no residual and no goodness of fit. It is a",
        "     shape imposed on three numbers, not tested against them.",
        "  5. Lahdeaho's dose was *achieved* intake, not assigned: patients who felt",
        "     worse ate less and stopped sooner. That biases the within-study dose",
        "     effect toward zero, which is the direction of the claim being made here.",
        "     Three of its ten moderate-dose patients withdrew early for symptoms.",
        "  6. Nothing here says a longer challenge is safe. It says it is more",
        "     detectable. Those are different questions and the second one is not a",
        "     licence to answer the first.",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
