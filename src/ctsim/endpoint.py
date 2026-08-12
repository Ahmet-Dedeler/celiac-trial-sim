"""Why the composite endpoint sometimes beats VH:CD and sometimes doesn't.

The literature contradicts itself on this and nobody has reconciled it.

    Syage et al. (CGH 2024) re-scored two trials on the same patients and the same
    biopsies and found the VCIEL composite beat both of its components in both:
    ALV003-1021 1.86 against 1.37 and 1.17, IMGX003 1.14 against 0.76 and 0.98.

    Takeda (UEG Week 2025, abstract MP739) reports of TAK-062 that "no endpoint
    outperformed Vh:Cd".

Both can be right, and the arithmetic says exactly when each one is.

VCIEL standardises each measure by its own SD and adds them with the sign flipped, so
the two injuries point the same way:

    VCIEL = (Vh:Cd - <Vh:Cd>)/sigma_v  -  (IEL - <IEL>)/sigma_i

Averaging two measurements beats either one only if both carry signal. If the standardized
change SDs are similar and the two measures are uncorrelated, the composite's standardized
effect is

    d_c = (d_v + d_i) / sqrt(2)

so the composite beats VH:CD alone exactly when

    d_i  >  (sqrt(2) - 1) * d_v   ~=  0.414 * d_v

IEL has to carry at least 41% of VH:CD's standardized signal to be worth adding. Below
that, the sqrt(2) the composite pays in noise costs more than the extra signal is worth.

The uncorrelated assumption is not an assumption. Syage measured it: R^2 between Vh:Cd and
IEL is 0.005-0.023, so |rho| is 0.07-0.15. Working backwards from their own reported effect
sizes gives an implied rho of -0.07 and +0.16. Two routes, same answer, near zero.

That threshold is what reconciles the two claims. A gluten challenge drives both endpoints,
so d_i/d_v lands well above 0.414 and the composite wins — which is the regime Syage
scored. A restoration design against SIGE is built so that the control arm does not
deteriorate, so neither endpoint has an injury signal to carry and the ratio collapses —
which is the regime Takeda scored. The two groups measured different designs and reported
the answers for their own.

The prediction, and it is testable rather than rhetorical: TAK-062's d_i/d_v should sit
below 0.414. Takeda measured the IELs and published them only qualitatively, so the one
number that would settle a public disagreement in the field is the one that is missing.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ctsim.published import Source

# Composite of two standardized, uncorrelated measures with similar change variability.
# Anything below this and adding IEL costs more noise than it contributes signal.
BREAK_EVEN_RATIO = math.sqrt(2.0) - 1.0


@dataclass(frozen=True)
class ArmPair:
    """One trial's placebo and drug arms on one endpoint."""

    placebo_delta: float
    placebo_sd: float
    placebo_n: int
    drug_delta: float
    drug_sd: float
    drug_n: int

    def cohens_d(self) -> float:
        num = ((self.placebo_n - 1) * self.placebo_sd**2
               + (self.drug_n - 1) * self.drug_sd**2)
        pooled = math.sqrt(num / (self.placebo_n + self.drug_n - 2))
        return abs(self.placebo_delta - self.drug_delta) / pooled


@dataclass(frozen=True)
class EndpointComparison:
    trial: str
    design: str
    vhcd: ArmPair
    iel: ArmPair
    source: Source

    @property
    def d_vhcd(self) -> float:
        return self.vhcd.cohens_d()

    @property
    def d_iel(self) -> float:
        return self.iel.cohens_d()

    @property
    def ratio(self) -> float:
        return self.d_iel / self.d_vhcd

    @property
    def d_composite(self) -> float:
        """Standardized effect of the composite, assuming rho = 0."""
        return (self.d_vhcd + self.d_iel) / math.sqrt(2.0)

    @property
    def composite_helps(self) -> bool:
        return self.ratio > BREAK_EVEN_RATIO

    @property
    def n_ratio(self) -> float:
        """Sample size on the composite as a fraction of sample size on VH:CD alone."""
        return (self.d_vhcd / self.d_composite) ** 2


# Only trials that posted a genuine mean and SD for the change on *both* endpoints in the
# same arms. That is a short list, and it is why this question has stayed open.
COMPARISONS = [
    EndpointComparison(
        trial="IMGX003 CeliacShield",
        design="prevention (2 g/day x 6 weeks)",
        vhcd=ArmPair(-0.35, 0.616, 22, -0.04, 0.466, 21),
        iel=ArmPair(24.75, 22.941, 22, 9.79, 15.741, 21),
        source=Source(
            citation="Murray JA et al. Gastroenterology 2022;163:1510-1521. "
                     "NCT03585478.",
            locator="Table 2, 'Histologic Efficacy Analysis - mITT Population', change "
                    "in Vh:Cd and change in IEL, mean (SD)",
            url="https://pmc.ncbi.nlm.nih.gov/articles/PMC9707643/",
            quote="Change in Vh:Cd: placebo -0.35 (0.616); IMGX003 -0.04 (0.466). "
                  "Change in IEL: placebo 24.75 (22.941); IMGX003 9.79 (15.741).",
        ),
    ),
    EndpointComparison(
        trial="TAK-101 Phase 2a",
        design="prevention (12->6 g/day x 14 days)",
        vhcd=ArmPair(-0.63, 0.657, 15, -0.18, 0.381, 13),
        iel=ArmPair(35.0, 18.756, 15, 28.62, 19.241, 13),
        source=Source(
            citation="Kelly CP et al. Gastroenterology 2021;161:66-80. Registry results "
                     "for NCT03738475.",
            locator="Outcome measures: change from baseline in Vh:Cd and in IEL count, "
                    "Baseline (Screening) to Day 29",
            url="https://clinicaltrials.gov/study/NCT03738475",
            quote="Change in Vh:Cd: placebo -0.63 (0.657); TIMP-GLIA -0.18 (0.381). "
                  "Change in IELs: placebo 35.00 (18.756); TIMP-GLIA 28.62 (19.241).",
        ),
    ),
]


# Syage's own numbers, to check this arithmetic against something published rather than
# only against itself. His effect sizes standardise by a baseline SD where the ones above
# use a pooled change SD, so the levels differ; the ratios are what must agree.
#   (trial, d_vhcd, d_iel, d_vciel)
SYAGE_REPORTED = [
    ("ALV003-1021", 1.37, 1.17, 1.86),
    ("IMGX003 CeliacShield", 0.76, 0.98, 1.14),
]

SYAGE_SOURCE = Source(
    citation="Syage JA et al. A Composite Morphometric Duodenal Biopsy Mucosal Scale "
             "for Celiac Disease Encompassing Both Morphology and Inflammation. "
             "Clin Gastroenterol Hepatol 2024;22:1238-1244.",
    locator="Results: effect sizes by endpoint; baseline correlation of Vh:Cd with IEL",
    url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12213069/",
    quote="The formulation of the VCIEL composite histologic scale was based on "
          "combining the Vh:Cd and IEL measurements for individual patients with equal "
          "weighting, by converting each scale to a fraction of their standard deviation "
          "and summing the results. ... low baseline correlation (R2 = 0.005-0.023).",
)

TAKEDA_SOURCE = Source(
    citation="Maxwell, Isola, Valimaki, Robert, Cheng, Zarei, Leffler. UEG Week 2025, "
             "abstract MP739 (TAK-062).",
    locator="Results",
    url="https://gutflix.eu/search/d/47fdfd28-a813-11f0-bca4-0242ac140006",
    quote="no endpoint outperformed Vh:Cd",
)

# What would settle it, and is not published.
NOT_AVAILABLE = {
    "TAK-062 (NCT05353985) IEL change by arm":
        "The trial behind the 'no endpoint outperformed Vh:Cd' claim measured IELs and "
        "reported them only qualitatively — no per-arm mean or SD in the registry, the "
        "abstract, or the protocol. So the one design that would test the threshold "
        "is the one whose numbers are missing. It is a restoration design against "
        "SIGE, which is exactly the regime the threshold predicts the composite should "
        "*lose* in.",
    "ALV003-1021 (Lahdeaho 2014) dispersion":
        "Syage scores it at 1.37/1.17/1.86 but the underlying paper publishes no SD, SE "
        "or CI for VH:CrD, so his numbers cannot be reproduced from source.",
}


def implied_correlation(d_v: float, d_i: float, d_c: float) -> float:
    """Back out the Vh:Cd-IEL correlation implied by a reported composite effect size.

    A cross-check on Syage rather than a new measurement: if the composite really is an
    equally weighted sum of two standardized measures, its effect size pins the
    correlation. It should land near the R^2 = 0.005-0.023 he reports separately, and it
    does. If these two ever disagreed badly it would mean the composite is not doing what
    its definition says.
    """
    mean_component = (d_v + d_i) / 2.0
    inflation = d_c / mean_component      # = sqrt(2 / (1 + rho))
    return 2.0 / inflation**2 - 1.0


def report() -> str:
    lines = [
        "ENDPOINT: when the VCIEL composite is worth using, and when it is not",
        "=" * 78,
        "  Syage says the composite beats both its components. Takeda says nothing beat",
        "  Vh:Cd. Both are right, in different designs, and the crossover is arithmetic.",
        "",
        "  A composite of two uncorrelated measures pays sqrt(2) in noise to gain the",
        "  second measure's signal, so it is worth adding IEL only when",
        "",
        f"      d_IEL  >  (sqrt(2) - 1) x d_VH:CD   =   {BREAK_EVEN_RATIO:.3f} x d_VH:CD",
        "",
        "-" * 78,
        "Trials that posted a mean and an SD for both endpoints in the same arms",
        "-" * 78,
        f"  {'trial':24}{'d VH:CD':>9}{'d IEL':>8}{'ratio':>8}"
        f"{'d VCIEL':>9}{'N vs Vh:Cd':>12}",
    ]
    for c in COMPARISONS:
        lines.append(
            f"  {c.trial[:22]:24}{c.d_vhcd:>9.3f}{c.d_iel:>8.3f}{c.ratio:>8.2f}"
            f"{c.d_composite:>9.3f}{c.n_ratio:>11.0%}"
        )
        lines.append(f"  {'':24}{c.design}")

    lines += [
        "",
        "  IMGX003 clears the threshold three times over, and the composite would have",
        "  cut its sample size by a third. That is the trial that missed on Vh:Cd at",
        "  p = 0.057 and hit on IEL at p = 0.018.",
        "",
        f"  TAK-101 sits on the threshold: {COMPARISONS[1].ratio:.2f} against"
        f" {BREAK_EVEN_RATIO:.3f}. The composite would have",
        "  done nothing for it. Worth stating plainly, because TAK-101 is the trial this",
        "  repo argues was underpowered, and the obvious fix is not the fix. It needed",
        "  patients, not a better endpoint.",
        "",
        "-" * 78,
        "Checking the arithmetic against Syage's own published effect sizes",
        "-" * 78,
        "  He standardises by a baseline SD and the rows above use a pooled change SD, so",
        "  the levels are not comparable. The ratios are.",
        "",
        f"  {'trial':24}{'d VH:CD':>9}{'d IEL':>8}{'ratio':>8}{'d VCIEL':>9}"
        f"{'implied rho':>13}",
    ]
    for trial, dv, di, dc in SYAGE_REPORTED:
        lines.append(
            f"  {trial[:22]:24}{dv:>9.2f}{di:>8.2f}{di / dv:>8.2f}{dc:>9.2f}"
            f"{implied_correlation(dv, di, dc):>13.2f}"
        )

    own = next(c for c in COMPARISONS if c.trial == "IMGX003 CeliacShield")
    syage_imgx = next(r for r in SYAGE_REPORTED if r[0] == "IMGX003 CeliacShield")
    lines += [
        "",
        f"  IMGX003 computed from Table 2 here: ratio {own.ratio:.2f}."
        f" Syage, independently: {syage_imgx[2] / syage_imgx[1]:.2f}.",
        "  The implied correlations sit either side of zero, against the R^2 of",
        "  0.005-0.023 he measures directly. Two routes, one answer: the two endpoints",
        "  are close to independent, which is the whole reason the composite gains.",
        "",
        "-" * 78,
        "What would settle the disagreement, and why it cannot be settled yet",
        "-" * 78,
    ]
    for what, why in NOT_AVAILABLE.items():
        lines.append(f"  {what}")
        for chunk in _wrap(why, 74):
            lines.append(f"    {chunk}")
        lines.append("")
    lines += [
        "  So the reconciliation offered here is a prediction, not a demonstration:",
        "  TAK-062's d_IEL/d_VH:CD should fall below 0.414, because a design built to",
        "  keep its control arm flat has removed the injury both endpoints measure.",
        "  Takeda has the number. It is one table.",
    ]
    return "\n".join(lines)


def _wrap(text: str, width: int) -> list[str]:
    words, line, out = text.split(), "", []
    for w in words:
        if len(line) + len(w) + 1 > width:
            out.append(line)
            line = w
        else:
            line = f"{line} {w}".strip()
    if line:
        out.append(line)
    return out


if __name__ == "__main__":
    print(report())
