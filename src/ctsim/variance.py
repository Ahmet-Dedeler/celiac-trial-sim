"""The noise is not a constant. It grows with the injury.

This repo used to quote one pooled SD for ΔVH:CD and apply it to every design. Once the
dataset was corrected and widened, a homogeneity test rejected that: the per-arm
variances are not draws from one common variance (Bartlett p < 0.001).

The pattern behind the rejection is simple and, in hindsight, obvious. Arms where the
mucosa barely moved are tight; arms that took a heavy gluten hit are spread out. Plot
SD against |mean change| across every published arm and it is close to a straight line:

    SD  ~=  floor  +  slope * |injury|

That makes sense mechanistically. The floor is what you measure when nothing happens —
biopsy siting, orientation, reading. The slope is patient heterogeneity in *response*:
if the average patient loses 1.5 of villous height, patients necessarily differ in how
much they lose, and that difference is proportional to the average.

It also changes a design conclusion this repo previously got wrong. Hitting patients
with more gluten raises the signal, but it raises the noise too, so the gain from a
stronger challenge saturates instead of compounding. The earlier claim — that a 10 g
challenge cuts a trial from 93 patients per arm to 15 — used one constant SD for both
and overstated the lever by roughly a factor of three.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import stats

from ctsim.model import EmpiricalSD, load_empirical


def all_arms(include_published: bool = True) -> list[EmpiricalSD]:
    """Every ΔVH:CD arm on the raw-ratio change scale, registry plus papers."""
    from ctsim import published

    rows = list(load_empirical())
    if include_published:
        rows += published.paper_rows()
    return rows


def _df_of(r: EmpiricalSD) -> int:
    return (r.n_arm - 1) if r.n_arm and r.n_arm > 1 else 1


@dataclass
class InjuryVarianceModel:
    floor: float          # SD when nothing happens to the mucosa
    slope: float          # extra SD per unit of injury
    floor_se: float
    slope_se: float
    r: float              # correlation of |injury| with SD
    n_arms: int
    n_trials: int

    def sd_at(self, injury: float) -> float:
        """Predicted between-patient SD of the change score at a given mean injury."""
        return self.floor + self.slope * abs(injury)

    def describe(self) -> str:
        return (
            f"SD = {self.floor:.3f} (+/-{self.floor_se:.3f})"
            f" + {self.slope:.3f} (+/-{self.slope_se:.3f}) * |injury|\n"
            f"  fitted on {self.n_arms} arms from {self.n_trials} trials; "
            f"r = {self.r:.3f}"
        )


def fit_injury_variance(rows: list[EmpiricalSD] | None = None) -> InjuryVarianceModel:
    """Weighted least-squares fit of arm SD against arm |mean change|.

    Weighted by sqrt(df): the sampling SE of an SD estimate is about sd/sqrt(2*df), so
    a 60-patient arm should not be outvoted by a 7-patient one.
    """
    rows = all_arms() if rows is None else rows
    x = np.array([abs(r.delta) for r in rows], dtype=float)
    y = np.array([r.sd for r in rows], dtype=float)
    w = np.sqrt(np.array([_df_of(r) for r in rows], dtype=float))

    slope, floor = np.polyfit(x, y, 1, w=w)

    # Standard errors on the two coefficients, from the weighted normal equations.
    X = np.column_stack([np.ones_like(x), x])
    W = np.diag(w**2)
    xtwx_inv = np.linalg.inv(X.T @ W @ X)
    resid = y - (floor + slope * x)
    dof = max(len(rows) - 2, 1)
    s2 = float((w**2 * resid**2).sum() / dof)
    se = np.sqrt(np.diag(xtwx_inv) * s2)

    return InjuryVarianceModel(
        floor=float(floor), slope=float(slope),
        floor_se=float(se[0]), slope_se=float(se[1]),
        r=float(np.corrcoef(x, y)[0, 1]),
        n_arms=len(rows), n_trials=len({r.nct_id for r in rows}),
    )


# ---------------------------------------------------------------------------
# What the model does to trial design
# ---------------------------------------------------------------------------

@dataclass
class ChallengeDesign:
    injury: float
    protection: float
    effect: float
    sd: float
    n_per_arm: int
    n_per_arm_constant_sd: int  # what the old constant-SD answer would have said


def required_n(injury: float, protection: float,
               model: InjuryVarianceModel | None = None,
               constant_sd: float | None = None,
               power: float = 0.80, alpha: float = 0.05) -> ChallengeDesign:
    """Patients per arm to detect a drug that blocks `protection` of `injury`.

    The point of returning both numbers is that the difference between them is the
    entire correction: a constant-SD calculation gets steadily more optimistic as the
    challenge gets harsher, because it credits the extra signal without charging for
    the extra noise.
    """
    model = fit_injury_variance() if model is None else model
    k = (stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(power)) ** 2
    effect = protection * abs(injury)
    sd = model.sd_at(injury)
    const = constant_sd if constant_sd is not None else model.sd_at(0.61)
    return ChallengeDesign(
        injury=abs(injury), protection=protection, effect=effect, sd=sd,
        n_per_arm=math.ceil(2 * k * sd**2 / effect**2),
        n_per_arm_constant_sd=math.ceil(2 * k * const**2 / effect**2),
    )


def report(rows: list[EmpiricalSD] | None = None) -> str:
    rows = all_arms() if rows is None else rows
    m = fit_injury_variance(rows)

    lines = [
        "Noise scales with injury — one pooled SD is the wrong summary",
        "-" * 78,
        "  " + m.describe().replace("\n", "\n  "),
        "",
        f"  {'|injury|':>9}{'SD obs':>9}{'SD fit':>9}{'n':>5}  trial / arm",
    ]
    for r in sorted(rows, key=lambda r: abs(r.delta)):
        lines.append(
            f"  {abs(r.delta):>9.3f}{r.sd:>9.3f}{m.sd_at(r.delta):>9.3f}"
            f"{r.n_arm or 0:>5}  {r.nct_id} {r.arm_label[:28]}"
        )

    lines += [
        "",
        "  The floor is what the assay costs you when nothing happens: biopsy siting,",
        "  orientation and reading. The slope is patient heterogeneity in response.",
        "",
        "-" * 78,
        "Consequence: a harsher gluten challenge helps, but it saturates",
        "-" * 78,
        "  patients per arm to detect a drug blocking half the injury:",
        "",
        (f"  {'challenge injury':>17}{'SD':>8}{'N/arm':>8}"
         f"{'N/arm if SD were constant':>28}"),
    ]
    for injury in (0.20, 0.61, 1.00, 1.53, 2.50):
        d = required_n(injury, 0.50, m)
        lines.append(
            f"  {d.injury:>17.2f}{d.sd:>8.3f}{d.n_per_arm:>8}"
            f"{d.n_per_arm_constant_sd:>28}"
        )
    lines += [
        "",
        "  The right-hand column is what this repo used to report. It is too optimistic",
        "  everywhere the challenge is harsher than the one it calibrated on, because a",
        "  bigger injury buys a bigger spread of injuries along with the bigger mean.",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
