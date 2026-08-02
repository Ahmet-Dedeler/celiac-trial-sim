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


# ---------------------------------------------------------------------------
# Holdout: does the model predict trials it has never seen?
# ---------------------------------------------------------------------------

# The two headline challenge trials. Fit without them; ask whether the injury
# model still recovers their placebo-arm SDs. If it does, the relationship is not
# just a within-dataset correlation.
HEADLINE_HOLDOUT = ("NCT06001177", "NCT03738475")


@dataclass
class HoldoutPrediction:
    nct_id: str
    arm_label: str
    injury: float
    sd_observed: float
    sd_predicted: float
    n_arm: int

    @property
    def abs_err(self) -> float:
        return abs(self.sd_observed - self.sd_predicted)


@dataclass
class HoldoutResult:
    held_out: tuple[str, ...]
    model: InjuryVarianceModel
    predictions: list[HoldoutPrediction]

    @property
    def mae(self) -> float:
        return float(np.mean([p.abs_err for p in self.predictions]))

    @property
    def max_abs_err(self) -> float:
        return max(p.abs_err for p in self.predictions)


def holdout_predict(
    held_out: tuple[str, ...] | list[str],
    rows: list[EmpiricalSD] | None = None,
) -> HoldoutResult:
    """Fit the injury model on every arm except `held_out` trials; predict those."""
    rows = all_arms() if rows is None else rows
    held = set(held_out)
    train = [r for r in rows if r.nct_id not in held]
    test = [r for r in rows if r.nct_id in held]
    if len(train) < 4:
        raise ValueError(f"holdout leaves only {len(train)} training arms")
    if not test:
        raise ValueError(f"no arms found for held-out trials {held_out}")
    model = fit_injury_variance(train)
    preds = [
        HoldoutPrediction(
            nct_id=r.nct_id, arm_label=r.arm_label, injury=abs(r.delta),
            sd_observed=r.sd, sd_predicted=model.sd_at(r.delta),
            n_arm=r.n_arm or 0,
        )
        for r in test
    ]
    return HoldoutResult(held_out=tuple(held_out), model=model, predictions=preds)


def leave_one_trial_out(rows: list[EmpiricalSD] | None = None) -> list[HoldoutResult]:
    """Leave-one-trial-out across every trial that contributes arms."""
    rows = all_arms() if rows is None else rows
    trials = sorted({r.nct_id for r in rows})
    return [holdout_predict((t,), rows) for t in trials]


# ---------------------------------------------------------------------------
# Biopsy averaging, using Takeda's measured variance shares
# ---------------------------------------------------------------------------

def sd_after_biopsy_averaging(
    sd: float,
    n_biopsies: int,
    *,
    n_baseline: int = 1,
    biopsy_share: float | None = None,
) -> float:
    """SD after averaging `n_biopsies` independent fragments instead of `n_baseline`.

    Uses Takeda's measured biopsy-level share (23% of total variance in TAK-062) so
    this is no longer a free knob. The non-biopsy share (patient + reader + residual)
    does not shrink with more fragments.
    """
    from ctsim.published import MEASURED_VARIANCE_SHARES

    if n_biopsies < 1 or n_baseline < 1:
        raise ValueError("biopsy counts must be >= 1")
    share = MEASURED_VARIANCE_SHARES["biopsy"] if biopsy_share is None else biopsy_share
    # σ²' / σ² = (1 - share) + share * (n_baseline / n_biopsies)
    factor = (1.0 - share) + share * (n_baseline / n_biopsies)
    return sd * math.sqrt(factor)


def biopsy_n_savings(
    injury: float,
    protection: float = 0.50,
    *,
    n_from: int = 1,
    n_to: int = 4,
    model: InjuryVarianceModel | None = None,
) -> tuple[int, int, float, float]:
    """(N at n_from biopsies, N at n_to, SD_from, SD_to) for a given challenge."""
    model = fit_injury_variance() if model is None else model
    sd0 = model.sd_at(injury)
    sd1 = sd_after_biopsy_averaging(sd0, n_to, n_baseline=n_from)
    k = (stats.norm.ppf(1 - 0.05 / 2) + stats.norm.ppf(0.80)) ** 2
    effect = protection * abs(injury)
    n_base = math.ceil(2 * k * sd0**2 / effect**2)
    n_avg = math.ceil(2 * k * sd1**2 / effect**2)
    return n_base, n_avg, sd0, sd1


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
        "",
        "-" * 78,
        "Holdout: fit without KAN-101 and TAK-101, predict their SDs",
        "-" * 78,
    ]
    ho = holdout_predict(HEADLINE_HOLDOUT, rows)
    lines.append(f"  training fit: {ho.model.describe().split(chr(10))[0]}")
    lines.append(f"  {'trial':<14}{'arm':<28}{'|inj|':>6}{'obs':>7}{'pred':>7}{'err':>7}")
    for p in sorted(ho.predictions, key=lambda p: p.injury):
        lines.append(
            f"  {p.nct_id:<14}{p.arm_label[:28]:<28}"
            f"{p.injury:>6.2f}{p.sd_observed:>7.3f}{p.sd_predicted:>7.3f}"
            f"{p.abs_err:>7.3f}"
        )
    lines.append(f"  MAE = {ho.mae:.3f}   max |err| = {ho.max_abs_err:.3f}")
    lines += [
        "",
        "-" * 78,
        "Biopsy lever (Takeda measured 23% of variance at the biopsy level)",
        "-" * 78,
        "  relative to a single-fragment assay, at injury 0.61 / 50% protection:",
    ]
    for k in (1, 2, 4, 8):
        n1, nk, sd1, sdk = biopsy_n_savings(0.61, 0.50, n_from=1, n_to=k, model=m)
        lines.append(
            f"  {k} biopsies: SD {sdk:.3f}  N/arm {nk}"
            + (f"  (saves {n1 - nk} vs single biopsy)" if k > 1 else "")
        )
    lines.append(
        "  Caveat: trials already take multiple fragments, so the gain vs current"
    )
    lines.append(
        "  practice is smaller than the gain vs a theoretical single-biopsy assay."
    )
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
