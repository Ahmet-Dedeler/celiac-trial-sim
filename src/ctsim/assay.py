"""Was there anything there to detect?

Power analysis asks whether a trial could resolve an effect of a given size. It takes
the effect size as given. For celiac gluten-challenge trials that is the wrong place to
stop, because the effect size is not a free parameter — it is capped by the design.

The logic of a challenge trial is: feed everyone gluten, let the control arm's mucosa
take the damage, and see whether the drug arm's does not. So

    * the control arm's observed deterioration is the *entire* signal available;
    * a drug preventing fraction f of that injury produces a between-arm difference of
      exactly f * |control change|;
    * f cannot exceed 1, because a drug cannot prevent more injury than gluten caused.

That gives a hard ceiling on the detectable effect, and therefore a number that matters
more than power:

    minimum detectable protection = MDE(80%) / |control arm change|

If that exceeds 1, no drug — not a perfect one — reaches 80% power in that trial. The
trial was not underpowered for its drug. It was unable to succeed.

This is the standard regulatory notion of *assay sensitivity*, applied to an endpoint
where the control-arm trajectory is public and nobody seems to have checked it.

Caveat, stated up front: the ceiling is |control change| only for a pure prevention
design. A drug could in principle push VH:CD *above* baseline by removing ongoing
gluten exposure and letting the mucosa heal, and the headroom for that is
(normal VH:CD - baseline VH:CD), not the challenge injury. None of these trials post a
baseline VH:CD, so that channel cannot be quantified from public data and is flagged
rather than assumed away.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from scipy import stats

from ctsim.model import LIT, EmpiricalSD, load_empirical, pooled_sd
from ctsim.simulate import mde


@dataclass
class AssaySensitivity:
    nct_id: str
    control_label: str
    n_per_arm: int
    control_delta: float
    control_delta_lo: float  # 95% CI on the control-arm change
    control_delta_hi: float
    mde_80: float
    sd_used: float
    min_protection: float          # using the point estimate of the injury
    min_protection_best_case: float  # using the largest injury the CI allows
    trial_has_control: bool = True

    @property
    def feasible(self) -> bool:
        """Could any drug, however good, have hit 80% power on this endpoint?"""
        return self.min_protection_best_case <= 1.0

    def describe(self) -> str:
        if self.min_protection > 1.0:
            verdict = "IMPOSSIBLE — a perfect drug still misses 80% power"
        elif self.min_protection > 0.75:
            verdict = f"needs {self.min_protection:.0%} protection — near-perfect drug only"
        else:
            verdict = f"needs {self.min_protection:.0%} protection"
        return f"{self.nct_id}: {verdict}"


def _injury_ci(delta: float, sd: float, n: int, conf: float = 0.95) -> tuple[float, float]:
    """CI on a single arm's mean change. The injury is an estimate too."""
    se = sd / math.sqrt(n)
    t = stats.t.ppf(1 - (1 - conf) / 2, df=max(n - 1, 1))
    return delta - t * se, delta + t * se


def assay_sensitivity(rows: list[EmpiricalSD] | None = None,
                      sd: float | None = None) -> list[AssaySensitivity]:
    """For each placebo-controlled trial, how good would the drug have had to be?

    `sd` overrides the noise estimate. Left unset, each trial is judged on **its own**
    measured SD rather than a common pooled one. That used to be the other way round,
    on the reasoning that a shared estimate stops a small noisy trial flattering itself.
    It has to be per-trial now: the noise is not a constant across designs — it scales
    with how hard the challenge hit (see `ctsim.variance`) — so a common SD would
    penalise gentle-challenge trials and flatter harsh-challenge ones.

    The trials' own numbers are also the more conservative choice here, and the headline
    comparison does not depend on the switch.
    """
    rows = load_empirical() if rows is None else rows

    by_trial: dict[str, list[EmpiricalSD]] = {}
    for r in rows:
        if r.n_arm:
            by_trial.setdefault(r.nct_id, []).append(r)

    out: list[AssaySensitivity] = []
    for nct, arms in sorted(by_trial.items()):
        controls = [a for a in arms if a.is_placebo]
        if not controls:
            # No placebo arm: a dose-comparison or methodology study. The protection
            # framing does not apply — there is no "untreated injury" to prevent.
            continue
        if arms[0].design != "prevention":
            # Restoration designs enrol patients who already have villous atrophy and
            # expect the drug to let it heal. A flat control arm there is not a missing
            # signal — the ceiling is the headroom to a normal mucosa, which needs a
            # baseline VH:CD nobody posts. Calling such a trial "impossible" on the
            # strength of its control arm not moving would be exactly backwards.
            continue
        ctrl = controls[0]
        n = min(a.n_arm for a in arms)
        sd_used = sd if sd is not None else pooled_sd(arms)
        m = mde(n, sd_used)
        lo, hi = _injury_ci(ctrl.delta, ctrl.sd, ctrl.n_arm or n)
        # Injury magnitude is the size of the control arm's deterioration. The best case
        # for the trial is the largest injury its CI still allows.
        injury = abs(ctrl.delta)
        injury_best = max(abs(lo), abs(hi))
        out.append(
            AssaySensitivity(
                nct_id=nct,
                control_label=ctrl.arm_label,
                n_per_arm=n,
                control_delta=ctrl.delta,
                control_delta_lo=lo,
                control_delta_hi=hi,
                mde_80=m,
                sd_used=sd_used,
                min_protection=m / injury if injury > 0 else float("inf"),
                min_protection_best_case=m / injury_best if injury_best > 0 else float("inf"),
            )
        )
    return out


# ---------------------------------------------------------------------------
# The design lever: how hard you challenge sets the ceiling
# ---------------------------------------------------------------------------

@dataclass
class ChallengeOption:
    source: str
    label: str
    injury: float
    min_protection: float
    n_for_half_protection: int


def challenge_dose_options(n_per_arm: int, sd: float | None = None,
                           rows: list[EmpiricalSD] | None = None) -> list[ChallengeOption]:
    """What the ceiling would be under each gluten-challenge protocol on record.

    The injury magnitude is not a fact of nature — it is chosen, when the sponsor picks
    a challenge dose and duration. Every *unmedicated* arm in the dataset is a
    measurement of what some protocol actually does to the mucosa, so they can be read
    across as design options for the next trial.

    Drug arms are excluded on purpose: their change reflects challenge minus whatever
    the drug did, which is the thing being measured, not a property of the protocol.

    Each protocol is costed at the SD that goes with *its own* injury, not a shared one.
    That matters: a harsher challenge raises the noise as well as the signal, so pricing
    every protocol at one SD makes hard challenges look better than they are.
    """
    from ctsim.variance import fit_injury_variance

    rows = load_empirical() if rows is None else rows
    model = fit_injury_variance()
    k = (stats.norm.ppf(0.975) + stats.norm.ppf(0.80)) ** 2

    out: list[ChallengeOption] = []
    for r in sorted(rows, key=lambda r: r.delta):
        if r.delta >= 0 or not r.is_unmedicated or r.design != "prevention":
            # Needs to be an unmedicated arm, in a challenge design, that actually moved.
            # A restoration trial's control arm is not measuring a challenge protocol.
            continue
        injury = abs(r.delta)
        sd_here = sd if sd is not None else model.sd_at(injury)
        half = injury / 2
        out.append(
            ChallengeOption(
                source=r.nct_id,
                label=r.arm_label,
                injury=injury,
                min_protection=mde(n_per_arm, sd_here) / injury,
                n_for_half_protection=math.ceil(2 * k * sd_here**2 / half**2),
            )
        )
    return out


# ---------------------------------------------------------------------------
# The benchmark: could each trial have seen the one effect anybody has seen?
# ---------------------------------------------------------------------------

@dataclass
class BenchmarkCheck:
    nct_id: str
    n_per_arm: int
    sd_used: float
    injury: float
    required_protection: float
    benchmark_protection: float
    detected: bool


def could_detect_benchmark(benchmark_protection: float,
                           rows: list[EmpiricalSD] | None = None,
                           sd: float | None = None) -> list[BenchmarkCheck]:
    """Which prevention trials could have detected a drug as good as the best one?

    "Underpowered" invites the reply that the effect might simply have been zero. This
    sidesteps that: take the largest protective effect anyone has actually demonstrated
    on this endpoint, hand it to each trial, and ask whether that trial would have found
    it. A trial that would have missed the field's best drug cannot be read as evidence
    that its own drug did not work.
    """
    rows = load_empirical() if rows is None else rows
    out: list[BenchmarkCheck] = []
    for a in assay_sensitivity(rows):
        out.append(
            BenchmarkCheck(
                nct_id=a.nct_id,
                n_per_arm=a.n_per_arm,
                sd_used=sd if sd is not None else pooled_sd(rows),
                injury=abs(a.control_delta),
                required_protection=a.min_protection,
                benchmark_protection=benchmark_protection,
                detected=benchmark_protection >= a.min_protection,
            )
        )
    return out


def report(rows: list[EmpiricalSD] | None = None) -> str:
    rows = load_empirical() if rows is None else rows
    meaningful = LIT["vhcd_clinically_significant"].value

    lines = [
        "Assay sensitivity: was there a detectable signal available at all?",
        "-" * 78,
        "  (each trial judged on its own measured SD; MDE at 80% power, two-sided)",
        "",
        (f"  {'trial':<14}{'n/arm':>6}{'own SD':>8}{'control d':>11}{'95% CI':>18}"
         f"{'MDE':>8}{'min protection':>16}"),
    ]
    for a in assay_sensitivity(rows):
        prot = ("impossible" if a.min_protection > 1
                else f"{a.min_protection:.0%}")
        lines.append(
            f"  {a.nct_id:<14}{a.n_per_arm:>6}{a.sd_used:>8.3f}{a.control_delta:>+11.3f}"
            f"{f'[{a.control_delta_lo:+.2f},{a.control_delta_hi:+.2f}]':>18}"
            f"{a.mde_80:>8.3f}{prot:>16}"
        )
        best = (f"{a.min_protection_best_case:.0%}"
                if a.min_protection_best_case <= 1 else
                f"{a.min_protection_best_case:.0%} — still impossible")
        lines.append(
            f"      taking the most favourable end of the control-arm CI: {best}"
        )

    lines += [
        "",
        "  A gluten challenge trial can only detect what the challenge causes. The",
        "  ceiling on the between-arm difference is the control arm's own injury.",
        "",
        "-" * 78,
        "The lever: challenge protocol sets the ceiling before N does",
        "-" * 78,
        "  observed injury under each challenge protocol on record, and what it would",
        "  take to exploit it (n=25/arm, the size KAN-101's Phase 2a actually ran):",
        "",
        (f"  {'source':<14}{'arm':<26}{'injury':>9}{'min protection':>16}"
         f"{'N/arm for 50%':>15}"),
    ]
    for c in challenge_dose_options(25, None, rows):
        prot = "impossible" if c.min_protection > 1 else f"{c.min_protection:.0%}"
        lines.append(
            f"  {c.source:<14}{c.label[:24]:<26}{c.injury:>9.3f}{prot:>16}"
            f"{c.n_for_half_protection:>15}"
        )
    lines += [
        "",
        f"  For reference the 'clinically meaningful' change is {meaningful:.2f}, so a",
        "  challenge that moves the control arm by less than that cannot produce a",
        "  clinically meaningful between-arm difference no matter how good the drug is.",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
