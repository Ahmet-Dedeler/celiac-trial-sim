"""How much should you trust the pooled SD?

The headline claim of this repo rests on a single number — the between-patient SD of
ΔVH:CD, pooled across a handful of arms. A point estimate with a prose caveat ("small
pool, but it looks stable") is not good enough to hang a power analysis on, because
required N scales with SD**2: a 20% error in the SD is a 44% error in the trial you
would tell someone to run.

So this module does three things the caveat cannot:

1. **Sampling interval.** A pooled variance from df degrees of freedom has a chi-square
   sampling distribution, which gives an exact CI on the SD *conditional on the arms
   sharing one variance*.

2. **Tests that conditional.** Bartlett's test asks whether the per-arm variances are
   mutually consistent. If they are not, pooling is the wrong operation and the CI in
   (1) is a fiction.

3. **Checks it isn't one trial's number.** Leave-one-trial-out, plus subsets by
   provenance (posted SD vs SE-derived) and by arm type (placebo only). If the answer
   moves a lot when one trial leaves, the pool is a single trial wearing a hat.

All three are reported. None of them is allowed to be silent.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from scipy import stats

from ctsim.model import EmpiricalSD, load_empirical, pooled_sd


def _df(rows: list[EmpiricalSD]) -> int:
    """Degrees of freedom backing a pooled variance, matching `pooled_sd`'s weights."""
    return sum((r.n_arm - 1) if r.n_arm and r.n_arm > 1 else 1 for r in rows)


# ---------------------------------------------------------------------------
# 1. Sampling interval on the pooled SD
# ---------------------------------------------------------------------------

@dataclass
class PooledEstimate:
    sd: float
    df: int
    n_arms: int
    n_trials: int
    ci_low: float
    ci_high: float
    conf: float

    def describe(self) -> str:
        return (
            f"pooled SD = {self.sd:.3f}  "
            f"({self.conf:.0%} CI {self.ci_low:.3f}-{self.ci_high:.3f}; "
            f"{self.n_arms} arms, {self.n_trials} trials, df={self.df})"
        )


def sd_confidence_interval(sd: float, df: int, conf: float = 0.95) -> tuple[float, float]:
    """Exact CI for a normal SD from a chi-square-distributed variance estimate.

    (df * s^2) / sigma^2 ~ chi2(df), so sigma lies in
    s * sqrt(df / chi2_{1-a/2}) .. s * sqrt(df / chi2_{a/2}).

    This is *sampling* uncertainty only. It assumes the pooled arms share one variance;
    `bartlett_homogeneity` is what checks that assumption.
    """
    a = (1 - conf) / 2
    lo = sd * math.sqrt(df / stats.chi2.ppf(1 - a, df))
    hi = sd * math.sqrt(df / stats.chi2.ppf(a, df))
    return lo, hi


def pooled_estimate(rows: list[EmpiricalSD] | None = None,
                    conf: float = 0.95) -> PooledEstimate:
    rows = load_empirical() if rows is None else rows
    sd = pooled_sd(rows)
    df = _df(rows)
    lo, hi = sd_confidence_interval(sd, df, conf)
    return PooledEstimate(
        sd=sd, df=df, n_arms=len(rows), n_trials=len({r.nct_id for r in rows}),
        ci_low=lo, ci_high=hi, conf=conf,
    )


# ---------------------------------------------------------------------------
# 2. Are the per-arm variances mutually consistent?
# ---------------------------------------------------------------------------

@dataclass
class Homogeneity:
    statistic: float
    df: int
    p_value: float
    min_sd: float
    max_sd: float

    @property
    def homogeneous(self) -> bool:
        return self.p_value >= 0.05

    def describe(self) -> str:
        verdict = (
            "consistent with a single common variance -> pooling is justified"
            if self.homogeneous else
            "HETEROGENEOUS -> a single pooled SD is the wrong summary"
        )
        return (
            f"Bartlett homogeneity of variances: T={self.statistic:.2f}, "
            f"df={self.df}, p={self.p_value:.3f}\n"
            f"  per-arm SD range {self.min_sd:.3f}-{self.max_sd:.3f}; {verdict}"
        )


def bartlett_homogeneity(rows: list[EmpiricalSD] | None = None) -> Homogeneity:
    """Bartlett's test that k arm variances come from one common variance.

    Bartlett is the right test here rather than Levene, because we only have posted
    summary statistics (s_i, n_i) — no patient-level values to take deviations from.
    It assumes normality, which for a change score across dozens of patients is the
    same assumption the trials' own t-tests already make.
    """
    rows = load_empirical() if rows is None else rows
    dfs = [(r.n_arm - 1) if r.n_arm and r.n_arm > 1 else 1 for r in rows]
    k = len(rows)
    N = sum(dfs)
    sp2 = sum(d * r.sd**2 for d, r in zip(dfs, rows)) / N
    numer = N * math.log(sp2) - sum(d * math.log(r.sd**2) for d, r in zip(dfs, rows))
    c = 1 + (sum(1 / d for d in dfs) - 1 / N) / (3 * (k - 1))
    T = numer / c
    return Homogeneity(
        statistic=T, df=k - 1,
        p_value=float(stats.chi2.sf(T, k - 1)),
        min_sd=min(r.sd for r in rows), max_sd=max(r.sd for r in rows),
    )


# ---------------------------------------------------------------------------
# 3. Is the pool really just one trial?
# ---------------------------------------------------------------------------

@dataclass
class Subset:
    label: str
    sd: float
    n_arms: int
    df: int


def leave_one_trial_out(rows: list[EmpiricalSD] | None = None) -> list[Subset]:
    """Recompute the pooled SD with each trial removed in turn."""
    rows = load_empirical() if rows is None else rows
    out: list[Subset] = []
    for nct in sorted({r.nct_id for r in rows}):
        rest = [r for r in rows if r.nct_id != nct]
        if not rest:
            continue
        out.append(Subset(f"without {nct}", pooled_sd(rest), len(rest), _df(rest)))
    return out


def provenance_subsets(rows: list[EmpiricalSD] | None = None) -> list[Subset]:
    """Split the pool along every axis that could be quietly driving the answer.

    `placebo/control arms only` is the most important one: an active arm's variance can
    be inflated by heterogeneous drug response, so the control arms are the cleaner
    estimate of pure endpoint noise. If those agree with the full pool, the SD is
    measuring the assay and not the pharmacology.
    """
    rows = load_empirical() if rows is None else rows
    control = [r for r in rows if "placebo" in r.arm_label.lower()]
    groups = [
        ("all arms", rows),
        ("placebo/control arms only", control),
        ("posted SD only (no SE conversion)", [r for r in rows if r.sd_source == "posted_sd"]),
        ("SE/CI-derived only", [r for r in rows if r.sd_source != "posted_sd"]),
    ]
    return [Subset(label, pooled_sd(sel), len(sel), _df(sel))
            for label, sel in groups if sel]


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def report(rows: list[EmpiricalSD] | None = None) -> str:
    rows = load_empirical() if rows is None else rows
    est = pooled_estimate(rows)
    hom = bartlett_homogeneity(rows)

    lines = [
        "How much should you trust the pooled SD?",
        "-" * 78,
        est.describe(),
        "",
        hom.describe(),
        "",
        "  per-arm detail (own CI shows how little a small arm constrains anything):",
    ]
    for r in sorted(rows, key=lambda r: (r.nct_id, r.arm_label)):
        n = r.n_arm or 2
        lo, hi = sd_confidence_interval(r.sd, max(n - 1, 1))
        lines.append(
            f"    {r.nct_id}  {r.arm_label[:26]:<28} n={n:<3} SD={r.sd:.3f} "
            f"[{lo:.3f}-{hi:.3f}]  ({r.sd_source})"
        )

    lines += ["", "  leave-one-trial-out:"]
    for s in leave_one_trial_out(rows):
        lines.append(f"    {s.label:<34} SD={s.sd:.3f}  ({s.n_arms} arms, df={s.df})")

    lines += ["", "  subsets by provenance and arm type:"]
    for s in provenance_subsets(rows):
        lines.append(f"    {s.label:<34} SD={s.sd:.3f}  ({s.n_arms} arms, df={s.df})")

    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
