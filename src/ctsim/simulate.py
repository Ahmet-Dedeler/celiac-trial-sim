"""Monte Carlo simulator and analytic power tools for celiac trial designs.

Core question this answers:

    Given the measurement noise that celiac histology endpoints actually carry,
    what effect size could a given trial design detect — and did the trials that
    "failed" ever have a realistic chance of succeeding?

The analytic path (`mde`, `required_n`) and the Monte Carlo path (`simulate_trial`)
are kept separate on purpose: the MC engine agrees with the closed form in the simple
case, which is the regression test, and then extends it to dropout, multi-biopsy
averaging and non-normal noise where no closed form exists.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import stats

from ctsim.model import (
    LIT,
    ancova_variance_ratio,
    baseline_correlation,
    decompose_vhcd,
    endpoint_comparison,
    load_empirical,
    pooled_sd,
    reader_share_range,
)


def _z(p: float) -> float:
    return float(stats.norm.ppf(p))


# ---------------------------------------------------------------------------
# Analytic
# ---------------------------------------------------------------------------

def mde(n_per_arm: int, sd: float, power: float = 0.80, alpha: float = 0.05) -> float:
    """Minimum detectable effect for a two-arm parallel trial, two-sided test."""
    return (_z(1 - alpha / 2) + _z(power)) * sd * math.sqrt(2.0 / n_per_arm)


def required_n(effect: float, sd: float, power: float = 0.80, alpha: float = 0.05) -> int:
    """Patients per arm needed to detect `effect` with the given power."""
    k = (_z(1 - alpha / 2) + _z(power)) ** 2
    return math.ceil(2 * k * (sd**2) / (effect**2))


def analytic_power(n_per_arm: int, effect: float, sd: float, alpha: float = 0.05) -> float:
    se = sd * math.sqrt(2.0 / n_per_arm)
    crit = _z(1 - alpha / 2)
    lam = effect / se
    return float(stats.norm.cdf(lam - crit) + stats.norm.cdf(-lam - crit))


# ---------------------------------------------------------------------------
# Monte Carlo
# ---------------------------------------------------------------------------

@dataclass
class SimResult:
    power: float
    mean_estimate: float
    se_estimate: float
    n_completers_per_arm: float
    n_sims: int

    def __str__(self) -> str:
        return (
            f"power={self.power:.3f}  est={self.mean_estimate:+.3f}"
            f"  se={self.se_estimate:.3f}  completers/arm={self.n_completers_per_arm:.1f}"
        )


def simulate_trial(
    n_per_arm: int,
    true_effect: float,
    total_sd: float,
    *,
    n_sims: int = 20_000,
    alpha: float = 0.05,
    dropout: float = 0.0,
    n_biopsies: int = 1,
    sampling_share: float = 0.0,
    seed: int = 0,
) -> SimResult:
    """Simulate a two-arm parallel trial on a continuous histology endpoint.

    Parameters
    ----------
    total_sd
        Between-patient SD of change-from-baseline, as observed in real trials.
    dropout
        Fraction of randomised patients lost before the endpoint biopsy.
    n_biopsies, sampling_share
        `sampling_share` is the fraction of *total variance* attributable to biopsy
        site/orientation sampling — i.e. the part that averaging more biopsy fragments
        per timepoint can reduce as 1/k. The rest (true biological change plus reader
        error) is not reducible by taking more fragments from the same patient.

        This is a sensitivity knob, not a measurement: the current public data cannot
        separate sampling from biological variance. Set it from an assumption and say so.
    """
    rng = np.random.default_rng(seed)

    var_total = total_sd**2
    var_sampling = sampling_share * var_total
    var_fixed = var_total - var_sampling
    # Averaging k independent fragments reduces only the sampling term.
    var_effective = var_fixed + var_sampling / max(n_biopsies, 1)
    sd_eff = math.sqrt(var_effective)

    n_complete = max(int(round(n_per_arm * (1 - dropout))), 2)

    ctrl = rng.normal(0.0, sd_eff, size=(n_sims, n_complete))
    trt = rng.normal(true_effect, sd_eff, size=(n_sims, n_complete))

    diff = trt.mean(axis=1) - ctrl.mean(axis=1)
    pooled_var = (ctrl.var(axis=1, ddof=1) + trt.var(axis=1, ddof=1)) / 2
    se = np.sqrt(2 * pooled_var / n_complete)
    t = diff / se
    crit = stats.t.ppf(1 - alpha / 2, df=2 * (n_complete - 1))
    reject = np.abs(t) > crit

    return SimResult(
        power=float(reject.mean()),
        mean_estimate=float(diff.mean()),
        se_estimate=float(se.mean()),
        n_completers_per_arm=float(n_complete),
        n_sims=n_sims,
    )


# ---------------------------------------------------------------------------
# Retrospective replay of real trials
# ---------------------------------------------------------------------------

@dataclass
class Replay:
    nct_id: str
    label: str
    n_per_arm: int
    sd_used: float
    mde_80: float
    ratio_to_meaningful: float
    power_for_meaningful: float
    power_for_half: float


def replay_trials(sd: float | None = None) -> list[Replay]:
    """For each trial with a usable ΔVH:CD arm, ask what it could have detected.

    `sd` defaults to the pooled empirical SD across all trials, so every design is
    judged against the same noise estimate rather than its own (which would let a
    small noisy trial flatter itself).
    """
    rows = load_empirical()
    sd = sd if sd is not None else pooled_sd(rows)
    meaningful = LIT["vhcd_clinically_significant"].value

    # One entry per trial, using the smallest arm N as the binding constraint.
    by_trial: dict[str, list] = {}
    for r in rows:
        if r.n_arm:
            by_trial.setdefault(r.nct_id, []).append(r)

    out: list[Replay] = []
    for nct, arms in sorted(by_trial.items()):
        n = min(a.n_arm for a in arms)
        m = mde(n, sd)
        out.append(
            Replay(
                nct_id=nct,
                label=" vs ".join(sorted({a.arm_label.split(":")[0][:22] for a in arms}))[:52],
                n_per_arm=n,
                sd_used=sd,
                mde_80=m,
                ratio_to_meaningful=m / meaningful,
                power_for_meaningful=analytic_power(n, meaningful, sd),
                power_for_half=analytic_power(n, meaningful / 2, sd),
            )
        )
    return out


def observed_sd_for(nct_id: str) -> tuple[float, int, str]:
    """The SD a trial actually got, on the scale its own analysis used."""
    from ctsim import published

    if nct_id.startswith("EudraCT"):
        sd, df = published.zed1227_residual_sd()
        return sd, df, "ancova_residual"
    arms = [r for r in load_empirical() if r.nct_id == nct_id]
    if not arms:
        return float("nan"), 0, "unknown"
    df = sum((a.n_arm - 1) for a in arms if a.n_arm and a.n_arm > 1)
    return pooled_sd(arms), df, arms[0].sd_scale


def planned_vs_actual() -> str:
    """Did the trials assume the right amount of noise?

    This is the question the rest of the repo backs into. A trial's sample size is a
    function of its assumed SD, and that assumption is written down in the protocol
    before anyone sees data. Comparing it to what the trial then measured turns
    "these trials were underpowered" into something falsifiable and attributable.
    """
    from ctsim import assay, published

    lines = [
        "Did the trials assume the right amount of noise?",
        "-" * 78,
        "  Every assumption below is quoted verbatim from a public protocol or SAP;",
        "  see ctsim.published for the exact sentence and URL.",
        "",
        (f"  {'trial':<20}{'target':>8}{'assumed SD':>12}{'observed':>10}"
         f"{'claimed':>9}{'actual power':>14}"),
    ]
    for a in published.ASSUMPTIONS:
        obs, _df, _scale = observed_sd_for(a.nct_id)
        actual = analytic_power(a.n_per_arm, a.target_effect, obs)
        lines.append(
            f"  {a.trial:<20}{a.target_effect:>8.2f}{a.assumed_sd:>12.2f}"
            f"{obs:>10.3f}{a.claimed_power:>9.0%}{actual:>14.0%}"
        )
    lines += [
        "",
        "  The one trial that overestimated its own noise is the one that worked.",
        "",
        "-" * 78,
        "Could each trial have detected the field's best drug?",
        "-" * 78,
    ]

    # ZED1227's placebo arm lost 0.61 of VH:CD; its 100 mg arm gave back 0.48 of that.
    delta, _lo, _hi = published.zed1227_control_injury()
    best_diff = 0.48  # ZED1227 100 mg vs placebo, NEJM 2021 Table 2
    benchmark = best_diff / abs(delta)
    lines.append(f"  Benchmark: ZED1227 100 mg prevented {best_diff:.2f} of a "
                 f"{abs(delta):.2f} injury = {benchmark:.0%} protection —")
    lines.append("  the largest protective effect demonstrated on this endpoint.")
    lines.append("")

    zed_sd, _df = published.zed1227_residual_sd()
    lines.append(f"  {'trial':<24}{'n/arm':>6}{'SD':>8}{'injury':>9}"
                 f"{'needs':>8}{'would it see 79%?':>20}")
    lines.append(f"  {'ZED1227 (CEC-3)':<24}{34:>6}{zed_sd:>8.3f}{abs(delta):>9.2f}"
                 f"{mde(34, zed_sd) / abs(delta):>7.0%}"
                 f"{'yes — it did':>20}")
    for c in assay.could_detect_benchmark(benchmark):
        own_sd, _d, _s = observed_sd_for(c.nct_id)
        needs = mde(c.n_per_arm, own_sd) / c.injury
        lines.append(
            f"  {c.nct_id:<24}{c.n_per_arm:>6}{own_sd:>8.3f}{c.injury:>9.2f}"
            f"{needs:>7.0%}{('yes' if benchmark >= needs else 'NO'):>20}"
        )
    lines += [
        "",
        "  KAN-101's placebo arm sustained exactly the same injury as ZED1227's (-0.61).",
        "  Judged on its own measured noise, it would have missed a drug as good as the",
        "  best one in the field. Its null result is not evidence about KAN-101.",
    ]
    return "\n".join(lines)


def report() -> str:
    # Imported here rather than at module scope: `assay` imports `mde` from this
    # module, so a top-level import would be circular.
    from ctsim import assay, uncertainty

    rows = load_empirical()
    est = uncertainty.pooled_estimate(rows)
    sd = est.sd
    meaningful = LIT["vhcd_clinically_significant"].value
    dec = decompose_vhcd(sd)

    lines: list[str] = []
    lines.append("=" * 78)
    lines.append("CELIAC TRIAL SIMULATOR — ΔVH:CD endpoint")
    lines.append("=" * 78)
    lines.append(f"pooled between-patient SD (empirical, {len(rows)} arms): {sd:.3f}"
                 f"  [95% CI {est.ci_low:.3f}-{est.ci_high:.3f}]")
    lines.append(f"clinically meaningful change:                        {meaningful:.2f}")
    lines.append(f"noise-to-signal ratio:                               {sd / meaningful:.2f}x")
    lines.append("")
    lines.append(dec.describe())
    lines.append("")
    lines.append("  reader share under every published estimate of reader error:")
    for label, reader_sd, share in reader_share_range(sd):
        lines.append(f"    {label:<62} sd={reader_sd:.3f}  {share:>5.1%}")
    lines.append("")
    lines.append(uncertainty.report(rows))
    lines.append("")

    lines.append("-" * 78)
    lines.append("Sample size required (80% power, alpha=0.05, two-sided)")
    lines.append("-" * 78)
    for eff, name in [
        (meaningful, "full clinically meaningful change (0.40)"),
        (meaningful * 0.75, "3/4 of meaningful (0.30)"),
        (meaningful / 2, "half of meaningful (0.20)"),
        (meaningful / 4, "quarter of meaningful (0.10)"),
    ]:
        n = required_n(eff, sd)
        lines.append(f"  detect {name:<42} {n:>5} per arm  ({2 * n:>5} total)")
    lines.append("")

    lines.append("-" * 78)
    lines.append("Retrospective replay: what could each trial actually detect?")
    lines.append("-" * 78)
    lines.append(f"  {'trial':<14}{'n/arm':>6}{'MDE(80%)':>10}{'vs 0.40':>9}"
                 f"{'power@0.40':>12}{'power@0.20':>12}")
    for r in replay_trials(sd):
        flag = "  <-- underpowered for a meaningful effect" if r.ratio_to_meaningful > 1 else ""
        lines.append(
            f"  {r.nct_id:<14}{r.n_per_arm:>6}{r.mde_80:>10.3f}"
            f"{r.ratio_to_meaningful:>8.2f}x{r.power_for_meaningful:>12.1%}"
            f"{r.power_for_half:>12.1%}{flag}"
        )
    lines.append("")
    lines.append("  the same trials judged against the ends of the SD confidence interval:")
    lines.append(f"  {'trial':<14}{'power@0.40 (SD low)':>22}{'(point)':>10}{'(SD high)':>12}")
    for r in replay_trials(sd):
        lines.append(
            f"  {r.nct_id:<14}"
            f"{analytic_power(r.n_per_arm, meaningful, est.ci_low):>22.1%}"
            f"{r.power_for_meaningful:>10.1%}"
            f"{analytic_power(r.n_per_arm, meaningful, est.ci_high):>12.1%}"
        )
    lines.append("")
    lines.append(assay.report(rows))
    lines.append("")
    lines.append(planned_vs_actual())
    lines.append("")

    lines.append("-" * 78)
    lines.append("Monte Carlo cross-check (should match the analytic column)")
    lines.append("-" * 78)
    for n in (25, 60, 100):
        s = simulate_trial(n, meaningful, sd, n_sims=40_000, seed=7)
        lines.append(
            f"  n={n:<4} true effect=0.40   MC power={s.power:.3f}   "
            f"analytic={analytic_power(n, meaningful, sd):.3f}"
        )
    lines.append("")

    lines.append("-" * 78)
    lines.append("Sensitivity: if X% of variance is biopsy sampling, averaging k fragments")
    lines.append("-" * 78)
    for share in (0.2, 0.4):
        row = [f"  sampling_share={share:.0%}: "]
        for k in (1, 2, 4, 8):
            s = simulate_trial(50, meaningful, sd, n_biopsies=k,
                               sampling_share=share, n_sims=20_000, seed=11)
            row.append(f"k={k}:{s.power:.2f}")
        lines.append("  ".join(row))
    lines.append("")
    lines.append("  (n=50/arm, true effect 0.40. sampling_share is an assumption, not a")
    lines.append("   measurement — public data cannot separate sampling from biology.)")
    lines.append("")

    lines.append("-" * 78)
    lines.append("Design levers, ranked by what they actually buy")
    lines.append("-" * 78)
    base_n = required_n(meaningful, sd)
    lines.append(f"  baseline: {base_n} per arm to detect 0.40 with a change-score analysis")
    lines.append("")

    corrs = baseline_correlation()
    if corrs:
        rhos = sorted(c.rho for c in corrs)
        lines.append("  1. Analyse by ANCOVA on baseline instead of a change score.")
        lines.append(f"     baseline-follow-up correlation recovered from posted SDs: "
                     f"rho = {rhos[0]:.2f}-{rhos[-1]:.2f}")
        for r in (rhos[0], rhos[-1]):
            ratio = ancova_variance_ratio(r)
            n = required_n(meaningful, sd * math.sqrt(ratio))
            lines.append(f"       rho={r:.2f}: {base_n} -> {n} per arm "
                         f"({1 - n / base_n:.0%} fewer patients)")
        lines.append("     Free — it is a line in the analysis plan — but it does not rescue"
                     " anything.")
        lines.append("")

    # Compare only protocols that actually produced measurable injury. Including the
    # arms where the challenge did nothing would give a true but useless ratio.
    opts = [o for o in assay.challenge_dose_options(25, sd, rows) if o.injury >= meaningful]
    if len(opts) >= 2:
        strongest = max(opts, key=lambda o: o.injury)
        weakest = min(opts, key=lambda o: o.injury)
        lines.append("  2. Use a gluten challenge that actually injures the mucosa.")
        for o in (weakest, strongest):
            lines.append(
                f"     {o.source} {o.label[:30]:<32} injury {o.injury:.3f}  -> "
                f"{o.n_for_half_protection} per arm to catch 50% protection"
            )
        lines.append(f"     Same drug, same endpoint, same analysis: "
                     f"{weakest.n_for_half_protection} patients per arm or "
                     f"{strongest.n_for_half_protection}.")
        lines.append("     Required N scales with 1/injury^2, so this dominates every other"
                     " lever.")
        lines.append("")
        skipped = [o for o in assay.challenge_dose_options(25, sd, rows)
                   if o.injury < meaningful]
        for o in skipped:
            lines.append(f"     (excluded: {o.source} {o.label[:26]} moved the mucosa only "
                         f"{o.injury:.3f} — no usable signal at any N)")
        if skipped:
            lines.append("")

    eps = endpoint_comparison()
    if eps:
        lines.append("  3. Switch endpoint to IEL density. (It does not help.)")
        for e in eps:
            lines.append(
                f"     {e.endpoint:<12} meaningful={e.meaningful:>7.3f}  "
                f"SD={e.pooled_sd:>7.3f}  standardized={e.standardized_effect:.3f}  "
                f"-> {required_n(e.meaningful, e.pooled_sd)} per arm"
            )
        lines.append("     Same population, same patients. The noise is in the biology and"
                     " the biopsy,")
        lines.append("     not in the choice of what to measure on the slide.")

    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
