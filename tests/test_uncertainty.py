"""Tests for the uncertainty layer.

The point of these is that the headline SD is a *statistic*, not a measurement, and
every claim downstream inherits its error bars. If the CI machinery silently breaks,
the repo goes back to asserting 0.740 as though it were exact.
"""

from __future__ import annotations

import math

import pytest

from ctsim.model import EmpiricalSD, load_empirical, pooled_sd
from ctsim.uncertainty import (
    bartlett_homogeneity,
    pooled_estimate,
    provenance_subsets,
    sd_confidence_interval,
)


def _arm(nct: str, n: int, sd: float, src: str = "posted_sd", placebo: bool = False):
    return EmpiricalSD(nct_id=nct, arm_label="placebo" if placebo else "drug",
                       n_arm=n, delta=-0.5, sd=sd, sd_source=src, time_frame="",
                       is_placebo=placebo)


# --- chi-square CI -----------------------------------------------------------

def test_ci_brackets_the_point_estimate():
    lo, hi = sd_confidence_interval(0.74, 177)
    assert lo < 0.74 < hi


def test_ci_narrows_as_df_grows():
    intervals = [sd_confidence_interval(0.74, df) for df in (10, 50, 200, 1000)]
    widths = [hi - lo for lo, hi in intervals]
    assert widths == sorted(widths, reverse=True)


def test_ci_matches_the_chi_square_definition():
    """Recompute the interval from first principles as an independent check."""
    from scipy import stats
    sd, df = 0.74, 177
    lo, hi = sd_confidence_interval(sd, df, conf=0.95)
    assert lo == pytest.approx(sd * math.sqrt(df / stats.chi2.ppf(0.975, df)))
    assert hi == pytest.approx(sd * math.sqrt(df / stats.chi2.ppf(0.025, df)))
    # And the interval must actually cover 95% of the sampling distribution.
    cover = stats.chi2.cdf(df * sd**2 / lo**2, df) - stats.chi2.cdf(df * sd**2 / hi**2, df)
    assert cover == pytest.approx(0.95, abs=1e-9)


# --- Bartlett ----------------------------------------------------------------

def test_bartlett_sees_identical_variances_as_homogeneous():
    rows = [_arm("A", 30, 0.7), _arm("B", 30, 0.7), _arm("C", 30, 0.7)]
    h = bartlett_homogeneity(rows)
    assert h.statistic == pytest.approx(0.0, abs=1e-9)
    assert h.p_value > 0.99
    assert h.homogeneous


def test_bartlett_flags_wildly_different_variances():
    rows = [_arm("A", 60, 0.2), _arm("B", 60, 1.5), _arm("C", 60, 0.3)]
    h = bartlett_homogeneity(rows)
    assert h.p_value < 0.001
    assert not h.homogeneous


def test_bartlett_agrees_with_scipy_on_generated_samples():
    """Our summary-statistic Bartlett must match scipy's sample-based one."""
    import numpy as np
    from scipy import stats

    rng = np.random.default_rng(0)
    samples = [rng.normal(0, s, size=n) for s, n in [(0.6, 40), (0.9, 35), (0.75, 50)]]
    rows = [_arm(f"T{i}", len(x), float(np.std(x, ddof=1))) for i, x in enumerate(samples)]

    ours = bartlett_homogeneity(rows)
    theirs = stats.bartlett(*samples)
    assert ours.statistic == pytest.approx(theirs.statistic, rel=1e-6)
    assert ours.p_value == pytest.approx(theirs.pvalue, rel=1e-6)


# --- robustness of the real pool ---------------------------------------------

def test_the_pool_is_heterogeneous_which_is_why_one_sd_is_not_quoted():
    """This test asserts the *failure* that motivates the injury-scaled model.

    An earlier version of this file asserted the opposite and passed, because the
    dataset then contained a bad SE->SD conversion that flattened the differences.
    With that fixed, the per-arm variances are demonstrably not draws from one common
    variance — so `ctsim.variance` exists, and no single pooled SD is quoted anywhere.
    """
    h = bartlett_homogeneity()
    assert not h.homogeneous, (
        f"variances now look homogeneous (p={h.p_value:.4f}); if that is real, the "
        "injury-scaled model may be unnecessary — check before deleting it"
    )


def test_the_injury_model_explains_the_heterogeneity():
    """The whole claim of ctsim.variance: injury accounts for the spread in SDs."""
    from ctsim.variance import all_arms, fit_injury_variance

    rows = all_arms()
    m = fit_injury_variance(rows)
    resid = [r.sd - m.sd_at(r.delta) for r in rows]
    raw = [r.sd for r in rows]
    spread_before = max(raw) - min(raw)
    spread_after = max(resid) - min(resid)
    assert spread_after < spread_before / 2, (
        f"injury does not explain the SD spread: {spread_before:.3f} -> {spread_after:.3f}"
    )
    assert m.r > 0.75, f"correlation of injury with SD is only {m.r:.2f}"
    assert m.slope > 2 * m.slope_se, "slope is not distinguishable from zero"


def test_the_noise_floor_is_well_determined():
    """The floor is the number a trial planner actually needs; pin it."""
    from ctsim.variance import fit_injury_variance

    m = fit_injury_variance()
    assert 0.30 < m.floor < 0.50
    assert m.floor_se < 0.06


def test_control_arms_agree_with_the_full_pool():
    """Noise must look the same in placebo arms as everywhere else.

    If it did not, the 'SD' would be measuring drug-response heterogeneity rather than
    endpoint noise, and using it as an assay property would be wrong.
    """
    subs = {s.label: s.sd for s in provenance_subsets()}
    assert abs(subs["placebo/control arms only"] - subs["all arms"]) < 0.10


def test_pooled_estimate_reports_the_df_it_actually_used():
    est = pooled_estimate()
    rows = load_empirical()
    assert est.df == sum((r.n_arm - 1) for r in rows if r.n_arm and r.n_arm > 1)
    assert est.sd == pytest.approx(pooled_sd(rows))
