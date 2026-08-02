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
    leave_one_trial_out,
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

def test_real_pool_is_homogeneous():
    """Regression guard on the claim that pooling these arms is legitimate.

    If new data breaks this, the right response is to stop quoting one pooled SD, not
    to relax the test.
    """
    h = bartlett_homogeneity()
    assert h.homogeneous, (
        f"per-arm variances are no longer mutually consistent (p={h.p_value:.4f}); "
        "a single pooled SD is no longer the right summary"
    )


def test_no_single_trial_drives_the_headline():
    """Dropping any one trial must not move the SD enough to change the conclusion."""
    full = pooled_estimate().sd
    for s in leave_one_trial_out():
        assert abs(s.sd - full) / full < 0.15, f"{s.label} moves the pooled SD too much"
        # The qualitative claim is SD > meaningful change (0.40). It must survive.
        assert s.sd > 0.40


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
