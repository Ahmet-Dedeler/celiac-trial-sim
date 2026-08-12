"""Tests for the VCIEL break-even rule.

The claim being pinned is that Syage and Takeda are not in conflict, and that the thing
deciding which of them applies is a single ratio with a derivable threshold.
"""

from __future__ import annotations

import math

import pytest

from ctsim.endpoint import (
    BREAK_EVEN_RATIO,
    COMPARISONS,
    NOT_AVAILABLE,
    SYAGE_REPORTED,
    implied_correlation,
)


def test_break_even_is_the_square_root_of_two_minus_one():
    """Not a fitted constant. It falls out of averaging two uncorrelated measures."""
    assert BREAK_EVEN_RATIO == pytest.approx(math.sqrt(2) - 1)
    assert BREAK_EVEN_RATIO == pytest.approx(0.414, abs=0.001)


def test_composite_beats_vhcd_exactly_at_the_threshold():
    """The rule has to be self-consistent: ratio above the threshold iff N improves."""
    for c in COMPARISONS:
        if c.ratio > BREAK_EVEN_RATIO:
            assert c.d_composite > c.d_vhcd and c.n_ratio < 1.0
        else:
            assert c.d_composite <= c.d_vhcd and c.n_ratio >= 1.0


def test_imgx003_is_the_case_where_the_composite_pays():
    c = next(x for x in COMPARISONS if x.trial.startswith("IMGX003"))
    assert c.ratio > 1.0, "IEL carried more signal than VH:CD in this trial"
    assert c.composite_helps
    assert c.n_ratio < 0.5, "the composite should cut required N by more than half"


def test_tak101_sits_on_the_threshold_and_the_composite_would_not_have_saved_it():
    """An honest negative, and it matters for this repo's own headline argument.

    TAK-101 is the trial argued elsewhere here to have been underpowered rather than
    ineffective. The tempting follow-on is that a better endpoint would have rescued it.
    It would not have: its IEL effect is right at break-even, so VCIEL buys it nothing.
    It needed patients.
    """
    c = next(x for x in COMPARISONS if x.trial.startswith("TAK-101"))
    assert c.ratio == pytest.approx(BREAK_EVEN_RATIO, abs=0.05)
    assert 0.95 < c.n_ratio < 1.10, "no meaningful sample-size change either way"


def test_arithmetic_reproduces_syage_independently():
    """Our pooled-change-SD route must agree with his baseline-SD route on the ratio.

    Levels differ because the standardisers differ; the ratio of the two endpoints'
    effects is the quantity that should survive, and it does to within about 5%.
    """
    ours = next(x for x in COMPARISONS if x.trial == "IMGX003 CeliacShield")
    _trial, dv, di, _dc = next(r for r in SYAGE_REPORTED
                               if r[0] == "IMGX003 CeliacShield")
    assert ours.ratio == pytest.approx(di / dv, rel=0.10)


def test_implied_correlation_matches_what_syage_measured_directly():
    """Syage reports R^2 = 0.005-0.023 between the two measures, so |rho| <= ~0.15.

    Backing rho out of his reported composite effect sizes has to land in the same
    place, or the composite is not the equally weighted sum its definition claims.
    """
    for _trial, dv, di, dc in SYAGE_REPORTED:
        rho = implied_correlation(dv, di, dc)
        assert abs(rho) < 0.25, f"implied rho {rho:.2f} is too far from independence"


def test_the_number_that_would_settle_it_is_recorded_as_missing():
    """The reconciliation is a prediction until TAK-062's IEL table exists.

    Recorded so it stays a stated gap rather than quietly becoming a claim.
    """
    assert any("TAK-062" in k for k in NOT_AVAILABLE)
    assert any("ALV003-1021" in k for k in NOT_AVAILABLE)


def test_every_comparison_is_sourced_to_a_table():
    for c in COMPARISONS:
        assert c.source.quote and c.source.url.startswith("http")
        # Both endpoints must come from the same arms, or the ratio is meaningless.
        assert c.vhcd.placebo_n == c.iel.placebo_n
        assert c.vhcd.drug_n == c.iel.drug_n
