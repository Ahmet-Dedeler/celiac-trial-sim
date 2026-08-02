"""Tests for hand-entered published data.

Everything in `ctsim.published` was typed in from a PDF, so it has a failure mode the
scraped dataset does not: a transcription error that is silently plausible. These tests
exist to catch that, mostly by checking the paper's numbers against each other.
"""

from __future__ import annotations

import math

import pytest

from ctsim.published import (
    ASSUMPTIONS,
    ZED1227_ARMS,
    ZED1227_LEVELS,
    zed1227_control_injury,
    zed1227_residual_sd,
    zed1227_rows,
)
from ctsim.simulate import analytic_power

# --- internal consistency of the transcribed ZED1227 table -------------------

def test_lsmean_differences_reconcile_with_the_published_contrasts():
    """Table 2 reports both per-arm changes and differences vs placebo.

    They must subtract. This catches a mistyped digit in any single cell, which is the
    realistic transcription failure.
    """
    changes = {label: delta for label, _n, delta, _lo, _hi, _p in ZED1227_ARMS}
    placebo = changes["Placebo"]
    published_differences = {
        "ZED1227 10 mg": 0.44,
        "ZED1227 50 mg": 0.49,
        "ZED1227 100 mg": 0.48,
    }
    for arm, expected in published_differences.items():
        assert changes[arm] - placebo == pytest.approx(expected, abs=0.005), arm


def test_confidence_intervals_bracket_their_point_estimates():
    for label, _n, delta, lo, hi, _p in ZED1227_ARMS:
        assert lo < delta < hi, label
        assert lo < hi, label


def test_levels_and_changes_point_the_same_way():
    """A negative LS-mean change must correspond to a fall from baseline to week 6."""
    changes = {label: delta for label, _n, delta, _lo, _hi, _p in ZED1227_ARMS}
    for label, _n, base, _bsd, wk6, _wsd in ZED1227_LEVELS:
        assert (wk6 - base < 0) == (changes[label] < 0), label
        # The raw and model-adjusted changes should not be wildly different.
        assert abs((wk6 - base) - changes[label]) < 0.30, label


def test_placebo_is_the_worst_arm():
    """The trial's whole claim is that the drug arms deteriorated less than placebo."""
    changes = {label: delta for label, _n, delta, _lo, _hi, _p in ZED1227_ARMS}
    placebo = changes.pop("Placebo")
    assert all(v > placebo for v in changes.values())


def test_analysed_n_never_exceeds_randomised_n():
    randomised = {"Placebo": 40, "ZED1227 10 mg": 41,
                  "ZED1227 50 mg": 41, "ZED1227 100 mg": 41}
    for label, n, _d, _lo, _hi, _p in ZED1227_ARMS:
        assert n <= randomised[label], label


# --- the derived residual SD -------------------------------------------------

def test_per_arm_residual_sds_agree_which_is_why_they_collapse_to_one():
    """All four CIs come from one model with a common residual variance.

    If they disagreed, treating them as one estimate would be wrong. Verifying they
    agree is what licenses the collapse — and stops anyone "fixing" it into four rows
    and quadrupling the apparent degrees of freedom.
    """
    z = 1.959963985
    per_arm = [((hi - lo) / (2 * z)) * math.sqrt(n)
               for _label, n, _d, lo, hi, _p in ZED1227_ARMS]
    assert max(per_arm) - min(per_arm) < 0.05, per_arm


def test_residual_sd_is_in_the_range_other_trials_report():
    sd, df = zed1227_residual_sd()
    assert 0.30 < sd < 1.10
    assert df == sum(n for _l, n, _d, _lo, _hi, _p in ZED1227_ARMS) - len(ZED1227_ARMS) - 1


def test_zed1227_row_carries_its_provenance_and_weight():
    (row,) = zed1227_rows()
    assert row.is_placebo and row.design == "prevention"
    assert row.sd_source == "from_lsmean_ci95"
    _sd, df = zed1227_residual_sd()
    assert row.n_arm - 1 == df, "row weight must reproduce the model's degrees of freedom"


def test_control_injury_matches_the_placebo_arm():
    delta, lo, hi = zed1227_control_injury()
    assert (delta, lo, hi) == (ZED1227_ARMS[0][2], ZED1227_ARMS[0][3], ZED1227_ARMS[0][4])
    assert delta < 0, "the placebo arm must have lost villous height"


# --- design assumptions ------------------------------------------------------

def test_every_assumption_carries_a_quote_and_a_url():
    assert ASSUMPTIONS
    for a in ASSUMPTIONS:
        assert a.source.url.startswith("http"), a.trial
        assert len(a.source.quote) > 40, a.trial
        assert a.source.citation, a.trial


def test_assumptions_reproduce_the_power_the_sponsors_claimed():
    """Recompute each trial's stated power from its own stated inputs.

    A mismatch means either a transcription error or that the sponsor used a design this
    two-sample formula does not describe — both worth knowing before quoting the number.
    """
    tolerances = {
        # Group-sequential with a triangular boundary; the sequential penalty costs
        # real power, so the fixed-sample formula reads high. Wide tolerance on purpose.
        "ZED1227 (CEC-3)": 0.12,
        "KAN-101 SynCeD": 0.03,
    }
    for a in ASSUMPTIONS:
        recomputed = analytic_power(a.n_per_arm, a.target_effect, a.assumed_sd)
        assert recomputed >= a.claimed_power - 0.02, (
            f"{a.trial}: recomputed {recomputed:.2f} below claimed {a.claimed_power:.2f}"
        )
        assert abs(recomputed - a.claimed_power) <= tolerances[a.trial], (
            f"{a.trial}: recomputed {recomputed:.2f} vs claimed {a.claimed_power:.2f}"
        )


def test_kan101_assumed_less_noise_than_it_measured():
    """The finding this repo now leads with, pinned so it cannot silently drift."""
    from ctsim.simulate import observed_sd_for

    a = next(x for x in ASSUMPTIONS if x.nct_id == "NCT06001177")
    observed, _df, _scale = observed_sd_for(a.nct_id)
    assert observed > a.assumed_sd * 1.2, "expected the assumption to be optimistic"
    assert analytic_power(a.n_per_arm, a.target_effect, observed) < a.claimed_power - 0.10


def test_zed1227_assumed_more_noise_than_it_measured():
    a = next(x for x in ASSUMPTIONS if x.nct_id.startswith("EudraCT"))
    observed, _df = zed1227_residual_sd()
    assert observed < a.assumed_sd, "expected the assumption to be conservative"


# --- paper-extracted arms ----------------------------------------------------

def test_every_paper_trial_has_a_source_with_a_quote():
    from ctsim.published import PAPER_ARM_SOURCES, PAPER_ARMS

    for _trial, nct, *_rest in PAPER_ARMS:
        src = PAPER_ARM_SOURCES.get(nct)
        assert src is not None, f"{nct} has no source"
        assert src.url.startswith("http") and len(src.quote) > 40


def test_paper_arm_values_are_physically_plausible():
    """VH:CD is a ratio of two lengths; a change score and its SD have hard limits."""
    from ctsim.published import PAPER_ARMS

    for trial, _nct, arm, n, delta, sd, _placebo, _design in PAPER_ARMS:
        assert n > 0, f"{trial} {arm}"
        assert abs(delta) < 2.0, f"{trial} {arm}: change of {delta} is not a change score"
        assert 0.0 < sd < 2.0, f"{trial} {arm}: SD {sd} out of range"


def test_paper_rows_merge_without_colliding_with_registry_rows():
    from ctsim.model import load_empirical
    from ctsim.published import paper_rows

    registry = {(r.nct_id, r.arm_label) for r in load_empirical()}
    for r in paper_rows():
        assert (r.nct_id, r.arm_label) not in registry, (
            f"{r.nct_id} {r.arm_label} would be counted twice"
        )
        assert r.sd_source == "paper_posted_sd"


def test_measured_variance_shares_are_shares():
    from ctsim.published import MEASURED_VARIANCE_SHARES, MEASURED_VARIANCE_SOURCE

    assert sum(MEASURED_VARIANCE_SHARES.values()) <= 1.0
    assert all(0 <= v <= 1 for v in MEASURED_VARIANCE_SHARES.values())
    # The direction is the load-bearing part: reader is the smallest term by far.
    assert MEASURED_VARIANCE_SHARES["reader"] < MEASURED_VARIANCE_SHARES["biopsy"]
    assert MEASURED_VARIANCE_SHARES["biopsy"] < MEASURED_VARIANCE_SHARES["patient"]
    assert MEASURED_VARIANCE_SOURCE.url.startswith("http")


def test_unavailable_data_is_recorded_rather_than_forgotten():
    from ctsim.published import NOT_AVAILABLE

    assert len(NOT_AVAILABLE) >= 4
    for trial, why in NOT_AVAILABLE.items():
        assert len(why) > 60, f"{trial}: reason is too thin to stop a re-search"
