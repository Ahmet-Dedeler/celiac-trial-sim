"""Tests for the assay-sensitivity layer.

These guard the strongest claim in the repo: that some of these trials could not have
detected their drug at any efficacy, because the control arm never sustained the injury
the drug was supposed to prevent.
"""

from __future__ import annotations

import math

import pytest

from ctsim.assay import assay_sensitivity, challenge_dose_options
from ctsim.model import (
    EmpiricalSD,
    ancova_variance_ratio,
    baseline_correlation,
    load_empirical,
    pooled_sd,
)
from ctsim.simulate import analytic_power, mde


def _arm(nct, label, n, delta, sd, placebo=False, mech="tolerance",
         design="prevention"):
    return EmpiricalSD(nct_id=nct, arm_label=label, n_arm=n, delta=delta, sd=sd,
                       sd_source="posted_sd", time_frame="", is_placebo=placebo,
                       mechanism=mech, design=design)


# --- the core identity -------------------------------------------------------

def test_min_protection_is_mde_over_injury():
    rows = [_arm("NCT1", "placebo", 40, -0.80, 0.74, placebo=True),
            _arm("NCT1", "drug", 40, -0.40, 0.74)]
    a = assay_sensitivity(rows, sd=0.74)[0]
    assert a.min_protection == pytest.approx(mde(40, 0.74) / 0.80)


def test_a_drug_at_exactly_the_min_protection_has_80pc_power():
    """The metric has to mean what it says: f_min * injury == the 80%-power effect."""
    rows = [_arm("NCT1", "placebo", 50, -1.20, 0.74, placebo=True),
            _arm("NCT1", "drug", 50, -0.60, 0.74)]
    a = assay_sensitivity(rows, sd=0.74)[0]
    effect_at_fmin = a.min_protection * abs(a.control_delta)
    assert analytic_power(a.n_per_arm, effect_at_fmin, 0.74) == pytest.approx(0.80, abs=0.005)


def test_a_flat_control_arm_makes_the_trial_infeasible():
    """No injury in the control arm means no prevention effect is available."""
    rows = [_arm("NCT1", "placebo", 60, -0.006, 0.77, placebo=True),
            _arm("NCT1", "drug", 60, -0.34, 0.76)]
    a = assay_sensitivity(rows, sd=0.74)[0]
    assert a.min_protection > 1.0
    assert not a.feasible


def test_a_large_injury_makes_a_small_trial_feasible():
    rows = [_arm("NCT1", "placebo", 25, -1.53, 0.74, placebo=True),
            _arm("NCT1", "drug", 25, -0.50, 0.74)]
    a = assay_sensitivity(rows, sd=0.74)[0]
    assert a.min_protection < 0.50
    assert a.feasible


def test_best_case_is_never_harsher_than_the_point_estimate():
    for a in assay_sensitivity():
        assert a.min_protection_best_case <= a.min_protection


def test_trials_without_a_control_arm_are_skipped():
    """A dose-comparison study has no 'untreated injury', so the metric must not fire."""
    rows = [_arm("NCT9", "3g gluten", 7, -0.06, 0.52,
                 mech="gluten_challenge_methodology"),
            _arm("NCT9", "10g gluten", 7, -1.53, 0.94,
                 mech="gluten_challenge_methodology")]
    assert assay_sensitivity(rows, sd=0.74) == []


# --- regression guards on the real trials ------------------------------------

def test_restoration_trials_are_excluded_from_the_protection_metric():
    """TAK-062's placebo arm barely moved (-0.006), which naively reads as "no signal
    was available". It is a restoration design — patients enrolled *with* atrophy and
    the drug's job was to let it heal, so a flat control arm is not a missing signal
    and the ceiling is healing headroom nobody posts. Applying the prevention metric
    would declare it impossible for the wrong reason."""
    assert not any(a.nct_id == "NCT05353985" for a in assay_sensitivity())
    rows = [r for r in load_empirical() if r.nct_id == "NCT05353985"]
    assert rows and all(r.design == "restoration" for r in rows)


def test_only_prevention_designs_get_a_protection_ceiling():
    for a in assay_sensitivity():
        rows = [r for r in load_empirical() if r.nct_id == a.nct_id]
        assert all(r.design == "prevention" for r in rows)


def test_kan101_needed_a_near_perfect_drug():
    a = next(x for x in assay_sensitivity() if x.nct_id == "NCT06001177")
    assert 0.85 < a.min_protection <= 1.05


# --- the challenge-protocol lever -------------------------------------------

def test_challenge_options_exclude_medicated_arms():
    opts = challenge_dose_options(25)
    assert opts, "expected at least one unmedicated challenge arm"
    labels = {o.label.lower() for o in opts}
    assert not any("kan-101" in x or "tak-062" in x for x in labels)


def test_stronger_challenge_needs_a_smaller_trial():
    """The design lever: injury magnitude buys more than sample size does."""
    opts = sorted(challenge_dose_options(25), key=lambda o: o.injury)
    ns = [o.n_for_half_protection for o in opts]
    assert ns == sorted(ns, reverse=True), "bigger injury must require fewer patients"


def test_required_n_scales_as_inverse_square_of_injury():
    opts = {o.label: o for o in challenge_dose_options(25)}
    a, b = sorted(opts.values(), key=lambda o: o.injury)[-2:]
    ratio_n = a.n_for_half_protection / b.n_for_half_protection
    ratio_injury_sq = (b.injury / a.injury) ** 2
    assert ratio_n == pytest.approx(ratio_injury_sq, rel=0.05)


# --- baseline correlation / ANCOVA ------------------------------------------

def test_ancova_ratio_bounds():
    assert ancova_variance_ratio(0.0) == pytest.approx(0.5)
    assert ancova_variance_ratio(1.0) == pytest.approx(1.0)
    ratios = [ancova_variance_ratio(r) for r in (0.0, 0.3, 0.6, 0.9)]
    assert ratios == sorted(ratios), "ANCOVA's advantage must shrink as rho rises"


def test_ancova_is_never_worse_than_a_change_score():
    assert all(ancova_variance_ratio(r) <= 1.0 for r in (0.0, 0.25, 0.5, 0.75, 1.0))


def test_recovered_rho_is_a_valid_correlation():
    corrs = baseline_correlation()
    assert corrs, "expected at least one arm posting both baseline and change SDs"
    for c in corrs:
        assert -1.0 <= c.rho <= 1.0, f"{c.nct_id} {c.arm_label}: rho={c.rho}"


def test_rho_recovery_inverts_the_change_score_identity():
    """Round-trip: build a change SD from a known rho, recover it."""
    sigma, rho = 0.80, 0.55
    sd_change = sigma * math.sqrt(2 * (1 - rho))
    assert 1 - sd_change**2 / (2 * sigma**2) == pytest.approx(rho)


# --- the bug this filter used to have ---------------------------------------

def test_positive_change_scores_survive_loading():
    """Baseline rows must be split off by their label, not by thresholding the value.

    A mucosal-healing arm has a positive change score; a value-based filter would drop
    it silently and bias the pool toward deteriorating arms.
    """
    rows = load_empirical()
    assert all("[baseline]" not in r.time_frame.lower() for r in rows)
    # Every retained row is a change score, so none should carry a raw baseline
    # magnitude (~2-3 in VH:CD units).
    assert all(abs(r.delta) < 2.0 for r in rows)
    assert pooled_sd(rows) == pytest.approx(0.740, abs=0.005)
