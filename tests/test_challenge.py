"""Tests for the dose-and-duration model of gluten-challenge injury.

Most of these pin claims rather than check arithmetic. The point of writing them down
is that the three live trials read out between now and 2027, and a claim that can drift
toward whatever the results turn out to be is not a claim.
"""

from __future__ import annotations

import math

import pytest

from ctsim.challenge import (
    CHALLENGE_ARMS,
    LAHDEAHO_IPD,
    SIGE_ARMS,
    THREE_GRAM_SERIES,
    VITAL_WHEAT_GLUTEN_PROTEIN_FRACTION,
    fit_duration_model,
    fit_long_challenge_dose_response,
    lahdeaho_ipd_fit,
    required_n,
    sige_injury_on_healed_mucosa,
)
from ctsim.variance import fit_injury_variance


# --- the dataset ------------------------------------------------------------

def test_every_arm_carries_a_quote_and_a_url():
    for a in CHALLENGE_ARMS:
        assert a.source.quote.strip(), f"{a.label} has no verbatim quote"
        assert a.source.url.startswith("http"), f"{a.label} has no retrievable source"
        assert a.source.locator.strip(), f"{a.label} says where but not exactly where"


def test_doses_are_all_on_the_protein_scale():
    """A dose stated in vital wheat gluten must be converted, not copied across.

    This is the trap that makes cross-trial dose comparisons wrong by ~30%. Sanofi
    states 250 mg of vital wheat gluten; the record here must not carry 250.
    """
    sanofi = next(s for s in SIGE_ARMS if "Sanofi" in s[0])
    assert sanofi[1] == pytest.approx(0.250 * VITAL_WHEAT_GLUTEN_PROTEIN_FRACTION)
    assert sanofi[1] < 0.250, "vital wheat gluten is heavier than the protein in it"


def test_tak101_dose_matches_its_own_protocol_schedule():
    """12 g/day for 3 days then 6 g/day for 11 days, averaged over the 14."""
    arm = next(a for a in CHALLENGE_ARMS if a.label.startswith("TAK-101"))
    assert arm.days == 14
    assert arm.cumulative_g == pytest.approx(12 * 3 + 6 * 11)
    assert arm.g_protein_per_day == pytest.approx(102.0 / 14.0)


def test_kan101_and_zed1227_agree_at_equal_cumulative_dose():
    """The coincidence that used to be misread as a replication at 3 g/day.

    Same cumulative exposure by two very different routes, same injury. Worth pinning
    because it is the single best piece of evidence that cumulative exposure, not daily
    dose, is what the mucosa integrates — and because misreading it caused a real error.
    """
    kan = next(a for a in CHALLENGE_ARMS if a.label.startswith("KAN-101"))
    zed = next(a for a in CHALLENGE_ARMS if a.label.startswith("ZED1227"))
    assert kan.g_protein_per_day == 9.0 and kan.days == 14
    assert zed.g_protein_per_day == 3.0 and zed.days == 42
    assert kan.cumulative_g == pytest.approx(zed.cumulative_g)
    assert kan.injury == pytest.approx(zed.injury, abs=0.01)


# --- duration is the lever --------------------------------------------------

def test_the_three_gram_series_is_three_independent_studies():
    arms = [a for a in CHALLENGE_ARMS if a.label in THREE_GRAM_SERIES]
    assert len(arms) == 3
    assert {a.days for a in arms} == {14, 42, 78}, "one dose, three durations"
    for a in arms:
        assert 2.9 <= a.g_protein_per_day <= 3.2, f"{a.label} is not a 3 g/day arm"


def test_injury_accumulates_with_duration_at_a_fixed_dose():
    m = fit_duration_model()
    assert m.per_day > 0.010, "the slope is the whole finding"
    assert m.per_day < 0.025
    # Ordering must hold even though the fit itself is only three points.
    assert m.injury_at(14) < m.injury_at(42) < m.injury_at(78)


def test_a_fortnight_at_three_grams_is_not_a_detectable_challenge():
    """The indictment of the field's default short challenge, stated as a number."""
    fortnight = required_n(14)
    six_weeks = required_n(42)
    assert fortnight.n_per_arm > 300
    assert six_weeks.n_per_arm < 100
    assert fortnight.n_per_arm > 5 * six_weeks.n_per_arm


def test_doubling_the_challenge_beats_doubling_nothing():
    six, twelve = required_n(42), required_n(84)
    assert twelve.n_per_arm < six.n_per_arm / 2, (
        "extending 6 weeks to 12 should more than halve N; if it stops doing so the "
        "variance model's slope has changed and the design advice needs revisiting"
    )


def test_dose_explains_almost_nothing_within_one_protocol():
    fit = lahdeaho_ipd_fit()
    assert fit.n == 21
    assert abs(fit.r_dose) < 0.35, "a fourfold dose range barely moves the outcome"
    assert fit.r2_joint < 0.15, "design explains little; the patient explains the rest"
    # And the residual spread is large enough to matter for powering.
    assert fit.injury_sd > 0.9


def test_lahdeaho_ipd_reproduces_its_published_summaries():
    """Guards the transcription of Table 2 against the paper's own reported means."""
    doses = [r[0] for r in LAHDEAHO_IPD]
    days = [r[1] for r in LAHDEAHO_IPD]
    assert len(LAHDEAHO_IPD) == 21, "21 completers reached the endpoint"
    assert min(doses) == pytest.approx(1.3) and max(doses) == pytest.approx(5.0)
    assert min(days) == 29 and max(days) == 103
    assert sum(doses) / len(doses) == pytest.approx(3.1, abs=0.05)


# --- SIGE -------------------------------------------------------------------

def test_sige_doses_sit_high_on_the_dose_response():
    """'Simulated inadvertent exposure' is not, by dose, an inadvertent exposure."""
    _imax, d50, _h = fit_long_challenge_dose_response()
    assert 0.02 < d50 < 0.20, "half-maximal dose in g/day of gluten protein"
    for label, g, _days, _entry, _stated in SIGE_ARMS:
        if g is None:
            continue
        assert g > 2 * d50, f"{label} is above twice the half-maximal dose"


def test_the_dose_response_passes_through_its_anchors():
    """Three parameters through three points: it must interpolate them exactly."""
    imax, d50, h = fit_long_challenge_dose_response()
    for label in ("Catassi 2007, 10 mg/day", "Catassi 2007, 50 mg/day",
                  "Lahdeaho 2011, 3.1 g/day mean"):
        a = next(x for x in CHALLENGE_ARMS if x.label == label)
        d = a.g_protein_per_day
        assert imax * d**h / (d50**h + d**h) == pytest.approx(a.injury, abs=1e-3)


def test_sige_prediction_only_applies_to_a_healed_mucosa():
    """The predicted injuries exceed what an atrophic entry criterion leaves available.

    That gap is the argument, not a bug: the dose-response comes entirely from healed
    patients, and every SIGE trial enrols patients under Vh:Cd 2.5. Pinned so the two
    populations never get quietly merged.
    """
    sanofi = next(s for s in SIGE_ARMS if "Sanofi" in s[0])
    implied = sige_injury_on_healed_mucosa(sanofi[1], sanofi[2])
    assert implied > 0.6, "on healed mucosa this is a substantial challenge"
    assert all("2.5" in s[3] for s in SIGE_ARMS), (
        "all three SIGE trials enrol atrophic patients; if one stops doing so it "
        "becomes the experiment that separates dose from floor effect"
    )


def test_flat_control_arm_costs_sample_size_if_sige_actually_injures():
    var = fit_injury_variance()
    from ctsim.simulate import mde

    flat = mde(34, var.sd_at(0.0))
    injured = mde(34, var.sd_at(0.4))
    assert injured > flat * 1.2, (
        "if SIGE injures, a restoration design needs a bigger difference than it "
        "planned for, and Sanofi's 34-vs-34 is the contrast that pays for it"
    )


# --- the composed calculator ------------------------------------------------

def test_required_n_chain_is_internally_consistent():
    d = required_n(42, protection=0.50)
    var = fit_injury_variance()
    assert d.sd == pytest.approx(var.sd_at(d.injury))
    # The MDE it settles on must be within reach of the effect it is chasing.
    assert d.mde_80 <= 0.50 * d.injury + 1e-9
    assert math.isfinite(d.protection_needed)


def test_harder_targets_cost_more_patients():
    easy = required_n(42, protection=0.80)
    hard = required_n(42, protection=0.30)
    assert hard.n_per_arm > easy.n_per_arm
