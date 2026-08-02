"""Tests for the prospective predictions.

These are unusual tests. They are not checking that code works — they are pinning claims
made in public before the trials reported, so that the claims cannot quietly drift toward
whatever the results turn out to be. If a prediction here changes, it should be because
someone deliberately changed it, in a commit, with a reason.
"""

from __future__ import annotations

import pytest

from ctsim.prospective import (
    AMLITELIMAB,
    BEST_RESTORATION_EFFECT,
    GLUTEN_DOSE_RESPONSE,
    LIVE_TRIALS,
    RESTORATION_BENCHMARKS,
    TEV_53408,
    ZED1227_CEC013,
    predict,
)

# --- the predictions themselves ---------------------------------------------

def test_tev53408_would_miss_the_best_known_drug():
    """Pinned prediction: reads out 2026-09."""
    p = predict(TEV_53408)
    assert p.trial.design == "prevention"
    assert p.min_protection is not None
    assert p.min_protection > 0.79, (
        "prediction is that TEV-53408 needs more protection than ZED1227 delivered"
    )
    assert 0.80 < p.min_protection < 0.95


def test_both_restoration_trials_need_more_than_has_ever_been_achieved():
    """Pinned prediction: reads out 2026-08 and 2027-06."""
    for trial in (AMLITELIMAB, ZED1227_CEC013):
        p = predict(trial)
        assert p.trial.design == "restoration"
        assert p.mde_80 > BEST_RESTORATION_EFFECT, (
            f"{trial.name}: MDE {p.mde_80:.3f} is no longer above the best observed "
            f"restoration effect {BEST_RESTORATION_EFFECT}"
        )


def test_larger_trials_resolve_smaller_effects():
    """Sanity: the ordering of MDEs must follow the ordering of sample sizes."""
    by_n = sorted(LIVE_TRIALS, key=lambda t: t.effective_n)
    mdes = [predict(t).mde_80 for t in by_n]
    prevention = [t.design == "prevention" for t in by_n]
    # Only compare within the restoration group, since the prevention arm uses a
    # different (injury-inflated) SD.
    rest = [m for m, isprev in zip(mdes, prevention) if not isprev]
    assert rest == sorted(rest, reverse=True)


# --- design facts that the predictions rest on ------------------------------

def test_gluten_dose_runs_opposite_to_baseline_damage():
    """The observation the audit turns on.

    The trial giving the largest gluten challenge enrols healed patients; the two giving
    the smallest enrol patients who already have atrophy.
    """
    biggest = max(LIVE_TRIALS, key=lambda t: t.gluten_mg_per_day)
    assert biggest.design == "prevention"
    assert all(t.design == "restoration"
               for t in LIVE_TRIALS if t is not biggest)


def test_dose_span_is_at_least_tenfold():
    doses = [t.gluten_mg_per_day for t in LIVE_TRIALS]
    assert max(doses) / min(doses) > 10


def test_sige_doses_sit_below_the_documented_injury_range():
    """SIGE doses are far below the 3 g/day behind every positive VH:CD result.

    This is not a criticism on its own — SIGE exists to suppress placebo improvement,
    not to injure. But it is why the protection-ceiling analysis does not apply to them.
    """
    sige = [t for t in LIVE_TRIALS if t.design == "restoration"]
    assert sige
    for t in sige:
        assert t.gluten_mg_per_day < 500
    prevention = [t for t in LIVE_TRIALS if t.design == "prevention"]
    assert all(t.gluten_mg_per_day >= 3000 for t in prevention)


def test_dose_response_is_monotonic():
    doses = [d for d, _v, _n in GLUTEN_DOSE_RESPONSE]
    injuries = [v for _d, v, _n in GLUTEN_DOSE_RESPONSE]
    assert doses == sorted(doses)
    assert injuries == sorted(injuries, reverse=True), "more gluten must mean more injury"


def test_only_one_trial_powers_on_the_histology_endpoint():
    powered = [t for t in LIVE_TRIALS if t.vhcd_is_powered]
    assert len(powered) == 1 and powered[0] is AMLITELIMAB, (
        "two of three run VH:CD without powering for it; if that changes, the audit's "
        "framing needs revisiting"
    )


# --- provenance -------------------------------------------------------------

def test_every_live_trial_carries_a_protocol_source():
    for t in LIVE_TRIALS:
        assert t.sources, t.name
        for s in t.sources:
            assert s.url.startswith("https://euclinicaltrials.eu"), (
                f"{t.name}: design parameters must come from the protocol on CTIS, not "
                "from ClinicalTrials.gov, which does not carry them"
            )
            assert len(s.quote) > 80


def test_restoration_benchmarks_include_the_negative_case():
    """The benchmark list must not quietly become a list of successes."""
    values = [v for v, _n in RESTORATION_BENCHMARKS]
    assert any(v < 0 for v in values), "TAK-062's wrong-direction result must stay in"
    assert max(values) >= BEST_RESTORATION_EFFECT


def test_effective_n_handles_unequal_allocation():
    """Dr Falk randomises 1.5:1, so per-arm N overstates what the contrast buys."""
    assert ZED1227_CEC013.n_per_arm == 72 and ZED1227_CEC013.n_control == 48
    assert ZED1227_CEC013.effective_n == pytest.approx(57.6, abs=0.1)
    assert ZED1227_CEC013.effective_n < ZED1227_CEC013.n_per_arm
