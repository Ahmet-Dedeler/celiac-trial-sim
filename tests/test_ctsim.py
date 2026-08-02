"""Regression tests. The important one is that Monte Carlo agrees with the closed form —
if those diverge, the simulator is wrong and every downstream claim is too.
"""

from __future__ import annotations

import math

import pytest

from ctsim.model import LIT, decompose_vhcd, load_empirical, pooled_sd
from ctsim.simulate import analytic_power, mde, required_n, simulate_trial


def test_pooled_sd_is_finite_and_plausible():
    rows = load_empirical()
    assert len(rows) >= 4, "expected several arms of empirical VH:CD data"
    sd = pooled_sd(rows)
    # Observed per-arm SDs across trials run 0.5-0.95; the pool must land in range.
    assert 0.4 < sd < 1.1


def test_mc_matches_analytic():
    """Monte Carlo power must track the normal-approximation closed form."""
    sd = 0.74
    for n in (25, 60, 100):
        a = analytic_power(n, 0.40, sd)
        m = simulate_trial(n, 0.40, sd, n_sims=40_000, seed=3).power
        # MC uses a t-test (slightly conservative vs the normal approximation).
        assert abs(a - m) < 0.03, f"n={n}: analytic={a:.3f} mc={m:.3f}"


def test_required_n_inverts_mde():
    sd = 0.74
    for n in (20, 50, 120):
        effect = mde(n, sd)
        assert required_n(effect, sd) == pytest.approx(n, abs=1)


def test_power_is_monotonic_in_n():
    sd = 0.74
    powers = [analytic_power(n, 0.4, sd) for n in (10, 25, 50, 100, 200)]
    assert powers == sorted(powers)


def test_variance_decomposition_conserves_variance():
    dec = decompose_vhcd(0.74)
    assert dec.reader_sd**2 + dec.residual_sd**2 == pytest.approx(0.74**2, rel=1e-9)
    assert dec.reader_share + dec.residual_share == pytest.approx(1.0, rel=1e-9)


def test_reader_error_is_a_minority_of_variance():
    """The headline qualitative claim: pathologist disagreement is not the main problem."""
    dec = decompose_vhcd(pooled_sd(load_empirical()))
    assert dec.reader_share < 0.25


def test_averaging_biopsies_only_helps_if_sampling_variance_exists():
    base = simulate_trial(50, 0.4, 0.74, n_biopsies=8, sampling_share=0.0,
                          n_sims=20_000, seed=5).power
    more = simulate_trial(50, 0.4, 0.74, n_biopsies=8, sampling_share=0.4,
                          n_sims=20_000, seed=5).power
    assert more > base + 0.03


def test_meaningful_change_constant_is_sourced():
    c = LIT["vhcd_clinically_significant"]
    assert c.value == 0.40 and "Taavela" in c.source


def test_null_effect_gives_nominal_type_i_error():
    r = simulate_trial(60, 0.0, 0.74, n_sims=40_000, seed=9)
    assert 0.035 < r.power < 0.065, f"type I error off nominal: {r.power}"
