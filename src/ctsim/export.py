"""Export computed results to JSON for the web page and for downstream reuse."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from ctsim import published, variance
from ctsim.assay import assay_sensitivity
from ctsim.model import LIT, decompose_vhcd
from ctsim.simulate import analytic_power, mde, observed_sd_for, required_n
from ctsim.variance import (
    HEADLINE_HOLDOUT,
    all_arms,
    biopsy_n_savings,
    fit_injury_variance,
    holdout_predict,
    required_n as injury_required_n,
)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "curated" / "results.json"

TRIAL_META = {
    "NCT03409796": ("Gluten challenge dose-finding", "3 g vs 10 g gluten, 14 days", "academic"),
    "NCT05353985": ("TAK-062 (zamaglutenase)", "Takeda — Phase 2, discontinued", "glutenase"),
    "NCT06001177": ("KAN-101 SynCeD", "Anokion — Phase 2a, program shelved", "tolerance"),
    "NCT03738475": ("TAK-101 Phase 2a", "COUR/Takeda — p=0.08 on VH:CD", "tolerance"),
    "NCT02637141": ("AMG 714", "Amgen/Celimmune — Phase 2a", "anti-IL-15"),
    "NCT02633020": ("AMG 714 (RCD-II)", "Amgen/Celimmune — Phase 2a", "anti-IL-15"),
    "NCT04424927": ("PRV-015 (ordesekimab)", "Provention/Sanofi — Phase 2b", "anti-IL-15"),
    "NCT03585478": ("IMGX003 CeliacShield", "Immunogenics — Phase 2", "glutenase"),
    "NCT01917630": ("CeliAction", "Alvine — Phase 2b latiglutenase", "glutenase"),
    "EudraCT2017-002241-30": ("ZED1227 (CEC-3)", "Dr Falk/Zedira — only histologic hit", "tg2_inhibitor"),
}


def _protection_table() -> list[dict]:
    """Headline delivered-vs-needed rows for the three challenge trials."""
    zed_delta, _lo, _hi = published.zed1227_control_injury()
    zed_sd, _ = published.zed1227_residual_sd()
    zed_delivered = 0.48 / abs(zed_delta)  # 100 mg arm vs placebo
    rows = [
        {
            "trial": "ZED1227 (CEC-3)",
            "nct_id": "EudraCT2017-002241-30",
            "n_per_arm": 34,
            "injury": abs(zed_delta),
            "sd": round(zed_sd, 4),
            "protection_delivered": round(zed_delivered, 3),
            "protection_needed": round(mde(34, zed_sd) / abs(zed_delta), 3),
            "reported": "p < 0.001",
            "hit": True,
        }
    ]
    # TAK-101 and KAN-101 from assay_sensitivity
    for nct, delivered, reported, hit in (
        ("NCT03738475", (0.63 - 0.18) / 0.63, "p = 0.08", False),
        ("NCT06001177", (0.61 - 0.85) / 0.61, "null", False),
    ):
        a = next(x for x in assay_sensitivity() if x.nct_id == nct)
        rows.append({
            "trial": TRIAL_META.get(nct, (nct, "", ""))[0],
            "nct_id": nct,
            "n_per_arm": a.n_per_arm,
            "injury": round(abs(a.control_delta), 4),
            "sd": round(a.sd_used, 4),
            "protection_delivered": round(delivered, 3),
            "protection_needed": round(a.min_protection, 3),
            "reported": reported,
            "hit": hit,
        })
    return rows


def build() -> dict:
    rows = all_arms()
    model = fit_injury_variance(rows)
    typical_injury = 0.61
    sd = model.sd_at(typical_injury)
    meaningful = LIT["vhcd_clinically_significant"].value
    dec = decompose_vhcd(sd)
    ho = holdout_predict(HEADLINE_HOLDOUT, rows)

    assumptions = []
    for a in published.ASSUMPTIONS:
        entry = {
            "trial": a.trial,
            "nct_id": a.nct_id,
            "target_effect": a.target_effect,
            "assumed_sd": a.assumed_sd,
            "n_per_arm": a.n_per_arm,
            "claimed_power": a.claimed_power,
            "endpoint_role": a.endpoint_role,
            "units": a.units,
            "source_url": a.source.url,
            "quote": a.source.quote,
        }
        if a.units == "vhcd_ratio":
            obs, _df, _scale = observed_sd_for(a.nct_id)
            entry["observed_sd"] = round(obs, 4)
            entry["actual_power"] = round(
                analytic_power(a.n_per_arm, a.target_effect, obs), 4
            )
        else:
            entry["observed_sd"] = None
            entry["actual_power"] = None
        assumptions.append(entry)

    challenge_table = []
    for injury in (0.20, 0.61, 1.00, 1.53, 2.50):
        d = injury_required_n(injury, 0.50, model)
        challenge_table.append({
            "injury": d.injury,
            "sd": round(d.sd, 4),
            "n_per_arm": d.n_per_arm,
            "n_per_arm_constant_sd": d.n_per_arm_constant_sd,
        })

    biopsy = []
    for k in (1, 2, 4, 8):
        n1, nk, _sd1, sdk = biopsy_n_savings(0.61, 0.50, n_from=1, n_to=k, model=model)
        biopsy.append({
            "n_biopsies": k,
            "sd": round(sdk, 4),
            "n_per_arm_50pct": nk,
            "saved_vs_single": n1 - nk,
        })

    # Replay each prevention trial at its own observed SD (not a global pool).
    replay = []
    for nct, n in (("NCT03738475", 13), ("NCT06001177", 25), ("NCT03409796", 7),
                   ("NCT03585478", 22)):
        obs, _df, _scale = observed_sd_for(nct)
        if obs != obs:  # NaN
            continue
        replay.append({
            "nct_id": nct,
            "label": nct,
            "trial": TRIAL_META.get(nct, ("", "", ""))[0],
            "sponsor_note": TRIAL_META.get(nct, ("", "", ""))[1],
            "n_per_arm": n,
            "sd_used": round(obs, 4),
            "mde_80": round(mde(n, obs), 4),
            "ratio_to_meaningful": round(mde(n, obs) / meaningful, 3),
            "power_for_meaningful": round(analytic_power(n, meaningful, obs), 4),
            "power_for_half": round(analytic_power(n, meaningful / 2, obs), 4),
        })

    return {
        "generated_by": "ctsim.export",
        "endpoint": "change from baseline in villous height : crypt depth ratio (VH:CD)",
        # Kept for back-compat with older page code; this is the model SD at a
        # typical challenge injury (0.61), NOT a pooled constant.
        "pooled_sd": round(sd, 4),
        "typical_injury": typical_injury,
        "n_arms": model.n_arms,
        "n_trials": model.n_trials,
        "clinically_meaningful_change": meaningful,
        "noise_to_signal": round(sd / meaningful, 3),
        "injury_variance_model": {
            "floor": round(model.floor, 4),
            "slope": round(model.slope, 4),
            "floor_se": round(model.floor_se, 4),
            "slope_se": round(model.slope_se, 4),
            "r": round(model.r, 4),
            "formula": f"SD = {model.floor:.3f} + {model.slope:.3f} × |injury|",
            "n_arms": model.n_arms,
            "n_trials": model.n_trials,
        },
        "variance_decomposition": {
            "total_sd": round(dec.total_sd, 4),
            "reader_sd": round(dec.reader_sd, 4),
            "residual_sd": round(dec.residual_sd, 4),
            "reader_share": round(dec.reader_share, 4),
            "residual_share": round(dec.residual_share, 4),
            "at_injury": typical_injury,
        },
        "measured_variance_shares": dict(published.MEASURED_VARIANCE_SHARES),
        "empirical_arms": [
            {
                "nct_id": r.nct_id,
                "trial": TRIAL_META.get(r.nct_id, ("", "", ""))[0],
                "arm": r.arm_label,
                "n": r.n_arm,
                "delta": r.delta,
                "sd": round(r.sd, 4),
                "sd_source": r.sd_source,
                "sd_fitted": round(model.sd_at(r.delta), 4),
            }
            for r in sorted(rows, key=lambda r: abs(r.delta))
        ],
        "challenge_dose_table": challenge_table,
        "biopsy_averaging": biopsy,
        "holdout": {
            "held_out": list(ho.held_out),
            "mae": round(ho.mae, 4),
            "max_abs_err": round(ho.max_abs_err, 4),
            "train_floor": round(ho.model.floor, 4),
            "train_slope": round(ho.model.slope, 4),
            "predictions": [
                {
                    "nct_id": p.nct_id,
                    "arm": p.arm_label,
                    "injury": round(p.injury, 4),
                    "sd_observed": round(p.sd_observed, 4),
                    "sd_predicted": round(p.sd_predicted, 4),
                    "abs_err": round(p.abs_err, 4),
                }
                for p in ho.predictions
            ],
        },
        "assumptions": assumptions,
        "not_powered_on_vhcd": [
            {"nct_id": nct, "quote": src.quote, "url": src.url}
            for nct, src in published.NOT_POWERED_ON_VHCD.items()
        ],
        "protection_table": _protection_table(),
        "required_n_per_arm": {
            str(round(eff, 2)): required_n(eff, sd)
            for eff in (0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50)
        },
        "replay": replay,
        "power_curve": [
            {
                "n_per_arm": n,
                "mde_80": round(mde(n, sd), 4),
                "power_at_0_40": round(analytic_power(n, 0.40, sd), 4),
                "power_at_0_20": round(analytic_power(n, 0.20, sd), 4),
            }
            for n in (10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400)
        ],
        "literature_constants": {
            k: {"value": c.value, "units": c.units, "source": c.source, "note": c.note}
            for k, c in LIT.items()
        },
    }


def main() -> None:
    data = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, indent=2))
    print(f"wrote {OUT}")
    m = data["injury_variance_model"]
    print(f"  model: {m['formula']}  r={m['r']}  arms={m['n_arms']}")
    print(f"  holdout MAE={data['holdout']['mae']}  "
          f"assumptions={len(data['assumptions'])}")


if __name__ == "__main__":
    main()
