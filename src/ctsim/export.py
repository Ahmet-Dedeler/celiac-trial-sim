"""Export computed results to JSON for the web page and for downstream reuse."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

from ctsim.model import LIT, decompose_vhcd, load_empirical, pooled_sd
from ctsim.simulate import analytic_power, mde, replay_trials, required_n

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "curated" / "results.json"

# Human-readable labels for the trials that surface on the site.
TRIAL_META = {
    "NCT03409796": ("Gluten challenge dose-finding", "3 g vs 10 g gluten, 14 days", "academic"),
    "NCT05353985": ("TAK-062 (zamaglutenase)", "Takeda — Phase 2, discontinued", "glutenase"),
    "NCT06001177": ("KAN-101 SynCeD", "Anokion — Phase 2a, program shelved", "tolerance"),
    "NCT02637141": ("AMG 714", "Amgen/Celimmune — Phase 2a", "anti-IL-15"),
    "NCT02633020": ("AMG 714 (RCD-II)", "Amgen/Celimmune — Phase 2a", "anti-IL-15"),
    "NCT04424927": ("PRV-015 (ordesekimab)", "Provention/Sanofi — Phase 2b", "anti-IL-15"),
}


def build() -> dict:
    rows = load_empirical()
    sd = pooled_sd(rows)
    meaningful = LIT["vhcd_clinically_significant"].value
    dec = decompose_vhcd(sd)

    return {
        "generated_by": "ctsim.export",
        "endpoint": "change from baseline in villous height : crypt depth ratio (VH:CD)",
        "pooled_sd": round(sd, 4),
        "n_arms": len(rows),
        "n_trials": len({r.nct_id for r in rows}),
        "clinically_meaningful_change": meaningful,
        "noise_to_signal": round(sd / meaningful, 3),
        "variance_decomposition": {
            "total_sd": round(dec.total_sd, 4),
            "reader_sd": round(dec.reader_sd, 4),
            "residual_sd": round(dec.residual_sd, 4),
            "reader_share": round(dec.reader_share, 4),
            "residual_share": round(dec.residual_share, 4),
        },
        "empirical_arms": [
            {
                "nct_id": r.nct_id,
                "trial": TRIAL_META.get(r.nct_id, ("", "", ""))[0],
                "arm": r.arm_label,
                "n": r.n_arm,
                "delta": r.delta,
                "sd": round(r.sd, 4),
                "sd_source": r.sd_source,
            }
            for r in rows
        ],
        "required_n_per_arm": {
            str(round(eff, 2)): required_n(eff, sd)
            for eff in (0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50)
        },
        "replay": [
            {
                **asdict(r),
                "trial": TRIAL_META.get(r.nct_id, ("", "", ""))[0],
                "sponsor_note": TRIAL_META.get(r.nct_id, ("", "", ""))[1],
                "mde_80": round(r.mde_80, 4),
                "ratio_to_meaningful": round(r.ratio_to_meaningful, 3),
                "power_for_meaningful": round(r.power_for_meaningful, 4),
                "power_for_half": round(r.power_for_half, 4),
                "sd_used": round(r.sd_used, 4),
            }
            for r in replay_trials(sd)
        ],
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
    print(f"  pooled_sd={data['pooled_sd']}  noise_to_signal={data['noise_to_signal']}x")


if __name__ == "__main__":
    main()
