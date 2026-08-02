"""Build the curated histology-endpoint dataset from cached ClinicalTrials.gov records.

Output: `data/curated/histology_endpoints.csv` — one row per (trial, outcome, arm),
with dispersion normalised to a **between-patient standard deviation** wherever the
posted dispersion allows it.

Provenance rules (important — this dataset is the input to every downstream claim):
  * `sd_source` records how the SD was obtained: as-posted, derived from SE, or
    derived from a 95% CI. Nothing is imputed.
  * Rows where dispersion cannot be converted keep `sd = None` and are excluded from
    pooling rather than guessed at.
  * `scale` distinguishes raw VH:CD ratio units from percent-change units. These are
    NOT pooled together.
"""

from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path

from ctsim.fetch import extract

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "ctgov_studies.json"
OUT = ROOT / "data" / "curated"

VHCD_PAT = re.compile(r"villous height|villus height|vh:cd|vh/cd", re.IGNORECASE)
IEL_PAT = re.compile(r"intraepithelial lymphocyte", re.IGNORECASE)

# z for a two-sided 95% CI
Z95 = 1.959963985


@dataclass
class Row:
    nct_id: str
    mechanism: str
    design: str  # prevention | restoration | unknown
    trial_n: int | None
    endpoint: str  # VHCD | IEL
    scale: str  # ratio | percent_change | cells_per_100
    outcome_type: str  # PRIMARY | SECONDARY
    param_type: str  # MEAN | LEAST_SQUARES_MEAN | ... — decides the SD's scale
    outcome_title: str
    time_frame: str
    arm_label: str
    is_placebo: bool
    n_arm: int | None
    value: float | None
    sd: float | None
    sd_source: str  # posted_sd | from_se | from_ci95 | unavailable
    posted_dispersion: str
    unit: str


def classify_endpoint(title: str) -> str | None:
    if VHCD_PAT.search(title):
        return "VHCD"
    if IEL_PAT.search(title):
        return "IEL"
    return None


def classify_scale(unit: str, title: str) -> str:
    u = (unit or "").lower()
    if "percent" in u or "percent change" in title.lower():
        return "percent_change"
    if "cell" in u or "iel per" in u:
        return "cells_per_100"
    if "ratio" in u or "unitless" in u:
        return "ratio"
    return u or "unknown"


def is_placebo(label: str) -> bool:
    return "placebo" in (label or "").lower()


def to_sd(value_spread: float | None, dispersion: str, n: int | None,
          lower: float | None, upper: float | None) -> tuple[float | None, str]:
    """Normalise a posted dispersion to a between-patient SD.

    SE -> SD requires the analysed N for that arm: SD = SE * sqrt(n).
    A 95% CI on a mean likewise implies SE = (upper-lower)/(2*z), then SD = SE*sqrt(n).

    Note: for LEAST_SQUARES_MEAN outcomes the posted SE comes from a model (often
    MMRM/ANCOVA) and includes covariate adjustment, so SD recovered this way is an
    approximation of the residual between-patient SD. Flagged via sd_source.
    """
    d = (dispersion or "").lower()
    if "standard deviation" in d and value_spread is not None:
        return value_spread, "posted_sd"
    if "standard error" in d and value_spread is not None and n:
        return value_spread * math.sqrt(n), "from_se"
    if "confidence interval" in d and lower is not None and upper is not None and n:
        se = (upper - lower) / (2 * Z95)
        return se * math.sqrt(n), "from_ci95"
    return None, "unavailable"


def build() -> list[Row]:
    studies = json.loads(RAW.read_text())
    records = [extract(s) for s in studies]
    rows: list[Row] = []

    for rec in records:
        for v in rec.outcome_values:
            ep = classify_endpoint(v.outcome_title)
            if ep is None or v.value is None:
                continue
            # Skip pure count/percentage-of-participants responder outcomes; those are
            # binary endpoints and need a different treatment than continuous SDs.
            if v.param_type in {"COUNT_OF_PARTICIPANTS", "NUMBER"} and "percent change" not in v.unit.lower():
                continue
            sd, src = to_sd(v.spread, v.dispersion_type, v.n_analyzed, v.lower, v.upper)
            rows.append(
                Row(
                    nct_id=rec.nct_id,
                    mechanism=rec.mechanism,
                    design=rec.design,
                    trial_n=rec.enrollment,
                    endpoint=ep,
                    scale=classify_scale(v.unit, v.outcome_title),
                    outcome_type=v.outcome_type,
                    param_type=v.param_type,
                    outcome_title=v.outcome_title,
                    time_frame=v.time_frame,
                    arm_label=v.arm_label,
                    is_placebo=is_placebo(v.arm_label),
                    n_arm=v.n_analyzed,
                    value=v.value,
                    sd=sd,
                    sd_source=src,
                    posted_dispersion=v.dispersion_type,
                    unit=v.unit,
                )
            )
    return rows


def main() -> None:
    rows = build()
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "histology_endpoints.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(rows[0]).keys()))
        w.writeheader()
        for r in rows:
            w.writerow(asdict(r))

    usable = [r for r in rows if r.sd is not None]
    print(f"wrote {len(rows)} rows -> {path}")
    print(f"  with usable SD: {len(usable)}")
    by_src: dict[str, int] = {}
    for r in rows:
        by_src[r.sd_source] = by_src.get(r.sd_source, 0) + 1
    print(f"  sd_source breakdown: {by_src}")


if __name__ == "__main__":
    main()
