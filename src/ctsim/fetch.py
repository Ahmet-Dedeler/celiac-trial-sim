"""Fetch celiac disease trial records from ClinicalTrials.gov API v2.

Pulls full protocol + results sections for a curated set of celiac drug trials,
caches raw JSON, and extracts the design/powering/outcome fields we need for the
variance model.

Everything here is provenance-tracked: every extracted number keeps its NCT ID and
the outcome-measure title it came from. Nothing is imputed at this stage.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import httpx

API = "https://clinicaltrials.gov/api/v2/studies"
RAW_DIR = Path(__file__).resolve().parents[2] / "data" / "raw"

# Curated set of celiac disease therapeutic trials.
# Grouped by mechanism so the provenance of the cohort is auditable.
CURATED: dict[str, list[str]] = {
    "glutenase": [
        "NCT01255696",  # ALV003 / latiglutenase, gluten challenge
        "NCT01917630",  # CeliAction, latiglutenase, n=500
        "NCT03701555",  # PvP001/002/003 healthy + CeD
        "NCT05353985",  # TAK-062 zamaglutenase Ph2  [RESULTS]
    ],
    "tight_junction": [
        "NCT00620451",  # larazotide, active CeD
        "NCT01396213",  # larazotide Ph2b gluten challenge
        "NCT03569007",  # larazotide Ph3 (discontinued)
    ],
    "tg2_inhibitor": [
        "NCT03766445",  # ZED1227 / TAK-227 proof of concept (NEJM 2021)
    ],
    "il15_axis": [
        "NCT02637141",  # AMG 714 Ph2 gluten challenge  [RESULTS]
        "NCT02633020",  # AMG 714 RCD-II  [RESULTS]
        "NCT04424927",  # PRV-015 Ph2b non-responsive  [RESULTS]
        "NCT06807463",  # TEV-53408 anti-IL-15 Ph2a
    ],
    "tolerance": [
        "NCT03644069",  # Nexvax2 Ph2
        "NCT04248855",  # KAN-101 ACeD Ph1
        "NCT05574010",  # KAN-101 ACeD-it Ph1b/2
        "NCT06001177",  # KAN-101 SynCeD Ph2a
    ],
    "other_immune": [
        "NCT06557772",  # amlitelimab anti-OX40L Ph2a/b
    ],
    "gluten_challenge_methodology": [
        "NCT03409796",  # 2-dose gluten challenge (IL-2 biomarker work)
        "NCT04614571",  # TCR-seq + transcriptional profiling under challenge
    ],
}

ALL_NCT: list[str] = [n for group in CURATED.values() for n in group]

MECHANISM_OF: dict[str, str] = {
    n: mech for mech, ids in CURATED.items() for n in ids
}


def fetch_studies(nct_ids: Iterable[str], batch: int = 20) -> list[dict[str, Any]]:
    """Fetch full study records (protocol + results) for the given NCT IDs."""
    ids = list(nct_ids)
    out: list[dict[str, Any]] = []
    with httpx.Client(timeout=60, headers={"accept": "application/json"}) as client:
        for i in range(0, len(ids), batch):
            chunk = ids[i : i + batch]
            r = client.get(
                API,
                params={
                    "filter.ids": ",".join(chunk),
                    "pageSize": batch,
                    # no `fields` param => full record including resultsSection
                },
            )
            r.raise_for_status()
            out.extend(r.json().get("studies", []))
            time.sleep(0.34)  # be polite to the API
    return out


def cache_raw(studies: list[dict[str, Any]]) -> Path:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    path = RAW_DIR / "ctgov_studies.json"
    path.write_text(json.dumps(studies, indent=2))
    return path


# --------------------------------------------------------------------------
# Extraction
# --------------------------------------------------------------------------


@dataclass
class OutcomeValue:
    """One reported value for one arm of one outcome measure."""

    nct_id: str
    outcome_title: str
    outcome_type: str  # PRIMARY / SECONDARY
    param_type: str  # MEAN, LEAST_SQUARES_MEAN, MEDIAN...
    dispersion_type: str  # STANDARD_DEVIATION, STANDARD_ERROR, 95% CI...
    unit: str
    arm_label: str
    n_analyzed: int | None
    value: float | None
    spread: float | None
    lower: float | None
    upper: float | None
    time_frame: str = ""


@dataclass
class TrialRecord:
    nct_id: str
    title: str
    mechanism: str
    sponsor: str
    phase: str
    status: str
    enrollment: int | None
    enrollment_type: str
    start_date: str
    completion_date: str
    has_results: bool
    primary_outcomes: list[str] = field(default_factory=list)
    arm_labels: list[str] = field(default_factory=list)
    outcome_values: list[OutcomeValue] = field(default_factory=list)


def _f(x: Any) -> float | None:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _i(x: Any) -> int | None:
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def extract_outcome_values(study: dict[str, Any], nct_id: str) -> list[OutcomeValue]:
    """Pull every reported per-arm value out of the results section.

    ClinicalTrials.gov nests these as:
      outcomeMeasures[] -> classes[] -> categories[] -> measurements[]
    with group ids resolved against outcomeMeasures[].groups[].
    """
    results = study.get("resultsSection") or {}
    om_module = results.get("outcomeMeasuresModule") or {}
    values: list[OutcomeValue] = []

    for om in om_module.get("outcomeMeasures", []) or []:
        groups = {g.get("id"): g.get("title", "") for g in om.get("groups", []) or []}
        # denoms[].counts[] keys the per-arm analysed N by `groupId` (not `id`).
        counts: dict[str, int | None] = {}
        for denom in om.get("denoms", []) or []:
            if denom.get("units", "").lower().startswith("particip"):
                for c in denom.get("counts", []) or []:
                    counts[c.get("groupId")] = _i(c.get("value"))
        title = om.get("title", "")
        otype = om.get("type", "")
        param = om.get("paramType", "")
        disp = om.get("dispersionType", "")
        unit = om.get("unitOfMeasure", "")
        tf = om.get("timeFrame", "")

        for cls in om.get("classes", []) or []:
            class_title = cls.get("title", "")
            for cat in cls.get("categories", []) or []:
                cat_title = cat.get("title", "")
                for m in cat.get("measurements", []) or []:
                    gid = m.get("groupId")
                    label_bits = [b for b in (class_title, cat_title) if b]
                    suffix = f" [{' / '.join(label_bits)}]" if label_bits else ""
                    values.append(
                        OutcomeValue(
                            nct_id=nct_id,
                            outcome_title=title + suffix,
                            outcome_type=otype,
                            param_type=param,
                            dispersion_type=disp,
                            unit=unit,
                            arm_label=groups.get(gid, gid or ""),
                            n_analyzed=counts.get(gid),
                            value=_f(m.get("value")),
                            spread=_f(m.get("spread")),
                            lower=_f(m.get("lowerLimit")),
                            upper=_f(m.get("upperLimit")),
                            time_frame=tf,
                        )
                    )
    return values


def extract(study: dict[str, Any]) -> TrialRecord:
    p = study["protocolSection"]
    ident = p["identificationModule"]
    design = p.get("designModule", {})
    status_mod = p.get("statusModule", {})
    nct = ident["nctId"]

    enroll = design.get("enrollmentInfo", {}) or {}
    outcomes = p.get("outcomesModule", {}) or {}
    arms = p.get("armsInterventionsModule", {}) or {}

    return TrialRecord(
        nct_id=nct,
        title=ident.get("officialTitle") or ident.get("briefTitle", ""),
        mechanism=MECHANISM_OF.get(nct, "unknown"),
        sponsor=(p.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {}) or {}).get(
            "name", ""
        ),
        phase="/".join(design.get("phases", []) or []),
        status=status_mod.get("overallStatus", ""),
        enrollment=_i(enroll.get("count")),
        enrollment_type=enroll.get("type", ""),
        start_date=(status_mod.get("startDateStruct", {}) or {}).get("date", ""),
        completion_date=(status_mod.get("completionDateStruct", {}) or {}).get("date", ""),
        has_results=bool(study.get("hasResults")),
        primary_outcomes=[o.get("measure", "") for o in outcomes.get("primaryOutcomes", []) or []],
        arm_labels=[a.get("label", "") for a in arms.get("armGroups", []) or []],
        outcome_values=extract_outcome_values(study, nct),
    )


def main() -> None:
    studies = fetch_studies(ALL_NCT)
    path = cache_raw(studies)
    records = [extract(s) for s in studies]
    found = {r.nct_id for r in records}
    missing = [n for n in ALL_NCT if n not in found]

    print(f"cached raw -> {path}")
    print(f"fetched {len(records)}/{len(ALL_NCT)} trials; missing: {missing or 'none'}")
    print(f"trials with posted results: {sum(r.has_results for r in records)}")
    print(f"total per-arm outcome values extracted: {sum(len(r.outcome_values) for r in records)}")


if __name__ == "__main__":
    main()
