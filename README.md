# Are celiac disease trials able to detect the drugs they're testing?

**Short answer: often, no.**

The endpoint that regulators want for celiac drug approval — change in the villous
height : crypt depth ratio (VH:CD) — carries a between-patient standard deviation of
**0.74**, measured from posted results of real trials. The change considered
clinically meaningful is **0.40**.

**The noise is 1.85× the effect size.**

At that noise level, detecting a full clinically meaningful improvement with 80% power
needs **~108 patients**. Detecting half of one needs **~430**. Celiac Phase 2 trials
routinely run 30–150 patients total.

This repo contains the data, the model, and a simulator so you can check that yourself.

```bash
uv run python -m ctsim.simulate     # full report
uv run pytest                       # 9 tests, incl. Monte Carlo vs closed form
```

Live version: **[ahmetdedeler.com/celiac](https://ahmetdedeler.com/celiac)**

---

## What we found

### 1. The endpoint noise is real, and it comes from the patients, not the pathologists

| Component | SD | Share of variance |
|---|---|---|
| Total between-patient (ΔVH:CD) | 0.740 | 100% |
| Reader error (single central reader) | 0.225 | **9.2%** |
| Residual: true biology + biopsy site + orientation | 0.705 | **90.8%** |

This matters commercially. The standard vendor pitch for celiac trial pathology is
*reader harmonization* — aligning pathologists before the trial starts. That addresses
roughly **9%** of the problem. The other 91% is patient heterogeneity and where in the
duodenum you put the forceps. Villous atrophy is patchy: 46.6% of patients show
different lesion grades at different duodenal sites and 16.9% vary within a single
biopsy fragment.

### 2. At least one "failed" trial was never able to succeed

Every trial judged against the same pooled noise estimate, so a small noisy trial
can't flatter itself:

| Trial | n/arm | Min. detectable effect (80%) | vs 0.40 | Power @0.40 | Power @0.20 |
|---|---|---|---|---|---|
| Gluten challenge dose-finding (NCT03409796) | 7 | 1.108 | 2.77× | 17.3% | 8.0% |
| **KAN-101 SynCeD** (NCT06001177) | 25 | **0.586** | **1.47×** | **48.1%** | 15.9% |
| TAK-062 (NCT05353985) | 59 | 0.382 | 0.95× | 83.6% | 31.2% |

**KAN-101's Phase 2a had a 48% chance of detecting a fully clinically meaningful
effect — a coin flip.** Anokion's program was shelved. Whether the drug worked is,
on this endpoint, not something that trial could have established.

TAK-062 is the contrast case: adequately powered at 84%, and it still showed the drug
arm doing *worse* than placebo (Δ −0.338 vs −0.006). That failure looks real.

### 3. Sample sizes nobody is running

| Effect to detect | Per arm | Total |
|---|---|---|
| 0.40 — full clinically meaningful change | 54 | 108 |
| 0.30 | 96 | 192 |
| 0.20 — half | 215 | 430 |
| 0.10 — quarter | 860 | 1720 |

---

## Data

All inputs are public. Nothing here is behind a data access committee.

- **Per-arm outcome values** scraped from ClinicalTrials.gov API v2 results sections
  for 19 celiac therapeutic trials (9 with posted results).
  → [`data/curated/histology_endpoints.csv`](data/curated/histology_endpoints.csv)
- **Computed results** → [`data/curated/results.json`](data/curated/results.json)
- **Literature constants** with citations live in
  [`src/ctsim/model.py`](src/ctsim/model.py) — measurement error from Taavela 2013
  (PLoS One), Marsh reproducibility from Corazza 2007 (CGH), patchiness from
  Bonamico 2010 (AJG).

### Provenance rules

Posted dispersions come in three flavours (SD, SE, 95% CI). Each is converted to a
between-patient SD and tagged with `sd_source` (`posted_sd` / `from_se` / `from_ci95`).
Rows that can't be converted are kept with `sd = None` and excluded from pooling
rather than guessed at. Raw-ratio and percent-change scales are never pooled together.

---

## Honest limitations

Read these before citing anything above.

1. **Small pool.** The headline SD rests on 6 arms across 3 trials. It is stable
   across them (0.52–0.94) and across very different designs, but it is not 30 trials.
2. **SE→SD conversion assumes the posted least-squares-mean standard errors reflect
   residual between-patient variance.** For MMRM/ANCOVA outcomes these are covariate-
   adjusted, so the recovered SD is an approximation. Flagged per-row.
3. **The variance decomposition uses a literature intra-observer estimate**, not a
   re-read study on these specific trials. The 9.2% reader share is an estimate with
   real uncertainty. It would take a genuinely large error to change the conclusion
   that reader disagreement is the minority term, but it is an estimate.
4. **Sampling vs biological variance cannot be separated from public data.** The
   simulator exposes this as an explicit `sampling_share` knob rather than pretending
   to know it. If you want the real number, the experiment is a multi-fragment
   re-read study, and it has not been done.
5. **This says nothing about whether any of these drugs work.** It says what these
   trials were capable of detecting. Underpowered is not the same as effective.

---

## Why this exists

A friend of the author has had celiac disease for 20 years. The dossier in
[`dossier/`](dossier/) is the background research: the 2026 therapeutic landscape, why
tolerance induction keeps failing, and where an outsider can contribute. This simulator
came out of that survey as the one genuinely open problem that needed no lab, no slides,
and no data access committee — just the observation that nobody had checked whether the
field's trials were arithmetically capable of finding what they were looking for.

Individual-patient-data placebo-response meta-analyses exist for ulcerative colitis,
functional dyspepsia, chronic constipation and osteoarthritis. Simulation-based power
analysis for noisy endpoints has been published for Alzheimer's. For celiac, neither
existed.

## Contributing

The most useful contributions, roughly in order:

1. **More trials.** Add NCT IDs to `CURATED` in `src/ctsim/fetch.py`. European trials
   (ZED1227/TAK-227 sits in EudraCT, not ClinicalTrials.gov) are the biggest gap.
2. **Published-paper extraction.** Several trials report VH:CD dispersion in the paper
   but not in the registry. Those numbers would roughly double the pool.
3. **Challenge the variance decomposition.** If you have re-read data, the 9.2% figure
   is the number to attack.
4. **IPD.** Takeda, Sanofi, Pfizer and Regeneron all share individual patient data
   through [Vivli](https://vivli.org). Nobody has requested the celiac trials.

Issues and PRs welcome. If you work on celiac clinically or run trials and think this
is wrong, please open an issue — being wrong in public and corrected quickly is the
point.

## License

MIT for code. Data is derived from public ClinicalTrials.gov records.
