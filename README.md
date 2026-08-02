# Could these celiac trials have detected their drugs?

Three gluten-challenge trials, one endpoint, near-identical injury to the placebo arm's
gut (−0.61, −0.61, −0.63 of villous height : crypt depth). Two of the drugs protected
about as well as each other. One is the only celiac drug ever to hit a histologic
primary endpoint. The other was written up as a miss.

| Trial | n/arm | Protection **delivered** | Protection **needed** for 80% power | Reported |
|---|---|---|---|---|
| **ZED1227** | 34 | **79%** | **54%** | p < 0.001 ✅ |
| **TAK-101** Ph2a | 13 | **71%** | **95%** | p = 0.08 ❌ |
| **KAN-101** SynCeD | 25 | −39% | 86% | null ❌ |

TAK-101 delivered essentially what ZED1227 delivered. At ZED1227's sample size it would
have had **92% power** instead of 55%. Its failure is a sample-size artifact, not a
result about the drug.

KAN-101's drug arm went the wrong way, so we can't say it worked — but we also can't say
it didn't. Handed a drug as good as ZED1227's, that trial still reports a null.

```bash
uv run python -m ctsim.simulate     # full report
uv run pytest                       # 66 tests
```

Live version: **[ahmetdedeler.com/celiac](https://ahmetdedeler.com/celiac)**

---

## 1. The noise is not a constant. It grows with the injury.

This repo used to quote one pooled standard deviation for ΔVH:CD and apply it to every
design. That was wrong twice over — once in the arithmetic (§4) and once in the concept.

Plot every published arm's SD against how much its mucosa actually moved:

**SD = 0.400 (±0.024) + 0.299 (±0.073) × |injury|** — 16 arms, 6 trials, r = 0.91

| |injury| | SD observed | SD fitted | n | arm |
|---|---|---|---|---|
| 0.006 | 0.417 | 0.402 | 60 | TAK-062 placebo |
| 0.06 | 0.516 | 0.418 | 7 | 3 g gluten × 14 d |
| 0.15 | 0.413 | 0.445 | 77 | latiglutenase 300 mg |
| 0.27 | 0.401 | 0.481 | 125 | CeliAction placebo |
| 0.35 | 0.616 | 0.504 | 22 | IMGX003 placebo |
| 0.61 | 0.614 | 0.582 | 25 | KAN-101 placebo |
| 0.63 | 0.657 | 0.588 | 15 | TAK-101 placebo |
| 0.85 | 0.707 | 0.654 | 25 | KAN-101 drug |
| 1.53 | 0.941 | 0.857 | 7 | 10 g gluten × 14 d |

Mechanistically this is not surprising. The **floor (0.40)** is what the assay costs you
when nothing happens: biopsy siting, orientation, reading. The **slope** is patient
heterogeneity in *response* — if the average patient loses 1.5 of villous height,
patients differ in how much they lose, and that spread scales with the mean.

A formal homogeneity test rejects a single pooled SD (Bartlett p < 0.001). The injury
model is what replaces it, and there is a test asserting the rejection, so nobody
reintroduces one average later.

### This weakens a design claim I published earlier

A harsher gluten challenge buys more signal — and more noise with it.

| Challenge injury | SD | N/arm for 50% protection | N/arm if SD were constant |
|---|---|---|---|
| 0.20 | 0.460 | 332 | 533 |
| 0.61 | 0.582 | 58 | 58 |
| 1.00 | 0.699 | 31 | 22 |
| **1.53** | 0.857 | **20** | **10** |
| 2.50 | 1.147 | 14 | 4 |

The right-hand column is the earlier version of this README. It claimed a 10 g/day
challenge cuts a trial from 93 patients per arm to 15. The real figure is **58 → 20**.
Still the single biggest lever, but it saturates: N falls roughly as 1/injury, not the
1/injury² that constant-SD arithmetic predicts.

## 2. The trials that guessed the noise wrong are the ones that failed

Sample size is a function of an assumed SD, written into the protocol before anyone sees
data. Those assumptions are public. Almost nobody checks them.

| Trial | Target effect | **Assumed SD** | **Observed SD** | Claimed power | Actual |
|---|---|---|---|---|---|
| **ZED1227** | 0.60 | **0.80** | **0.481** | 80% | ~100% |
| **KAN-101 SynCeD** | 0.50 | **0.50** | **0.662** | 94% | **78%** |
| **IMGX003 CeliacShield** | 0.40 | **0.45** | **0.548** | 86% | **73%** |
| AMG 714 | 40 pp (%-chg) | 36 | n/a — different scale | 89% | — |

ZED1227 assumed 66% more noise than it found and cleared its endpoint with room to
spare. KAN-101 and IMGX003 both assumed less than they measured, and both missed.
AMG 714 powered on %-change in VH:CD (SD = 36) — a different scale from the ratio
units above, so it cannot be scored in the observed column, but the SAP is public.

Three more trials ran VH:CD and never powered on it at all: **TAK-101** sized itself
on IFN-γ SFUs, **TAK-062** on a CDSD symptom score (Cohen's d = 0.55), **PRV-015** on
CeD PRO abdominal symptoms. Their histology results were exploratory by design.

Every row is quoted verbatim with a URL in
[`src/ctsim/published.py`](src/ctsim/published.py).

## 3. What was available to detect

A gluten-challenge trial can only detect what the challenge causes. If the control arm
loses 0.61, a drug preventing fraction *f* produces a between-arm difference of exactly
*f* × 0.61, and *f* cannot exceed 1. That caps the effect before sample size enters —
the standard regulatory notion of assay sensitivity, applied to an endpoint whose
control-arm trajectory is public.

Each trial is judged on **its own** measured SD, since noise is not a constant across
designs.

## 4. Corrections to this repo's own earlier claims

Everything below was wrong in a previous version and is fixed here.

**The headline SD was 0.740. It is ~0.52, and even that shouldn't be quoted alone.**
The pipeline converted *least-squares-mean* standard errors to SDs via SE×√n. That is
invalid — an LS-mean SE carries model terms that cancel in the between-arm contrast. On
TAK-062, which supplied two-thirds of the pooled degrees of freedom, it inflated the SD
from 0.418 to 0.775. The registry's own posted p-value settles it:

| Route | Implied SE of difference | Test statistic | p |
|---|---|---|---|
| SE×√n on arm LS means | 0.141 | 2.35 | 0.019 |
| Posted ANCOVA contrast | 0.077 | 4.32 | **0.000015** |
| *What the registry posts* | | | **< 0.001** ✓ |

Only the contrast reproduces it. LS-mean rows now derive their SD from the posted
contrast or drop out entirely. Where trials post genuine means with genuine SDs the two
routes agree (KAN-101: 0.654 vs 0.662), so this doesn't bend the data — it only bites
where the shortcut was invalid.

**Reader error was 9.2%, then 4.6%, and is really a range.** Taavela's "error margins"
are Bland–Altman repeatability coefficients — the SD of a *paired difference*, not of one
read. A change score is also a paired difference, so the extra √2 was a double count.
The paper's own limits of agreement confirm it (span 0.891 vs 2×1.96×0.227 = 0.890).
Since only the injury term grows, the reader share depends on the challenge: **3–16%**
under Taavela 2013, up to 9–44% under the harsher 2021 H&E figures.

Takeda then measured it directly (UEG Week 2025, abstract MP739):

> "52% of the variability was at the patient level and 23% at the biopsy level. Only 1%
> of the variability was due to reader effects."

That settles the direction and retires the old "sampling can't be separated from biology"
limitation — the answer is 23%. Relative to a single-fragment assay, averaging 4
biopsies at a typical challenge (injury 0.61) cuts required N from **58 → 48** per arm
for 50% protection; 8 biopsies get you to 46. It saturates fast because patient-level
variance (52%) does not shrink with more fragments. Caveat: trials already take multiple
fragments, so the gain vs current practice is smaller than vs a theoretical single
biopsy. Their 1% reader share is after averaging readers, so it is a floor, not a
like-for-like comparison with Taavela's single-reader figures.

**"0.40 is the clinically meaningful change."** It isn't. Taavela derived it from their
own reader error — *"a cautious new cut-off value of 0.4 ... a clinically relevant
difference between measurements"* — a reading-reproducibility floor for one patient,
rounded up. Tampere consensus adopted >0.4 citing only that paper, graded **D**. FDA's
2022 draft guidance names Marsh-Oberhuber and never mentions VH:CD; EMA has no celiac
guideline. Trials don't target it either: ZED1227 powered on 0.6, KAN-101 on 0.50.

**"TAK-062's failure looks real."** Withdrawn. Its VH:CD endpoint was *secondary* and it
powered on a symptom score. It is also a restoration design, where a flat control arm is
not a missing signal, so the protection ceiling does not apply — trials are now
classified prevention vs restoration and scored accordingly.

**Citations.** PMID 24098545 is a tamoxifen meta-analysis; Taavela 2013 is **24146832**.
The 30% IEL threshold originates in Pollock 1992 (*Ann Clin Biochem*), not Taavela.
`NCT03766445` does not exist. ZED1227 is EudraCT **2017-002241-30**, not 2018-002603-14.

## 5. Two levers that still don't save you

**ANCOVA instead of a change score.** Baseline–follow-up correlation recovered from
posted SDs: ρ = 0.22–0.73, cutting required N from 34 to 21–29 per arm. Free, and not
enough. (Arms where the implied ρ goes negative are excluded and flagged: that means the
follow-up spread exceeded baseline — the injury model again — not a real anticorrelation.)

**Switching to IEL density.** Required N depends on the standardized effect, so a noisier
endpoint with a proportionally larger threshold costs nothing. In the challenge
population VH:CD needs 46/arm and IEL needs 56. IEL is worse. Takeda reached the same
conclusion independently: *"no endpoint outperformed Vh:Cd."*

## 6. Data and provenance

- **22 celiac trials** from the ClinicalTrials.gov API v2, 11 with posted results
  → [`data/curated/histology_endpoints.csv`](data/curated/histology_endpoints.csv)
- **Papers, hand-entered** with a verbatim quote and URL each, in
  [`src/ctsim/published.py`](src/ctsim/published.py): ZED1227 (EudraCT-only, no NCT
  number exists), IMGX003/CeliacShield, CeliAction.
- **Literature constants** with citations in [`src/ctsim/model.py`](src/ctsim/model.py).

Rules the numbers depend on:

**Dispersion conversion.** SD / SE / CI each convert differently, and least-squares means
convert only through a posted contrast. Every row is tagged `sd_source`. Rows that cannot
be converted validly keep `sd = None` and drop out — losing data is the correct outcome
when the alternative is inventing it.

**One row per variance estimate.** ZED1227's four arms all imply the same SD
(0.475–0.487) because they come from one model with a common residual variance. Entering
four rows would quadruple the apparent degrees of freedom; it is stored once at the
model's real df (137).

**Recorded absences.** `NOT_AVAILABLE` lists what was searched for and genuinely isn't
published — ALV003-1021's dispersion (nowhere, in paper or three supplements), TAK-062's
raw SD, AMG 714's raw ratio units (both 2019 papers closed-access), Nexvax2's within-arm
change SDs — so the next person doesn't spend an afternoon rediscovering the gap.

## 7. The model survives holdout

Fit without KAN-101 and TAK-101 — the two headline trials — and ask it to predict their
arm SDs from injury alone:

| Held-out arm | \|injury\| | SD observed | SD predicted | \|err\| |
|---|---|---|---|---|
| TAK-101 drug | 0.18 | 0.381 | 0.450 | 0.069 |
| KAN-101 placebo | 0.61 | 0.614 | 0.539 | 0.075 |
| TAK-101 placebo | 0.63 | 0.657 | 0.543 | 0.114 |
| KAN-101 drug | 0.85 | 0.707 | 0.589 | 0.118 |

MAE = 0.094. The slope stays positive (0.207) and the correlation holds (r = 0.90) on
the training set alone. Leave-one-trial-out across all six trials never produces MAE
above 0.30. This is the check that would have caught a coincidence.

## 8. Honest limitations

1. **The injury model is 16 arms from 6 trials.** r = 0.91 and the slope is four SEs from
   zero, but it is a straight line through a modest cloud, and the harshest-challenge
   point (injury 1.53) rests on n = 7.
2. **Scale is confounded with population.** Change-score SDs and ANCOVA-residual SDs are
   different quantities, and the trials supplying each are also the trials with different
   designs. Both are tagged; neither is silently averaged into the other.
3. **The protection ceiling applies only to prevention designs.** Restoration trials'
   ceiling is headroom to a normal mucosa, which needs a baseline VH:CD none of them post.
4. **Reader share is a range, not a number** — 3–16% under one study, up to 44% under
   another, 1% as measured by Takeda with averaged readers.
5. **None of this says whether any of these drugs work.** It says what these trials could
   have detected. Underpowered is not the same as ineffective — which is the entire point
   of the TAK-101 row at the top.

## Why this exists

A friend of the author has had celiac disease for 20 years. The dossier in
[`dossier/`](dossier/) is the background research: the 2026 therapeutic landscape, why
tolerance induction keeps failing, and where an outsider can contribute. This came out of
that survey as the one open problem needing no lab, no slides, and no data access
committee — just the observation that nobody had checked whether the field's trials were
arithmetically capable of finding what they were looking for.

## Contributing

1. **More arms for the injury model.** Any trial reporting a mean change in VH:CD with a
   genuine SD extends the fit. Add to `PAPER_ARMS` with a verbatim quote. AMG 714's raw
   ratio units (both 2019 Lancet Gastro papers are closed-access) would help most.
2. **Challenge the variance decomposition.** Takeda's 52/23/1 split is one abstract with
   no methods detail. If you have re-read data, that is the number to attack.
3. **IPD.** Takeda, Sanofi, Pfizer and Regeneron all share individual patient data
   through [Vivli](https://vivli.org). Nobody has requested the celiac trials.
4. **CeliAction / larazotide sample-size sections.** No SAP is posted on CT.gov for
   those; if you have the protocol PDF, the planned-vs-observed table wants another row.

Issues and PRs welcome. If you work on celiac clinically or run trials and think this is
wrong, please open an issue — being wrong in public and corrected quickly is the point.

## License

MIT for code. Data is derived from public ClinicalTrials.gov records and published papers.
