# Could these celiac trials have detected their drugs?

Two celiac trials ran a gluten challenge against the same endpoint. In both, the placebo
arm lost **exactly 0.61** of villous height : crypt depth ratio — the same injury, the
same disease, the same measurement.

**ZED1227** needed its drug to prevent **54%** of that injury to reach 80% power. Its
drug prevented 72–79%. It worked, and remains the first and so far only celiac drug to
hit a histologic primary endpoint.

**KAN-101's Phase 2a** needed its drug to prevent **86%**. Handed a drug exactly as good
as ZED1227's, it would still have reported a null.

The difference was not the biology. It was 25 patients per arm instead of 34, and a
protocol that assumed the endpoint was quieter than it is.

```bash
uv run python -m ctsim.simulate     # full report
uv run pytest                       # 51 tests
```

Live version: **[ahmetdedeler.com/celiac](https://ahmetdedeler.com/celiac)**

---

## 1. The trials that guessed the noise wrong are the ones that failed

Sample size is a function of an assumed standard deviation, written into the protocol
before anyone sees data. Those assumptions are public. Almost nobody checks them against
what the trial then measured.

| Trial | Target effect | **Assumed SD** | **Observed SD** | Claimed power | Actual power |
|---|---|---|---|---|---|
| **ZED1227** (CEC-3) | 0.60 | **0.80** | **0.481** | 80% | ~100% |
| **KAN-101 SynCeD** | 0.50 | **0.50** | **0.662** | 94% | **78%** |

ZED1227 assumed 66% more noise than it found and cleared its endpoint with room to
spare. KAN-101 assumed 25% less noise than it found, and the 94% power in its statistical
analysis plan was really 78% — or 67% against the pooled estimate below.

Both assumptions are quoted verbatim with URLs in
[`src/ctsim/published.py`](src/ctsim/published.py). KAN-101's is from the SAP posted on
ClinicalTrials.gov; ZED1227's appears identically in the NEJM paper and the protocol.

Neither sponsor cited a source for its assumed SD. Neither number appears to have come
from data.

## 2. What the noise actually is

Pooled between-patient SD of ΔVH:CD, from posted per-arm results:

**0.740, 95% CI 0.670–0.826** (6 arms, 3 trials, df = 177)

The README this replaces called that a small pool and moved on. It survives being
attacked properly:

| Check | Result |
|---|---|
| Bartlett test of variance homogeneity | T = 3.87, df = 5, **p = 0.57** — pooling is justified |
| Leave-one-trial-out | 0.683 – 0.767 |
| Placebo/control arms only | 0.732 |
| Posted SDs only, no SE conversion | 0.683 |
| **ZED1227, entirely independent, not in the pool** | **0.481** (ANCOVA-residual scale) |

Every subset lands between 0.68 and 0.77. ZED1227 sits lower, but on a different scale —
see the provenance section.

### Reader error is the small term, and probably smaller than previously stated

| Component | SD | Share of variance |
|---|---|---|
| Total between-patient (ΔVH:CD) | 0.740 | 100% |
| Reader (intra-observer, single central reader) | 0.159 | **4.6%** |
| Residual: biology + biopsy site + orientation | 0.723 | **95.4%** |

An earlier version of this repo put the reader share at 9.2%. That was a double count.
Taavela's "error margins" are Bland–Altman repeatability coefficients — twice the SD of
the *paired difference* between two reads, not twice the SD of one read. The paper's own
limits of agreement settle it: reported as −0.516 to +0.375, a span of 0.891, and
2 × 1.96 × 0.227 = 0.890. A change-from-baseline score is itself a difference of two
reads, so it inherits that SD directly, with no further √2.

Reader error is not a single number, though, and quoting only the friendliest estimate
would be cheating. The same group, same SOP, reported considerably worse reproducibility
eight years later:

| Estimate | Reader SD | Share |
|---|---|---|
| Taavela 2013, intraobserver | 0.159 | 4.6% |
| Taavela 2021, intraobserver APOA4 | 0.194 | 6.9% |
| Taavela 2021, intraobserver H&E | 0.264 | 12.7% |
| Taavela 2021, **inter**observer H&E | 0.508 | 47.2% |

The conclusion — that reader disagreement is the minority term — holds across every
*intra*observer estimate, which is the right comparator for a trial with one blinded
central reader. It fails if that assumption fails. That is the number to attack.

Also note what the reader-error study did *not* capture: it re-read **the same paraffin
blocks**. Zero biopsy-site variability, zero between-endoscopy variability. It is a floor
on reading error, not on measuring a change across two endoscopies.

## 3. What was available to detect

A gluten-challenge trial can only detect what the challenge causes. If the control arm
loses 0.61, a drug preventing fraction *f* of that produces a between-arm difference of
exactly *f* × 0.61, and *f* cannot exceed 1. That caps the detectable effect before
sample size enters.

| Trial | n/arm | Control injury | MDE (80%) | Protection needed |
|---|---|---|---|---|
| ZED1227 | 34 | −0.61 | 0.327 | **54%** — and it delivered 72–79% |
| KAN-101 SynCeD | 25 | −0.61 | 0.525 | **86%** |
| KAN-101, judged on the pooled SD | 25 | −0.61 | 0.586 | **96%** |

This is the standard regulatory notion of assay sensitivity, applied to an endpoint where
the control-arm trajectory is public.

**The challenge protocol is the dominant design lever**, because required N scales with
1/injury². The dataset contains direct measurements of what different protocols do:

| Protocol | Control-arm injury | N/arm to detect 50% protection |
|---|---|---|
| 10 g/day gluten, 14 days (NCT03409796) | **−1.53** | **15** |
| KAN-101 SynCeD 2-week challenge | −0.61 | 93 |
| 3 g/day gluten, 14 days (NCT03409796) | −0.06 | no usable signal at any N |

A 3 g/day challenge for two weeks does not measurably flatten villi. A 10 g/day challenge
does, and makes a 15-patient arm sufficient where a weak challenge needs 93.

### Two levers that do not save you

**ANCOVA instead of a change score.** The baseline-to-follow-up correlation is recoverable
from trials that post both a baseline SD and a change SD: ρ = 0.44–0.73. ANCOVA cuts the
required N by 13–28%. It is free — a line in the analysis plan — and it is not enough.

**Switching to IEL density.** Required N depends on the standardized effect, so a noisier
endpoint with a proportionally larger meaningful change costs nothing. In the challenge
population: VH:CD needs 46/arm, IEL density needs 67/arm. IEL is *worse*. The noise is in
the biology and the biopsy, not in the choice of what to measure on the slide.

## 4. What "0.40" actually is

Most of the field, and the first version of this repo, treats 0.40 as the clinically
meaningful change in VH:CD. It is not, and the distinction matters because it is used as
a power target.

Taavela derived it from their own reader error:

> "twice the standard deviation was only 0.318 in the intraobserver Bland-Altman
> analysis. A cautious new cut-off value of 0.4 could thus be assigned to represent a
> clinically relevant difference **between measurements**."

It is a reading-reproducibility floor for one patient's paired biopsies, rounded up. The
value it replaced (0.5) was itself a lab convention. The Tampere consensus (*Gut* 2018)
adopted >0.4 citing only this paper and graded it **D**. There is no anchor-based MCID
study for VH:CD.

Three corrections follow, all of which apply to this repo's own earlier claims:

1. **FDA has never named VH:CD.** The 2022 draft guidance on celiac drug development
   asks for "histology using a clinically accepted scale (e.g., Marsh-Oberhuber
   classification)" and specifies no effect size. The previous claim that VH:CD is "the
   endpoint regulators want" was wrong. EMA has no celiac guideline at all.
2. **"Noise is 1.85× the effect size" oversold it.** A between-patient SD and a
   within-patient reading threshold are different quantities. SD/δ is still exactly what
   sets sample size, so the arithmetic downstream stands — but the phrase implied the
   endpoint is unusable, which does not follow.
3. **Trials do not target 0.40.** ZED1227 powered on 0.6, KAN-101 on 0.50. Judging a
   trial against a threshold it never adopted is the same error, pointed the other way —
   so each trial is now also scored against its own registered target.

## 5. Data and provenance

- **20 celiac trials** pulled from the ClinicalTrials.gov API v2, 10 with posted results.
  → [`data/curated/histology_endpoints.csv`](data/curated/histology_endpoints.csv)
- **ZED1227** is hand-entered from the NEJM paper — it is registered only in EudraCT
  (2017-002241-30) and has no NCT number, so no API reaches it. Every value carries a
  verbatim quote and URL in [`src/ctsim/published.py`](src/ctsim/published.py).
- **Literature constants** with citations in [`src/ctsim/model.py`](src/ctsim/model.py).

Three provenance rules the numbers depend on:

**Dispersion conversion.** Posted dispersions come as SD, SE or 95% CI. Each is converted
to a between-patient SD and tagged `sd_source`. Rows that cannot be converted keep
`sd = None` and are excluded rather than guessed at.

**Variance scale.** A raw change score has variance 2σ²(1−ρ); an ANCOVA least-squares mean
has σ²(1−ρ²). These are not the same quantity, and the dataset now tags which is which
(`param_type` → `sd_scale`). Pooling by scale: change-score arms give 0.683 (df 60),
ANCOVA-residual arms give 0.768 (df 117). They cannot be reconciled through ρ alone
because scale is confounded with population here — the ANCOVA arms are also the only
restoration-design arms. Both are reported; neither is silently averaged into the other.

**ZED1227 is one row, not four.** All four of its arms imply nearly the same SD
(0.475–0.487) because they come from one model with a common residual variance. Entering
them as four arms would quadruple the apparent degrees of freedom. It is stored at the
model's real df (137).

## 6. Honest limitations

1. **Small pool.** Six arms, three trials for the headline SD, plus ZED1227 independently.
   It now survives a homogeneity test and leave-one-trial-out, which is more than the
   previous version could say — but it is not thirty trials.
2. **Reader share is a range, not a number.** 4.6% under Taavela 2013, 12.7% under the
   2021 H&E figure, 47% if the single-central-reader assumption fails. The qualitative
   claim survives all the intraobserver estimates.
3. **Sampling vs biological variance cannot be separated from public data.** The
   simulator exposes this as an explicit `sampling_share` knob rather than pretending to
   know it. The experiment that would settle it is a multi-fragment re-read study.
4. **The protection-ceiling analysis applies only to prevention designs.** Trials
   enrolling patients with active atrophy (TAK-062, PRV-015) expect healing, so their
   ceiling is the headroom to a normal mucosa, which needs a baseline VH:CD none of them
   post. Those trials are excluded from that table rather than scored wrongly.
   In particular the earlier claim that TAK-062's failure "looks real" is withdrawn: its
   VH:CD endpoint was **secondary**, and it powered itself on a symptom score, not on
   histology at all.
5. **SE→SD conversion assumes posted LS-mean standard errors reflect residual
   between-patient variance.** For MMRM/ANCOVA outcomes these are covariate-adjusted.
   Flagged per row.
6. **None of this says whether any of these drugs work.** It says what these trials were
   capable of detecting. Underpowered is not the same as ineffective — which is exactly
   why KAN-101's null result should not be read as evidence about KAN-101.

## Why this exists

A friend of the author has had celiac disease for 20 years. The dossier in
[`dossier/`](dossier/) is the background research: the 2026 therapeutic landscape, why
tolerance induction keeps failing, and where an outsider can contribute. This came out of
that survey as the one open problem needing no lab, no slides, and no data access
committee — just the observation that nobody had checked whether the field's trials were
arithmetically capable of finding what they were looking for.

## Contributing

1. **More trials.** Latiglutenase/IMGX003, TAK-101 and the ALV003 gluten-challenge study
   all report VH:CD in papers but not the registry. Add them to
   [`src/ctsim/published.py`](src/ctsim/published.py) with a verbatim quote.
2. **More design assumptions.** The planned-vs-observed table is the most useful thing
   here and currently has two rows. Every protocol PDF on ClinicalTrials.gov with a
   sample-size section is another row.
3. **Challenge the variance decomposition.** If you have re-read data, the reader share
   is the number to attack.
4. **IPD.** Takeda, Sanofi, Pfizer and Regeneron all share individual patient data
   through [Vivli](https://vivli.org). Nobody has requested the celiac trials.

Issues and PRs welcome. If you work on celiac clinically or run trials and think this is
wrong, please open an issue — being wrong in public and corrected quickly is the point.

## License

MIT for code. Data is derived from public ClinicalTrials.gov records and published papers.
