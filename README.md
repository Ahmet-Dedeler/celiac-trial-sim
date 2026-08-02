# Could these celiac trials have detected their drugs?

## TL;DR

Three trials fed patients gluten and measured how badly it wrecked the lining of their
gut. Two of the drugs worked about as well as each other: ZED1227 prevented 79% of the
damage, TAK-101 prevented 71%. ZED1227 is the only celiac drug that has ever passed this
test. TAK-101 was written up as a failure.

The difference is that TAK-101 only had 13 patients per arm. That small, the drug had to
prevent 95% of the damage for the trial to have a fair shot at showing anything, and no
celiac drug has ever come close to 95%. Give TAK-101 the 34 patients per arm that
ZED1227 had and it passes comfortably. KAN-101 had the same problem at 25 per arm: it
needed 86%, so even a drug as good as ZED1227's would have come back looking like nothing.

These trials keep coming out too small because they assume the measurement noise is a
fixed number. It isn't. The harder the gluten hits, the more patients differ from one
another, so the noise grows right along with the damage
(**SD = 0.400 + 0.299 × injury**, from 16 arms across 6 trials). Get that number wrong
when you write the protocol and you have decided the outcome before anyone enrolls.

None of this shows the drugs work. It shows we don't know. A trial too small to see
anything is not evidence that there was nothing to see, and that is how the field has
been reading these.

I got several things wrong in earlier versions of this. Those are in §4, with the
arithmetic that shows why.

| Trial | n/arm | Protection **delivered** | Protection **needed** for 80% power | Reported |
|---|---|---|---|---|
| **ZED1227** | 34 | **79%** | **54%** | p < 0.001 ✅ |
| **TAK-101** Ph2a | 13 | **71%** | **95%** | p = 0.08 ❌ |
| **KAN-101** SynCeD | 25 | −39% | 86% | null ❌ |

KAN-101's drug arm went the wrong way, so we can't say it worked. We also can't say it
didn't.

```bash
uv run python -m ctsim.simulate     # full report
uv run pytest                       # 81 tests
```

Live version: **[ahmetdedeler.com/celiac](https://ahmetdedeler.com/celiac)**

---

## 1. The noise is not a constant. It grows with the injury.

I used to quote one pooled standard deviation for ΔVH:CD and apply it to every design.
That was wrong twice: once in the arithmetic (§4) and once in the idea itself.

Plot every published arm's SD against how much its mucosa actually moved:

**SD = 0.400 (±0.024) + 0.299 (±0.073) × |injury|** (16 arms, 6 trials, r = 0.91)

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

This shouldn't be surprising. The **floor (0.40)** is what the assay costs you when
nothing happens: biopsy siting, orientation, reading. The **slope** is patient
heterogeneity in *response*. If the average patient loses 1.5 of villous height, patients
differ in how much they lose, and that spread scales with the mean.

A formal homogeneity test rejects a single pooled SD (Bartlett p < 0.001). The injury
model replaces it. There's a test asserting the rejection too, so nobody quietly
reintroduces one average later (including me).

### This weakens a design claim I published earlier

A harsher gluten challenge buys more signal, and more noise with it.

| Challenge injury | SD | N/arm for 50% protection | N/arm if SD were constant |
|---|---|---|---|
| 0.20 | 0.460 | 332 | 533 |
| 0.61 | 0.582 | 58 | 58 |
| 1.00 | 0.699 | 31 | 22 |
| **1.53** | 0.857 | **20** | **10** |
| 2.50 | 1.147 | 14 | 4 |

The right-hand column is what I said before. I claimed a 10 g/day challenge cuts a trial
from 93 patients per arm to 15. The real figure is **58 → 20**. Still the biggest lever
here, but it saturates: N falls roughly as 1/injury, not the 1/injury² that constant-SD
arithmetic predicts.

## 2. The trials that guessed the noise wrong are the ones that failed

Sample size is a function of an assumed SD, written into the protocol before anyone sees
data. Those assumptions are public. Almost nobody checks them.

| Trial | Target effect | **Assumed SD** | **Observed SD** | Claimed power | Actual |
|---|---|---|---|---|---|
| **ZED1227** | 0.60 | **0.80** | **0.481** | 80% | ~100% |
| **KAN-101 SynCeD** | 0.50 | **0.50** | **0.662** | 94% | **78%** |
| **IMGX003 CeliacShield** | 0.40 | **0.45** | **0.548** | 86% | **73%** |
| AMG 714 | 40 pp (%-chg) | 36 | n/a (different scale) | 89% | n/a |

ZED1227 assumed 66% more noise than it found and cleared its endpoint with room to spare.
KAN-101 and IMGX003 both assumed less than they measured, and both missed. The one trial
that was pessimistic about its own noise is the one that worked.

AMG 714 powered on %-change in VH:CD (SD = 36), a different scale from the ratio units
above, so it can't be scored in the observed column. Its SAP is public anyway.

Three more trials ran VH:CD and never powered on it at all: **TAK-101** sized itself on
IFN-γ SFUs, **TAK-062** on a CDSD symptom score (Cohen's d = 0.55), **PRV-015** on CeD PRO
abdominal symptoms. Their histology results were exploratory by design.

Every row is quoted verbatim with a URL in
[`src/ctsim/published.py`](src/ctsim/published.py).

## 3. What was available to detect

A gluten-challenge trial can only detect what the challenge causes. If the control arm
loses 0.61, a drug preventing fraction *f* produces a between-arm difference of exactly
*f* × 0.61, and *f* can't exceed 1. That caps the effect before sample size even enters.
It's the standard regulatory idea of assay sensitivity, applied to an endpoint whose
control-arm trajectory happens to be public.

Each trial is judged on **its own** measured SD, since noise isn't constant across
designs.

## 4. Things I got wrong

Everything below was wrong in a previous version and is fixed here.

**The headline SD was 0.740. It's ~0.52, and even that shouldn't be quoted on its own.**
My pipeline converted *least-squares-mean* standard errors to SDs via SE×√n. That's
invalid, because an LS-mean SE carries model terms that cancel in the between-arm
contrast. On TAK-062, which supplied two-thirds of the pooled degrees of freedom, it
inflated the SD from 0.418 to 0.775. The registry's own posted p-value settles it:

| Route | Implied SE of difference | Test statistic | p |
|---|---|---|---|
| SE×√n on arm LS means | 0.141 | 2.35 | 0.019 |
| Posted ANCOVA contrast | 0.077 | 4.32 | **0.000015** |
| *What the registry posts* | | | **< 0.001** ✓ |

Only the contrast reproduces it. LS-mean rows now derive their SD from the posted
contrast or drop out entirely. Where trials post genuine means with genuine SDs the two
routes agree (KAN-101: 0.654 vs 0.662), so this doesn't bend the data. It only bites
where the shortcut was invalid.

**Reader error was 9.2%, then 4.6%, and is really a range.** Taavela's "error margins"
are Bland–Altman repeatability coefficients: the SD of a *paired difference*, not of one
read. A change score is also a paired difference, so my extra √2 was counting the same
thing twice. The paper's own limits of agreement confirm it (span 0.891 vs
2×1.96×0.227 = 0.890). Since only the injury term grows, the reader share depends on the
challenge: **3–16%** under Taavela 2013, up to 9–44% under the harsher 2021 H&E figures.

Takeda then measured it directly (UEG Week 2025, abstract MP739):

> "52% of the variability was at the patient level and 23% at the biopsy level. Only 1%
> of the variability was due to reader effects."

That settles the direction, and it retires my old "sampling can't be separated from
biology" limitation. The answer is 23%. Against a single-fragment assay, averaging 4
biopsies at a typical challenge (injury 0.61) cuts required N from **58 → 48** per arm for
50% protection, and 8 biopsies get you to 46. It saturates fast because patient-level
variance (52%) doesn't shrink no matter how many fragments you take. Two caveats: trials
already take multiple fragments, so the gain over current practice is smaller than the
gain over a theoretical single biopsy, and Takeda's 1% reader share is measured after
averaging readers, so it's a floor rather than a like-for-like comparison with Taavela's
single-reader numbers.

**"0.40 is the clinically meaningful change."** It isn't. Taavela derived it from their
own reader error (*"a cautious new cut-off value of 0.4 ... a clinically relevant
difference between measurements"*), which is a reading-reproducibility floor for one
patient, rounded up. Tampere consensus adopted >0.4 citing only that paper, graded **D**.
FDA's 2022 draft guidance names Marsh-Oberhuber and never mentions VH:CD, and EMA has no
celiac guideline at all. The trials don't target it either: ZED1227 powered on 0.6,
KAN-101 on 0.50.

**"TAK-062's failure looks real."** Withdrawn. Its VH:CD endpoint was *secondary* and it
powered on a symptom score. It's also a restoration design, where a flat control arm
isn't a missing signal, so the protection ceiling doesn't apply. Trials are now classified
prevention vs restoration and scored accordingly.

**Citations.** PMID 24098545 is a tamoxifen meta-analysis; Taavela 2013 is **24146832**.
The 30% IEL threshold originates in Pollock 1992 (*Ann Clin Biochem*), not Taavela.
`NCT03766445` does not exist. ZED1227 is EudraCT **2017-002241-30**, not 2018-002603-14.

## 4b. Predictions about three trials that have not reported yet

Everything above is retrospective, which is the easy direction. A model that explains
finished trials has had every chance to be fitted to them. So I also want to say, before
the fact, what three running trials can resolve. They read out between **September 2026
and mid-2027**.

```bash
uv run python -m ctsim.prospective
```

None of the load-bearing design detail is on ClinicalTrials.gov. It comes from the
sponsors' own protocols on the EU Clinical Trials Information System, and two of them were
partly redacted and then recovered elsewhere in the same package. Sanofi's gluten dose is
blacked out in the protocol but printed in the patient information sheet. Dr Falk's is
blacked out in protocol v9.0 but readable in the tracked-changes v8.0.

| | Teva TEV-53408 | Sanofi amlitelimab | Dr Falk ZED1227 |
|---|---|---|---|
| Gluten | **3 g/day** | ~250 mg/day | ~214 mg/day |
| Challenge window | 6 wk | 12 wk (of a 24-wk period) | 15 wk |
| Entry | Vh:Cd **≥2.0** (healed) | Vh:Cd **<2.5** (atrophic) | Vh:Cd **≤2.5** (atrophic) |
| n (contrast) | 20 v 20 | 34 v 34 | 72 v 48 |
| VH:CD powered? | **No** (precision, SD 0.65) | Yes, assumptions redacted | **No**, rides on a symptom score |
| MDE (80%) | 0.52 | 0.27 | 0.21 |

**The gluten dose spans 14-fold, and it runs opposite to the damage at entry.** The trial
delivering the largest challenge is the one enrolling healed patients. The two delivering
the smallest are enrolling patients who already have atrophy.

**Teva** is the only one built to manufacture a histologic signal (3 g/day is the dose
behind every positive VH:CD result on record, going into a mucosa with room to fall), and
it's also the only one that openly declines to power for the endpoint. Its exposure is
just arithmetic. At 20 evaluable per arm it needs **85% protection**, and the best ZED1227
arm delivered 80%. Its own ±0.40 precision target spans everything from no effect to
complete protection.

That's the sharpest claim here, so it's worth saying how little reverses it. **The verdict
flips at 23 evaluable per arm.** The protocol targets 40 evaluable (20 per arm) but permits
up to 48 randomised (24 per arm), and the registry already lists 50 enrolled. So the honest
form is conditional: at ~20 evaluable per arm TEV-53408 could not have seen a best-in-class
drug, at 23+ it could. That's the first number to check at readout, before any p-value.

**Sanofi and Dr Falk** need healing of **0.27** and **0.21**. The largest drug-minus-placebo
VH:CD difference any restoration-design celiac trial has ever produced is **+0.14**, in
ZED1227's own CEC-004 in 397 patients (0.26 vs 0.12), and that trial still missed its
primary endpoint. So both need more than anyone has managed: 1.9× and 1.5×.

That's the prediction, and it's falsifiable. It's pinned in
[`tests/test_prospective.py`](tests/test_prospective.py) so it can't drift toward whatever
the results turn out to be.

Every protocol number above is quoted in
[`src/ctsim/prospective.py`](src/ctsim/prospective.py) with the CTIS document id it came
from. CTIS serves these PDFs through a signed-URL handshake instead of a plain link, which
is probably most of the reason nobody reads this corpus:

```bash
curl -s https://euclinicaltrials.eu/ctis-public-api/documents/<ct-number>/<uuid>/download
```

## 5. Levers

**ANCOVA instead of a change score.** Baseline–follow-up correlation recovered from
posted SDs: ρ = 0.22–0.73, which cuts required N from 34 to 21–29 per arm. Free, and not
enough. (Arms where the implied ρ goes negative are excluded and flagged. That means the
follow-up spread exceeded baseline, which is the injury model again, not a real
anticorrelation.)

**Endpoint choice. This one is real, and I got it wrong before.**

I said IEL density is the worse endpoint (VH:CD 46/arm vs IEL 56) and cited Takeda's
*"no endpoint outperformed Vh:Cd"* as agreement. But that comparison asks a narrow
question: how many patients to detect a fixed *threshold* change, where the IEL threshold
is a 30% relative change converted through a single posted baseline. That isn't the
question a trial designer actually has.

Syage et al. asked the better one: the standardized effect each endpoint actually
achieved on the same patients and the same biopsies ([*Clin Gastroenterol Hepatol*
2024;22:1238, PMC12213069](https://pmc.ncbi.nlm.nih.gov/articles/PMC12213069/)):

| Trial | ΔVh:Cd | ΔIEL | **ΔVCIEL** (composite) |
|---|---|---|---|
| ALV003-1021 | 1.37 (p=0.038) | 1.17 (p=0.005) | **1.86 (p=0.004)** |
| IMGX003 CeliacShield | 0.76 (**p=0.057**) | 0.98 (**p=0.018**) | **1.14 (p=0.007)** |

Read the IMGX003 row again. **The same trial missed on VH:CD at p=0.057 and hit on IEL at
p=0.018 and on the composite at p=0.007.** Same patients, same slides, same gluten
challenge. Which way it went was decided by which number came off the microscope. IMGX003
sits in my dataset as a miss, and on the endpoint that best captured its effect it wasn't
one.

So endpoint choice isn't a dead lever. It's probably the cheapest live one, since the
composite costs nothing beyond measurements trials already take, and it beat both of its
own components in both trials.

Two caveats worth keeping. Syage's effect size is the *observed* drug effect over a
baseline SD, so it's partly a function of how well each drug worked rather than a pure
assay property, and a well-powered head-to-head on a fixed target would be much better
evidence. And **the literature contradicts itself here.** Takeda's MP739 reports *"no
endpoint outperformed Vh:Cd"* on TAK-062, directly against Syage on two other trials.
Nobody has reconciled that. It's an open question, not a settled one.

## 6. Data and provenance

- **22 celiac trials** from the ClinicalTrials.gov API v2, 11 with posted results
  → [`data/curated/histology_endpoints.csv`](data/curated/histology_endpoints.csv)
- **Papers, hand-entered** with a verbatim quote and URL each, in
  [`src/ctsim/published.py`](src/ctsim/published.py): ZED1227 (EudraCT-only, no NCT
  number exists), IMGX003/CeliacShield, CeliAction.
- **Literature constants** with citations in [`src/ctsim/model.py`](src/ctsim/model.py).

Rules the numbers depend on:

**Dispersion conversion.** SD / SE / CI each convert differently, and least-squares means
only convert through a posted contrast. Every row is tagged `sd_source`. Rows that can't
be converted validly keep `sd = None` and drop out. Losing data is the right outcome when
the alternative is inventing it.

**One row per variance estimate.** ZED1227's four arms all imply the same SD (0.475–0.487)
because they come from one model with a common residual variance. Entering four rows would
quadruple the apparent degrees of freedom, so it's stored once at the model's real df
(137).

**Recorded absences.** `NOT_AVAILABLE` lists what I looked for and genuinely isn't
published: ALV003-1021's dispersion (nowhere, not in the paper or three supplements),
TAK-062's raw SD, AMG 714's raw ratio units (both 2019 papers closed-access), Nexvax2's
within-arm change SDs. So the next person doesn't spend an afternoon rediscovering the
same gap.

## 7. The model survives holdout

Fit without KAN-101 and TAK-101, the two headline trials, then ask it to predict their
arm SDs from injury alone:

| Held-out arm | \|injury\| | SD observed | SD predicted | \|err\| |
|---|---|---|---|---|
| TAK-101 drug | 0.18 | 0.381 | 0.450 | 0.069 |
| KAN-101 placebo | 0.61 | 0.614 | 0.539 | 0.075 |
| TAK-101 placebo | 0.63 | 0.657 | 0.543 | 0.114 |
| KAN-101 drug | 0.85 | 0.707 | 0.589 | 0.118 |

MAE = 0.094. The slope stays positive (0.207) and the correlation holds (r = 0.90) on the
training set alone. Leave-one-trial-out across all six trials never produces MAE above
0.30. If the relationship were a coincidence, this is where it would have shown up.

## 8. Honest limitations

1. **The injury model is 16 arms from 6 trials.** r = 0.91 and the slope is four SEs from
   zero, but it's still a straight line through a modest cloud, and the harshest-challenge
   point (injury 1.53) rests on n = 7.
2. **Scale is confounded with population.** Change-score SDs and ANCOVA-residual SDs are
   different quantities, and the trials supplying each are also the trials with different
   designs. Both are tagged, and neither is silently averaged into the other.
3. **The protection ceiling only applies to prevention designs.** For restoration trials
   the ceiling is headroom to a normal mucosa, which needs a baseline VH:CD that none of
   them post.
4. **Reader share is a range, not a number**: 3–16% under one study, up to 44% under
   another, 1% as measured by Takeda with averaged readers.
5. **None of this says whether any of these drugs work.** It says what these trials could
   have detected. Underpowered isn't the same as ineffective, which is the whole point of
   the TAK-101 row at the top.

## Why this exists

The dossier in [`dossier/`](dossier/) is the background research: the 2026 therapeutic
landscape, why tolerance induction keeps failing, and where an outsider can actually
contribute. This came out of that survey as the one open problem that needs no lab, no
slides, and no data access committee. Nobody had checked whether these trials were
arithmetically capable of finding what they were looking for.

## Contributing

1. **More arms for the injury model.** Any trial reporting a mean change in VH:CD with a
   genuine SD extends the fit. Add it to `PAPER_ARMS` with a verbatim quote. AMG 714's raw
   ratio units (both 2019 Lancet Gastro papers are closed-access) would help most.
2. **Challenge the variance decomposition.** Takeda's 52/23/1 split is one abstract with
   no methods detail. If you have re-read data, that's the number to attack.
3. **IPD.** Takeda, Sanofi, Pfizer and Regeneron all share individual patient data through
   [Vivli](https://vivli.org). Nobody has requested the celiac trials.
4. **CeliAction / larazotide sample-size sections.** No SAP is posted on CT.gov for those.
   If you have the protocol PDF, the planned-vs-observed table wants another row.

Issues and PRs welcome. If you work on celiac clinically or run trials and think this is
wrong, please open an issue. I'd rather be wrong in public and corrected fast.

## License

MIT for code. Data is derived from public ClinicalTrials.gov records and published papers.
