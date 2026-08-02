# Software Gap Analysis — What Already Exists, and What Doesn't

Survey run 2026-08-02, before committing compute. Half the value here is the **ruled-out** list.

---

## Ruled out — already built, often several times over

| Idea | What already exists | Verdict |
|---|---|---|
| Gluten epitope / immunopeptidome database | **ProPepper** (37,914 in-silico-digested prolamin peptides, 833 epitopes); **Frontiers 2022** CD-peptide DB (1,041 causative peptides); published **TG2-deamidation + multi-enzyme digestion + NetMHCIIpan 4.3** framework; Tye-Din quantitative epitope mapping | Dead. Rebuilding = arriving where the field is |
| Celiac single-cell atlas | **Cell Reports 2025**: 203,555 cells, 21 active CeD + 11 control duodenal samples. **Nat Immunol 2025**: immune–epithelial–stromal ecosystem. **Nature 2024**: integrated GI atlas, 1.6M cells, 271 donors, 137 cell types | Dead. Extremely well covered |
| AI histology classifier for celiac | Wei et al. 2019 (F1 93.5%, AUC 0.993); **CeliacNet** (AUC >0.96 all Marsh classes); UVA celiac-vs-EE work; medRxiv 2023 + Sci Rep 2024 cell/tissue classification; **medRxiv Mar 2025** interpretable segmentation of crypts/villi/IELs/enterocytes *estimating villus-to-crypt ratio*; **MeasureNet** (arXiv Dec 2024). Commercial: **PathAI** AIM-CD3 deployed in trials | Crowded *and* blocked — every dataset (UVA 162 slides, Radboudumc 1,031 slides, a 1,364-patient set) is private. No public celiac WSI benchmark exists. Can't be fixed from a laptop |
| Patient registry | **iCureCeliac** — 17,000+ participants, 90+ countries, largest global registry. **Go Beyond Celiac** + MyGreenSpace app | Exists |
| Trial matching | **iQualifyCeliac** eligibility screening | Exists |
| Consumer GF apps | Saturated | Irrelevant to a cure |

---

## The gap

**Celiac drug trials are run on primary endpoints whose measurement error is comparable to — possibly larger than — the effect size they are powered to detect. Nobody has quantified this.**

### The numbers

From the PLOS One morphometry validation study, using *expert quantitative morphometry* (the good method, not Marsh eyeballing):

**VH:CD ratio**
- Interobserver limits of agreement: **−0.516 to +0.375**; SD 0.227; **error margin (2SD) = 0.454**
- Intraobserver error margin (2SD): 0.318
- **Threshold for a clinically significant change: 0.4**

**IEL density**
- Intraobserver error margin: **34.2%** (CD3⁺ paraffin), 33.0% (CD3⁺ frozen), **53.2% (H&E)**
- **Threshold for a clinically significant change: 30%**

So on the two endpoints FDA cares about, reader error alone sits at or above the effect size. On H&E-based IEL counting, the error is nearly **2×** the meaningful change.

### And that's only one of three variance components

**Reader variance** is the one people talk about. Two more are larger and almost never modelled:

1. **Biological sampling variance.** Villous atrophy is patchy: **46.6% of patients have different lesions at different duodenal sites**, and **16.9% have variable lesions within a single biopsy**. Severity increases proximal→distal, so which fragment you grab changes the number.
2. **Orientation/sectioning variance.** VH:CD requires well-oriented sections through the villus–crypt axis. Poorly oriented fragments are unusable or biased — this is why there's a published paper using an **ApoA4 immunostain just to define where the villus ends and the crypt begins.**

### Meanwhile, how trials are actually powered

The AMG 714 trial: **~63 subjects, 21/arm, 80% power**, based on a single assumed **SD of 36** for baseline-to-week-12 %-change in VH:CD. One lumped SD. No decomposition into reader / sampling / orientation / biological components.

Compare recent trial sizes: TEV-CeD2 n≈48. TAK-101 Phase 2 n=102. FB102 Phase 1b n=32 (24 active / 8 placebo). Nexvax2, latiglutenase, larazotide, TAK-062 — all "failed."

### The commercial vendor doesn't close this either

**PathAI** — the leading digital pathology vendor in this space — deploys **AIM-CD3 for IEL quantification but has no AI model for Vh:Cd**, the actual primary endpoint. For VH:CD they rely on human readers plus a "reader harmonization program," and publish **no quantitative data on variability reduction** — only a *simulation* of the effect of harmonization, with no baseline error rates. Both products are Research Use Only.

So the field's best commercial tooling does not quantify the noise on its own primary endpoint.

### Confirmed absence

Individual-patient-data placebo-response meta-analyses exist for **ulcerative colitis, functional dyspepsia, chronic constipation, and osteoarthritis**. Simulation-based power analysis for noisy endpoints has been published for **Alzheimer's**. For celiac — the disease where placebo/Hawthorne response is the leading suspect in drug failure — **neither exists.**

---

## The project: an open, calibrated celiac trial simulator

Answer the question the field cannot currently answer: **given the real measurement noise, was any of these trials ever capable of detecting the effect it was looking for?**

### Inputs — all public, all free

- Published reproducibility studies for VH:CD, IEL, Marsh (kappa **0.35**), VCIEL
- Patchiness/sampling-variability studies
- **Full protocols and Statistical Analysis Plans are posted on clinicaltrials.gov as PDFs** — AMG 714 (NCT02637141), PRV-015 (NCT04424927), TAK-062 (NCT05353985), TEV-53408 (NCT06807463), NCT03409796, and more. These contain the exact powering assumptions
- Published results for every completed trial: effect sizes, SDs, dropout, placebo arms

### Build

1. **Variance decomposition.** Pool published estimates into a hierarchical model separating reader / section / biopsy-site / patient / true-treatment variance. Nobody has done this for celiac; it's the core scientific contribution.
2. **Monte Carlo trial simulator.** Parameters: true effect size, N, duration, endpoint (VH:CD / IEL / VCIEL / PRO / blood IL-2), reader design (single central vs multi-reader), biopsies per timepoint, dropout, and a **Hawthorne/adherence-drift term** for the placebo arm.
3. **Retrospective replay.** Re-run every failed celiac trial in the simulator. Output per trial: *the probability this design would have detected an effect of the size the drug plausibly had.* If that number is 0.3, the trial didn't fail — the design did.
4. **Prospective design tool.** Given a target effect, return the required N per endpoint and per design. Show the tradeoff explicitly: how many patients does VH:CD cost you versus a blood-based readout with better reproducibility.

### Validation

- **Sanity check:** simulate trials with the published assumptions and confirm the simulator reproduces the observed placebo-arm variance in trials where it's reported. If it doesn't, the variance model is wrong.
- **Prediction check:** fit on pre-2022 trials, predict outcome distributions of post-2022 trials (TAK-101 Ph2, FB102 Ph1b), compare.
- Publish the failure modes, not just the headline.

### Honest caveats

- **Single central blinded readers remove *between*-reader variance.** So the crude "0.454 > 0.4" comparison overstates the case for trials that use one central reader. Correcting this properly *is* the project, not an objection to it — sampling and orientation variance remain regardless, and are plausibly larger.
- Published reproducibility estimates are few and heterogeneous. The variance model will have wide credible intervals. Report them honestly; wide intervals that still cross the effect size is itself the finding.
- This is a *critique-plus-tool*, not a drug. Its value is in redirecting the next trial and flagging shelved assets that may have been mis-killed.

### Why it's worth doing

- **Zero gated inputs.** No lab, no institution, no data access committee, no slides, no EGA.
- **Runs on the laptop.** Monte Carlo on an M5 Max. No Azure needed at all.
- **It's the missing piece under everything else in this repo** — it quantifies Wall 4 in [02-why-no-cure.md](02-why-no-cure.md), it prices the value of a better endpoint (which is the argument for [05-project-tcr-atlas.md](05-project-tcr-atlas.md) and for blood IL-2), and it produces a number that sponsors, regulators, and patient foundations all care about.
- **Commercially interesting:** if a shelved drug was underpowered rather than ineffective, that's an asset sitting on a shelf with human safety data. TAK-101 and latiglutenase are the first two to check.

---

## Sources

- [Validation of Morphometric Analyses of Small-Intestinal Biopsy Readouts — PLOS One](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0076163)
- [Comparison of the Interobserver Reproducibility With Different Histologic Criteria Used in Celiac Disease — CGH](https://www.cghjournal.org/article/S1542-3565(07)00327-8/fulltext)
- [How patchy is patchy villous atrophy? — Am J Gastroenterol](https://pubmed.ncbi.nlm.nih.gov/20372112/)
- [A Composite Morphometric Duodenal Biopsy Mucosal Scale (VCIEL) — CGH](https://www.sciencedirect.com/science/article/pii/S1542356523009163)
- [Apolipoprotein A4 Defines the Villus-Crypt Border — Front Immunol](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.713854/full)
- [FDA Draft Guidance: Celiac Disease — Developing Drugs (2022)](https://www.fda.gov/media/157682/download)
- [AMG 714 Statistical Analysis Plan (NCT02637141)](https://cdn.clinicaltrials.gov/large-docs/41/NCT02637141/SAP_001.pdf)
- [PathAI end-to-end workflow for celiac disease clinical trials](https://www.pathai.com/blog/precision-at-every-step-pathais-end-to-end-workflow-for-celiac-disease-clinical-trials)
- [ProPepper database](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4597978/)
- [Immune–epithelial–stromal networks in celiac small intestine — Nature Immunology 2025](https://www.nature.com/articles/s41590-025-02146-2)
- [Single-cell integration reveals metaplasia in inflammatory gut diseases — Nature 2024](https://www.nature.com/articles/s41586-024-07571-1)
- [Interpretable Machine Learning based Detection of Coeliac Disease — medRxiv 2025](https://www.medrxiv.org/content/10.1101/2025.03.11.25323763v1.full)
- [iCureCeliac registry](https://celiac.org/icureceliac/)
