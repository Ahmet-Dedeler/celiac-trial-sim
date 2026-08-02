# Why There's No Cure — The Three Walls

Celiac is the best-characterized autoimmune disease we have. We know the antigen, the HLA
restriction, the enzyme that makes the neo-epitope, the T cells, and the effector mechanism.
And we still can't cure it. Here is exactly where it breaks.

---

## Wall 1 — The memory T cells are permanent

Gluten-specific CD4+ T cells in HLA-DQ2.5/DQ8 individuals **do not decay**. Sollid's group
tracked clonotypes in the same patients across 16–28 years on a strict gluten-free diet and
found **up to 53% clonotype overlap** between samples decades apart. The same clones sit in
blood and gut, in a narrow, stable phenotype that has been described as "T cell scarring."

Implications, which the field mostly refuses to internalize:

- A gluten-free diet does not deplete the pathogenic repertoire. It only starves it of antigen.
- Any true cure requires **eradication or genuine re-education** of clones that have survived
  20+ years of antigen deprivation. Peptide dosing was never plausibly going to do that —
  these are exactly the cells most resistant to deletion-by-antigen.
- This is *unlike* type 1 diabetes (target tissue destroyed, so tolerance arrives too late) and
  *unlike* MS (unclear antigen). Celiac is the one where the target tissue fully regenerates.
  **If you cleared the clones, the patient would actually be cured.** That's the prize.

## Wall 2 — The antigen is a repertoire, not a peptide

There is no single "gluten epitope." There are dozens of immunodominant and subdominant
DQ2.5- and DQ8-restricted epitopes spread across α-, γ-, and ω-gliadins, low- and
high-molecular-weight glutenins, plus barley hordeins and rye secalins. The canonical shortlist
(DQ2.5-glia-α1a, -α2, -ω1, -ω2) covers the dominant response — not the whole thing.

- **TG2 deamidation manufactures new epitopes.** Glutamine → glutamate at specific positions
  creates the negative charge that HLA-DQ2.5's P4/P6/P7 pockets want. So the antigen space is
  the *post-translationally modified* proteome, not the genomic one.
- **Epitope spreading refills what you tolerize.** Nexvax2 tolerized three peptides. TAK-101
  used whole gliadin but delivered it systemically. Neither covers the space.
- This is also why CRISPR wheat is at 97.7% and not 100%: knock out the D-genome α-gliadins and
  the A/B genome and the glutenins are still there. **97.7% is not safe. 20 ppm is the legal
  threshold and celiac T cells respond well below it.**

## Wall 3 — The tissue is licensed independently of the T cells

The epithelium in celiac isn't a passive victim. IL-15 (and IFN-λ) put enterocytes and
intraepithelial lymphocytes into a cytotoxic, NKG2D-driven state where IELs kill enterocytes.
This arm runs semi-autonomously from the CD4 response.

Two consequences that explain almost every trial result of the last five years:

1. **You can protect the villi without touching the T cells.** That's exactly what FB102
   (anti-CD122) and the anti-IL-15 antibodies do — best histology data in the field, no
   tolerance induced. Chronic drug, works, not a cure.
2. **You probably cannot induce tolerance in an IL-15-high gut.** IL-15 is the classic signal
   that converts a tolerogenic antigen encounter into a priming one — it renders T cells
   resistant to Treg suppression and to TGF-β. Every tolerance trial dosed patients whose gut
   was inflamed or about to be inflamed by a gluten challenge. **They were trying to teach
   tolerance in the one microenvironment that forbids it.**

---

## The synthesis nobody is testing

> Quiet the IL-15 arm **first**, then tolerize into the permissive window, then withdraw.

This is the type-1-diabetes playbook (anti-CD3 + antigen beats antigen alone) applied to the
one disease where the antigen is actually known and the target tissue actually regenerates.

Every single celiac trial to date has been **monotherapy**. There is no registered combination
trial of an IL-15/CD122 blocker with a tolerogenic agent. The pieces now exist and are in
different companies' hands:

- IL-15 arm: FB102 (argenx), TEV-53408 (Teva), ordesekimab (Sanofi)
- Tolerance arm: TPM502 (Topas), VTP-1000 (Barinthus), TAK-101 (shelved, licensable)
- Antigen-load arm: TAK-227 (Takeda), glutenases

The combination is scientifically the most defensible untested hypothesis in the field, and
it's blocked by **business structure, not biology**. That is an unusually tractable kind of
obstacle.

---

## Wall 4 (the meta-wall) — the field can't measure success

Not biology, but it has killed more drugs here than biology has.

- **Symptom endpoints are poisoned by Hawthorne bias.** Patients enrolled in celiac trials get
  *better* at gluten avoidance than they were at baseline. The placebo arm improves. The drug
  signal drowns. TAK-062 and latiglutenase both plausibly died of this.
- **Histology requires endoscopy**, so you get 2 timepoints, tiny n, and no dose-response
  resolution.
- **Symptoms and histology don't correlate well** in celiac — you can be asymptomatic with
  villous atrophy and vice versa. So the two endpoints regularly disagree, and FDA wants both.
- FDA requires **52-week** trials for chronic dosing with validated PROs *and* histology.
- Children are excluded from nearly all trials, so the first approved drug will be adults-only.

**The fix already exists and is underused:** gluten-stimulated whole-blood IL-2 release. Plasma
IL-2 rises within ~4 hours of gluten ingestion, is the earliest and most sensitive marker of
gluten-specific CD4 activation known, and works in patients on a gluten-free diet. It turns a
scoped-biopsy readout into a blood draw. Barinthus already used it as the Phase 1 pharmacology
readout. **A cheap, decentralized, high-frequency IL-2-based endpoint would unblock the entire
field's trial-design problem.**

---

## Wall 5 — Money

- NIH celiac funding: roughly **$9M/year** (FY2021), up from ~$3M. Crohn's — a quarter of the
  prevalence — has historically drawn several times more.
- ~1% of the population, ~70 million people worldwide, highest mortality ratio among the
  GI diseases in one comparative analysis.
- The commercial counter-argument that suppresses investment: *"the diet works, so the bar for
  a drug is high and payers won't reimburse."* This is why programs get killed for "economic
  reasons" (KAN-101) even with positive data.
- Consequence: **few investigators → few grants → few reviewers with expertise → fewer grants.**
  A self-reinforcing funding trap.

---

## Sources

- [Disease-driving CD4+ T cell clonotypes persist for decades in celiac disease — JCI](https://www.jci.org/articles/view/98819)
- [Gluten-Free Diet Induces Rapid Changes in Phenotype and Survival Properties of Gluten-Specific T Cells — Gastroenterology 2024](https://www.gastrojournal.org/article/S0016-5085(24)00351-2/fulltext)
- [Therapeutic and Diagnostic Implications of T Cell Scarring in Celiac Disease and Beyond](https://www.sciencedirect.com/science/article/pii/S1471491419301261)
- [Update 2020: nomenclature and listing of celiac disease-relevant gluten epitopes recognized by CD4+ T cells](https://link.springer.com/article/10.1007/s00251-019-01141-w)
- [A Mouse Model of Celiac Disease (DQ8-Dd-villin-IL-15tg) — Abadie/Jabri](https://pmc.ncbi.nlm.nih.gov/articles/PMC10203737/)
- [Breakthrough mouse model of celiac disease — UChicago Medicine](https://www.uchicagomedicine.org/forefront/gastrointestinal-articles/2020/february/celiac-mouse-model)
- [Plasma IL-2 and Symptoms Response after Acute Gluten Exposure — PubMed](https://pubmed.ncbi.nlm.nih.gov/34797778/)
- [Whole blood IL-2 release test for rare circulating gluten-specific T cells — PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8119810/)
- [Disparities in NIH Funding Between GI Diseases (CDF)](https://celiac.org/wp-content/uploads/2017/09/NIH-Disparities-in-Disease-Funding.pdf)
