# Project 1: The Gluten-Specific TCR Atlas + Clone-Burden Score

**Status:** spec, not started. Written 2026-08-02.
**Supersedes** P1 in [03-attack-plan.md](03-attack-plan.md) as the recommended first project.

---

## Why not the epitope atlas (revision)

Prior-art check on the immunopeptidome idea came back crowded:

- **ProPepper** already holds 37,914 in-silico-digested prolamin peptides and 833 epitopes across
  cereal grain protein families.
- A **Frontiers 2022 sequence-searchable database** of CD-associated peptides/proteins exists for
  novel-food risk assessment — 1,041 causative peptides, 76 representative proteins.
- An **integrated in-silico framework with simulated TG2 deamidation + multi-enzyme GI digestion,
  scored with NetMHCIIpan 4.3**, has already been applied to gliadin-derived peptides.
- Tye-Din et al. published comprehensive quantitative T-cell epitope mapping of gluten years ago.

The deamidation-aware enumeration is therefore *done*. What remains unowned there is narrow —
the **minimal epitope cover-set number** and the cross-grain falsification test — and both are
now a ~1-day analysis layered on ProPepper, not a headline project. Keep it as a side quest.

---

## The pick

**Build the instrument that measures whether someone has been cured.**

A cure for celiac means one thing at the cellular level: the gluten-specific CD4+ memory
clonotypes are *gone*. Sollid's group showed those clones persist 20–30 years on a strict
gluten-free diet, with up to 53% clonotype overlap across samples taken decades apart. So:

- Every tolerance therapy that "worked" so far moved IL-2/IFN-γ — cytokine output, i.e.
  **suppression**. None of them showed clone clearance, because nobody measured clone clearance.
- There is currently **no cheap, scalable way to ask "did the clones actually go away?"** from a
  blood draw. It requires HLA-DQ2.5:gluten tetramer staining in a specialist lab.

If a bulk TCRβ repertoire from ordinary blood can be scored for gluten-specific clone burden,
the entire field gets the endpoint it lacks — and it's the endpoint that distinguishes a *cure*
from a *drug you take forever*.

## Why celiac is uniquely suited to this (and most diseases are not)

TCR-specificity prediction usually fails because repertoires are private — your TCRs aren't
mine. Celiac breaks that pattern on every axis:

1. **The HLA is fixed.** ~90% of patients are HLA-DQ2.5 (DQA1\*05:01/DQB1\*02:01). One restriction
   element instead of a population-wide MHC zoo.
2. **The antigen set is small and shared.** Four immunodominant epitopes carry the bulk of the
   response: DQ2.5-glia-α1a, -α2, -ω1, -ω2.
3. **Consequently the repertoire is unusually public** — median **37%** of a given patient's
   gluten-specific clonotypes are public, shared across individuals. That is an extraordinary
   number for antigen-specific TCR work and it is what makes a cross-patient classifier viable.
4. **Strong known motif structure** — the conserved CDR3β R-motif on TRAV26-1/TRBV7-2 for
   DQ2.5-glia-α2 (10.5% of clonotypes), plus 16 further annotated paired CDR3α:CDR3β motifs.

## Ground truth that exists

| Resource | Content | Access |
|---|---|---|
| Dahal-Koirala / Risnes et al., Front Immunol 2021 (Oslo) | **3,122 tetramer-sorted clonotypes, 63 patients**, 4 epitopes; 2,764 clonotypes / 34 patients used for motif discovery; 325 public sequences (145 α, 102 β, 78 paired); 17 annotated motifs | Supplementary Excel **open**; sequences EGA + SRA |
| SRA **SRP102399**, **SRP102402** | Raw repertoire sequencing | **Open** — needs MiXCR-class processing |
| EGA **EGAS00001003245**, **EGAS00001005047** | Deposited sequences | **Controlled access** — DAC application needed |
| Author compilation | 2,918 clonotypes, nt + aa | "on reasonable request" — **email Oslo** |
| Sci Rep 2017 deep TCRβ blood + gut | Gluten-challenge timepoints, celiac vs control | Open / immuneACCESS |
| PLOS One 2021 | scRNA-seq of gluten-specific T cells | GEO |
| VDJdb / McPAS-TCR | Curated celiac entries | Open |
| bioRxiv Dec 2025 | "Shared TCRs in peripheral blood offer robust celiac disease [detection]" | Open preprint — **read first, it may be adjacent to this** |

**The bottleneck is acquisition and normalization, not modeling.** That is precisely the shape
of work a parallel agent fleet is good at.

---

## Plan

### Phase 0 — Acquisition sprint (parallel, ~6h)
One agent per source. Pull supplementary tables, SRA runs, GEO series, immuneACCESS projects,
VDJdb/McPAS celiac slices. Normalize everything to one schema:
`{v_gene, j_gene, cdr3_aa, cdr3_nt, chain, epitope, patient_id, tissue, study, sorting_method, confidence}`.
Use IMGT gene nomenclature throughout; V/J naming inconsistency across studies is the single
biggest silent corrupter in TCR meta-analysis.

In parallel and immediately: **email the Oslo group** for the 2,918-clonotype compilation, and
**start the EGA DAC application** (free, slow — begin day 1 so it lands while other work runs).

### Phase 1 — The atlas (~6h)
Deduplicate, resolve conflicts, annotate with published motifs, ship as a versioned open
dataset + a paper-quality description. **This artifact alone does not currently exist in clean
form and is worth publishing on its own**, independent of whether the classifier works.

### Phase 2 — Classifier (~12h)
- Baselines first: TCRdist3 nearest-neighbour, GIANA/ClusTCR clustering, exact public-clonotype
  matching. Do not skip these — for public repertoires, simple matching is a brutally strong
  baseline and half the published deep models fail to beat it.
- Then embeddings: a TCR language model (TCR-BERT class) or ESM-derived CDR3 representation,
  fine-tuned for gluten-specificity and epitope assignment.
- Output: **P(gluten-specific)** per clonotype, and an aggregate **clone-burden score** per
  repertoire, weighted by clone frequency.

### Phase 3 — Validation, and this is the part that matters (~12h)
The project is worthless without these. Three independent tests, pre-registered before running:

1. **Hold out by patient, never by clone.** Splitting by clonotype leaks public clones across
   train/test and manufactures a fake AUC. This is the standard failure mode in this subfield.
2. **Negative control: HLA-DQ2.5⁺ healthy donors.** They carry naive gluten-reactive precursors
   but no expanded memory clones. Score should be low-but-nonzero. If healthy DQ2.5⁺ people
   score like patients, the classifier is reading HLA background, not disease.
3. **Positive control: gluten challenge dynamics.** Gluten-specific cells peak in blood around
   day 6 of challenge. The score must rise on challenge and fall after. The Sci Rep 2017 dataset
   has these timepoints.

If any of the three fails, publish that. A negative result here is genuinely useful to the field
and costs the same two days.

### Phase 4 — The ask
Take the validated score to the people running trials — argenx (FB102), Teva, Topas, Barinthus,
and Oslo — as a proposed secondary endpoint: *does your drug clear clones, or only silence them?*
No sponsor currently reports this, and the answer determines whether any of these is a cure.

---

## Compute

Almost none. TCR models are small; the atlas is tens of thousands of rows. An M5 Max / 64GB
handles Phases 0–2 locally. Budget **$50–150** of Azure for one A100-day if fine-tuning a TCR
language model, and only then. Do not spend credits to feel productive — the bottleneck here is
data provenance and evaluation discipline, not FLOPs.

## Failure modes to watch

- **Leakage via public clones.** Addressed by patient-level splits. Check twice.
- **The Dec 2025 preprint may already do the diagnostic version.** Read it first. If so, pivot
  hard to the *therapeutic-monitoring* framing — clone clearance as a cure endpoint — which is
  a different and unclaimed question.
- **Adaptive Biotechnologies** has done TCR-based celiac classification work. Same response:
  diagnosis is not the target; longitudinal clearance is.
- **EGA access may not arrive.** Everything above is designed to work on open data alone; EGA is
  upside, not a dependency.
