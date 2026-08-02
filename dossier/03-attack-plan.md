# Attack Plan — What a Small Technical Team Can Actually Do

Framing: you cannot run a clinical trial, do wet-lab work, or dose a patient. What you *can*
do is build the missing shared substrate, fund the missing experiment, and force the
combination nobody is testing. Ranked by **leverage × accessibility**.

---

## P1 — Build the Gluten Immunopeptidome Atlas
**Compute-only. Weeks, not years. Genuinely unowned.**

**The gap:** the canonical epitope list (Sollid et al. 2020 nomenclature) is a *curated
literature list* — a few dozen experimentally confirmed epitopes. It is not a systematic map.
Nobody has computed, across all sequenced *Triticum / Hordeum / Secale* cultivar proteomes,
the complete set of TG2-deamidatable 9-mers that bind HLA-DQ2.5 and DQ8. Every downstream
program is flying blind because of it:

- Tolerance therapies pick 3–4 peptides and hope. (Nexvax2 picked 3. It failed.)
- CRISPR-wheat programs claim "97.7% gluten reduction" with no principled celiac-safety score.
- DONQ52 had to be characterized empirically against 25+ pHLA complexes.

**What you build:**
1. Ingest every gliadin/glutenin/hordein/secalin sequence across cultivars (UniProt, wheat
   pangenome projects, IWGSC).
2. Apply the **TG2 deamidation rule** (QXP motifs, position-specific Q→E) combinatorially to
   generate the real post-translational antigen space — this is the step everyone skips.
3. Predict HLA-DQ2.5 / DQ2.2 / DQ8 binding across all resulting 9-mer cores. Class II binding
   prediction is the weak link — use an ensemble (NetMHCIIpan-class + structure-based rescoring
   with AlphaFold3/Boltz-class pMHC modeling) rather than trusting any single predictor.
4. Rescore for **TCR accessibility**, not just binding — a peptide that binds DQ2.5 but presents
   a flat surface isn't pathogenic. This is where modern structure prediction earns its keep and
   where the literature is thinnest.
5. Calibrate the whole thing against the ~40 experimentally confirmed epitopes as a held-out set.

**Deliverables:**
- A public **"celiac safety score" for any grain protein sequence or cultivar** — the thing
  wheat breeders and CRISPR groups currently do not have.
- A ranked **minimal epitope cover set**: how many peptides must a tolerance therapy include to
  cover ≥95% of the pathogenic repertoire? If the answer is 30 rather than 4, that single number
  reframes the entire tolerance field and explains a decade of failures.
- Per-cultivar ranking: which existing wheat landraces are already closest to safe.

**Why this is the right first move:** it's the only item on this list that is fully accessible
to a software team, produces a real scientific artifact in weeks, is publishable, and buys you
credibility with every lab in P4/P5. It also has a genuine chance of being *the* insight —
"the tolerance trials were dosing a fraction of the antigen space" is a testable, high-value,
currently-unproven claim.

**Risk:** class II epitope prediction is genuinely hard and historically unreliable. If your
held-out calibration on the 40 known epitopes is poor, say so publicly and pivot to the
structure-based subset rather than shipping a confident-looking wrong atlas.

---

## P2 — Fix the endpoint problem: a decentralized IL-2 gluten-response platform
**Software + a lab partnership. This is the "build a company" option.**

**The gap:** the field's drugs keep dying on endpoints, not on biology (see
[02-why-no-cure.md](02-why-no-cure.md), Wall 4). Meanwhile a validated blood readout exists —
gluten-stimulated whole-blood IL-2 release, plasma IL-2 peaking ~4h post-exposure — and it is
being used as a *diagnostic* rather than as the trial infrastructure it could be.

**What you build:** the celiac equivalent of a CGM-plus-registry.
- Standardized at-home or draw-center whole-blood gluten-stimulation kit → central lab IL-2.
- Longitudinal cohort: serology (tTG-IgA, DGP), symptom PRO (CeD-GSRS/CeD-PRO), diet logs,
  ideally urine/stool gluten immunogenic peptides for objective adherence.
- **The killer feature: objective adherence measurement.** Hawthorne bias kills celiac trials
  because nobody can tell "the drug worked" from "the patient got better at the diet." GIP
  testing + IL-2 separates those two for the first time.
- Sell/licence the cohort + endpoint package to the sponsors: argenx, Teva, Takeda, Topas,
  Barinthus all need exactly this and none of them owns it.

**Why it's good for you specifically:** it's a data/software business with a real moat, it
directly unblocks every drug program in the field, and it doesn't require you to be a
biologist. It's also the fastest route to being *inside* the field rather than outside it.

**Risk:** regulated-diagnostic territory (CLIA/LDT, and EU IVDR if you go there), and cohort
recruitment is slow. Partner with an existing patient org (CDF, Beyond Celiac, national
associations) rather than building recruitment from zero.

---

## P3 — Force the combination experiment
**A document + ~$200–500k. Highest science-per-dollar on this list.**

**The hypothesis** (from Wall 3): tolerance induction fails because it is administered into an
IL-15-licensed gut. Sequence it instead — IL-15/CD122 blockade first, tolerogenic antigen into
the permissive window, then withdraw both and challenge.

**The experiment:** the DQ8-Dd-villin-IL-15tg mouse (Abadie & Jabri, Chicago) is the only model
that develops real gluten- and HLA-dependent villous atrophy and heals on a gluten-free diet.
Four arms: vehicle / anti-IL-15 alone / tolerogenic nanoparticle alone / sequenced combination.
Primary readout: villous height:crypt depth after gluten re-challenge **after drug withdrawal**
— i.e. durable tolerance, not suppression.

**Your role:** write the rationale dossier, fund the study, and get it run in a lab that already
has the model. This is one grant-sized cheque and one well-argued 10-page document. If it works
in the mouse, it is directly translatable — every component already has human safety data.

**Risk:** the mouse is low-throughput and closely held; you need Chicago (or a licensee)
to cooperate. Approach as a funder-collaborator, not as an outsider with a theory.

---

## P4 — Protein design: the two open molecular problems
**Design is compute; validation needs a wet-lab partner.**

**(a) A glutenase that actually finishes the job.** TAK-062 failed and latiglutenase failed
because residual immunogenic peptides at ppm levels still trigger T cells — and stomach pH
inactivates most proteases. A Barcelona group published a designed molecule active at **pH 2**
in 2026, which was the field's hard blocker. The remaining problem is a well-posed engineering
one: proline-specific endopeptidase activity, acid-stable, gastric-residence-time kinetics,
degrading the *specific* epitope set from P1 below the T-cell activation threshold in a real
meal matrix. RFdiffusion/ProteinMPNN + directed evolution is the correct modern toolkit and
almost nobody is pointing it at this.

**(b) A better pHLA-DQ2.5 blocker.** DONQ52 (Chugai) proved the concept: an antibody that
recognizes 25+ distinct gluten:HLA-DQ2.5 complexes, neutralizes gluten-specific T-cell
activation, and leaves systemic immunity alone. It's a large IgG in Phase 1/2. The open space is
**smaller, cheaper, orally-viable blockers** — cyclic peptides, minibinders — designed against
the P1 epitope cover set. Design is now genuinely tractable; the bottleneck is a partner lab
with pHLA reagents and gluten-specific T-cell clones (Oslo, Leiden, Melbourne).

---

## P5 — The actual cure: delete the clones
**Highest ceiling, needs a real biotech. This is the 10-year answer.**

Two published routes, both orphaned:

**(a) Engineered gluten-specific Tregs** — Porret et al., *Sci Transl Med* 2025. Orthotopic
TCR replacement in human Tregs using DQ2.5-glia-α1a/α2-specific TCRs; suppresses effector T
cells in vitro and in vivo, and shows **bystander suppression** across epitopes (which partly
routes around Wall 2). No sponsor. Engineered-Treg platforms are entering the clinic in other
indications, so the manufacturing path exists.

**(b) pMHC-targeted deletion (the "inverted CAR")** — build a cytotoxic cell displaying
HLA-DQ2.5:deamidated-gliadin as *bait* so that any T cell recognizing it gets killed. This is
the CAAR-T logic proven in pemphigus (DSG3-CAAR-T), applied to a peptide-MHC ligand. **Nobody
has published this for celiac.** It is the only approach that would actually clear the
decades-persistent clonotypes, and celiac is arguably the best possible test case because the
target tissue fully regenerates and the antigen is exogenous and avoidable.

Why celiac is the *right* disease to try this in and why nobody has: the safety bar for a cell
therapy in a disease "managed by diet" is brutal. That's a real objection, not a stupid one.
The counter is refractory celiac disease type II — a pre-lymphoma condition with genuinely poor
prognosis where an aggressive cell therapy is clearly justified, and where you'd establish the
platform before moving to ordinary celiac. **That's the wedge indication.**

---

## P6 — Capital and coordination
- The **CDF Impact Fund** (launched Mar 2026, $15M target phase one, with Triple G Ventures)
  invests across detection, prevention, and cures. This is a live, named door.
- The structural fix is a **focused research organization**: celiac has the funding trap
  (few investigators → few grants → few expert reviewers → few investigators) that FROs exist
  to break. ~$25–50M over 5 years aimed *specifically* at the memory-T-cell eradication problem
  would be more than the entire field's annual NIH allocation and would target the one question
  nobody is funded to ask.
- A **prize** for the first demonstrated durable clearance of gluten-specific clonotypes in a
  human (measured by TCR-seq + IL-2 release, both now standard) would be cheap and would
  redirect attention to the actual endpoint of cure rather than symptom management.

---

## What "five strong engineers for three weeks" should actually do

Not "cure celiac." Produce one field-moving artifact:

| Week | Work |
|---|---|
| 1 | Assemble the corpus: all prolamin sequences across cultivars; the ~40 experimentally confirmed epitopes as ground truth; TG2 deamidation rules from the enzymology literature |
| 2 | Build the pipeline: deamidation enumeration → DQ2.5/DQ2.2/DQ8 binding ensemble → structural rescoring → TCR-accessibility filter. Calibrate on held-out known epitopes; report honest AUC |
| 3 | Ship: public repo, a per-cultivar celiac-safety score, and **the minimal epitope cover-set number** with confidence intervals. Preprint it. Send it to Sollid (Oslo), Jabri/Abadie (Chicago), the CSIC Córdoba wheat group, and Topas/Barinthus |

If the cover-set number comes back much larger than 4, you have a concrete, publishable
explanation for why every tolerance trial has failed — and a specification for the one that
wouldn't. That's a genuine contribution, and it costs three weeks and some GPU time.

---

## What not to do

- **Don't build another gluten-detection app or GF restaurant finder.** Solved, crowded, and
  irrelevant to a cure.
- **Don't self-experiment or run anything on the person you know.** Undiagnosed persistent
  villous atrophy after 20 years is a real risk that needs a gastroenterologist, not a project.
- **Don't try to do biology remotely.** Every wet-lab item here needs a real lab. Fund one,
  join one, or partner with one — but don't pretend a simulation is an experiment.
- **Don't trust any single class-II binding predictor.** This subfield has a long history of
  confident, wrong predictions. Ensemble, calibrate, and publish the failure modes.
