# NK developmental origin from gene-therapy clonal tracking (Six et al. 2020) — re-analysis

**Question.** Do NK cells share a clonal (HSPC) origin more with T cells (shared *lymphoid*
progenitor) or with myeloid cells (semi-distinct, myeloid-leaning NK ontogeny)? This is the
NK developmental-origin question the mtDNA modality could not answer (the ReDeeM adult-marrow
probe was a clean negative — clone resolution too coarse, common ancestor too deep to leave a
shared private mark on surviving leaves). Lentiviral integration-site (IS) clonal tracking in
HSPC gene-therapy patients answers it directly and at population-clone resolution.

**Why this data can answer it (and mtDNA could not).** Each transduced HSPC carries a unique
lentiviral integration site (`posid`) inherited by all its progeny. FACS-sorting the mature
lineages and sequencing IS per gate tells you which HSPC clones contributed to which lineage.
A clone shared between the sorted NK gate and the sorted T gate = one HSPC ancestor produced
both — a *direct* readout of shared progenitor output, caught at the source rather than inferred
between two surviving leaves. This sidesteps the mtscATAC coverage wall entirely.

**Data.** Six E et al., *Blood* 2020, "Clonal tracking in gene therapy patients reveals a
diversity of human hematopoietic differentiation programs" (doi 10.1182/blood.2019002350).
Processed IS table deposited directly in github.com/BushmanLab/HSC_diversity
(`data/intSites.mergedSamples.collapsed.csv.gz`, 404,335 IS rows × 28 cols; raw reads SRA
SRP139090, not needed). Six patients (4 Wiskott-Aldrich LV-GT: WAS2/4/5/7; two
β-hemoglobinopathy: b0/bE, bS/bS), five FACS-sorted lineages
(GRANULOCYTES, MONOCYTES, TCELLS, BCELLS, NKCELLS), timepoints to M78. No realignment,
no IS calling — pure re-analysis of the deposited table.

**Method.** Per patient, define a clone as a unique integration site (`posid`) and take its
presence/absence in each sorted lineage. For each partner lineage L, count NK∩L shared clones
and compare to the count expected if NK and L clones were drawn independently from that
patient's clone universe (**hypergeometric fold-enrichment = observed / expected**, with the
hypergeometric survival-function p-value). Fold-enrichment is size-controlled: it corrects for
the very different clone counts per gate, which is the coverage-matched control the ReDeeM probe
showed is mandatory. Headline contrast = NK–Monocyte vs NK–T fold, paired across the 6 patients
(Wilcoxon signed-rank).

## Result: NK clonal output leans myeloid, not T-lymphoid

Per-patient fold-enrichment of NK clone-sharing (obs/expected; 1.0 = no enrichment):

| patient | NK–T | NK–B | NK–Mono | NK–Gran | NK leans |
|---|---:|---:|---:|---:|---|
| WAS2  | 0.33 | 0.71 | 0.77 | 0.62 | myeloid |
| WAS4  | 0.63 | 2.12 | 4.53 | 3.02 | myeloid |
| WAS5  | 0.32 | 0.50 | 0.59 | 0.51 | myeloid |
| WAS7  | 0.40 | 1.23 | 1.92 | 1.15 | myeloid |
| b0/bE | 0.65 | 0.47 | 0.59 | 0.52 | T/lymphoid (marginal) |
| bS/bS | 1.62 | 1.00 | 1.75 | 1.21 | myeloid |
| **median** | **0.51** | **0.86** | **1.26** | **0.89** | |

**NK–Monocyte fold exceeds NK–T fold in 5 of 6 patients** (Wilcoxon signed-rank p = 0.0625,
the minimum attainable at n=6 with one dissenter; NK–myeloid-mean vs NK–T 4/6, p = 0.16).
NK–T sharing sits **below** chance in 5/6 patients (median fold 0.51) — NK clones are
*depleted* of T-cell co-membership — while NK–Monocyte sits at or above chance (median 1.26).
The lone exception, b0/bE, is marginal (NK–T 0.65 vs NK–Mono 0.59, both below 1.0) rather than
a genuine T-leaning case.

## Interpretation

- **NK is clonally closer to myeloid than to T at the HSPC level.** This reproduces and extends
  Biasco/Aiuti 2016 (Cell Stem Cell, TIGET/Milan cohort) — which first placed NK clonally
  nearer myeloid than other lymphoid — in an independent human cohort with an explicit
  size-controlled statistic. It is also *consistent in spirit* with the Dunbar rhesus-macaque
  barcoding work, but note the difference: the macaque studies mainly show NK (especially the
  CD16+ subset) is clonally *distinct/self-sustaining* — maintained by its own HSPC clones,
  separate from T and B — which is not the same as "aligned with myeloid." Both results point
  away from a simple shared T/NK lymphoid progenitor; only Six 2020 shows the specific
  myeloid co-occurrence. Do not state the macaque data as showing myeloid alignment.
- **The human evidence is genuinely MIXED — do not present myeloid-leaning as settled.**
  Scala 2021 (Nat Commun, doi 10.1038/s41467-021-21834-9) pulls the *opposite* direction: in
  patients who had lost the tag from myeloid and B cells, the tag persisted in *both* T and NK
  for 15 years — direct evidence of a long-term *lymphoid* progenitor sustaining T+NK together.
  This analysis (bulk of the clonal repertoire leans myeloid) and Scala (a minority long-term
  lymphoid progenitor feeds T+NK) can coexist — a rare lymphoid-restricted T+NK progenitor and
  a dominant myeloid-proximal NK output are not mutually exclusive — but the honest summary is
  that the two human clonal-tracking results disagree in emphasis, so NK origin is unsettled,
  not resolved in favour of myeloid.
- **Directional, not decisive.** n=6 patients, p=0.0625 — this is a consistent trend, not a
  significance-certified result. The disease context (WAS + β-hemoglobinopathy gene therapy,
  post-transplant reconstitution) may not reflect steady-state ontogeny, and IS clone-sharing
  is population-level (which HSPC clones feed which gate), not a per-cell lineage tree.

## Caveats

1. **Population-level, not per-cell.** "Shared clone" = the same HSPC integration-site clone
   contributed to both sorted gates; it is not a per-cell parent-daughter link.
2. **Transplant/gene-therapy setting**, not steady-state or thymic development. Reconstitution
   dynamics after conditioning may bias lineage output.
3. **Sort purity.** Cross-lineage contamination between FACS gates would inflate apparent
   sharing; the paper reports per-timepoint contamination matrices (`data/crossOverReports.tsv`)
   — a sensitivity analysis correcting for these is the natural next step if this is pursued.
4. **n=6 patients.** The direction is consistent (5/6) but the paired test is at its resolution
   floor; more IS-tracked cohorts with a sorted NK gate would firm it up. Candidate cohorts (by
   their actual groups — do not conflate): Scala 2021 (Naldini/Aiuti, TIGET/Milan); the broader
   TIGET WAS and ADA-SCID series (Biasco/Aiuti, Milan); and the Paris cohorts of Six/Cavazzana
   (Institut Imagine/Necker), of which the present study is one. NB: Six 2020 is a **Paris**
   cohort, NOT TIGET — TIGET is the Milan institute behind Biasco 2016 and Scala 2021.

## Bottom line

The gene-therapy clonal-tracking data give a direct, if underpowered, readout on the origin
question that mtDNA could not: **in these 6 patients, NK cells share HSPC clones with myeloid
cells more than with T cells** (5/6 patients, size-controlled, p≈0.06). This argues against NK
sharing a *simple* clonal origin with T cells — consistent with Biasco 2016 (NK clonally near
myeloid) and with the macaque finding that NK is clonally *distinct* from T/B. It stands in
tension with Scala 2021, which found a lymphoid progenitor sustaining T+NK together. **The
honest verdict is that human NK developmental origin is unsettled / mixed, with this analysis
adding one myeloid-leaning data point — not that the question is resolved.**
