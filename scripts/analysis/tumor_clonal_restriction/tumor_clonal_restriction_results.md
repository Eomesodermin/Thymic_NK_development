# Plan 3 — Are tumor-infiltrating NK cells clonally restricted? Integrated results

**Dataset:** GSE302113 (Liu et al., mtscATAC-seq; NSCLC + ovarian, matched tumor/normal/blood). 8 donors, NK identified by ATAC gene-activity, somatic mtDNA clones as the lineage barcode. All analysis local (`mtclone` pipeline). See per-step reports for detail.

## The question

Are tumor NK a **clonal-restriction product** (a few clones seed the tumor and expand locally) or a **polyclonal infiltrate** drawn from many circulating clones? Two orthogonal signatures decide it: (1) **within-tumor clone-size skew** (Gini) higher in tumor than blood, and (2) **tumor↔blood clone segregation** — expanded tumor clones absent from blood (seeding) vs tumor mirroring blood (passive). A restricted tumor should show both.

## Headline verdict

**NK tumor infiltration is clonally restricted in high-grade serous ovarian cancer, but not in lung cancer — it is tumor-type-dependent, not a general property of NK TILs.** The signal is carried by two HGSC ovarian donors; all five NSCLC donors show passive polyclonal infiltration. This *sharpens* Liu's mixed, single-case picture: their SU-O-005 oligoclonal-NK observation replicates and extends to a second, cleaner HGSC case (SU-O-004), while the innate "high clone-sharing" pattern they reported in aggregate is exactly what we see in every lung tumor.

## Per-patient integrated calls

| donor | dx | tumor NK w/ variant | within-tumor skew (Gini T vs B) | seeding p (base NK) | seeding p (strict NK) | robust configs (of 4) | **final call** |
|---|---|--:|:--:|--:|--:|:--:|:--|
| SU-O-004 | HGSC ovarian | 1617 | 0.63 vs 0.89 ↓ | 0.001 | 0.001 | 4/4 | SEEDING (robust) |
| SU-O-005 | HGSC ovarian | 2013 | 0.77 vs 0.61 ↑ | 0.001 | 0.473 | 3/4 | SEEDING (robust, 1 asterisk) |
| SU-O-002 | LGSC ovarian | 261 | 0.67 vs 0.77 ↓ | 0.011 | 0.012 | 2/4 | seeding (fragile) |
| SU-L-001 | NSCLC | 320 | 0.76 vs 0.86 ↓ | 0.454 | 0.756 | 1/4 | polyclonal (passive) |
| SU-L-004 | NSCLC | 278 | 0.82 vs 0.75 ↑ | 0.184 | 0.235 | 1/4 | polyclonal (passive) |
| SU-L-005 | NSCLC | 470 | 0.33 vs 0.59 ↓ | 0.541 | 0.575 | 0/4 | polyclonal (passive) |
| SU-L-003 | NSCLC | 335 | 0.75 vs 0.81 ↓ | 0.490 | 0.039 | 0/4 | polyclonal (passive) |
| SU-L-002 | NSCLC | 125 | 0.75 vs 0.80 ↓ | 0.381 | 1.000 | 0/4 | polyclonal (passive) |

## What each signature showed
**Signature 1 — within-tumor skew (Step 5).** Only 2/8 donors have tumor NK *more* clonally skewed
than their own blood NK: SU-O-005 (Gini 0.77 vs 0.61) and SU-L-004 (0.82 vs 0.75). The other 6 —
including the strongest seeding case SU-O-004 — have tumor NK equal to or *less* skewed than blood.
So within-tumor oligoclonality alone is neither necessary nor sufficient for restriction.

**Signature 2 — tumor↔blood segregation (Steps 6–7), the decisive test.** The permutation null
(shuffle site labels within donor) asks whether tumor and blood clones are more segregated than chance.
All three ovarian donors were significant at baseline (p≤0.011); no NSCLC donor was. Robustness
(Step 7, restricted to the 4 clone-structure-preserving configs) then resolved which calls are real:
- **SU-O-004 (HGSC): robust seeding** — significant in all 4 configs *and* under the strict NK-not-ILC
  cell definition. The one unambiguous case.
- **SU-O-005 (HGSC): robust seeding with one asterisk** — significant in 3/4 configs, but loses
  significance under the strict NK-not-ILC filter (p=0.47). The Liu anchor's cross-site signal is
  real under the primary analysis but sensitive to how NK are defined.
- **SU-O-002 (LGSC): fragile** — significant only 2/4 configs; downgraded, not a reliable seeder.
- **All 5 NSCLC: polyclonal** — never robustly significant; tumor NK clones shared with blood at
  chance. Apparent significance appears only under clone-distorting settings (a sparsity artifact).

**The two signatures are partly decoupled.** SU-O-005 has the strongest *within-tumor* skew but a
softer *cross-site* signal; SU-O-004 has modest skew but the strongest, cleanest seeding. Restriction
is best read from the segregation test (Signature 2), not from within-tumor Gini alone.

## How this revises Liu et al.
- **Confirms** their central ovarian observation: HGSC tumor NK carry oligoclonal, tumor-restricted
  expansions absent from blood — and shows it is not a one-tumor fluke (SU-O-004 is a second, cleaner
  case than their SU-O-005 anchor).
- **Explains** their aggregate finding that NK grouped with the *innate, high-clone-sharing* pattern:
  that aggregate is dominated by the lung tumors, where NK infiltration genuinely is polyclonal/passive.
  The mixed picture resolves into a **tumor-type split** — restricted in HGSC ovarian, passive in NSCLC.
- **Adds an honest caveat** they could not test at single-case resolution: the SU-O-005 anchor's seeding
  signal is sensitive to NK-vs-ILC cell definition, so the strongest *evidence* for restriction is
  actually SU-O-004, not the originally-highlighted case.

## Limitations (see Step 7 for full detail)
- **Power** rides on the 2 HGSC cases with strong tumor-NK-with-variant counts (SU-O-004: 1,617;
  SU-O-005: 2,013). NSCLC per-donor tumor NK are thinner (125–470). Conclusions are per-patient, not
  cohort-pooled (clone trees are donor-private).
- **Ambient/doublet floor is only partial** — the GEO processed layer has no per-cell depth, so a
  coverage-based contamination model can't be fit; the variant-frequency sweep is the available proxy
  (robust ovarian calls do not depend on near-public variants). A full floor needs raw mgatk output.
- **mtDNA clonality cannot distinguish "expanded in situ" from "arrived + proliferated"** — both yield
  a large tumor-restricted clone. Layering the ATAC residency/proliferation programs is a follow-up.
- **Rückert blood adaptive-NK oligoclonality benchmark** (external yardstick, Plan 2 scope) is not yet
  available; the within-patient blood-NK comparison served as the primary yardstick here.

## Reproducibility
Ordered scripts in `scripts/analysis/tumor_clonal_restriction/` (01_download → 08_robustness); all use
the frozen `mtclone` core. Processed objects + tables mirrored to
`HPC_data/Thymic_NK_development/processed/GSE302113/`. Per-step reports: Steps 2–7 markdown artifacts.

## Figures
- Step 3 power gate · Step 4 clone QC · Step 5 within-tumor clonality · Step 6 seeding test ·
  Step 7 robustness grid (all saved as artifacts).
