# Plan 3 — Tumor NK clonal-restriction analysis (GSE302113)

Ordered, self-contained pipeline. All scripts use the frozen `mtclone` core
(`../../mtclone/`) and run locally in the `mtclone` conda env. Data live under
`~/Documents/HPC_data/Thymic_NK_development/{raw,processed}/GSE302113/`.

## Question
Are tumor-infiltrating NK cells clonally restricted (few clones seed + expand locally) or a
polyclonal passive infiltrate? mtDNA somatic clones = lineage barcode.

## Scripts (run in order)
| # | script | does | key output |
|---|---|---|---|
| 01 | download_fragments.py | fetch 36 ATAC fragment files (~35.5 GB) from GEO | raw/GSE302113/*_fragments.tsv.gz |
| 02 | ingest_gse302113.py | per-cell heteroplasmy TSV -> per-donor AnnData (8 donors) | processed/GSE302113/<donor>.h5ad |
| 03 | gene_activity_from_fragments.py | ArchR-style gene activity from fragments (31-gene panel) | geneactivity/<GSM>.geneactivity.h5ad |
| 04 | nk_label_and_power_gate.py | NK ID (mtclone.classify.label_nk atac) + power gate | nk_labelled/<donor>.nk_labelled.h5ad, nk_power_*.csv |
| 05 | clone_inference.py | mtclone.clones.call_clones per donor, joint across sites | clones/<donor>.{nk,nk_not_ilc}.clones.h5ad |
| 06 | within_tumor_clonality.py | Gini + normalized Shannon per site (Signature 1) | within_tumor_clonality.csv |
| 07 | tumor_blood_overlap.py | clone sharing + seeding permutation null (Signature 2) | tumor_blood_{sharing,permutation}.csv |
| 08 | robustness.py | seeding robustness across variant/caller settings | robustness_{grid,pvalue_matrix}.csv |

## Headline result
NK tumor infiltration is clonally restricted in **HGSC ovarian** (SU-O-004 robust; SU-O-005 robust
with an NK-definition asterisk) but **passive/polyclonal in all 5 NSCLC** donors — tumor-type-dependent.
See `tumor_clonal_restriction_results.md` (integrated writeup) and per-step markdown reports.

## Notes
- NK identification requires ATAC fragments (steps 01+03); the mtDNA heteroplasmy layer (step 02) is
  small and fragment-free.
- Two NK sets carried throughout: base NK and strict `nk_not_ilc` (ILC-excluded) as a sensitivity check.
- Step 07's stringent variant filters (max_pb<=0.1) and aggressive edge cutoffs DISTORT clone structure
  and produce spurious significance — the verdict uses only the 4 structure-preserving configs.
