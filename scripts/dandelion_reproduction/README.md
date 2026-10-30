# Dandelion NK/ILC-from-DN reproduction + driver dissection

Reproduces Suo et al. (Dandelion, *Nat Biotechnol* 2023, 10.1038/s41587-023-01734-7)
Fig 5b/c — the pseudotime trajectory in a TRBJ V(D)J feature space that places ILC/NK
cells on an aborted T-development branch — and then dissects whether that branch is
**TCR-footprint-driven** or **transcriptome-driven**.

Full write-up: `results/dandelion_reproduction/dandelion_reproduction_memo.md`.

## Pipeline order

| step | script | what |
|---|---|---|
| 1 | `cellranger vdj --chain TR` (run on an HPC cluster, not included) | align the 64 αβTCR libraries (E-MTAB-11388) on `all_contig`; re-download and re-align truncated FASTQs |
| 1b | `sc-dandelion` preprocessing (HPC, not included) | IgBLAST re-annotation to 64 AIRR TSVs |
| 2 | `01_subset_trajectory_cells.py` | subset the Suo lymphoid GEX to the recipe cell types (111,706 cells), keep raw counts |
| 3 | `03_build_trajectory.py` | join VDJ and GEX, check contigs, scVI, neighbourhoods, TRBJ pseudobulk, Palantir |
| 4 | `04_driver_dissection.py` | (A) GEX-only control, (B) TRBJ permutation null (n = 200), (C) productive-only |
| 5 | `05_plot_fig5.py` | Fig 5b (neighbourhood feature UMAP) and Fig 5c (per-cell branch probabilities) |

There is no `02_*`: the alignment and preprocessing in steps 1 and 1b are the second stage, and
`01_*` is the GEX subset that the VDJ output is joined onto in `03_*`.

## Inputs and outputs

All paths are under `$NK_DATA_DIR` (default `./data`; see `scripts/mtclone/paths.py`).

- `raw/suo_fetal/lymphoid.h5ad`: the CELLxGENE Suo lymphoid object (GEX only; CELLxGENE strips TCR).
- `processed/dandelion_reproduction/dandelion_tsv/`: the 64 Dandelion AIRR TSVs from step 1b.
- `processed/dandelion_reproduction/trajectory_gex_subset.h5ad`: step 2 output.
- `processed/dandelion_reproduction/trajectory_adata_scvi_nhoods.h5ad`: scVI and neighbourhood
  checkpoint, so Palantir can be re-run without the ~25-minute scVI step.
- `processed/dandelion_reproduction/trajectory_{pseudobulk,percell}.h5ad`: trajectory results.
- `dandelion_reproduction/suo_abTCR_vdj_manifest.json`: the library manifest (accession, path and
  sample for each of the 64 libraries).
- `results/dandelion_reproduction/`: tables (`branch_prob_by_celltype.csv`,
  `driver_dissection_summary.csv`) and the memo. Figures are written to `$NK_RESULTS_DIR`
  (default `./results`) and the report copies are in `report_assets/tcr_fig*.png`.
