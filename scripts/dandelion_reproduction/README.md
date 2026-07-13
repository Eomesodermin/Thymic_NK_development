# Dandelion NK/ILC-from-DN reproduction + driver dissection

Reproduces Suo et al. (Dandelion, *Nat Biotechnol* 2023, 10.1038/s41587-023-01734-7)
Fig 5b/c — the pseudotime trajectory in a TRBJ V(D)J feature space that places ILC/NK
cells on an aborted T-development branch — and then dissects whether that branch is
**TCR-footprint-driven** or **transcriptome-driven**.

Full write-up: `results/dandelion_reproduction/dandelion_reproduction_memo.md`.

## Pipeline order

| step | script | where | what |
|---|---|---|---|
| 1 | *(HPC)* `cr_vdj_suo_one.sh` + `submit_suo_all.sh` | Marvin | align 64 αβTCR libraries (E-MTAB-11388), `cellranger vdj --chain TR` on `all_contig` |
| 1b | *(HPC)* `cr_vdj_suo_fix*.sh`, `resubmit_suo_failed.sh` | Marvin | re-download truncated FASTQs (size-verify vs Content-Length) + realign |
| 1c | *(HPC)* `dandelion_preprocess_all.sh` | Marvin | sc-dandelion IgBLAST reannotation → 64 AIRR TSVs |
| 2 | `01_subset_trajectory_cells.py` | local | subset Suo lymphoid GEX to recipe cell types (111,706 cells), keep raw counts |
| 3 | `03_build_trajectory.py` | local | join VDJ↔GEX → check_contigs → scVI → neighborhoods → TRBJ pseudobulk → Palantir |
| 4 | `04_driver_dissection.py` | local | (A) GEX-only control, (B) TRBJ permutation null n=200, (C) productive-only |
| 5 | `05_plot_fig5.py` | local | Fig 5b (nhood feature UMAP) + Fig 5c (per-cell branch probs) |

The HPC scripts live in `HPC_workflows/thymic_nk_development/dandelion_reproduction/`
(mirrored here for traceability; they ran on Marvin under
`/lustre/scratch/data/dcorvino_hpc-thymic_nk_development/`). There is no local `02_*`:
step 2's alignment/preprocessing is the HPC block above; `01_*` is the GEX subset that
step 2's VDJ is joined onto in `03_*`.

## Inputs / outputs (under `HPC_data/Thymic_NK_development/`)

- `raw/suo_fetal/lymphoid.h5ad` — CELLxGENE Suo lymphoid object (GEX only; CZI strips TCR).
- `processed/dandelion_reproduction/dandelion_tsv/` — 64 Dandelion AIRR TSVs (from HPC).
- `processed/dandelion_reproduction/trajectory_gex_subset.h5ad` — step 2 output.
- `processed/dandelion_reproduction/trajectory_adata_scvi_nhoods.h5ad` — checkpoint (scVI +
  neighborhoods) so Palantir can be re-run without the ~25-min scVI prefix.
- `processed/dandelion_reproduction/trajectory_{pseudobulk,percell}.h5ad` — trajectory results.
- `results/dandelion_reproduction/` — figures (`fig5b_*`, `fig5c_*`, `fig_driver_dissection.png`),
  tables (`branch_prob_by_celltype.csv`, `driver_dissection_summary.csv`), and the memo.

VDJ library manifest: `data/dandelion_reproduction/suo_abTCR_vdj_manifest.{csv,json}`.

## Environment

Conda env `ddl-traj`: scanpy 1.10.4, dandelion 0.5.7, scvi-tools 1.3.3, palantir 1.4.4,
pertpy 0.10.0. Run the compute scripts standalone with
`NUMBA_CACHE_DIR=/tmp/numba_cache_ddltraj NUMBA_DISABLE_JIT=1 python <script>.py`
(numba `/dev/fd/3` locator bug otherwise breaks the scanpy import).

## Key result (one line)

The branch reproduces in direction (NK prob_NK/ILC 0.62 vs T ~0.48) but is statistically
fragile (permutation null p=0.29) and **not TCR-footprint-driven in this re-alignment** — the
V-less non-productive TRB relics the paper leans on do not appear in a standard
`cellranger vdj --chain TR` re-derivation.
