# Plan 2 — CD56-bright / CD56-dim NK: one lineage or two?

Somatic-mtDNA clonal lineage tracing to test whether CD56-bright and CD56-dim NK cells are a
single lineage (continuum / shared progenitors) or two independent developmental pathways.
Two datasets, deliberately complementary:

- **Rückert** (GSE197008/197037) — sorted **blood** NK, ASAP-seq (mtDNA via mgatk) + surface
  protein (CD56/CD16 gate bright/dim directly). High-n, carries the statistical weight.
- **ReDeeM** (Weng GSE219014; Figshare 24418966 + 23290004) — **bone-marrow** BMMC, UMI-consensus
  mtDNA + paired RNA (Multiome). Progenitor-proximal complement; 4 usable donors.

## Run order

| # | Script | Does |
|---|--------|------|
| 1a | `export_ruckert_mgatk.R` | mgatk allele counts → per-cell heteroplasmy + coverage (strand-concordance QC, NA-drop) |
| 1b | `ingest_ruckert.py` | build `.h5ad` per mtASAP sample; donor demux (hashtags/donor.csv); 21-plex ADT → obsm |
| 2a | `export_redeem_seurat.R` | export STD.CellType + RNA markers + wnn.umap from annotated Seurat |
| 2b | `ingest_redeem.py` | join Consensus.final heteroplasmy × Seurat NK annotation via ATAC↔RNA barcode translation |
| 3 | `step3_classify_power.py` | subset NK, classify bright/dim (Rückert protein gate, ReDeeM RNA), power gate (≥1 informative variant) |
| 4 | `step4_clone_inference.py` | donor-private clone inference (ReDeeM graph, Rückert variant_group) |
| 5 | `step5_sharing_test.py` | **the core test:** bright/dim clone sharing vs within-donor permutation null |
| 6 | `step6_clonal_partition.py` | converse: does the clonal partition map onto the bright/dim axis? (η² vs null) |
| 7 | `step7_robustness.py` | QC/caller sweep + doublet floor + evidence-asymmetry + tissue caveats |

Steps 3–7 read/write `HPC_data/Thymic_NK_development/processed/plan2_step{3..7}/`. Ingestion
adapters live in `scripts/mtclone/` (`io.read_mtscatac_heteroplasmy`, `io.read_redeem_consensus`);
the two dataset-agnostic helpers (mgatk strand-concordance, Multiome barcode translation) are also
published as the personal skill `sc-mtdna-clone-ingest`.

## Verdict

**Continuum, not two lineages** — see `development_brightdim_results.md`. Bright and dim NK share
somatic-mtDNA clones exactly as much as a one-well-mixed-lineage null predicts, replicated across
blood + marrow, two assays, and two label modalities, with demonstrated power to have detected
segregation.
