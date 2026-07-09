# `mtclone` — shared mtDNA lineage-tracing pipeline

Assay-agnostic toolkit for somatic-mtDNA clonal analysis of single cells. Built once, used by
two questions: **NK development** (are CD56-bright and CD56-dim NK one lineage or two?) and
**tumor clonal restriction** (are tumor-infiltrating NK a few expanded clones or a polyclonal
infiltrate?). Both reduce to the same operations on a cell × mtDNA-variant matrix.

## Install

```bash
conda env create -f environment.yml
conda activate mtclone
# run from the scripts/ dir (parent of mtclone/) so `import mtclone` resolves,
# or add scripts/ to PYTHONPATH.
```

## Pipeline

```
io       ingest heterogeneous deposits -> one AnnData schema (see schema.md)
qc       variant QC + informative-variant selection + binarization
clones   infer clones (shared-variant graph -> communities), per donor
metrics  clone-size distribution, Gini/Shannon, cross-group sharing, permutation null
classify NK identification + bright/dim (marker or data-driven) labelling
```

**Design fact:** both anchor deposits ship *already-called* per-cell mtDNA matrices
(GSE302113 heteroplasmy TSVs; GSE219014 ReDeeM `.mtx`), so the core engine runs on a called
matrix. Re-calling from raw reads (mgatk/redeemR) is the optional, deferred `recall/` module.

## Worked example

```python
import mtclone
from mtclone import io, qc, clones, metrics, classify

# 1. ingest (mtscATAC / GSE302113 style)
a = io.read_mtscatac_heteroplasmy(
        "GSM9096509.cell_heteroplasmic_df.tsv.gz",
        "GSM9096509.variant_stats.tsv.gz",
        dataset="GSE302113", sample="GSM9096509",
        donor="SU-L-003", tissue="PBMC", site="blood")
io.validate_schema(a)                       # enforces schema.md invariants

# 2. QC + informative-variant selection (+ binarize)
sel = qc.select_informative_variants(a, min_cell_coverage=10,
                                     max_pseudobulk_het=0.01,   # <1% pseudobulk window
                                     min_cells_detected=5,
                                     binarize_threshold=0.07)

# 3. clones (per donor — clones are donor-private)
called = clones.call_clones_per_donor(sel, method="graph",
                                      min_shared_variants=1, edge_weight_cutoff=0.5,
                                      min_clone_size=2)

# 4. metrics
print(metrics.clonality_summary(called))            # Gini, Shannon, max clone frac
share = metrics.clone_sharing_matrix(called, group_key="site")   # tumor vs blood etc.

# development question: are two subsets clonally mixed?
res  = metrics.between_vs_within_sharing(called, "bright", "dim", group_key="nk_subset")
null = metrics.permutation_null(
          called,
          lambda ad_: metrics.between_vs_within_sharing(ad_, "bright","dim",
                                                         group_key="nk_subset")["frac_mixed_clones"],
          group_key="nk_subset", stratify_by="donor", n=1000)
# low observed vs null  -> two lineages;  observed ~ null -> continuum
```

## Thresholds (defaults + rationale)

| parameter | default | rationale |
|---|---|---|
| `min_cell_coverage` | 10 | mtscATAC field standard ≥10× mito coverage for confident calls (Lareau 2021). Only applied as an absolute floor when `coverage_source=='depth_file'`; otherwise coverage is a relative proxy and only zero-coverage cells are dropped (a warning is issued). |
| `max_pseudobulk_het` | 0.01 | drop variants >1% pseudobulk heteroplasmy — germline/homoplasmic/ubiquitous, not clone-informative (Liu et al. conservative threshold). Loosen toward 0.9 only for synthetic/structural tests. |
| `min_cells_detected` | 5 | a variant seen in <5 cells can't define a clone and inflates noise. |
| `binarize_threshold` | 0.07 | heteroplasmy ≥0.07 → alt call (Liu default); exact fractions are unreliable at low mito coverage, so the graph works on presence/absence. Sweep with `qc.sweep_binarization`. |
| `edge_weight_cutoff` | 0.5 | min Jaccard (shared/union of alt variants) for a cell–cell edge before community detection. |
| `min_clone_size` | 2 | singletons carry no lineage information. |

## Clone inference

`method='graph'` (default): cell–cell graph weighted by Jaccard of shared alt variants →
Leiden communities (falls back to connected components if `leidenalg` is unavailable).
`method='variant_group'`: cells sharing their rarest common variant form a clone — an
independent cross-check. Validated on synthetic data: graph caller recovers known clone
structure at ARI > 0.9 (`tests/test_synthetic.py`).

## Testing

```bash
# unit tests (synthetic clone structure — no downloads)
PYTHONPATH=.. python -m pytest tests/ -q

# end-to-end gate on one real sample per dataset (needs the testdata files)
PYTHONPATH=.. python tests/run_end_to_end.py --testdata <dir> --out <dir>
```

## Known data-availability note (feeds Plan 2)

GSE219014's ReDeeM `*_mtDNA_Variants_matrix.mtx.gz` ships on GEO **without** companion
`barcodes`/`features` files, and its dimensions match neither axis of the paired RNA h5. Cell
barcodes and variant identities are therefore **not recoverable from GEO alone** — real
ReDeeM analysis must source them from the authors' redeemR repository / Zenodo deposit. The
`read_redeem_mtx` adapter is correct; the inputs it needs (barcodes, features, ideally depth)
must be obtained before Plan 2 can use this dataset for biology.

## `recall/` (deferred)

Optional consistency re-calling (mgatk for mtscATAC BAMs; redeemR for ReDeeM raw genotypes).
Not on the critical path — published calls are the primary input. First task when built: an
availability check for the required raw inputs (BAMs usually SRA/EGA controlled-access;
ReDeeM raw-genotype tables' presence for the human atlas is unconfirmed).
