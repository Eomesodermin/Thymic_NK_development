"""mtclone — shared mtDNA lineage-tracing pipeline.

A small, assay-agnostic toolkit for somatic-mtDNA clonal analysis of single cells.
Both scientific questions it serves (NK development bright-vs-dim; tumor clonal
restriction) reduce to the same operations on a cell x variant heteroplasmy matrix:

    io       -> ingest heterogeneous deposits into one AnnData schema (see schema.md)
    qc       -> variant QC + informative-variant selection + binarization
    clones   -> infer clones (shared-variant graph -> communities)
    metrics  -> clone-size distribution, Gini/Shannon, cross-group sharing, permutation null
    classify -> NK identification + bright/dim (or data-driven 2-cluster) labelling

Design note: both anchor deposits (GSE302113 mtscATAC, GSE219014 ReDeeM) ship *already
called* per-cell mtDNA matrices, so the core engine operates on a called matrix. Re-calling
from raw reads (mgatk/redeemR) lives in the optional `recall/` submodule and is not on the
critical path.
"""

SCHEMA_VERSION = "1.0"

REQUIRED_OBS = [
    "cell_id", "donor", "sample", "tissue", "site", "assay", "dataset",
    "coverage", "coverage_source",
]
REQUIRED_VAR = [
    "variant_id", "pos", "ref", "alt", "pseudobulk_heteroplasmy", "n_cells_detected",
]

__version__ = "0.1.0"

# Lazy imports so `import mtclone` is cheap and submodules can be developed independently.
from . import io, qc, clones, metrics, classify  # noqa: E402,F401

__all__ = ["io", "qc", "clones", "metrics", "classify",
           "SCHEMA_VERSION", "REQUIRED_OBS", "REQUIRED_VAR", "__version__"]
