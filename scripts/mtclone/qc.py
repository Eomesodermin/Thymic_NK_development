"""qc.py — variant QC, informative-variant selection, binarization.

All thresholds are parameters with documented defaults drawn from the Liu et al. mtscATAC
workflow (and standard mgatk practice). Nothing here is dataset-specific; the caller passes
the thresholds appropriate to the assay.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import anndata as ad

# ---- documented defaults (see README / Liu et al. methods) -------------------------------
DEFAULTS = dict(
    min_cell_coverage=10.0,      # per-cell mtDNA coverage floor (x). mtscATAC standard >=10x.
    max_pseudobulk_het=0.01,     # drop variants with >1% pseudobulk heteroplasmy
                                 #   (germline/homoplasmic/ubiquitous — not clone-informative).
    min_cells_detected=5,        # variant must be seen in >=5 cells to be usable.
    min_strand_correlation=None, # optional strand-balance floor (e.g. 0.65); None = off.
    binarize_threshold=0.07,     # heteroplasmy >= 0.07 -> alt call (Liu default).
)


def compute_qc(adata: ad.AnnData) -> ad.AnnData:
    """Attach/refresh per-variant pseudobulk stats. Idempotent."""
    het = adata.X
    adata.var["pseudobulk_heteroplasmy"] = np.asarray(het.mean(axis=0)).ravel().astype(np.float32)
    adata.var["n_cells_detected"] = np.asarray((het > 0).sum(axis=0)).ravel().astype(int)
    adata.obs["n_variants_detected"] = np.asarray((het > 0).sum(axis=1)).ravel().astype(int)
    return adata


def select_informative_variants(
    adata: ad.AnnData,
    *,
    min_cell_coverage: float = DEFAULTS["min_cell_coverage"],
    max_pseudobulk_het: float = DEFAULTS["max_pseudobulk_het"],
    min_cells_detected: int = DEFAULTS["min_cells_detected"],
    min_strand_correlation=DEFAULTS["min_strand_correlation"],
    binarize_threshold: float = DEFAULTS["binarize_threshold"],
    inplace: bool = False,
) -> ad.AnnData:
    """Filter cells by coverage, select clone-informative variants, and binarize.

    Steps:
      1. drop cells below `min_cell_coverage` (skipped for derived/relative coverage — see note);
      2. select variants with pseudobulk heteroplasmy in (0, max_pseudobulk_het],
         detected in >= min_cells_detected cells, passing optional strand filter;
      3. write layers['binary'] = (heteroplasmy >= binarize_threshold).

    Returns a filtered copy (or mutates in place). Records applied thresholds in
    uns['thresholds'].

    NOTE on coverage: when obs['coverage_source'] != 'depth_file', `coverage` is a *relative*
    proxy (variant-count based), not an absolute x-coverage, so the absolute
    `min_cell_coverage` floor is not comparable. In that case the floor is applied as a
    *percentile-safe* minimum of >=1 (drop only zero-coverage cells) and a warning is issued.
    """
    a = adata if inplace else adata.copy()
    a = compute_qc(a)

    # ---- 1. cell coverage filter
    cov = a.obs["coverage"].values.astype(float)
    src = a.obs["coverage_source"].iloc[0] if "coverage_source" in a.obs else "unknown"
    if src == "depth_file":
        cell_mask = cov >= min_cell_coverage
    else:
        import warnings
        warnings.warn(
            f"coverage_source='{src}': coverage is a relative proxy, not absolute x; "
            f"applying floor as coverage>=1 instead of >={min_cell_coverage}."
        )
        cell_mask = cov >= 1
    a = a[cell_mask].copy()
    a = compute_qc(a)

    # ---- 2. informative-variant selection
    pv = a.var["pseudobulk_heteroplasmy"].values
    nd = a.var["n_cells_detected"].values
    var_mask = (pv > 0) & (pv <= max_pseudobulk_het) & (nd >= min_cells_detected)
    if min_strand_correlation is not None and "strand_correlation" in a.var:
        sc = a.var["strand_correlation"].values
        var_mask &= np.where(np.isnan(sc), True, sc >= min_strand_correlation)
    a = a[:, var_mask].copy()

    # ---- 3. binarize
    binarize_heteroplasmy(a, threshold=binarize_threshold, inplace=True)

    a.uns["thresholds"] = dict(
        min_cell_coverage=min_cell_coverage, max_pseudobulk_het=max_pseudobulk_het,
        min_cells_detected=min_cells_detected, min_strand_correlation=min_strand_correlation,
        binarize_threshold=binarize_threshold,
        n_cells_kept=int(a.n_obs), n_variants_kept=int(a.n_vars),
    )
    return a


def binarize_heteroplasmy(adata: ad.AnnData, *, threshold: float = 0.07,
                          inplace: bool = True) -> ad.AnnData:
    """layers['binary'] = (X >= threshold) as uint8 sparse."""
    a = adata if inplace else adata.copy()
    X = a.X.tocsr()
    b = X.copy()
    b.data = (b.data >= threshold).astype(np.uint8)
    b.eliminate_zeros()
    a.layers["binary"] = b
    return a


def sweep_binarization(adata: ad.AnnData, thresholds=(0.03, 0.05, 0.07, 0.10, 0.20)) -> dict:
    """Return {threshold: (n_alt_calls, mean_variants_per_cell)} for sensitivity analysis."""
    out = {}
    X = adata.X.tocsr()
    for t in thresholds:
        n_alt = int((X.data >= t).sum())
        vpc = np.asarray((X >= t).sum(axis=1)).ravel().mean()
        out[float(t)] = (n_alt, float(vpc))
    return out


# ----------------------------------------------------------------------------- QC plots
def qc_plots(adata: ad.AnnData, path=None):
    """Coverage histogram, pseudobulk-heteroplasmy distribution, variants-per-cell.

    Returns the matplotlib Figure. Saves to `path` if given.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    a = compute_qc(adata)
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6))

    axes[0].hist(a.obs["coverage"].values, bins=50, color="#3b6ea5")
    axes[0].set(title="Per-cell mtDNA coverage", xlabel="coverage (x or proxy)", ylabel="cells")

    pv = a.var["pseudobulk_heteroplasmy"].values
    pv = pv[pv > 0]
    axes[1].hist(np.log10(pv + 1e-6), bins=50, color="#a53b5a")
    axes[1].set(title="Pseudobulk heteroplasmy", xlabel="log10 pseudobulk het", ylabel="variants")

    axes[2].hist(a.obs["n_variants_detected"].values, bins=50, color="#4a8c5a")
    axes[2].set(title="Variants detected per cell", xlabel="n variants", ylabel="cells")

    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight")
    return fig
