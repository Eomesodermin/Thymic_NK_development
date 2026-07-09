"""classify.py — NK identification + bright/dim (or data-driven) subset labelling.

These operate on a paired-modality object: an expression AnnData (RNA counts) or an
ATAC gene-activity AnnData, aligned to the mtDNA object by cell_id. They write labels onto
obs which the mtDNA object then inherits by a join on cell_id.

The development question explicitly allows "two clonal lineages that MAY map to bright/dim"
rather than requiring a priori labels, so both a marker-threshold labeller and a data-driven
2-cluster option are provided.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import anndata as ad

# ---- marker panels -----------------------------------------------------------------------
NK_MARKERS_POS = ["NCAM1", "NCR1", "GNLY", "KLRD1", "NKG7", "KLRF1"]
NK_MARKERS_NEG = ["CD3D", "CD3E", "CD3G", "CD19", "MS4A1", "CD14"]
# ATAC gene-activity criterion (Liu): NK = NCR1/GNLY accessible, CD3D/CD8A closed.
NK_ATAC_POS = ["NCR1", "GNLY", "KLRD1"]
NK_ATAC_NEG = ["CD3D", "CD8A"]

BRIGHT_MARKERS = ["NCAM1", "GZMK", "SELL", "XCL1", "IL7R"]   # CD56-bright, cytokine-type
DIM_MARKERS = ["FCGR3A", "PRF1", "FGFBP2", "GZMB", "CX3CR1"] # CD56-dim/CD16+, cytotoxic


def _score(adata: ad.AnnData, genes) -> np.ndarray:
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return np.zeros(adata.n_obs)
    X = adata[:, present].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    # standardize each gene then average -> a simple signature score
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1
    return ((X - mu) / sd).mean(1)


def label_nk(adata: ad.AnnData, *, modality="rna", pos_thresh=0.0, neg_thresh=0.25,
             key_added="cell_type") -> ad.AnnData:
    """Label cells NK vs other by marker score. Writes obs[key_added] in {'NK','other'}.

    modality='rna'  -> NK_MARKERS_POS high & NK_MARKERS_NEG low.
    modality='atac' -> gene-activity NK_ATAC_POS open & NK_ATAC_NEG closed (Liu criterion).
    """
    pos_genes = NK_ATAC_POS if modality == "atac" else NK_MARKERS_POS
    neg_genes = NK_ATAC_NEG if modality == "atac" else NK_MARKERS_NEG
    pos = _score(adata, pos_genes)
    neg = _score(adata, neg_genes)
    is_nk = (pos > pos_thresh) & (neg < neg_thresh)
    adata.obs[key_added] = np.where(is_nk, "NK", "other")
    adata.obs["nk_pos_score"] = pos
    adata.obs["nk_neg_score"] = neg
    return adata


def label_bright_dim(adata: ad.AnnData, *, nk_mask_key="cell_type", nk_value="NK",
                     method="marker", key_added="nk_subset") -> ad.AnnData:
    """Label NK cells bright vs dim.

    method='marker'      : bright if BRIGHT score > DIM score, else dim.
    method='data_driven' : 2-means on the (DIM - BRIGHT) axis among NK, then orient the
                           cluster with higher FCGR3A as 'dim'. Supports the honest
                           "distinct subset vs remainder" framing when markers are ambiguous.
    Non-NK cells get 'unassigned'.
    """
    sub = adata.obs.get(nk_mask_key)
    nk = np.ones(adata.n_obs, bool) if sub is None else (adata.obs[nk_mask_key] == nk_value).values
    out = np.array(["unassigned"] * adata.n_obs, dtype=object)

    if nk.sum() == 0:
        adata.obs[key_added] = out
        return adata

    a_nk = adata[nk]
    bright = _score(a_nk, BRIGHT_MARKERS)
    dim = _score(a_nk, DIM_MARKERS)

    if method == "marker":
        lab = np.where(bright >= dim, "bright", "dim")
    elif method == "data_driven":
        from sklearn.cluster import KMeans
        axis = (dim - bright).reshape(-1, 1)
        km = KMeans(n_clusters=2, n_init=10, random_state=0).fit(axis)
        cl = km.labels_
        # orient: cluster with higher mean FCGR3A -> dim
        fcg = _score(a_nk, ["FCGR3A"])
        mean0, mean1 = fcg[cl == 0].mean() if (cl == 0).any() else -np.inf, \
                       fcg[cl == 1].mean() if (cl == 1).any() else -np.inf
        dim_cluster = 0 if mean0 >= mean1 else 1
        lab = np.where(cl == dim_cluster, "dim", "bright")
    else:
        raise ValueError(f"unknown method {method!r}")

    out[nk] = lab
    adata.obs[key_added] = out
    adata.obs.loc[adata.obs_names[nk], "nk_bright_score"] = bright
    adata.obs.loc[adata.obs_names[nk], "nk_dim_score"] = dim
    return adata


def transfer_labels(mt_adata: ad.AnnData, expr_adata: ad.AnnData,
                    cols=("cell_type", "nk_subset"), on="cell_id") -> ad.AnnData:
    """Join classification labels from the expression object onto the mtDNA object by cell_id."""
    right = expr_adata.obs.set_index(on)[[c for c in cols if c in expr_adata.obs.columns]]
    left_key = mt_adata.obs[on] if on in mt_adata.obs else mt_adata.obs_names.to_series()
    joined = left_key.map(lambda cid: right.loc[cid].to_dict() if cid in right.index else {})
    for c in cols:
        if c in right.columns:
            mt_adata.obs[c] = [d.get(c, "unassigned") if isinstance(d, dict) else "unassigned"
                               for d in joined]
    return mt_adata


def labelled_umap(adata: ad.AnnData, color=("cell_type", "nk_subset"), path=None):
    """Quick QC UMAP coloured by labels. Computes a UMAP if absent. Returns the Figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    a = adata.copy()
    if "X_umap" not in a.obsm:
        if "X_pca" not in a.obsm:
            sc.pp.pca(a, n_comps=min(30, a.n_vars - 1))
        sc.pp.neighbors(a, n_neighbors=15)
        sc.tl.umap(a)
    cols = [c for c in color if c in a.obs.columns]
    fig, axes = plt.subplots(1, len(cols), figsize=(6 * len(cols), 5), squeeze=False)
    for ax, c in zip(axes[0], cols):
        cats = pd.Categorical(a.obs[c].astype(str))
        for cat in cats.categories:
            m = cats == cat
            ax.scatter(a.obsm["X_umap"][m, 0], a.obsm["X_umap"][m, 1], s=4, label=str(cat))
        ax.set(title=c, xticks=[], yticks=[])
        ax.legend(markerscale=3, fontsize=7, loc="best")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=150, bbox_inches="tight")
    return fig
