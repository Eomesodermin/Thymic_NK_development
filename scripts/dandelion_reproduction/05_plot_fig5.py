#!/usr/bin/env python
"""
Step 5 output: reproduce Dandelion Fig 5b/c.

Fig 5b = UMAP of the neighborhood V(D)J (TRBJ) feature space, coloured by
         (i) cell type / lineage, (ii) pseudotime.
Fig 5c = per-cell branch probabilities (T vs NK/ILC) on the cell UMAP, plus
         the pseudotime overlay.

Inputs : processed/dandelion_reproduction/trajectory_pseudobulk.h5ad  (nhood-level, Fig 5b)
         processed/dandelion_reproduction/trajectory_percell.h5ad     (cell-level, Fig 5c)
Outputs: fig5b_vdj_feature_umap.png, fig5c_branch_probabilities.png
Env: ddl-traj (NUMBA_CACHE_DIR set).
"""
import warnings; warnings.filterwarnings("ignore")
import os, numpy as np, pandas as pd, scanpy as sc
import matplotlib as mpl, matplotlib.pyplot as plt

BASE = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/dandelion_reproduction"
OUT  = os.path.expanduser("~/Documents/Github/Eomesodermin/Thymic_NK_development/results/dandelion_reproduction")
os.makedirs(OUT, exist_ok=True)

# --- minimal publication rc (figure-style rules; ddl-traj has no skill kernel) ---
mpl.rcParams.update({
    "figure.dpi":150, "savefig.dpi":300, "font.size":8, "axes.titlesize":8,
    "axes.labelsize":8, "xtick.labelsize":6, "ytick.labelsize":6, "legend.fontsize":6,
    "axes.spines.top":False, "axes.spines.right":False, "figure.facecolor":"white",
})

def umap_panel(ax, adata, color, title, cmap=None, categorical=False, umap_key="X_umap"):
    xy = adata.obsm[umap_key]
    if categorical:
        cats = adata.obs[color].astype("category")
        codes = cats.cat.codes
        pal = plt.cm.tab20(np.linspace(0,1,len(cats.cat.categories)))
        for i,c in enumerate(cats.cat.categories):
            m = cats==c
            ax.scatter(xy[m,0], xy[m,1], s=3, color=pal[i], label=str(c), linewidths=0)
    else:
        vals = adata.obs[color].values if color in adata.obs else np.asarray(adata[:,color].X).ravel()
        sca = ax.scatter(xy[:,0], xy[:,1], s=3, c=vals, cmap=cmap or "coolwarm", linewidths=0)
        plt.colorbar(sca, ax=ax, fraction=0.046, pad=0.02)
    ax.set_title(title); ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")

def main():
    pb = sc.read_h5ad(f"{BASE}/trajectory_pseudobulk.h5ad")
    bdata = sc.read_h5ad(f"{BASE}/trajectory_percell.h5ad")

    # ---- Fig 5b: neighborhood V(D)J feature space ----
    ct_col = "celltype_annotation" if "celltype_annotation" in pb.obs else "anno"
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
    if "X_umap" not in pb.obsm and "X_pca" in pb.obsm:
        sc.pp.neighbors(pb, use_rep="X_pca"); sc.tl.umap(pb)
    umap_panel(axes[0], pb, ct_col, "Fig 5b: V(D)J feature space (cell type)", categorical=True)
    axes[0].legend(loc="center left", bbox_to_anchor=(1.0,0.5), frameon=False, markerscale=2, ncol=1)
    if "pseudotime" in pb.obs:
        umap_panel(axes[1], pb, "pseudotime", "pseudotime", cmap="viridis")
    prob_cols = [c for c in pb.obs.columns if c.startswith("prob_")]
    if prob_cols:
        umap_panel(axes[2], pb, prob_cols[0], prob_cols[0], cmap="coolwarm")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig5b_vdj_feature_umap.png", bbox_inches="tight")
    print("wrote fig5b_vdj_feature_umap.png")

    # ---- Fig 5c: per-cell branch probabilities ----
    pcols = [c for c in bdata.obs.columns if c.startswith("prob_")]
    n = 2 + len(pcols)
    fig, axes = plt.subplots(1, n, figsize=(3.6*n, 3.4))
    umap_panel(axes[0], bdata, "celltype_annotation", "Fig 5c: cell type", categorical=True)
    axes[0].legend(loc="center left", bbox_to_anchor=(1.0,0.5), frameon=False, markerscale=2, fontsize=5)
    if "pseudotime" in bdata.obs:
        umap_panel(axes[1], bdata, "pseudotime", "pseudotime", cmap="viridis")
    for k,pc in enumerate(pcols):
        umap_panel(axes[2+k], bdata, pc, pc, cmap="coolwarm")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig5c_branch_probabilities.png", bbox_inches="tight")
    print("wrote fig5c_branch_probabilities.png")

if __name__ == "__main__":
    main()
