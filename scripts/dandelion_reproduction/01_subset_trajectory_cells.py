#!/usr/bin/env python
"""
Step 2 (part 1): subset the Suo et al. fetal lymphoid atlas to the trajectory recipe
cell types, keeping raw counts for scVI. Reproduces the cell selection of
Suo/Dandelion Fig 5 (before the TRBJ-with-VDJ filter, which is applied after VDJ join).

Recipe cell types (Dandelion Methods, Fig 5b/c):
  DN(early)_T, DN(P)_T, DN(Q)_T, DP(P)_T, DP(Q)_T, ILC2, ILC3, CYCLING_ILC, NK, CYCLING_NK
ABT(ENTRY) is also carried as the branch-point compartment the paper emphasises.

Input : raw/suo_fetal/lymphoid.h5ad      (CELLxGENE Suo collection b1a879f6..., DOI 10.1126/science.abo0510)
Output: processed/dandelion_reproduction/trajectory_gex_subset.h5ad
"""
import warnings; warnings.filterwarnings("ignore")
import scanpy as sc, numpy as np, anndata as ad, os

RAW = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/raw/suo_fetal/lymphoid.h5ad"
OUT_DIR = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/dandelion_reproduction"
os.makedirs(OUT_DIR, exist_ok=True)

RECIPE = ["DN(early)_T","DN(P)_T","DN(Q)_T","DP(P)_T","DP(Q)_T",
          "ILC2","ILC3","CYCLING_ILC","NK","CYCLING_NK"]
BRANCHPOINT = ["ABT(ENTRY)"]

adata = sc.read_h5ad(RAW)
print("full lymphoid:", adata.shape)
# CELLxGENE stores raw counts in .raw or layers; X may be normalised. Check.
print("X dtype/max:", adata.X.dtype, adata.X[:100].max())
keep = adata.obs["celltype_annotation"].isin(RECIPE + BRANCHPOINT)
sub = adata[keep].copy()
print("subset:", sub.shape)
print(sub.obs["celltype_annotation"].value_counts())

# ensure raw counts available for scVI: CELLxGENE convention = adata.raw.X is raw counts
if sub.raw is not None:
    print("raw present:", sub.raw.shape, "raw max:", sub.raw.X[:100].max())

sub.write_h5ad(os.path.join(OUT_DIR, "trajectory_gex_subset.h5ad"))
print("wrote trajectory_gex_subset.h5ad")
