#!/usr/bin/env python
"""
Step 6 (the novel question): is the Dandelion NK/ILC-from-DN branch TCR-DRIVEN or
TRANSCRIPTOME-DRIVEN?

The paper builds the trajectory in a V(D)J feature space (TRBJ usage per neighborhood)
and concludes ILC/NK cells sit on an aborted T-development trajectory. But the branch
could be an artifact of transcriptional similarity alone (NK/ILC and DN thymocytes share
progenitor programs) rather than being carried by the TCR-rearrangement footprint.

Three dissection tests:

  (A) GEX-ONLY control: rebuild the trajectory on the scVI (transcriptome) feature space
      with NO TRBJ. If the NK/ILC branch survives unchanged, the branch is transcriptome-
      driven and the TCR feature space adds nothing. If it collapses/degrades, TCR carries it.

  (B) TRBJ ABLATION / PERMUTATION: shuffle the TRBJ-usage vectors across neighborhoods
      (destroying the real TCR structure while preserving marginal composition), re-run
      Palantir, and measure how much the T-vs-NK/ILC branch probabilities degrade
      (correlation of per-cell branch prob real vs permuted; n permutations -> null).

  (C) PRODUCTIVE-only vs +NON-PRODUCTIVE TRBJ: rebuild the TRBJ feature space using
      (i) only productive TRBJ, (ii) productive + non-productive. If the aborted
      (non-productive) contigs specifically sharpen the NK/ILC branch, that is direct
      evidence the *aborted* rearrangement — not just any TCR signal — carries the branch.

Inputs : processed/dandelion_reproduction/trajectory_percell.h5ad  (from 03_build_trajectory)
         processed/dandelion_reproduction/trajectory_pseudobulk.h5ad
Outputs: driver_dissection_summary.csv + figures
Env: ddl-traj (run with NUMBA_CACHE_DIR set).
"""
import warnings; warnings.filterwarnings("ignore")
import os, numpy as np, pandas as pd, scanpy as sc
import dandelion as ddl, palantir

BASE = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/dandelion_reproduction"
CKPT = f"{BASE}/trajectory_adata_scvi_nhoods.h5ad"   # scVI + nhoods (from 03_build_trajectory)
CT_COL = "celltype_annotation"
CD34_ENS = "ENSG00000174059"
N_PERM = 200
T_TYPES = {"DN(early)_T","DN(P)_T","DN(Q)_T","DP(P)_T","DP(Q)_T","ABT(ENTRY)"}
NKILC_TYPES = {"NK","CYCLING_NK","ILC2","ILC3","CYCLING_ILC"}
def lineage_of(ct):
    if ct in T_TYPES: return "T"
    if ct in NKILC_TYPES: return "NK_ILC"
    return "other"

# ----------------------------------------------------------------------------
# Shared machinery. root + terminal neighborhoods are properties of the CELL graph
# (CD34 for the root, lineage composition for the two tips) and are held FIXED across
# all feature spaces so that the only thing that changes between tests is the feature
# matrix the pseudotime is computed on.
# ----------------------------------------------------------------------------
def branch_from_pb(pbX, pb_names, root_idx, term_idx, term_labels):
    """Run Palantir on an arbitrary neighborhood feature matrix pbX (n_nhoods x n_feat).
    Returns a DataFrame of branch probabilities (columns = term_labels), indexed by nhood."""
    import anndata as ad
    pb = ad.AnnData(np.asarray(pbX, dtype=float))
    pb.obs_names = pb_names
    sc.pp.pca(pb, n_comps=min(20, pb.shape[1]-1))
    terminal_states = pd.Series(term_labels, index=pb.obs_names[term_idx])
    proj = pd.DataFrame(pb.obsm["X_pca"], index=pb.obs_names)
    dm = palantir.utils.run_diffusion_maps(proj, n_components=5)
    ms = palantir.utils.determine_multiscale_space(dm); ms.index = ms.index.astype(str)
    pr = palantir.core.run_palantir(ms, pb.obs_names[root_idx], num_waypoints=500,
                                    terminal_states=terminal_states.index)
    bp = pr.branch_probs.copy()
    bp.columns = [terminal_states[c] for c in bp.columns]
    return bp

def project_to_cells(bp_nhood, nhoods):
    """Project neighborhood branch probs to cells: mean over the nhoods each cell belongs to."""
    # nhoods: cells x n_nhoods binary; bp_nhood aligned to nhood columns
    W = nhoods                                   # scipy sparse cells x nhoods
    denom = np.asarray(W.sum(1)).ravel() + 1e-9
    out = {}
    for col in bp_nhood.columns:
        v = bp_nhood[col].values
        out[col] = np.asarray(W @ v).ravel() / denom
    return pd.DataFrame(out)

def separation(cellprob_nk, lineage):
    """Branch-separation metric: mean NK/ILC-branch prob in NK/ILC cells minus in T cells.
    Positive + large => NK/ILC cells are pushed onto their own branch tip."""
    m_nk = lineage.values == "NK_ILC"; m_t = lineage.values == "T"
    return float(np.nanmean(cellprob_nk[m_nk]) - np.nanmean(cellprob_nk[m_t]))

def main():
    adata = sc.read_h5ad(CKPT)
    adata.obs["lineage"] = adata.obs[CT_COL].map(lineage_of).astype("category")
    nh = adata.obsm["nhoods"]; nh_size = np.asarray(nh.sum(0)).ravel() + 1e-9
    n_nhood = nh.shape[1]
    pb_names = pd.Index([f"nh{i}" for i in range(n_nhood)])

    # fixed root (highest-CD34 nhood) and terminal tips (max frac_T, max frac_NK/ILC)
    col = adata[:, CD34_ENS].X
    cd34 = np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()
    root_idx = int(np.argmax((nh.T @ cd34) / nh_size))
    is_T = (adata.obs["lineage"].values == "T").astype(float)
    is_NK = (adata.obs["lineage"].values == "NK_ILC").astype(float)
    frac_T = (nh.T @ is_T) / nh_size; frac_NK = (nh.T @ is_NK) / nh_size
    t_idx = int(np.argmax(frac_T)); nk_idx = int(np.argmax(frac_NK))
    if t_idx == nk_idx: nk_idx = int(np.argsort(frac_NK)[-2])
    term_idx = [t_idx, nk_idx]; term_labels = ["T", "NK_ILC"]
    print(f"root nhood={root_idx} (CD34); tips T={t_idx} NK_ILC={nk_idx}")

    # ---------- REAL TRBJ feature space ----------
    print("=== REAL: TRBJ feature space ===")
    pb_real = ddl.tl.vdj_pseudobulk(adata, pbs=nh, obs_to_take=[CT_COL, "lineage"],
                                    extract_cols=["j_call_abT_VDJ_main"])
    X_trbj = pb_real.X.toarray() if hasattr(pb_real.X, "toarray") else np.asarray(pb_real.X)
    bp_real = branch_from_pb(X_trbj, pb_names, root_idx, term_idx, term_labels)
    cell_real = project_to_cells(bp_real, nh)
    sep_real = separation(cell_real["NK_ILC"].values, adata.obs["lineage"])
    print(f"  REAL branch separation (NK/ILC prob: NK cells - T cells) = {sep_real:.4f}")

    results = {"sep_real_TRBJ": sep_real, "n_nhood": int(n_nhood),
               "n_NKILC_cells": int(is_NK.sum()), "n_T_cells": int(is_T.sum())}

    # ---------- (A) GEX-only control: feature = scVI transcriptome mean per nhood ----------
    print("=== (A) GEX-only (scVI transcriptome) control ===")
    scvi = adata.obsm["X_scvi"]
    X_gex = (nh.T @ scvi) / nh_size[:, None]        # n_nhood x 20 scVI means
    bp_gex = branch_from_pb(X_gex, pb_names, root_idx, term_idx, term_labels)
    cell_gex = project_to_cells(bp_gex, nh)
    sep_gex = separation(cell_gex["NK_ILC"].values, adata.obs["lineage"])
    results["sep_A_GEXonly"] = sep_gex
    print(f"  GEX-only branch separation = {sep_gex:.4f}")

    # ---------- (B) TRBJ permutation null ----------
    # Shuffle the TRBJ feature vectors ACROSS neighborhoods (each nhood gets a random nhood's
    # TRBJ profile), destroying the real TRBJ<->nhood correspondence while keeping the marginal
    # feature distribution. Root/tips stay fixed. If separation collapses to ~0, the real TRBJ
    # structure is what carries the branch (TCR-driven). Empirical p = P(perm sep >= real sep).
    print(f"=== (B) TRBJ permutation null (n={N_PERM}) ===")
    rng = np.random.default_rng(0); seps = []
    for p in range(N_PERM):
        Xp = X_trbj[rng.permutation(n_nhood)]
        try:
            bp_p = branch_from_pb(Xp, pb_names, root_idx, term_idx, term_labels)
            cell_p = project_to_cells(bp_p, nh)
            seps.append(separation(cell_p["NK_ILC"].values, adata.obs["lineage"]))
        except Exception:
            pass
    seps = np.array([s for s in seps if np.isfinite(s)])
    pval = float((np.sum(seps >= sep_real) + 1) / (len(seps) + 1))
    results["perm_null_mean_sep"] = float(seps.mean())
    results["perm_null_sd_sep"] = float(seps.std())
    results["perm_n"] = int(len(seps))
    results["perm_pvalue"] = pval
    print(f"  permuted separation: mean {seps.mean():.4f} sd {seps.std():.4f} (n={len(seps)})")
    print(f"  real={sep_real:.4f} -> empirical p = {pval:.4f}")

    # ---------- (C) productive-only vs +non-productive TRBJ ----------
    # Rebuild the TRBJ feature space counting ONLY productive-beta cells (null the j_call of
    # cells whose VDJ/beta contig is non-productive). If separation DROPS vs REAL (which
    # includes non-productive), the non-productive/aborted TRB specifically carries the branch.
    print("=== (C) productive-only TRBJ ===")
    ad_prod = adata.copy()
    prodmask = ad_prod.obs.get("productive_abT_VDJ", pd.Series(index=ad_prod.obs_names)).astype(str)
    is_nonprod = ~prodmask.str.upper().str.startswith("T")
    ad_prod.obs.loc[is_nonprod, "j_call_abT_VDJ_main"] = np.nan
    print(f"  cells nulled (non-productive beta): {int(is_nonprod.sum())} of {ad_prod.n_obs}")
    pb_prod = ddl.tl.vdj_pseudobulk(ad_prod, pbs=nh, obs_to_take=[CT_COL, "lineage"],
                                    extract_cols=["j_call_abT_VDJ_main"])
    X_prod = pb_prod.X.toarray() if hasattr(pb_prod.X, "toarray") else np.asarray(pb_prod.X)
    bp_prod = branch_from_pb(X_prod, pb_names, root_idx, term_idx, term_labels)
    cell_prod = project_to_cells(bp_prod, nh)
    sep_prod = separation(cell_prod["NK_ILC"].values, adata.obs["lineage"])
    results["sep_C_productive_only"] = sep_prod
    print(f"  productive-only branch separation = {sep_prod:.4f} (vs real {sep_real:.4f})")

    df = pd.DataFrame([results])
    df.to_csv(f"{BASE}/driver_dissection_summary.csv", index=False)
    # also save the null distribution for the figure
    np.save(f"{BASE}/driver_perm_null.npy", seps)
    print("WROTE driver_dissection_summary.csv +  driver_perm_null.npy")
    print(df.T)

if __name__ == "__main__":
    main()
