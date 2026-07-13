#!/usr/bin/env python
"""
Steps 2-5: assemble the Suo trajectory object from Dandelion-processed VDJ + GEX,
then reproduce the Dandelion Fig 5b/c NK/ILC-from-DN pseudotime + branch probabilities.

Pipeline (exact API from zktuong/dandelion docs/notebooks/8-pseudobulk-trajectory.ipynb
and src/dandelion/tools/_trajectory.py; paper Methods 10.1038/s41587-023-01734-7 Fig 5):

  1. Read per-library Dandelion AIRR TSVs (all_contig_dandelion.tsv, NON-productive retained)
     -> add per-lib cell prefix matching Suo obs index (<GEX_id>-<barcode>) -> ddl.concat
  2. ddl.pp.check_contigs(vdj, adata, library_type="tr-ab")  -> chain_status
  3. ddl.tl.find_clones(vdj)
  4. ddl.tl.transfer(adata, vdj)  -> writes v_call_abT_VDJ / productive_abT_VDJ / chain_status obs cols
  5. ddl.tl.setup_vdj_pseudobulk(adata, mode="abT")  -> keep TRBJ-bearing productive-pair cells
  6. scVI on raw counts (batch=donor) -> X_scvi   [CZI object lacks stored latent]
  7. sc.pp.neighbors(use_rep=X_scvi, k=50) -> Milo make_nhoods -> umap
  8. ddl.tl.vdj_pseudobulk(pbs=nhoods, obs_to_take=celltype) -> V(D)J feature space (TRBJ usage)
  9. Palantir: diffusion maps (5 PCs), root=highest-CD34 nhood, terminals=UMAP1 extremes (T vs NK/ILC)
 10. pseudotime_transfer -> project_pseudotime_to_cell -> Fig 5b (nhood UMAP) + Fig 5c (per-cell)

Env: ddl-traj. Run with NUMBA_CACHE_DIR set (see run wrapper).
"""
import warnings; warnings.filterwarnings("ignore")
import os, glob, sys
import scanpy as sc, numpy as np, pandas as pd, anndata as ad
import dandelion as ddl
import palantir

BASE = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development"
GEX_SUBSET = f"{BASE}/processed/dandelion_reproduction/trajectory_gex_subset.h5ad"
VDJ_DIR    = f"{BASE}/processed/dandelion_reproduction/dandelion_tsv"   # per-lib all_contig_dandelion.tsv
MANIFEST   = "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/data/dandelion_reproduction/suo_abTCR_vdj_manifest.json"
OUT_DIR    = f"{BASE}/processed/dandelion_reproduction"
CT_COL     = "celltype_annotation"

# ---- lineage grouping: T-lineage vs NK/ILC-lineage for terminal states ----
T_TYPES   = {"DN(early)_T","DN(P)_T","DN(Q)_T","DP(P)_T","DP(Q)_T","ABT(ENTRY)"}
NKILC_TYPES = {"NK","CYCLING_NK","ILC2","ILC3","CYCLING_ILC"}

def lineage_of(ct):
    if ct in T_TYPES: return "T"
    if ct in NKILC_TYPES: return "NK_ILC"
    return "other"

def make_nhoods(adata, prop=0.1, seed=0):
    """Milo neighborhood sampling (verbatim from pertpy.tools._milo.Milo.make_nhoods, sans R).
    Samples vertices on the kNN connectivity graph, refines each to the cell nearest the
    neighborhood median in the embedding, and writes a binary cells x nhoods matrix to
    adata.obsm['nhoods']. Requires sc.pp.neighbors() to have been run (uses obsp connectivities
    and uns['neighbors']['params']['use_rep'])."""
    import random
    from sklearn.metrics import euclidean_distances
    use_rep = adata.uns["neighbors"]["params"].get("use_rep", "X_pca")
    knn_graph = adata.obsp["connectivities"].copy()
    X_dimred = adata.obsm[use_rep]
    n_ixs = int(np.round(adata.n_obs * prop))
    knn_graph[knn_graph != 0] = 1
    random.seed(seed)
    rv = sorted(random.sample(range(adata.n_obs), k=n_ixs))
    ixs_nn = knn_graph[rv, :]
    nzr = ixs_nn.nonzero()[0]; nzc = ixs_nn.nonzero()[1]
    refined = np.empty(len(rv))
    for i in range(len(rv)):
        pts = X_dimred[nzc[nzr == i], :]
        nh_pos = np.median(pts, 0).reshape(-1, 1)
        nn_ixs = nzc[nzr == i]
        dists = euclidean_distances(pts, nh_pos.T)
        refined[i] = nn_ixs[dists.argmin()]
    refined = np.unique(refined.astype("int")); refined.sort()
    adata.obsm["nhoods"] = knn_graph[:, refined]
    adata.obs["nhood_ixs_refined"] = adata.obs_names.isin(adata.obs_names[refined]).astype(int)
    adata.uns["nhood_neighbors_key"] = None

def main():
    import json
    manifest = json.load(open(MANIFEST))
    lane2gex = {m["lane"]: m["gex"] for m in manifest}

    # ---- 1. read + concat VDJ ----
    # Suo obs index = <GEX_id>-<barcode> (no trailing -1). Dandelion barcodes are <barcode>-1.
    # read_10x_airr(prefix=<GEX_id>, remove_trailing_hyphen_number=True) -> cell_id = <GEX_id>_<barcode>
    # NOTE: Suo uses "-" between gex and barcode; dandelion default sep is "_". Set sep="-".
    print("=== reading Dandelion VDJ TSVs ===")
    vdj_list = []
    for tsv in sorted(glob.glob(f"{VDJ_DIR}/*_dandelion.tsv")):
        base = os.path.basename(tsv)
        lane = None
        for L in lane2gex:
            if L in base:
                lane = L; break
        if lane is None:
            print("  skip (no lane match):", base); continue
        gex = lane2gex[lane]
        try:
            v = ddl.read_10x_airr(tsv, prefix=gex, sep="-", remove_trailing_hyphen_number=True)
        except ValueError as e:
            # near-empty libraries (a handful of contigs) drop all-empty required columns
            # (e.g. junction_aa) on read and fail metadata init. They contribute negligible
            # cells; skip them and log.
            print(f"  SKIP {lane} -> {gex}: malformed/near-empty ({str(e)[:40]})")
            continue
        vdj_list.append(v)
        print(f"  {lane} -> {gex}: {v.data.shape[0]} contigs")
    print(f"libraries read OK: {len(vdj_list)}")
    vdj = ddl.concat(vdj_list)
    print("concat VDJ:", vdj.data.shape)

    # ---- load GEX subset ----
    adata = sc.read_h5ad(GEX_SUBSET)
    print("GEX subset:", adata.shape)
    # verify barcode overlap before proceeding
    vdj_bc = set(vdj.data["cell_id"].unique())
    gex_bc = set(adata.obs_names)
    ov = len(vdj_bc & gex_bc)
    print(f"barcode overlap: VDJ cells {len(vdj_bc)} | GEX cells {len(gex_bc)} | overlap {ov}")
    assert ov > 0, "no barcode overlap - check prefix/sep harmonisation"

    # ---- 2-4. check_contigs, find_clones, transfer ----
    # CRITICAL (paper Methods, "nonproductive recombination as a fossil record"): the NK/ILC
    # branch rests on NON-PRODUCTIVE TRB contigs (most without a V gene). check_contigs defaults
    # to productive_only=True, which would DISCARD exactly those cells. Set productive_only=False
    # so the non-productive TRB fossils are retained for the TRBJ feature space.
    print("=== check_contigs (tr-ab, productive_only=False) ===")
    vdj, adata = ddl.pp.check_contigs(vdj, adata, library_type="tr-ab", productive_only=False)
    ddl.tl.find_clones(vdj)
    ddl.tl.transfer(adata, vdj)
    print("obs cols after transfer:", [c for c in adata.obs.columns if "abT" in c or c=="chain_status"])
    print("chain_status:", {k:int(v) for k,v in adata.obs["chain_status"].value_counts().items()})

    # ---- 5. setup_vdj_pseudobulk ----
    # Paper: "We used TRBJ frequency to construct a V(D)J feature space because all T/ILC/NK
    # cells express TRBJ." The NK/ILC branch is carried by NON-PRODUCTIVE TRB (fossils, mostly
    # V-gene-less). So the feature space must be agnostic to productivity: keep both
    # productive_vdj=False and productive_vj=False. allowed_chain_status still restricts to
    # cells with a real (VDJ/beta) rearrangement rather than No_contig.
    print("=== setup_vdj_pseudobulk(mode=abT, productive_vdj=False, productive_vj=False) ===")
    adata_full = adata.copy()   # keep the full matched object for the GEX-only control graph
    adata = ddl.tl.setup_vdj_pseudobulk(adata, mode="abT",
                                        productive_vdj=False, productive_vj=False)
    print("cells with TRBJ-bearing abTCR (productive OR non-productive beta):", adata.shape)
    print("  celltype breakdown:", {k:int(v) for k,v in adata.obs[CT_COL].value_counts().items()})

    # ---- 6. scVI ----
    # GEX subset .X is log-normalised; scVI needs raw counts, which live in .raw (integer
    # counts, verified). Rebuild a counts-layer AnnData restricted to the retained cells.
    print("=== scVI ===")
    import scvi
    raw = adata.raw.to_adata()                    # integer counts, same var (Ensembl IDs)
    raw = raw[adata.obs_names].copy()             # restrict to TRBJ-bearing cells kept above
    raw.obs = adata.obs.copy()
    raw.layers["counts"] = raw.X.copy()
    sc.pp.highly_variable_genes(raw, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    hv = raw[:, raw.var.highly_variable].copy()
    scvi.model.SCVI.setup_anndata(hv, layer="counts", batch_key="donor_id")
    m = scvi.model.SCVI(hv, n_latent=20)
    m.train(max_epochs=200, accelerator="cpu", early_stopping=True)
    adata.obsm["X_scvi"] = m.get_latent_representation()

    # ---- 7. neighbors + Milo + umap ----
    # pertpy's pt.tl.Milo() imports rpy2/R at construction (for the edgeR differential-abundance
    # GLM), which we don't need. make_nhoods itself is pure Python (sample vertices on the kNN
    # graph, refine to the median-profile nearest cell, build a binary cells x nhoods matrix).
    # Reimplemented verbatim from pertpy.tools._milo.Milo.make_nhoods to avoid the R dependency.
    print("=== neighbors + make_nhoods ===")
    sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=50)
    make_nhoods(adata, prop=0.1, seed=0)
    print(f"  neighborhoods: {adata.obsm['nhoods'].shape[1]}")
    sc.tl.umap(adata)
    # checkpoint: everything above (read TSVs + check_contigs + scVI ~25 min) is expensive;
    # save so downstream Palantir iterations can resume cheaply.
    os.makedirs(OUT_DIR, exist_ok=True)
    adata.write_h5ad(f"{OUT_DIR}/trajectory_adata_scvi_nhoods.h5ad")
    print("CHECKPOINT: trajectory_adata_scvi_nhoods.h5ad")

    # ---- 8. vdj_pseudobulk (V(D)J feature space) ----
    print("=== vdj_pseudobulk ===")
    adata.obs["lineage"] = adata.obs[CT_COL].map(lineage_of).astype("category")
    # Paper Fig 5b feature space = TRBJ usage (J gene of the beta/VDJ chain).
    # extract_cols default is [v/j_call_abT_VDJ_main, v/j_call_abT_VJ_main] (full alpha+beta feature
    # space). We compute BOTH: the full space (matches the notebook demo) and a TRBJ-only variant.
    pb = ddl.tl.vdj_pseudobulk(adata, pbs=adata.obsm["nhoods"], obs_to_take=[CT_COL,"lineage"])
    sc.tl.pca(pb)
    # TRBJ-only feature space (paper Fig 5b): restrict to the VDJ J-gene call
    pb_trbj = ddl.tl.vdj_pseudobulk(adata, pbs=adata.obsm["nhoods"], obs_to_take=[CT_COL,"lineage"],
                                    extract_cols=["j_call_abT_VDJ_main"])
    sc.tl.pca(pb_trbj)

    # ---- 9. Palantir ----
    print("=== Palantir ===")
    # root = neighborhood with highest mean CD34 (progenitor). var_names are Ensembl IDs;
    # CD34 = ENSG00000174059. Fall back to feature_name lookup, then PC0 argmax.
    CD34_ENS = "ENSG00000174059"
    cd34_id = None
    if CD34_ENS in adata.var_names:
        cd34_id = CD34_ENS
    elif "feature_name" in adata.var.columns:
        hit = adata.var_names[adata.var["feature_name"].astype(str) == "CD34"]
        cd34_id = hit[0] if len(hit) else None
    if cd34_id is not None:
        col = adata[:, cd34_id].X
        cd34 = np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()
        nhood_cd34 = adata.obsm["nhoods"].T @ cd34 / (np.asarray(adata.obsm["nhoods"].sum(0)).ravel() + 1e-9)
        rootcell = int(np.argmax(nhood_cd34))
        print(f"  root = highest-CD34 nhood (idx {rootcell})")
    else:
        rootcell = int(np.argmax(pb.obsm["X_pca"][:,0]))
        print(f"  CD34 not found -> root = PC0 argmax (idx {rootcell})")
    # terminal states = the two branch tips we want to contrast (Fig 5c): the mature-T tip and
    # the NK/ILC tip. Score each neighborhood by its lineage composition (fraction of member
    # cells that are T vs NK/ILC, from the binary nhoods matrix) and take the argmax nhood for
    # each. This yields UNIQUE terminal labels ("T", "NK_ILC") whatever the PC geometry.
    nh = adata.obsm["nhoods"]                                   # cells x nhoods (binary)
    nh_size = np.asarray(nh.sum(0)).ravel() + 1e-9
    is_T = (adata.obs["lineage"].values == "T").astype(float)
    is_NKILC = (adata.obs["lineage"].values == "NK_ILC").astype(float)
    frac_T = (nh.T @ is_T) / nh_size
    frac_NK = (nh.T @ is_NKILC) / nh_size
    # pb rows are the neighborhoods in the same column order as nhoods; align by position
    t_idx = int(np.argmax(frac_T)); nk_idx = int(np.argmax(frac_NK))
    if t_idx == nk_idx:                                         # degenerate; fall back to 2nd-best NK
        nk_idx = int(np.argsort(frac_NK)[-2])
    term_idx = [t_idx, nk_idx]
    terminal_states = pd.Series(["T", "NK_ILC"], index=pb.obs_names[term_idx])
    print(f"  terminal T nhood frac_T={frac_T[t_idx]:.2f}; NK_ILC nhood frac_NK={frac_NK[nk_idx]:.2f}")
    pca_proj = pd.DataFrame(pb.obsm["X_pca"], index=pb.obs_names)
    dm_res = palantir.utils.run_diffusion_maps(pca_proj, n_components=5)
    ms = palantir.utils.determine_multiscale_space(dm_res); ms.index = ms.index.astype(str)
    pr = palantir.core.run_palantir(ms, pb.obs_names[rootcell], num_waypoints=500,
                                    terminal_states=terminal_states.index)
    # relabel branch-prob columns from cell-id -> lineage label (unique, so no collision)
    pr.branch_probs.columns = [terminal_states[c] for c in pr.branch_probs.columns]

    # ---- 10. transfer + project ----
    pb = ddl.tl.pseudotime_transfer(pb, pr)
    bdata = ddl.tl.project_pseudotime_to_cell(adata, pb, terminal_states.values)

    os.makedirs(OUT_DIR, exist_ok=True)
    pb.write_h5ad(f"{OUT_DIR}/trajectory_pseudobulk.h5ad")
    bdata.write_h5ad(f"{OUT_DIR}/trajectory_percell.h5ad")
    print("WROTE trajectory_pseudobulk.h5ad + trajectory_percell.h5ad")
    print("terminal states:", dict(terminal_states))
    print("bdata obs:", [c for c in bdata.obs.columns if "prob" in c or c=="pseudotime"])

    # ==================================================================
    # DRIVER DISSECTION (A): GEX-only control trajectory.
    # Build a pseudobulk whose FEATURE SPACE is the scVI transcriptome mean per
    # neighborhood (NOT TRBJ usage). Run the identical Palantir recipe. If the
    # NK/ILC branch survives here as well, the branch is transcriptome-driven and
    # the TRBJ feature space is not what carries it. Save for 04_driver_dissection.
    # ==================================================================
    print("=== DRIVER (A): GEX-only (scVI) pseudobulk control ===")
    nh = adata.obsm["nhoods"]                          # (cells x nhoods) sparse
    w  = np.asarray(nh.sum(0)).ravel() + 1e-9
    scvi_mean = (nh.T @ adata.obsm["X_scvi"]) / w[:, None]   # (nhoods x 20) mean scVI per nhood
    pb_gex = pb.copy()
    pb_gex.obsm["X_scvi_mean"] = scvi_mean
    from sklearn.decomposition import PCA
    npc = min(10, scvi_mean.shape[1])
    pb_gex.obsm["X_pca"] = PCA(n_components=npc).fit_transform(scvi_mean)
    pc1g = pb_gex.obsm["X_pca"][:,1]
    term_idx_g = [int(np.argmax(pc1g)), int(np.argmin(pc1g))]
    term_labels_g = [pb_gex.obs["lineage"].iloc[i] for i in term_idx_g]
    terminal_states_g = pd.Series(term_labels_g, index=pb_gex.obs_names[term_idx_g])
    proj_g = pd.DataFrame(pb_gex.obsm["X_pca"], index=pb_gex.obs_names)
    dm_g = palantir.utils.run_diffusion_maps(proj_g, n_components=5)
    ms_g = palantir.utils.determine_multiscale_space(dm_g); ms_g.index = ms_g.index.astype(str)
    pr_g = palantir.core.run_palantir(ms_g, pb_gex.obs_names[rootcell], num_waypoints=500,
                                      terminal_states=terminal_states_g.index)
    pr_g.branch_probs.columns = terminal_states_g[pr_g.branch_probs.columns]
    pb_gex = ddl.tl.pseudotime_transfer(pb_gex, pr_g)
    pb_gex.write_h5ad(f"{OUT_DIR}/trajectory_pseudobulk_gexonly.h5ad")
    print("WROTE trajectory_pseudobulk_gexonly.h5ad (GEX-only control)")
    print("GEX-only terminal states:", dict(terminal_states_g))

if __name__ == "__main__":
    main()
