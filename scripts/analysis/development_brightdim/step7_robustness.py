#!/usr/bin/env python
"""step7_robustness.py — Plan 2 Step 7: robustness + honest limits.

Re-runs the Step-5 core statistic (pooled per-arm frac_mixed_clones vs within-donor permutation
null) under varied ingestion/caller parameters, plus a doublet/ambient floor check. The question
throughout: does the CONTINUUM conclusion (observed mixing ~ null) survive?

(1) Parameter sweep — binarization threshold in {0.03,0.05,0.10}, coverage floor in {3,5,10}x,
    each arm using its Step-4 caller (ReDeeM graph, Ruckert variant_group). Clones are re-derived
    per donor from the Step-3 labelled NK objects under each setting.
(2) Caller cross-check — run each arm under the OTHER caller.
(3) Doublet floor — mixed clones vs pure clones: variants/cell ratio (doublets ~2x variants);
    and drop-top-5%-vpc reanalysis (if mixing were doublet-driven, observed would collapse to null).

Outputs (processed/plan2_step7/): plan2_step7_sweep.csv, plan2_step7_doublet_check.csv,
plan2_step7_robustness.png
"""
import sys, os, glob, warnings
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad
warnings.simplefilter("ignore")
sys.path.insert(0, "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts")
from mtclone import qc, clones

S3 = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/plan2_step3"
S4 = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/plan2_step4"
OUT = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/plan2_step7"
os.makedirs(OUT, exist_ok=True)

def frac_mixed_fast(code, ib, K, min_size=2):
    b = np.bincount(code[ib], minlength=K); d = np.bincount(code[~ib], minlength=K)
    keep = (b + d) >= min_size
    return np.nan if keep.sum() == 0 else float(((b >= 1) & (d >= 1))[keep].mean())

def perm_null(code, ib, dcode, K, n=500, seed=0, min_size=2):
    rng = np.random.default_rng(seed); obs = frac_mixed_fast(code, ib, K, min_size)
    idx = {s: np.where(dcode == s)[0] for s in np.unique(dcode)}
    null = np.empty(n)
    for k in range(n):
        p = ib.copy()
        for _, ix in idx.items():
            p[ix] = rng.permutation(ib[ix])
        null[k] = frac_mixed_fast(code, p, K, min_size)
    null = null[~np.isnan(null)]
    return obs, float(null.mean()), float((np.sum(null <= obs) + 1) / (null.size + 1)), \
        float((np.sum(null >= obs) + 1) / (null.size + 1))

def load_objs():
    objs = {}
    for fp in sorted(glob.glob(f"{S3}/*.labelled.h5ad")):
        name = os.path.basename(fp).replace(".labelled.h5ad", "")
        a = ad.read_h5ad(fp)
        a = a[(a.obs.cell_type == "NK") & (a.obs.nk_subset.isin(["bright", "dim"]))].copy()
        if a.n_obs > 0:
            objs[name] = a
    return objs

def clones_for_unit(sub, binz, cov, method, ckw):
    b = sub.copy()
    b.X = sp.csr_matrix(np.nan_to_num(b.X.toarray() if sp.issparse(b.X) else np.asarray(b.X)))
    f = qc.select_informative_variants(b, max_pseudobulk_het=0.90, min_cells_detected=5,
                                       min_cell_coverage=cov, binarize_threshold=binz)
    if f.n_vars == 0 or f.n_obs < 2:
        return None
    clones.call_clones(f, method=method, min_clone_size=2, **ckw)
    return f.obs[["nk_subset", "clone_id", "donor"]].copy()

def arm_stat(objs, names, binz, cov, method, ckw, nperm=500):
    frames = []
    for name in names:
        a = objs[name]
        for d in a.obs.donor.astype(str).unique():
            o = clones_for_unit(a[a.obs.donor.astype(str) == d], binz, cov, method, ckw)
            if o is None:
                continue
            o = o[o.clone_id >= 0].copy()
            if o.empty:
                continue
            o["du"] = f"{name}.{d}"; o["ck"] = o["du"] + "::" + o["clone_id"].astype(str)
            frames.append(o)
    if not frames:
        return None
    oa = pd.concat(frames, ignore_index=True)
    code = pd.factorize(oa["ck"])[0]; K = int(code.max()) + 1
    return perm_null(code, (oa.nk_subset.values == "bright"), pd.factorize(oa["du"])[0], K, n=nperm)

GRAPH = dict(min_shared_variants=1, edge_weight_cutoff=0.5)

def run():
    objs = load_objs()
    red = [n for n in objs if not n.startswith("mtASAP")]
    ruck = [n for n in objs if n.startswith("mtASAP")]
    # (1)+(2) sweep + caller cross-check
    rows = []
    grid = [(0.03, 5.0), (0.05, 5.0), (0.10, 5.0), (0.05, 10.0), (0.05, 3.0)]
    for binz, cov in grid:
        for arm, names, m, ckw in [("ReDeeM", red, "graph", GRAPH),
                                    ("Ruckert", ruck, "variant_group", {})]:
            r = arm_stat(objs, names, binz, cov, m, ckw)
            if r:
                rows.append(dict(arm=arm, binarize=binz, cov_floor=cov, caller=m,
                                 observed=r[0], null_mean=r[1], p_less=r[2], p_greater=r[3]))
    sw = pd.DataFrame(rows); sw.to_csv(f"{OUT}/plan2_step7_sweep.csv", index=False)
    # (3) doublet floor
    checks = []
    for fp in sorted(glob.glob(f"{S4}/*.nkclones.h5ad")):
        name = os.path.basename(fp).replace(".nkclones.h5ad", "")
        arm = "Ruckert" if name.startswith("mtASAP") else "ReDeeM"
        a = ad.read_h5ad(fp); a = a[(a.obs.clone_id >= 0) & (a.obs.nk_subset.isin(["bright", "dim"]))].copy()
        if a.n_obs < 20:
            continue
        X = a.layers["binary"].tocsr() if "binary" in a.layers else (a.X >= 0.05)
        vpc = np.asarray((X > 0).sum(1)).ravel()
        dfc = pd.DataFrame({"cid": a.obs.clone_id.values, "br": (a.obs.nk_subset.values == "bright"),
                            "dm": (a.obs.nk_subset.values == "dim"), "vpc": vpc})
        g = dfc.groupby("cid").agg(nb=("br", "sum"), nd=("dm", "sum"), mvpc=("vpc", "mean"), sz=("cid", "size"))
        g = g[g.sz >= 2]; mixed = (g.nb >= 1) & (g.nd >= 1)
        checks.append(dict(object=name, arm=arm,
                           vpc_mixed=float(g.mvpc[mixed].mean()) if mixed.any() else np.nan,
                           vpc_pure=float(g.mvpc[~mixed].mean()) if (~mixed).any() else np.nan))
    cc = pd.DataFrame(checks)
    cc["vpc_ratio_mixed_over_pure"] = cc.vpc_mixed / cc.vpc_pure
    cc.to_csv(f"{OUT}/plan2_step7_doublet_check.csv", index=False)
    return sw, cc

if __name__ == "__main__":
    sw, cc = run()
    print("=== parameter + caller sweep ===")
    print(sw.to_string(index=False))
    print("\n=== doublet floor (variants/cell mixed vs pure) ===")
    print(cc.to_string(index=False))
