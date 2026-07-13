"""Step 7 — robustness + honest limits (Plan 3).

Re-runs the Step-4 clone call and the Step-6 seeding permutation across a grid that varies ONE
knob at a time from the Step-6 baseline (max_pseudobulk_het=0.90, binarize=0.07, edge=0.5,
min_shared=1). Reports whether each donor's seeding CALL (p_less<0.05) is stable.

Axes:
  A. max_pseudobulk_het (near-public/ambient FLOOR — the key test): 0.01, 0.10, 0.50, [0.90 base]
     Low = drop high-frequency/near-public variants (aggressive contamination clean).
  B. binarize_threshold: 0.05, [0.07 base], 0.10
  C. edge_weight_cutoff (clone-caller stringency): 0.30, [0.50 base], 0.70

Baseline row is recomputed here too (self-contained), so the table is internally comparable.

Usage:
    python robustness_step7.py --nk-dir <nk_labelled> --out <dir> [--n-perm 500]
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

sys.path.insert(0, "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts")
from mtclone import qc, clones, metrics


def clean(a):
    X = a.X
    if sp.issparse(X): X = X.toarray()
    a.X = sp.csr_matrix(np.nan_to_num(np.asarray(X, np.float32), nan=0.0))
    return a


def call_and_seed(a_nk, *, max_pb, binz, edge, min_shared, n_perm):
    """Select variants -> call clones -> seeding permutation on tumor+blood. Returns dict."""
    sub = clean(a_nk.copy())
    Xd = sub.X.toarray()
    sub.var["pseudobulk_heteroplasmy"] = Xd.mean(0)
    sub.var["n_cells_detected"] = (Xd >= 0.01).sum(0)
    try:
        sel = qc.select_informative_variants(sub, min_cell_coverage=10, max_pseudobulk_het=max_pb,
                                              min_cells_detected=5, binarize_threshold=binz)
    except Exception as e:
        return dict(error=str(e))
    if sel.n_vars == 0:
        return dict(n_vars=0, p_less=np.nan)
    clones.call_clones(sel, method="graph", min_shared_variants=min_shared,
                       edge_weight_cutoff=edge, min_clone_size=2)
    if not {"tumor", "blood"} <= set(sel.obs["site"].unique()):
        return dict(n_vars=int(sel.n_vars), p_less=np.nan)
    tb = sel[sel.obs["site"].isin(["tumor", "blood"])].copy()
    stat = lambda a: metrics.between_vs_within_sharing(a, "tumor", "blood", group_key="site",
                                                       min_clone_size=2)["frac_mixed_clones"]
    pn = metrics.permutation_null(tb, stat, group_key="site", stratify_by="donor", n=n_perm, seed=0)
    return dict(n_vars=int(sel.n_vars), n_clones=int(np.unique(sel.obs["clone_id"][sel.obs["clone_id"]>=0]).size),
                observed=pn["observed"], null_mean=pn["null_mean"], p_less=pn["p_less"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nk-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-perm", type=int, default=500)
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    BASE = dict(max_pb=0.90, binz=0.07, edge=0.50, min_shared=1)
    configs = [("baseline", BASE)]
    for v in [0.01, 0.10, 0.50]:      configs.append((f"max_pb={v}", {**BASE, "max_pb": v}))
    for v in [0.05, 0.10]:            configs.append((f"binarize={v}", {**BASE, "binz": v}))
    for v in [0.30, 0.70]:            configs.append((f"edge={v}", {**BASE, "edge": v}))

    donors = {}
    for h in sorted(Path(args.nk_dir).glob("*.nk_labelled.h5ad")):
        d = h.stem.replace(".nk_labelled", "")
        a = ad.read_h5ad(h)
        donors[d] = a[a.obs["is_nk"].values].copy()

    rows = []
    for cfg_name, cfg in configs:
        for d, a_nk in donors.items():
            r = call_and_seed(a_nk, n_perm=args.n_perm, **cfg)
            r.update(config=cfg_name, donor=d)
            rows.append(r)
        print(f"[{cfg_name}] done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "robustness_seeding_grid.csv", index=False)

    # stability summary: seeding call per donor per config
    df["seeding_sig"] = df["p_less"] < 0.05
    piv = df.pivot_table(index="donor", columns="config", values="p_less")
    piv.to_csv(out / "robustness_pvalue_matrix.csv")
    print("\n=== seeding p_less across configs ===")
    cols = ["baseline","max_pb=0.01","max_pb=0.1","max_pb=0.5","binarize=0.05","binarize=0.1","edge=0.3","edge=0.7"]
    cols = [c for c in cols if c in piv.columns]
    print(piv[cols].round(3).to_string())


if __name__ == "__main__":
    main()
