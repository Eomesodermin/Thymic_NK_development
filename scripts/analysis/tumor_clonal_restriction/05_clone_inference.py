"""Step 4 — clone inference per donor across sites (Plan 3).

For each donor: subset to NK cells (base is_nk), clean the heteroplasmy matrix (outer-join concat
leaves NaN where a variant was absent from a site's cells -> 0), select variants informative
WITHIN that donor's NK cells (so clone structure is defined by NK-carried variants), binarize,
then call clones jointly across all of the donor's sites (tumor+normal+blood+met share one clonal
tree — this is what makes cross-site overlap measurable). Donors are held separate by construction
(one object per donor). Reports clone-size distribution per site + clone QC.

Also emits a NK-not-ILC (strict) clone call in parallel for the Step-5/6 sensitivity check.

Usage:
    python clone_inference.py --nk-dir <nk_labelled> --out <dir>
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

sys.path.insert(0, "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts")
from mtclone import qc, clones


def clean_matrix(a):
    """NaN->0, ensure CSR sparse float32 (concat left X dense with NaN)."""
    X = a.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.nan_to_num(np.asarray(X, dtype=np.float32), nan=0.0)
    a.X = sp.csr_matrix(X)
    return a


def call_for_mask(a_donor, mask, *, tag, min_shared=1, edge=0.5, min_clone=2, resolution=1.0):
    """Select informative variants within the masked (NK) subset, call clones, return sub + qc."""
    sub = a_donor[mask].copy()
    if sub.n_obs < 3:
        return None, dict(tag=tag, n_cells=int(sub.n_obs), skipped="too_few_cells")
    sub = clean_matrix(sub)
    # recompute pseudobulk / detection on the subset so variant selection is NK-specific
    Xd = sub.X.toarray()
    sub.var["pseudobulk_heteroplasmy"] = Xd.mean(0)
    sub.var["n_cells_detected"] = (Xd >= 0.01).sum(0)
    try:
        sel = qc.select_informative_variants(
            sub, min_cell_coverage=10, max_pseudobulk_het=0.90,
            min_cells_detected=5, binarize_threshold=0.07)
    except Exception as e:
        return None, dict(tag=tag, n_cells=int(sub.n_obs), error=str(e))
    if sel.n_vars == 0:
        return None, dict(tag=tag, n_cells=int(sub.n_obs), n_vars=0, skipped="no_informative_variants")
    clones.call_clones(sel, method="graph", min_shared_variants=min_shared,
                        edge_weight_cutoff=edge, min_clone_size=min_clone, resolution=resolution)
    q = clones.clone_qc(sel)
    q.update(tag=tag, n_vars=int(sel.n_vars))
    return sel, q


def clone_size_by_site(sel):
    """Per-site clone-size distribution for assigned cells."""
    df = sel.obs[sel.obs["clone_id"] >= 0]
    rows = []
    for site in sorted(df["site"].unique()):
        d = df[df["site"] == site]
        sizes = d["clone_id"].value_counts()
        rows.append(dict(site=site, n_cells=int(len(d)), n_clones=int(d["clone_id"].nunique()),
                         max_clone=int(sizes.max()) if len(sizes) else 0,
                         top_clone_frac=float(sizes.max() / len(d)) if len(d) else 0.0))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nk-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    qc_rows, size_rows = [], []
    for h in sorted(Path(args.nk_dir).glob("*.nk_labelled.h5ad")):
        donor = h.stem.replace(".nk_labelled", "")
        a = ad.read_h5ad(h)
        for tag, key in [("nk", "is_nk"), ("nk_not_ilc", "nk_not_ilc")]:
            mask = a.obs[key].values.astype(bool)
            sel, q = call_for_mask(a, mask, tag=tag)
            q["donor"] = donor
            qc_rows.append(q)
            if sel is not None:
                sel.write_h5ad(out / f"{donor}.{tag}.clones.h5ad")
                sd = clone_size_by_site(sel); sd["donor"] = donor; sd["tag"] = tag
                size_rows.append(sd)
                print(f"{donor} [{tag}]: {q['n_clones']} clones, "
                      f"{q['frac_assigned']:.2f} assigned, {q['n_vars']} vars", flush=True)
            else:
                print(f"{donor} [{tag}]: {q.get('skipped') or q.get('error')}", flush=True)

    qdf = pd.DataFrame(qc_rows)
    qdf.to_csv(out / "clone_qc.csv", index=False)
    if size_rows:
        pd.concat(size_rows, ignore_index=True).to_csv(out / "clone_size_by_site.csv", index=False)
    print("\n=== clone QC (base NK) ===")
    print(qdf[qdf["tag"] == "nk"].to_string(index=False))


if __name__ == "__main__":
    main()
