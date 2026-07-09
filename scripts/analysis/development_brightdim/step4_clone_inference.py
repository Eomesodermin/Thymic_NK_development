#!/usr/bin/env python
"""step4_clone_inference.py — Plan 2 Step 4: NK clone inference (donor-private).

Clones are DONOR-PRIVATE (a somatic mtDNA variant arose in one individual), so inference is
ALWAYS per donor — never pool donors. For each donor:
  1. subset that donor's NK cells (Step-3 labelled objects carry cell_type + nk_subset);
  2. re-select variants informative WITHIN this donor's NK (Step-3 QC pooled NK across a
     Rückert sample's 2-4 donors; here we tighten to the single donor so clone edges rest on
     variants that actually vary among this donor's NK);
  3. mtclone.clones.call_clones (graph/Leiden: min_shared_variants=1, jaccard>=0.5,
     min_clone_size=2) — same settings as Plan 3.

Reports clone-size distribution, singleton/unassigned fraction, and — the number that decides
whether Step 5 can run — how many multi-cell clones contain >=2 CLASSIFIABLE (bright|dim) NK.

Outputs (processed/plan2_step4/):
  plan2_step4_clone_qc.csv        — per donor: n NK, n clones, frac assigned, max/median clone,
                                     n multi-cell clones, n informative-for-sharing clones
  <object>.<donor>.nkclones.h5ad  — per-donor NK object with obs['clone_id'] (for Step 5)
"""
import sys, os, glob, warnings
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp

MTCLONE = "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts"
sys.path.insert(0, MTCLONE)
from mtclone import qc, clones

STEP3 = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/plan2_step3"
OUT = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/plan2_step4"
os.makedirs(OUT, exist_ok=True)

QC_KW = dict(max_pseudobulk_het=0.90, min_cells_detected=5,
             min_cell_coverage=5.0, binarize_threshold=0.05)

# Clone definition is matched to the data type (confirmed with user, Step 4):
#   * ReDeeM  = UMI-consensus mtDNA -> clean multi-variant lineages -> graph/Leiden
#     (min_shared_variants=1, jaccard>=0.5), same as Plan 3.
#   * Rückert = raw mtASAP-seq ~13x -> too many false variants for confident pairwise
#     Jaccard clustering (<1% of cells get an edge at jaccard>=0.5). Use the standard
#     mtscATAC clone definition (Ludwig/Lareau): a clone = cells sharing a specific rare
#     somatic variant -> mtclone `variant_group`.
CLONE_KW = {
    "graph":         dict(method="graph", min_shared_variants=1, edge_weight_cutoff=0.5, min_clone_size=2),
    "variant_group": dict(method="variant_group", min_clone_size=2),
}
# object prefix -> arm/method. ReDeeM objects: Young*/Old*; Rückert: mtASAP*.
def _method_for(obj):
    return "variant_group" if obj.startswith("mtASAP") else "graph"


def infer_one_donor(nk, obj, donor):
    """nk: AnnData of one donor's NK cells (heteroplasmy in X). Returns (qc_dict, labelled nk)."""
    # clean X: densify absent entries as 0 (NaN-safe), keep as sparse
    if sp.issparse(nk.X):
        nk.X = sp.csr_matrix(np.nan_to_num(nk.X.toarray()))
    else:
        nk.X = sp.csr_matrix(np.nan_to_num(np.asarray(nk.X)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        f = qc.select_informative_variants(nk, **QC_KW)
    if f.n_vars == 0 or f.n_obs < 2:
        return None, None
    method = _method_for(obj)
    clones.call_clones(f, **CLONE_KW[method])
    q = clones.clone_qc(f)
    q["clone_method"] = method
    # sharing-informative clones: multi-cell clones with >=2 classifiable (bright|dim) NK
    lab = f.obs["clone_id"].values
    sub = f.obs["nk_subset"].values
    classif = np.isin(sub, ["bright", "dim"])
    n_multi = n_share = 0
    for cid in np.unique(lab[lab >= 0]):
        m = lab == cid
        if m.sum() >= 2:
            n_multi += 1
            if classif[m].sum() >= 2:
                n_share += 1
    q.update(object=obj, donor=donor, n_informative_variants=int(f.n_vars),
             n_multicell_clones=int(n_multi), n_sharing_informative_clones=int(n_share),
             n_bright=int((sub == "bright").sum()), n_dim=int((sub == "dim").sum()))
    return q, f


def run():
    rows = []
    for fp in sorted(glob.glob(f"{STEP3}/*.labelled.h5ad")):
        obj = os.path.basename(fp).replace(".labelled.h5ad", "")
        a = ad.read_h5ad(fp)
        nk_all = a[a.obs["cell_type"] == "NK"].copy()
        if nk_all.n_obs == 0:
            continue
        for donor in nk_all.obs["donor"].astype(str).unique():
            nk = nk_all[nk_all.obs["donor"].astype(str) == donor].copy()
            q, f = infer_one_donor(nk, obj, donor)
            if q is None:
                rows.append(dict(object=obj, donor=donor, n_cells=int(nk.n_obs),
                                 n_clones=0, frac_assigned=0.0, note="too few informative variants/cells"))
                continue
            rows.append(q)
            f.write_h5ad(f"{OUT}/{obj}.{donor}.nkclones.h5ad")
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/plan2_step4_clone_qc.csv", index=False)
    return df


if __name__ == "__main__":
    df = run()
    cols = ["object","donor","n_cells","n_clones","frac_assigned","max_clone_size",
            "median_clone_size","n_multicell_clones","n_sharing_informative_clones",
            "n_bright","n_dim","n_informative_variants"]
    cols = [c for c in cols if c in df.columns]
    print(df[cols].to_string(index=False))
