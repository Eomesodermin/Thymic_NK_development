#!/usr/bin/env python
"""step3_classify_power.py — Plan 2 Step 3: NK id + bright/dim classification + power gate.

Two arms, two classification strategies (the protein-gate is the reason Rückert is primary):
  * RÜCKERT (GSE197008, blood, sorted NK): every cell is NK by sort. Classify bright vs dim on
    SURFACE PROTEIN (ADT) — the gold standard: CD56-bright = CD56-hi / CD16-lo, CD56-dim =
    CD16-hi. We CLR-normalize the ADT panel, build the bright axis = z(CD56) - z(CD16), and
    split data-driven (2-component GMM), orienting the higher-CD56 component as bright.
  * ReDeeM (GSE219014, marrow): NK identity is the authors' STD.CellType == 'NK'. Classify
    bright vs dim on RNA markers via mtclone.classify.label_bright_dim(method='data_driven').

Then the GO/NO-GO POWER GATE (the binding constraint of the whole plan): per donor, run
mtclone.qc.select_informative_variants and count how many bright / dim NK carry >=1 informative
(clone-usable) mtDNA variant. Clones are donor-private, so QC is per donor.

Outputs (to processed/plan2_step3/):
  plan2_step3_counts.csv     — per (dataset, donor, subset): n NK, n with >=1 informative variant,
                               median informative variants/cell, median coverage
  plan2_step3_power.md       — human-readable power assessment + go/no-go per donor
  plan2_step3_<arm>_umap.png — labelled NK bright/dim embedding per arm
Also writes labelled per-donor .h5ad (nk_subset + informative-variant layer) to processed/plan2_step3/.
"""
import sys, os, glob, warnings
import numpy as np
import pandas as pd
import anndata as ad

MTCLONE = "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts"
sys.path.insert(0, MTCLONE)
from mtclone import qc, classify

RUCK = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/GSE197008"
REDEEM = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/GSE219014"
OUT = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/plan2_step3"
os.makedirs(OUT, exist_ok=True)

# QC thresholds for informative-variant selection (clone-informative = somatic, sub-homoplasmic)
QC_KW = dict(max_pseudobulk_het=0.90,      # drop germline/homoplasmic (present in ~all cells)
             min_cells_detected=5,          # variant seen in >=5 cells within donor
             min_cell_coverage=5.0,          # real depth floor (both arms have depth_file coverage)
             binarize_threshold=0.05)        # a cell "carries" a variant at >=5% heteroplasmy


def _clr(M):
    X = np.asarray(M, dtype=float) + 1.0
    gm = np.exp(np.log(X).mean(1, keepdims=True))
    return np.log(X / gm)


def classify_ruckert_protein(a, method="canonical", bright_q=0.75, dim_q=0.25):
    """bright/dim on CLR-normalized ADT surface protein. Returns (labels, bright_axis).

    method='canonical' (PRIMARY): the flow-cytometry definition — CD56-bright = CD56-HIGH AND
        CD16-LOW. Gate on CD56 > `bright_q` quantile AND CD16 < `dim_q` quantile → bright;
        the rest → dim. This recovers the true ~5-10% blood bright fraction (a forced 2-way
        split does not — it cuts at the median). Quantiles are per-sample (batch-robust).
    method='data_driven' (SECONDARY / sensitivity): 2-component GMM on z(CD56)-z(CD16); the
        higher-CD56 component = bright. Kept for the "distinct subset vs remainder" framing.
    """
    names = list(a.uns["protein_names"])
    pc = pd.DataFrame(_clr(a.obsm["protein"]), columns=names, index=a.obs_names)
    z = lambda s: (s - s.mean()) / (s.std() or 1.0)
    axis = (z(pc["CD56"]) - z(pc["CD16"])).values                 # high = bright-like
    if method == "canonical":
        is_bright = (pc["CD56"].values > pc["CD56"].quantile(bright_q)) & \
                    (pc["CD16"].values < pc["CD16"].quantile(dim_q))
        lab = np.where(is_bright, "bright", "dim")
    elif method == "data_driven":
        from sklearn.mixture import GaussianMixture
        gm = GaussianMixture(n_components=2, n_init=5, random_state=0).fit(axis.reshape(-1, 1))
        cl = gm.predict(axis.reshape(-1, 1))
        m0 = pc["CD56"].values[cl == 0].mean() if (cl == 0).any() else -np.inf
        m1 = pc["CD56"].values[cl == 1].mean() if (cl == 1).any() else -np.inf
        lab = np.where(cl == (0 if m0 >= m1 else 1), "bright", "dim")
    else:
        raise ValueError(method)
    return lab, axis


def power_gate(a_donor):
    """select informative variants for one donor; return (n_cells, per-cell informative count)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            f = qc.select_informative_variants(a_donor, **QC_KW)
        except Exception as e:
            return None, None, str(e)
    if f.n_vars == 0 or f.n_obs == 0:
        return f, np.zeros(a_donor.n_obs, int), None
    B = f.layers["binary"]
    per_cell = np.asarray((B > 0).sum(1)).ravel()
    # map back to donor cells by cell_id (QC may have dropped low-coverage cells)
    s = pd.Series(per_cell, index=f.obs_names)
    return f, s, None


def run_ruckert():
    rows, embed = [], []
    for fp in sorted(glob.glob(f"{RUCK}/mtASAP*.h5ad")):
        a = ad.read_h5ad(fp)
        samp = os.path.basename(fp).replace(".h5ad", "")
        a.obs["cell_type"] = "NK"                    # sorted
        lab, axis = classify_ruckert_protein(a)
        a.obs["nk_subset"] = lab
        a.obs["bright_axis"] = axis
        for donor in a.obs["donor"].unique():
            sub = a[a.obs["donor"] == donor].copy()
            f, per_cell, err = power_gate(sub)
            n_inf_var = int(f.n_vars) if f is not None else 0
            for subset in ["bright", "dim"]:
                cells = sub.obs_names[sub.obs["nk_subset"] == subset]
                pc = per_cell.reindex(cells).fillna(0) if per_cell is not None else pd.Series(0, index=cells)
                rows.append(dict(dataset="Ruckert", arm="blood", sample=samp, donor=donor,
                                 subset=subset, n_nk=len(cells),
                                 n_with_informative=int((pc >= 1).sum()),
                                 med_inf_var_per_cell=float(np.median(pc)) if len(pc) else 0,
                                 med_coverage=round(float(np.nanmedian(sub.obs.loc[cells, "coverage"])), 1),
                                 n_informative_variants=n_inf_var))
        embed.append((samp, a.obs["bright_axis"].values, a.obs["nk_subset"].values,
                      _clr(a.obsm["protein"])[:, list(a.uns["protein_names"]).index("CD16")]))
        a.write_h5ad(f"{OUT}/{samp}.labelled.h5ad")
    return pd.DataFrame(rows), embed


def run_redeem():
    rows = []
    for fp in sorted(glob.glob(f"{REDEEM}/*.BMMC.h5ad")):
        a = ad.read_h5ad(fp)
        samp = os.path.basename(fp).replace(".BMMC.h5ad", "")
        donor = a.obs["donor"].iloc[0]
        a.obs["cell_type"] = np.where(a.obs["STD.CellType"] == "NK", "NK", "other")
        # build an expression AnnData from the RNA marker panel for classify._score
        mk = ad.AnnData(X=np.asarray(a.obsm["rna_markers"], dtype=np.float32),
                        obs=a.obs[["cell_type"]].copy())
        mk.var_names = list(a.uns["rna_marker_names"])
        classify.label_bright_dim(mk, nk_mask_key="cell_type", nk_value="NK",
                                  method="data_driven", key_added="nk_subset")
        a.obs["nk_subset"] = mk.obs["nk_subset"].values
        nk = a[a.obs["cell_type"] == "NK"].copy()
        if nk.n_obs == 0:
            rows.append(dict(dataset="ReDeeM", arm="marrow", sample=samp, donor=donor,
                             subset="(none)", n_nk=0, n_with_informative=0,
                             med_inf_var_per_cell=0, med_coverage=0, n_informative_variants=0))
            continue
        f, per_cell, err = power_gate(nk)
        n_inf_var = int(f.n_vars) if f is not None else 0
        for subset in ["bright", "dim"]:
            cells = nk.obs_names[nk.obs["nk_subset"] == subset]
            pc = per_cell.reindex(cells).fillna(0) if per_cell is not None else pd.Series(0, index=cells)
            rows.append(dict(dataset="ReDeeM", arm="marrow", sample=samp, donor=donor,
                             subset=subset, n_nk=len(cells),
                             n_with_informative=int((pc >= 1).sum()),
                             med_inf_var_per_cell=float(np.median(pc)) if len(pc) else 0,
                             med_coverage=round(float(np.nanmedian(nk.obs.loc[cells, "coverage"])), 1),
                             n_informative_variants=n_inf_var))
        a.write_h5ad(f"{OUT}/{samp}.labelled.h5ad")
    return pd.DataFrame(rows)


if __name__ == "__main__":
    dr, embed = run_ruckert()
    dd = run_redeem()
    df = pd.concat([dr, dd], ignore_index=True)
    df.to_csv(f"{OUT}/plan2_step3_counts.csv", index=False)
    print(df.to_string(index=False))
    import pickle
    pickle.dump(embed, open(f"{OUT}/_ruckert_embed.pkl", "wb"))
