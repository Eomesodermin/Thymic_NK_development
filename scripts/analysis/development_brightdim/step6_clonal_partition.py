#!/usr/bin/env python
"""step6_clonal_partition.py — Plan 2 Step 6: does the clonal partition map onto bright/dim?

The converse of Step 5. Step 5 started from the a priori bright/dim labels and asked whether they
share clones. Step 6 starts from the *clonal partition* (which cells share lineages, ignoring
labels) and asks whether the major lineages carry ANY bright/dim transcriptional signal — the
"ideally maps to bright/dim" step, and a direct test of the user's "see two clonal lineages"
framing.

Axis score (continuous, per cell):
  * ReDeeM  : z-scored RNA bright markers (SELL,XCL1,IL7R,GZMK,NCAM1) minus dim markers
              (FCGR3A,PRF1,FGFBP2,GZMB,CX3CR1).  + = bright-like, - = dim-like.
  * Rückert : z-scored CD56 - CD16 surface protein (bright = CD56hi/CD16lo).

Statistic: eta^2 = fraction of the axis-score variance explained by the clonal partition
(clones with >=10 cells). Null: permute clone labels among these cells, 1000x. eta^2 >> null =>
lineages track bright/dim (two-lineage-like); eta^2 ~ null and small => lineages are
transcriptionally mixed (continuum).

NOTE: Rückert variant_group clones are mostly size 2-4, so 6 of 8 donor-units have <2 clones
reaching >=10 cells and cannot support the large-lineage eta^2 test (NaN). The 2 exceptions
(mtASAP2 CMVpos3 n=14, CMVpos4 n=18) ARE tested and give eta^2 indistinguishable from null
(p=0.29, 0.42) — corroborating the ReDeeM result. Rückert's arm otherwise rests on Step 5.
NOTE: ENKP/ILCP (2024 Nat Immunol) signatures are NOT computed — the exported 20-gene panel shares
only IL7R with those progenitor sets; recovering them needs a full-transcriptome re-export from the
Seurat objects (flagged as a limitation, not run).

Outputs (processed/plan2_step6/): plan2_step6_eta2.csv, plan2_step6_clonal_partition.png
"""
import sys, os, glob
import numpy as np, pandas as pd, anndata as ad

sys.path.insert(0, "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts")
S4 = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/plan2_step4"
OUT = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/plan2_step6"
os.makedirs(OUT, exist_ok=True)

BRIGHT = ["SELL", "XCL1", "IL7R", "GZMK", "NCAM1"]
DIM = ["FCGR3A", "PRF1", "FGFBP2", "GZMB", "CX3CR1"]

def redeem_score(a):
    M = pd.DataFrame(a.obsm["rna_markers"], index=a.obs_names, columns=list(a.uns["rna_marker_names"]))
    Z = (M - M.mean()) / M.std(ddof=0).replace(0, 1)
    return (Z[BRIGHT].mean(1) - Z[DIM].mean(1)).values

def ruckert_score(a):
    P = pd.DataFrame(a.obsm["protein"], index=a.obs_names, columns=list(a.uns["protein_names"]))
    Z = (P - P.mean()) / P.std(ddof=0).replace(0, 1)
    return (Z["CD56"] - Z["CD16"]).values

def eta2_and_null(sc, cid, min_size=10, nperm=1000, seed=0):
    vc = pd.Series(cid).value_counts(); big = vc[vc >= min_size].index.tolist()
    if len(big) < 2:
        return np.nan, np.nan, np.nan, len(big)
    mask = np.isin(cid, big); s = sc[mask]; c = cid[mask]
    def eta2(s, c):
        grand = s.mean()
        ssb = sum((c == k).sum() * (s[c == k].mean() - grand) ** 2 for k in np.unique(c))
        return ssb / ((s - grand) ** 2).sum()
    obs = eta2(s, c)
    rng = np.random.default_rng(seed)
    null = np.array([eta2(s, rng.permutation(c)) for _ in range(nperm)])
    p = (np.sum(null >= obs) + 1) / (nperm + 1)
    return obs, float(null.mean()), float(p), len(big)

def run():
    rows = []
    for fp in sorted(glob.glob(f"{S4}/*.nkclones.h5ad")):
        b = os.path.basename(fp).replace(".nkclones.h5ad", ""); obj, donor = b.rsplit(".", 1)
        arm = "Ruckert" if obj.startswith("mtASAP") else "ReDeeM"
        a = ad.read_h5ad(fp); a = a[a.obs["clone_id"] >= 0].copy()
        if a.n_obs < 20:
            continue
        sc = ruckert_score(a) if arm == "Ruckert" else redeem_score(a)
        obs, nm, p, nbig = eta2_and_null(sc, a.obs["clone_id"].values)
        rows.append(dict(object=obj, donor=donor, arm=arm, n_clones_ge10=nbig,
                         eta2_obs=obs, eta2_null=nm, p_perm=p))
    return pd.DataFrame(rows)

if __name__ == "__main__":
    df = run()
    df.to_csv(f"{OUT}/plan2_step6_eta2.csv", index=False)
    print(df.to_string(index=False))
