#!/usr/bin/env python
"""step5_sharing_test.py — Plan 2 Step 5: the core one-lineage-vs-two test.

Question: do bright and dim NK share somatic-mtDNA clones (one well-mixed lineage / continuum)
or are their clones disjoint (two independent developmental pathways)?

Statistic (per donor): `frac_mixed_clones` from mtclone.metrics.between_vs_within_sharing —
of clones containing >=1 bright or dim NK, the fraction containing BOTH. High => shared
lineages (continuum); ~0 => disjoint (two lineages).

Null: shuffle bright/dim labels WITHIN donor, holding clone structure fixed, 1000x. This holds
clone count/size and the (rare) bright base-rate fixed, so it asks precisely: given this clone
structure and this bright fraction, how much bright/dim mixing would one well-mixed lineage
produce by chance? Observed << null => two lineages; observed ~ null => continuum.
  - p_less  = frac(null <= observed): SMALL => observed sharing is unusually LOW => two lineages.
  - p_greater = frac(null >= observed): SMALL => unusually HIGH sharing.

Aggregation: clones are donor-private, so within each ARM we concatenate donor-units, make
clone_id globally unique (donor-prefixed), and run one within-donor-stratified permutation over
the pooled statistic (mean of per-donor frac_mixed, weighted by testable clones). Arms use
different clone definitions (Step 4), so we test WITHIN arm and compare, never pool across arms.

Outputs (processed/plan2_step5/):
  plan2_step5_perdonor.csv    — per donor-unit: composition + observed/null/p
  plan2_step5_pooled.csv      — per arm: pooled observed, null mean/sd, p, CI
  plan2_step5_sharing.png     — headline: composition + observed-vs-null per arm
"""
import sys, os, glob, warnings
import numpy as np, pandas as pd, anndata as ad

MT = "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts"
sys.path.insert(0, MT)
from mtclone import metrics

S4 = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/plan2_step4"
OUT = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/plan2_step5"
os.makedirs(OUT, exist_ok=True)
NPERM = 1000

def _arm(obj):
    return "Ruckert" if obj.startswith("mtASAP") else "ReDeeM"

# ---- vectorized statistic + within-donor permutation (validated == metrics reference; the
#      metrics set-lambda groupby recomputed 1000x is too slow on ~1,600-clone Ruckert units).
def frac_mixed_fast(clone_code, is_bright, K, min_size=2):
    """frac of testable clones (size>=min_size) that contain BOTH a bright and a dim NK."""
    b = np.bincount(clone_code[is_bright], minlength=K)
    d = np.bincount(clone_code[~is_bright], minlength=K)
    keep = (b + d) >= min_size
    if keep.sum() == 0:
        return np.nan
    mixed = (b >= 1) & (d >= 1)
    return float(mixed[keep].mean())

def perm_test_fast(clone_code, is_bright, donor_code, K, n=NPERM, seed=0, min_size=2):
    rng = np.random.default_rng(seed)
    obs = frac_mixed_fast(clone_code, is_bright, K, min_size)
    idx_by_d = {s: np.where(donor_code == s)[0] for s in np.unique(donor_code)}
    lab = is_bright.copy()
    null = np.empty(n)
    for k in range(n):
        perm = lab.copy()
        for _, idx in idx_by_d.items():
            perm[idx] = rng.permutation(lab[idx])   # shuffle bright/dim WITHIN donor
        null[k] = frac_mixed_fast(clone_code, perm, K, min_size)
    null = null[~np.isnan(null)]
    return dict(observed=obs, null_mean=float(null.mean()), null_sd=float(null.std()),
                p_less=float((np.sum(null <= obs) + 1) / (null.size + 1)),
                p_greater=float((np.sum(null >= obs) + 1) / (null.size + 1)),
                n_null=int(null.size), null=null)

def load_units():
    """Return list of (object, donor, arm, AnnData of classifiable NK with clone_id)."""
    units = []
    for fp in sorted(glob.glob(f"{S4}/*.nkclones.h5ad")):
        base = os.path.basename(fp).replace(".nkclones.h5ad", "")
        obj, donor = base.rsplit(".", 1)
        a = ad.read_h5ad(fp)
        a = a[a.obs["nk_subset"].isin(["bright", "dim"])].copy()   # only classifiable NK
        a = a[a.obs["clone_id"] >= 0].copy()                        # only clone-assigned
        if a.n_obs < 2:
            continue
        units.append((obj, donor, _arm(obj), a))
    return units

def per_donor(units):
    rows = []
    for obj, donor, arm, a in units:
        code, uniq = pd.factorize(a.obs["clone_id"].values); K = len(uniq)
        ib = (a.obs["nk_subset"].values == "bright")
        dc = np.zeros(a.n_obs, dtype=int)   # single donor per object
        r = perm_test_fast(code, ib, dc, K, n=NPERM, seed=0)
        b = np.bincount(code[ib], minlength=K); d = np.bincount(code[~ib], minlength=K)
        keep = (b + d) >= 2
        rows.append(dict(
            object=obj, donor=donor, arm=arm,
            n_bright=int(ib.sum()), n_dim=int((~ib).sum()),
            n_clones_total=int(keep.sum()),
            n_mixed=int(((b >= 1) & (d >= 1))[keep].sum()),
            n_bright_only=int(((b >= 1) & (d == 0))[keep].sum()),
            n_dim_only=int(((d >= 1) & (b == 0))[keep].sum()),
            obs_frac_mixed=r["observed"], null_mean=r["null_mean"], null_sd=r["null_sd"],
            p_less=r["p_less"], p_greater=r["p_greater"],
        ))
    return pd.DataFrame(rows)

def pooled_by_arm(units):
    """Concatenate within arm, donor-prefix clone_id, run one within-donor permutation."""
    rows = []
    null_by_arm = {}
    for arm in ["ReDeeM", "Ruckert"]:
        sub = [(obj, donor, a) for obj, donor, ar, a in units if ar == arm]
        if not sub:
            continue
        parts = []
        for obj, donor, a in sub:
            b = a.obs[["nk_subset", "clone_id"]].copy()
            b["donor_unit"] = f"{obj}.{donor}"
            b["clone_key_str"] = b["donor_unit"] + "::" + b["clone_id"].astype(str)
            parts.append(b)
        # donor-prefixed clone strings -> integer codes so clones stay donor-private; the
        # permutation shuffles labels WITHIN donor_unit, so pooled mixing is within-donor.
        obs_all = pd.concat(parts, ignore_index=True)
        code = pd.factorize(obs_all["clone_key_str"])[0]; K = int(code.max()) + 1
        ib = (obs_all["nk_subset"].values == "bright")
        dcode = pd.factorize(obs_all["donor_unit"])[0]
        r = perm_test_fast(code, ib, dcode, K, n=NPERM, seed=0)
        null = r["null"]; null_by_arm[arm] = null
        ci = (np.percentile(null, 2.5), np.percentile(null, 97.5)) if null.size else (np.nan, np.nan)
        rows.append(dict(arm=arm, n_donor_units=len(sub),
                         observed=r["observed"], null_mean=r["null_mean"], null_sd=r["null_sd"],
                         null_ci_lo=ci[0], null_ci_hi=ci[1],
                         p_less=r["p_less"], p_greater=r["p_greater"], n_null=r["n_null"]))
    return pd.DataFrame(rows), null_by_arm

if __name__ == "__main__":
    units = load_units()
    pd_df = per_donor(units)
    pd_df.to_csv(f"{OUT}/plan2_step5_perdonor.csv", index=False)
    pool_df, nulls = pooled_by_arm(units)
    pool_df.to_csv(f"{OUT}/plan2_step5_pooled.csv", index=False)
    np.savez(f"{OUT}/plan2_step5_nulls.npz", **{k: v for k, v in nulls.items()})
    print("=== per donor-unit ===")
    print(pd_df[["object","donor","arm","n_clones_total","n_mixed","obs_frac_mixed",
                 "null_mean","p_less","p_greater"]].to_string(index=False))
    print("\n=== pooled by arm ===")
    print(pool_df.to_string(index=False))
