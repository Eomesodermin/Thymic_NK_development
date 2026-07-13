"""Step 6 — tumor-vs-blood NK clone overlap (the seeding test), Plan 3.

For each donor with tumor+blood NK, using Step-4 clone assignments:
  1. directional clone sharing (clone_sharing_matrix): fraction of TUMOR clones also in blood,
     and reciprocal fraction of BLOOD clones also in tumor; plus tumor+met+blood where present.
  2. tumor-restricted vs shared EXPANDED clones: of clones with >=`expanded_min` tumor cells,
     how many are tumor-restricted (absent from blood) vs shared with blood.
  3. permutation null (permutation_null): shuffle site labels WITHIN the donor, recompute the
     tumor<->blood mixing statistic (between_vs_within_sharing frac_mixed_clones). Observed mixing
     BELOW null => clones more site-segregated than chance => local restriction/seeding.
     Observed ~ null => tumor mirrors blood => passive polyclonal infiltration.

Run for base-NK and strict NK-not-ILC. n=1000 permutations, seed=0.

Usage:
    python tumor_blood_overlap.py --clone-dir <clones> --out <dir> [--n-perm 1000]
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad

sys.path.insert(0, "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts")
from mtclone import metrics


def restricted_vs_shared(sel, *, expanded_min=5, clone_key="clone_id"):
    """Of clones with >= expanded_min tumor cells, count tumor-restricted vs shared-with-blood."""
    obs = sel.obs[[ "site", clone_key]].copy()
    obs = obs[obs[clone_key] >= 0]
    tum = obs[obs["site"] == "tumor"]
    blo_clones = set(obs[obs["site"] == "blood"][clone_key].unique())
    tum_sizes = tum.groupby(clone_key).size()
    expanded = tum_sizes[tum_sizes >= expanded_min].index
    n_exp = len(expanded)
    shared = sum(1 for c in expanded if c in blo_clones)
    restricted = n_exp - shared
    return dict(n_expanded_tumor_clones=int(n_exp),
                n_tumor_restricted=int(restricted),
                n_shared_with_blood=int(shared),
                frac_tumor_restricted=float(restricted / n_exp) if n_exp else np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clone-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--expanded-min", type=int, default=5)
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    share_rows, perm_rows = [], []
    for tag in ["nk", "nk_not_ilc"]:
        for h in sorted(Path(args.clone_dir).glob(f"*.{tag}.clones.h5ad")):
            donor = h.stem.replace(f".{tag}.clones", "")
            sel = ad.read_h5ad(h)
            sites = set(sel.obs["site"].unique())
            if not {"tumor", "blood"} <= sites:
                continue

            # 1+2 directional sharing + restricted/shared expanded
            M = metrics.clone_sharing_matrix(sel, group_key="site", min_clone_size=2)
            rs = restricted_vs_shared(sel, expanded_min=args.expanded_min)
            row = dict(donor=donor, tag=tag,
                       frac_tumor_clones_in_blood=float(M.loc["tumor", "blood"]),
                       frac_blood_clones_in_tumor=float(M.loc["blood", "tumor"]))
            if "met" in sites:
                row["frac_tumor_clones_in_met"] = float(M.loc["tumor", "met"])
                row["frac_met_clones_in_blood"] = float(M.loc["met", "blood"])
            row.update(rs)
            share_rows.append(row)

            # 3 permutation null on tumor<->blood mixing (restrict to tumor+blood cells)
            tb = sel[sel.obs["site"].isin(["tumor", "blood"])].copy()
            stat = lambda a: metrics.between_vs_within_sharing(
                a, "tumor", "blood", group_key="site", min_clone_size=2)["frac_mixed_clones"]
            pn = metrics.permutation_null(tb, stat, group_key="site",
                                          stratify_by="donor", n=args.n_perm, seed=0)
            perm_rows.append(dict(donor=donor, tag=tag, observed_frac_mixed=pn["observed"],
                                  null_mean=pn["null_mean"], null_sd=pn["null_sd"],
                                  p_less=pn["p_less"], p_greater=pn["p_greater"], n_null=pn["n_null"]))
            print(f"{donor} [{tag}]: tumor->blood share {row['frac_tumor_clones_in_blood']:.2f}, "
                  f"restricted {rs['n_tumor_restricted']}/{rs['n_expanded_tumor_clones']}, "
                  f"mixing obs {pn['observed']:.3f} vs null {pn['null_mean']:.3f} "
                  f"(p_less={pn['p_less']:.3f})", flush=True)

    sdf = pd.DataFrame(share_rows); pdf = pd.DataFrame(perm_rows)
    sdf.to_csv(out / "tumor_blood_sharing.csv", index=False)
    pdf.to_csv(out / "tumor_blood_permutation.csv", index=False)
    print("\n=== sharing (base NK) ===")
    print(sdf[sdf["tag"]=="nk"].to_string(index=False))
    print("\n=== permutation null (base NK) ===")
    print(pdf[pdf["tag"]=="nk"].to_string(index=False))


if __name__ == "__main__":
    main()
