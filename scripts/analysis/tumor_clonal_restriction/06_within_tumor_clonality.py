"""Step 5 — within-tumor NK clonality (Plan 3).

For each donor with tumor NK, compute clone-size distribution + oligoclonality indices (Gini,
normalized Shannon / Pielou evenness, max-clone fraction) per SITE, using the Step-4 clone
assignments. Compare tumor-NK clonality to the same patient's blood-NK and normal-NK (the
within-patient reference — restriction should show as HIGHER skew in tumor). Run for both the
base-NK and strict NK-not-ILC clone sets (sensitivity).

External oligoclonality benchmark (Ruckert blood adaptive-NK, GSE197008/197037) is Plan-2 scope
and not yet produced; the within-patient blood-NK comparison is used as the primary yardstick
here and the external benchmark is left as a documented follow-up.

Usage:
    python within_tumor_clonality.py --clone-dir <clones> --out <dir>
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad

sys.path.insert(0, "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts")
from mtclone import metrics


def per_site_clonality(sel, tag, donor):
    rows = []
    for site in sorted(sel.obs["site"].unique()):
        s = metrics.clonality_summary(sel, group=("site", site))
        s.update(donor=donor, tag=tag, site=site)
        rows.append(s)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clone-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    rows = []
    for tag in ["nk", "nk_not_ilc"]:
        for h in sorted(Path(args.clone_dir).glob(f"*.{tag}.clones.h5ad")):
            donor = h.stem.replace(f".{tag}.clones", "")
            sel = ad.read_h5ad(h)
            rows += per_site_clonality(sel, tag, donor)
    df = pd.DataFrame(rows)
    df.to_csv(out / "within_tumor_clonality.csv", index=False)

    # tumor-vs-blood contrast per donor (base NK)
    base = df[df["tag"] == "nk"]
    piv_g = base.pivot_table(index="donor", columns="site", values="gini")
    piv_e = base.pivot_table(index="donor", columns="site", values="normalized_shannon")
    contrast = pd.DataFrame({
        "gini_tumor": piv_g.get("tumor"), "gini_blood": piv_g.get("blood"),
        "gini_normal": piv_g.get("normal"),
        "even_tumor": piv_e.get("tumor"), "even_blood": piv_e.get("blood"),
    })
    contrast["tumor_more_skewed_than_blood"] = contrast["gini_tumor"] > contrast["gini_blood"]
    contrast.to_csv(out / "tumor_vs_blood_clonality_contrast.csv")
    print("=== per-site clonality (base NK) ===")
    print(base[["donor","site","n_clones","n_cells","gini","normalized_shannon","max_clone_frac"]]
          .to_string(index=False))
    print("\n=== tumor vs blood (base NK) ===")
    print(contrast.to_string())


if __name__ == "__main__":
    main()
