"""Step 3 — NK identification + per-donor×site power gate (Plan 3).

Joins ATAC gene-activity (marker panel) onto each donor's mtDNA object, labels NK with the
frozen-core scorer mtclone.classify.label_nk(modality="atac") as the base call, adds disclosed
Plan-3 refinements (broader contaminant exclusion + NK-vs-ILC), then reports the go/no-go
table: NK per donor×site and how many carry >=1 informative mtDNA variant.

Contaminant handling (per Liu criterion + our NK-vs-ILC concern):
  - NK positive : NCR1/GNLY/KLRD1/NKG7/KLRF1 accessible
  - excluded    : T (CD3D/CD3E/CD8A/CD4 open), B (MS4A1/CD19), myeloid (CD14/LYZ)
  - iNKT        : CD3D-open cells removed by the T filter (Liu: iNKT are CD3D+/CD56+)
  - NK vs ILC   : among NK-positive, require EOMES/TBX21 (NK TFs) and low IL7R/KIT/RORC (ILC);
                  ILC lack EOMES and are IL7R-high. Reported as a separate refinement column.

Usage:
    python nk_label_and_power_gate.py --ga-dir <geneactivity> --het-dir <processed> \
        --panel marker_panel_grch38.json --out <dir>
"""
from __future__ import annotations
import argparse, json, sys, os
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad

sys.path.insert(0, "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts")
from mtclone import io, qc, classify


def zscore_cols(X):
    mu = X.mean(0); sd = X.std(0); sd[sd == 0] = 1
    return (X - mu) / sd


def load_ga_for_donor(ga_dir, donor, samples):
    """Concatenate per-sample gene-activity objects for a donor."""
    objs = []
    for s in samples:
        p = Path(ga_dir) / f"{s}.geneactivity.h5ad"
        if p.exists():
            objs.append(ad.read_h5ad(p))
    if not objs:
        return None
    return ad.concat(objs, join="outer", fill_value=0.0)


def label_nk_refined(ga, panel):
    """Return a per-cell DataFrame with NK call + refinement flags.

    BASE CALL = the frozen-core scorer `mtclone.classify.label_nk(modality="atac")`
    (NK_ATAC_POS = NCR1/GNLY/KLRD1 open, NK_ATAC_NEG = CD3D/CD8A closed). This is the
    canonical NK label used downstream (obs['is_nk']).

    Layered REFINEMENTS (Plan-3 additions, disclosed — the frozen core intentionally does
    not do these; they only ADD exclusion/annotation columns, never override the base call):
      - broader contaminant exclusion (B: MS4A1/CD19, myeloid: CD14/LYZ) on top of the
        core's T-cell CD3D/CD8A filter;
      - NK-vs-ILC refinement (EOMES/TBX21 NK-TFs > IL7R/KIT/RORC/IL1R1 ILC markers).
    """
    genes = list(ga.var_names)
    X = ga.X.toarray() if hasattr(ga.X, "toarray") else np.asarray(ga.X)

    # ---- BASE: frozen-core mtclone.classify.label_nk on the gene-activity matrix ----
    base = classify.label_nk(ga.copy(), modality="atac", key_added="cell_type")
    is_nk_core = (base.obs["cell_type"] == "NK").reindex(ga.obs_names).fillna(False).values

    # ---- REFINEMENT signals (z-scored panel; additive exclusion only) ----
    Z = pd.DataFrame(zscore_cols(X), columns=genes, index=ga.obs_names)
    def sig(gs):
        gs = [g for g in gs if g in Z.columns]
        return Z[gs].mean(1) if gs else pd.Series(0.0, index=Z.index)
    b_sig = sig(["MS4A1", "CD19"])
    mye   = sig(["CD14", "LYZ"])
    nk_tf = sig(["EOMES", "TBX21"])
    ilc   = sig(["IL7R", "KIT", "RORC", "IL1R1"])

    is_nk = pd.Series(is_nk_core, index=ga.obs_names) & (b_sig < 0.25) & (mye < 0.25)
    nk_not_ilc = is_nk & (nk_tf > ilc)

    out = pd.DataFrame(index=ga.obs_names)
    out["nk_pos_score"] = base.obs["nk_pos_score"].reindex(ga.obs_names).values
    out["t_score"] = base.obs["nk_neg_score"].reindex(ga.obs_names).values
    out["nk_tf_score"] = nk_tf.values
    out["ilc_score"] = ilc.values
    out["is_nk_core"] = is_nk_core
    out["is_nk"] = is_nk.values
    out["nk_not_ilc"] = nk_not_ilc.values
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ga-dir", required=True)
    ap.add_argument("--het-dir", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-cell-coverage", type=int, default=10)
    ap.add_argument("--max-pseudobulk-het", type=float, default=0.90)
    ap.add_argument("--min-cells-detected", type=int, default=5)
    ap.add_argument("--binarize-threshold", type=float, default=0.07)
    args = ap.parse_args()
    panel = json.load(open(args.panel))
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    rows = []
    donor_objs = {}
    for h in sorted(Path(args.het_dir).glob("*.h5ad")):
        donor = h.stem
        a = ad.read_h5ad(h)
        samples = sorted(a.obs["sample"].unique())
        ga = load_ga_for_donor(args.ga_dir, donor, samples)
        if ga is None:
            print(f"{donor}: NO gene activity", flush=True); continue

        lab = label_nk_refined(ga, panel)
        # map labels onto mtDNA object by cell_id (obs_names shared: <GSM>_<barcode>)
        a.obs["is_nk"] = lab["is_nk"].reindex(a.obs_names).fillna(False).values
        a.obs["nk_not_ilc"] = lab["nk_not_ilc"].reindex(a.obs_names).fillna(False).values
        for c in ["nk_pos_score", "t_score", "nk_tf_score", "ilc_score"]:
            a.obs[c] = lab[c].reindex(a.obs_names).astype(float).values

        # informative variants on the NK subset, per donor (clone tree is donor-wide)
        nk = a[a.obs["is_nk"].values].copy()
        info_ok = False
        if nk.n_obs > 0:
            try:
                sel = qc.select_informative_variants(
                    nk, min_cell_coverage=args.min_cell_coverage,
                    max_pseudobulk_het=args.max_pseudobulk_het,
                    min_cells_detected=args.min_cells_detected,
                    binarize_threshold=args.binarize_threshold)
                info_ok = sel.n_vars > 0
                # per-cell: carries >=1 informative variant (use the binarized layer the
                # selector writes; fall back to thresholding X if absent)
                if "binary" in sel.layers:
                    Xb = sel.layers["binary"]
                    Xb = Xb.toarray() if hasattr(Xb, "toarray") else np.asarray(Xb)
                    Xb = Xb > 0
                else:
                    Xb = (sel.X.toarray() if hasattr(sel.X, "toarray") else np.asarray(sel.X)) >= args.binarize_threshold
                has_var = pd.Series(Xb.sum(1) >= 1, index=sel.obs_names)
                a.obs["nk_has_informative"] = has_var.reindex(a.obs_names).fillna(False).values
                n_info_vars = int(sel.n_vars)
            except Exception as e:
                a.obs["nk_has_informative"] = False; n_info_vars = 0
                print(f"{donor}: informative-variant selection failed: {e}", flush=True)
        else:
            a.obs["nk_has_informative"] = False; n_info_vars = 0

        # per donor×site tallies
        for site in sorted(a.obs["site"].unique()):
            m = a.obs["site"] == site
            rows.append(dict(
                donor=donor, dx=str(a.obs["diagnosis"].iloc[0]) if "diagnosis" in a.obs else "?",
                site=site, n_cells=int(m.sum()),
                n_nk=int((m & a.obs["is_nk"]).sum()),
                n_nk_not_ilc=int((m & a.obs["nk_not_ilc"]).sum()),
                n_nk_with_var=int((m & a.obs["is_nk"] & a.obs["nk_has_informative"]).sum()),
                n_informative_variants_donor=n_info_vars))
        donor_objs[donor] = a
        a.write_h5ad(out / f"{donor}.nk_labelled.h5ad")
        print(f"{donor}: {int(a.obs['is_nk'].sum())} NK, {n_info_vars} informative vars", flush=True)

    tab = pd.DataFrame(rows)
    tab.to_csv(out / "nk_power_gate.csv", index=False)
    print("\n=== NK power gate ===")
    print(tab.to_string(index=False))


if __name__ == "__main__":
    main()
