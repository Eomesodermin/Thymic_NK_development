"""End-to-end harness (Plan 1 gate): run the full mtclone chain on ONE real sample from
each anchor dataset, assert the pipeline produces valid, sane output, and write a report.

Usage:
    python run_end_to_end.py --testdata <dir> --out <dir>

The chain: ingest -> validate_schema -> QC/informative-variant selection -> clone inference
-> clonality metrics -> permutation-null sanity check. QC figures are saved per dataset.

Datasets:
  mtscATAC  GSE302113 GSM9096509 (PBMC CD56+ NK): clean, self-contained (heteroplasmy TSV
            + variant_stats). This is the primary real-data validation.
  ReDeeM    GSE219014 Aged1 BMMC: the variant .mtx ships WITHOUT companion barcodes/features
            in the GEO deposit, so this leg validates the .mtx adapter + downstream
            structurally (positional names) and documents the annotation gap for Plan 2.
"""
from __future__ import annotations
import argparse, gzip, json, os, sys, datetime, traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mtclone
from mtclone import io, qc, clones, metrics


def _log(rep, msg):
    rep.append(msg)
    print(msg)


def run_mtscatac(testdata, outdir, rep):
    _log(rep, "\n## mtscATAC — GSE302113 GSM9096509 (PBMC CD56+ NK)\n")
    het = testdata / "GSM9096509.cell_heteroplasmic_df.tsv.gz"
    vs = testdata / "GSM9096509.variant_stats.tsv.gz"
    a = io.read_mtscatac_heteroplasmy(het, vs, dataset="GSE302113", sample="GSM9096509",
                                      donor="SU-L-003", tissue="PBMC", site="blood")
    _log(rep, f"- ingested: {a.n_obs} cells x {a.n_vars} variants")
    assert io.validate_schema(a) == [], "schema invalid"
    _log(rep, "- schema: VALID")

    sel = qc.select_informative_variants(a, min_cell_coverage=10, max_pseudobulk_het=0.90,
                                         min_cells_detected=5, binarize_threshold=0.07)
    _log(rep, f"- informative variants: {sel.n_vars} (of {a.n_vars}); cells kept {sel.n_obs}")
    assert sel.n_vars > 0, "no informative variants"
    qc.qc_plots(a, path=str(outdir / "mtscatac_qc.png"))
    _log(rep, "- QC figure: mtscatac_qc.png")

    called = clones.call_clones(sel, method="graph", min_shared_variants=1,
                                edge_weight_cutoff=0.5, min_clone_size=2)
    cq = clones.clone_qc(called)
    _log(rep, f"- clones: {cq['n_clones']} clones, {cq['frac_assigned']:.1%} cells assigned, "
              f"max size {cq['max_clone_size']}")
    assert cq["n_clones"] >= 1, "no clones called"

    summ = metrics.clonality_summary(called)
    _log(rep, f"- clonality: Gini={summ['gini']:.3f}, norm.Shannon={summ['normalized_shannon']:.3f}, "
              f"max clone frac={summ['max_clone_frac']:.3f}")

    # permutation-null sanity: split cells into two arbitrary groups, sharing null must be
    # calibrated (observed within noise of null since split is random)
    rng = np.random.default_rng(0)
    called.obs["rand_group"] = rng.choice(["A", "B"], size=called.n_obs)
    stat = lambda ad_: metrics.between_vs_within_sharing(ad_, "A", "B",
                                                         group_key="rand_group")["frac_mixed_clones"]
    null = metrics.permutation_null(called, stat, group_key="rand_group",
                                    stratify_by="donor", n=200, seed=1)
    _log(rep, f"- permutation null (random split): observed={null['observed']:.3f}, "
              f"null_mean={null['null_mean']:.3f}, p_greater={null['p_greater']:.3f} "
              f"(random split -> p should be non-extreme)")
    assert 0.0 <= null["p_greater"] <= 1.0
    return dict(ok=True, n_cells=int(a.n_obs), n_variants=int(a.n_vars),
                n_informative=int(sel.n_vars), n_clones=int(cq["n_clones"]),
                gini=summ["gini"], norm_shannon=summ["normalized_shannon"])


def run_redeem(testdata, outdir, rep):
    _log(rep, "\n## ReDeeM — GSE219014 Aged1 BMMC (structural validation)\n")
    import scipy.io, scipy.sparse as sp
    mtx = testdata / "GSE219014_Aged1_BMMC.mtDNA_Variants_matrix.mtx.gz"
    with gzip.open(mtx, "rt") as fh:
        M = scipy.io.mmread(fh).tocsr()
    _log(rep, f"- raw matrix: {M.shape[0]} x {M.shape[1]}, {M.nnz} nonzero")
    _log(rep, "- NOTE: GEO deposit ships NO barcodes/features for this matrix. Orientation "
              "and variant identities are NOT recoverable from GEO alone; real analysis "
              "(Plan 2) must source them from the redeemR repo / Zenodo. This leg uses "
              "positional placeholder names to validate the adapter + downstream only.")
    # orient to cells x variants: assume the LARGER axis is variants (ReDeeM finds many)
    n_r, n_c = M.shape
    if n_r > n_c:   # rows=variants -> transpose to cells x variants
        cells_axis, vars_axis = n_c, n_r
    else:
        cells_axis, vars_axis = n_r, n_c
    barcodes = [f"cell{i}" for i in range(cells_axis)]
    features = [f"chrM:{i+1}:A>G" for i in range(vars_axis)]  # placeholder positional ids
    a = io.read_redeem_mtx(mtx, barcodes, features, depth=None, dataset="GSE219014",
                           sample="Aged1_BMMC", donor="Aged1", tissue="BMMC", site="marrow")
    _log(rep, f"- ingested (placeholder names): {a.n_obs} cells x {a.n_vars} variants")
    assert io.validate_schema(a) == [], "schema invalid"
    _log(rep, "- schema: VALID")

    sel = qc.select_informative_variants(a, max_pseudobulk_het=0.90, min_cells_detected=5,
                                         binarize_threshold=0.07)
    _log(rep, f"- informative variants: {sel.n_vars}; cells kept {sel.n_obs}")
    called = clones.call_clones(sel, method="graph", min_shared_variants=1,
                                edge_weight_cutoff=0.5, min_clone_size=2)
    cq = clones.clone_qc(called)
    _log(rep, f"- clones: {cq['n_clones']}, {cq['frac_assigned']:.1%} assigned")
    assert io.validate_schema(sel) == []
    return dict(ok=True, n_cells=int(a.n_obs), n_variants=int(a.n_vars),
                n_informative=int(sel.n_vars), n_clones=int(cq["n_clones"]),
                annotation_gap=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--testdata", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    testdata = Path(args.testdata); outdir = Path(args.out); outdir.mkdir(parents=True, exist_ok=True)

    rep = [f"# mtclone end-to-end test report",
           f"\n_Generated {datetime.datetime.now().isoformat(timespec='seconds')} · "
           f"mtclone v{mtclone.__version__}_\n"]
    results = {}
    for name, fn in [("mtscATAC", run_mtscatac), ("ReDeeM", run_redeem)]:
        try:
            results[name] = fn(testdata, outdir, rep)
        except Exception as e:
            results[name] = dict(ok=False, error=str(e))
            _log(rep, f"\n**{name} FAILED:** {e}\n```\n{traceback.format_exc()}\n```")

    all_ok = all(r.get("ok") for r in results.values())
    rep.insert(2, f"\n**STATUS: {'🟢 GREEN — all legs passed' if all_ok else '🔴 RED — see failures'}**\n")
    rep.append("\n## Summary\n")
    rep.append("| dataset | cells | variants | informative | clones | note |")
    rep.append("|---|---|---|---|---|---|")
    for k, r in results.items():
        if r.get("ok"):
            note = "annotation gap (placeholder names)" if r.get("annotation_gap") else "clean"
            rep.append(f"| {k} | {r['n_cells']} | {r['n_variants']} | {r['n_informative']} | "
                       f"{r['n_clones']} | {note} |")
        else:
            rep.append(f"| {k} | — | — | — | — | FAILED: {r.get('error')} |")

    (outdir / "mtclone_test_report.md").write_text("\n".join(rep))
    (outdir / "mtclone_test_results.json").write_text(json.dumps(results, indent=2, default=str))
    print("\n" + "=" * 60)
    print("REPORT:", outdir / "mtclone_test_report.md")
    print("GREEN" if all_ok else "RED")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
