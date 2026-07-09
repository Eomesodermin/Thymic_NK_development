"""Unit tests on synthetic matrices with KNOWN clone structure.

The goal is to prove each stage does what it claims on data whose answer we control:
 - disjoint clone sets  -> zero cross-group sharing, permutation p at the expected tail;
 - one merged clone set  -> high cross-group sharing.
Run: pytest -q  (from scripts/mtclone/, or `pytest scripts/mtclone`).
"""
import os, sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mtclone
from mtclone import io, qc, clones, metrics, classify


# ------------------------------------------------------------------ synthetic builders
def _build(n_clones, cells_per_clone, variants_per_clone, n_donors=1,
           het=0.3, subset_assignment=None, seed=0):
    """Build a schema-conforming AnnData where each clone carries a private variant block.

    Returns (adata, truth_clone_labels). subset_assignment: optional dict clone->group to
    stamp obs['group'] (for sharing tests).
    """
    rng = np.random.default_rng(seed)
    n_cells = n_clones * cells_per_clone
    n_vars = n_clones * variants_per_clone
    X = np.zeros((n_cells, n_vars), dtype=np.float32)
    truth = np.empty(n_cells, dtype=int)
    groups = np.empty(n_cells, dtype=object)
    donors = np.empty(n_cells, dtype=object)
    ci = 0
    for c in range(n_clones):
        vblock = slice(c * variants_per_clone, (c + 1) * variants_per_clone)
        for _ in range(cells_per_clone):
            X[ci, vblock] = het * (rng.random(variants_per_clone) > 0.1)  # mostly-present
            truth[ci] = c
            groups[ci] = (subset_assignment or {}).get(c, "A")
            donors[ci] = f"D{c % n_donors}"
            ci += 1

    var = pd.DataFrame({
        "variant_id": [f"chrM:{1000+i}:A>G" for i in range(n_vars)],
        "pos": [1000 + i for i in range(n_vars)],
        "ref": "A", "alt": "G",
    })
    obs = pd.DataFrame({
        "cell_id": [f"cell{i}" for i in range(n_cells)],
        "donor": donors, "sample": "s1", "tissue": "t", "site": "s",
        "group": groups,
    })
    a = io._finalize(X, obs, var, assay="mtscATAC", dataset="SYNTH",
                     coverage=np.full(n_cells, 50.0, np.float32), coverage_source="depth_file")
    a.obs["group"] = groups
    a.obs["truth"] = truth
    return a, truth


# ------------------------------------------------------------------ schema / io
def test_schema_validates():
    a, _ = _build(3, 10, 4)
    assert io.validate_schema(a, strict=True) == []


def test_schema_catches_bad_range():
    a, _ = _build(2, 5, 3)
    a.X = sp.csr_matrix(a.X.toarray() * 5.0)  # values >1
    with pytest.raises(ValueError):
        io.validate_schema(a, strict=True)


def test_variant_id_parsing():
    for s, exp in [("3243A>G", (3243, "A", "G")), ("chrM:3243:A:G", (3243, "A", "G")),
                   ("m.3243A>G", (3243, "A", "G")), ("3243_A_G", (3243, "A", "G"))]:
        assert io.parse_variant_string(s) == exp


def test_validator_catches_varname_mismatch():
    a, _ = _build(2, 5, 3)
    a.var["variant_id"] = ["x"] * a.n_vars  # break var_names==variant_id
    with pytest.raises(ValueError):
        io.validate_schema(a, strict=True)


# ------------------------------------------------------------------ qc
def test_informative_selection_drops_ubiquitous():
    a, _ = _build(3, 10, 4, het=0.3)
    # inject a ubiquitous (germline-like) variant present in all cells at high het
    Xd = a.X.toarray()
    ubiq = np.full((a.n_obs, 1), 0.9, np.float32)
    import numpy as _np
    Xd = _np.hstack([Xd, ubiq])
    newvar = pd.concat([a.var, pd.DataFrame({
        "variant_id": ["chrM:9999:A>G"], "pos": [9999], "ref": ["A"], "alt": ["G"]},
        index=["chrM:9999:A>G"])])
    a2 = io._finalize(Xd, a.obs.copy(), newvar.reset_index(drop=True),
                      assay="mtscATAC", dataset="SYNTH",
                      coverage=a.obs["coverage"].values, coverage_source="depth_file")
    sel = qc.select_informative_variants(a2, max_pseudobulk_het=0.5, min_cells_detected=2)
    assert "chrM:9999:A>G" not in set(sel.var_names)   # ubiquitous dropped
    assert sel.n_vars > 0
    assert "binary" in sel.layers


# ------------------------------------------------------------------ clones
def test_graph_clones_recover_truth():
    a, truth = _build(4, 12, 5, het=0.3)
    sel = qc.select_informative_variants(a, max_pseudobulk_het=0.9, min_cells_detected=2)
    called = clones.call_clones(sel, method="graph", min_shared_variants=2,
                                edge_weight_cutoff=0.3, min_clone_size=2)
    from sklearn.metrics import adjusted_rand_score
    m = called.obs["clone_id"].values >= 0
    ari = adjusted_rand_score(truth[m], called.obs["clone_id"].values[m])
    assert ari > 0.9, f"ARI too low: {ari}"


def test_variant_group_caller_runs():
    a, truth = _build(3, 10, 4)
    sel = qc.select_informative_variants(a, max_pseudobulk_het=0.9, min_cells_detected=2)
    called = clones.call_clones(sel, method="variant_group", min_clone_size=2)
    assert called.obs["clone_id"].max() >= 0


# ------------------------------------------------------------------ metrics
def test_disjoint_groups_zero_sharing():
    # two groups on DISJOINT clones -> mixed fraction 0
    assign = {0: "A", 1: "A", 2: "B", 3: "B"}
    a, truth = _build(4, 12, 5, subset_assignment=assign)
    sel = qc.select_informative_variants(a, max_pseudobulk_het=0.9, min_cells_detected=2)
    called = clones.call_clones(sel, method="graph", min_shared_variants=2,
                                edge_weight_cutoff=0.3, min_clone_size=2)
    res = metrics.between_vs_within_sharing(called, "A", "B", group_key="group")
    assert res["frac_mixed_clones"] == 0.0


def test_merged_group_high_sharing():
    # both groups distributed WITHIN every clone -> high mixed fraction
    a, truth = _build(4, 20, 5)
    # assign group per-cell alternating within each clone
    grp = np.array(["A" if i % 2 == 0 else "B" for i in range(a.n_obs)], dtype=object)
    a.obs["group"] = grp
    sel = qc.select_informative_variants(a, max_pseudobulk_het=0.9, min_cells_detected=2)
    called = clones.call_clones(sel, method="graph", min_shared_variants=2,
                                edge_weight_cutoff=0.3, min_clone_size=2)
    res = metrics.between_vs_within_sharing(called, "A", "B", group_key="group")
    assert res["frac_mixed_clones"] > 0.8


def test_gini_shannon_bounds():
    assert metrics.gini([10, 10, 10]) < 0.05          # even
    assert metrics.gini([100, 1, 1]) > 0.5            # skewed
    assert abs(metrics.normalized_shannon([5, 5, 5, 5]) - 1.0) < 1e-6  # max even


def test_permutation_null_calibrated():
    # disjoint groups: observed mixed=0, null (shuffling within donor) should also be ~0
    # -> use a case with real structure: groups disjoint by clone but same donor so shuffle mixes
    assign = {0: "A", 1: "A", 2: "B", 3: "B"}
    a, truth = _build(4, 15, 5, n_donors=1, subset_assignment=assign)
    sel = qc.select_informative_variants(a, max_pseudobulk_het=0.9, min_cells_detected=2)
    called = clones.call_clones(sel, method="graph", min_shared_variants=2,
                                edge_weight_cutoff=0.3, min_clone_size=2)
    stat = lambda ad_: metrics.between_vs_within_sharing(
        ad_, "A", "B", group_key="group")["frac_mixed_clones"]
    null = metrics.permutation_null(called, stat, group_key="group",
                                    stratify_by="donor", n=200, seed=1)
    # observed sharing (0) should be LOWER than shuffled null (which mixes A/B across clones)
    assert null["observed"] <= null["null_mean"]
    assert null["p_less"] < 0.10, f"expected significant depletion, got {null}"


# ------------------------------------------------------------------ classify (marker logic)
def test_bright_dim_marker_split():
    rng = np.random.default_rng(0)
    # FCGR3A is already in DIM_MARKERS; keep gene list unique
    genes = list(dict.fromkeys(classify.BRIGHT_MARKERS + classify.DIM_MARKERS + ["FCGR3A"]))
    n = 40
    X = rng.normal(size=(n, len(genes))).astype(np.float32)
    # first 20 cells: high bright markers; last 20: high dim markers + FCGR3A
    for gi, g in enumerate(genes):
        if g in classify.BRIGHT_MARKERS:
            X[:20, gi] += 3
        if g in classify.DIM_MARKERS or g == "FCGR3A":
            X[20:, gi] += 3
    a = ad.AnnData(X=X, var=pd.DataFrame(index=genes),
                   obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]))
    a.obs["cell_type"] = "NK"
    classify.label_bright_dim(a, method="marker")
    lab = a.obs["nk_subset"].values
    assert (lab[:20] == "bright").mean() > 0.8
    assert (lab[20:] == "dim").mean() > 0.8
