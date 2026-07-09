"""Unit test: read_redeem_consensus on a synthetic Consensus.final fixture.

Builds a tiny RawGenotypes.Sensitive.StrandBalance + QualifiedTotalCts with KNOWN
heteroplasmy so we can assert exact fractions (Freq/Depth), correct orientation
(cells x variants), canonical variant ids, real per-cell coverage, and schema validity.
"""
import os, sys
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mtclone import io


# 14-col layout: UMI Cell Pos Variants Call Ref FamSize GT_Cts CSS DB_Cts SG_Cts Plus Minus Depth
def _row(umi, cell, pos, variant, depth):
    ref, alt = variant.split("_")[1], variant.split("_")[2]
    return [umi, cell, pos, variant, alt, ref, 3, 3, 1, 2, 1, 2, 1, depth]


def _write_fixture(tmp):
    d = tmp / "Young1.T1.BMMC.Consensus.final"
    d.mkdir(parents=True)
    rows = []
    # cellA: variant 3243_A_G with 3 alt-UMIs at depth 30  -> hetero = 3/30 = 0.1
    for u in range(3):
        rows.append(_row(f"u{u}", "cellA", 3243, "3243_A_G", 30))
    # cellA: variant 5000_C_T with 5 alt-UMIs at depth 10   -> hetero = 5/10 = 0.5
    for u in range(5):
        rows.append(_row(f"v{u}", "cellA", 5000, "5000_C_T", 10))
    # cellB: variant 3243_A_G with 6 alt-UMIs at depth 12    -> hetero = 6/12 = 0.5
    for u in range(6):
        rows.append(_row(f"w{u}", "cellB", 3243, "3243_A_G", 12))
    # cellB: variant 8000_G_A with 20 alt-UMIs at depth 10   -> raw 2.0, clipped to 1.0
    for u in range(20):
        rows.append(_row(f"x{u}", "cellB", 8000, "8000_G_A", 10))
    pd.DataFrame(rows).to_csv(d / "RawGenotypes.Sensitive.StrandBalance",
                              sep="\t", header=False, index=False)
    # QualifiedTotalCts: cell, pos, coverage
    qtc = [["cellA", 3243, 30], ["cellA", 5000, 10],
           ["cellB", 3243, 12], ["cellB", 8000, 10]]
    pd.DataFrame(qtc).to_csv(d / "QualifiedTotalCts", sep="\t", header=False, index=False)
    return d


def test_redeem_consensus_heteroplasmy_exact(tmp_path):
    d = _write_fixture(tmp_path)
    a = io.read_redeem_consensus(d, thr="S", donor="Young1", sample="Young1_BMMC",
                                 tissue="BMMC", site="marrow")
    # orientation + shape: 2 cells x 3 variants
    assert a.shape == (2, 3)
    assert io.validate_schema(a, strict=True) == []
    # canonical variant ids
    assert set(a.var_names) == {"chrM:3243:A>G", "chrM:5000:C>T", "chrM:8000:G>A"}
    X = a.to_df()
    ca = "Young1_BMMC_cellA"; cb = "Young1_BMMC_cellB"
    assert np.isclose(X.loc[ca, "chrM:3243:A>G"], 0.1, atol=1e-6)
    assert np.isclose(X.loc[ca, "chrM:5000:C>T"], 0.5, atol=1e-6)
    assert np.isclose(X.loc[cb, "chrM:3243:A>G"], 0.5, atol=1e-6)
    assert np.isclose(X.loc[cb, "chrM:8000:G>A"], 1.0, atol=1e-6)   # clipped from 2.0
    assert X.loc[ca, "chrM:8000:G>A"] == 0.0                        # absent -> 0


def test_redeem_consensus_real_coverage(tmp_path):
    d = _write_fixture(tmp_path)
    a = io.read_redeem_consensus(d, thr="S", donor="Young1", sample="Young1_BMMC")
    assert (a.obs["coverage_source"] == "depth_file").all()
    # cellA mean coverage = mean(30, 10) = 20 ; cellB = mean(12, 10) = 11
    cov = a.obs.set_index("cell_id")["coverage"]
    assert np.isclose(cov["Young1_BMMC_cellA"], 20.0, atol=1e-4)
    assert np.isclose(cov["Young1_BMMC_cellB"], 11.0, atol=1e-4)


def test_redeem_consensus_whitelist(tmp_path):
    d = _write_fixture(tmp_path)
    a = io.read_redeem_consensus(d, thr="S", cell_whitelist=["cellA"])
    assert a.n_obs == 1
    assert a.obs["cell_id"].iloc[0].endswith("cellA")


def test_redeem_consensus_missing_thr(tmp_path):
    d = _write_fixture(tmp_path)
    with pytest.raises(FileNotFoundError):
        io.read_redeem_consensus(d, thr="VS")   # Specific file not written


def test_redeem_consensus_barcode_translation(tmp_path):
    # ATAC-space folder barcodes cellA/cellB -> RNA-space rnaA/rnaB via map.
    d = _write_fixture(tmp_path)
    bmap = {"cellA": "rnaA", "cellB": "rnaB", "other": "rnaX"}
    a = io.read_redeem_consensus(d, thr="S", sample="S", barcode_map=bmap,
                                 barcode_space="atac")
    # returned cell_ids should be in RNA space
    assert set(x.split("_", 1)[1] for x in a.obs["cell_id"]) == {"rnaA", "rnaB"}
    # heteroplasmy preserved through translation
    X = a.to_df()
    ra = "S_rnaA"
    assert np.isclose(X.loc[ra, "chrM:3243:A>G"], 0.1, atol=1e-6)
    # coverage still real (translated barcodes)
    assert (a.obs["coverage_source"] == "depth_file").all()
    cov = a.obs.set_index("cell_id")["coverage"]
    assert np.isclose(cov["S_rnaA"], 20.0, atol=1e-4)


def test_redeem_consensus_barcode_translation_df(tmp_path):
    # map as a DataFrame with atac/rna columns
    d = _write_fixture(tmp_path)
    bm = pd.DataFrame({"atac": ["cellA", "cellB"], "rna": ["rnaA", "rnaB"]})
    a = io.read_redeem_consensus(d, thr="S", sample="S", barcode_map=bm)
    assert a.n_obs == 2
    assert set(x.split("_", 1)[1] for x in a.obs["cell_id"]) == {"rnaA", "rnaB"}
