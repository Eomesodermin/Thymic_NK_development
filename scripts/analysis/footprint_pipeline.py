"""
footprint_pipeline.py — per-cell TCR footprint counting from cellranger all-contig tables.

Core reusable unit for the NK developmental-origins project. Applies identically to:
  - E-MTAB-12524 (warm-up / NMD control)
  - Domínguez Conde atlas NK cells (after cellranger vdj re-alignment on Marvin)

KEY LESSON (validated on E-MTAB-12524, 2026-07-07):
  cellranger's is_cell / high_confidence filters REMOVE the non-productive/partial
  contigs that ARE the developmental footprint. Counting from the raw all_contig
  table: 634 nonprod_only cells (6.3%). Counting after is_cell+high_confidence: 0.
  => ALWAYS operate on the raw all-contig table; filters are optional sensitivity knobs.

Definitions (locked):
  has_productive    : cell has >=1 productive TCR contig
  has_nonproductive : cell has >=1 non-productive (out-of-frame / PTC) TCR contig
  nonprod_only      : non-productive present AND no productive  <-- the aborted-relic target
"""
import pandas as pd
import numpy as np
import os
import glob

TCR_CHAINS = ("TRA", "TRB", "TRD", "TRG")

# NOTE (E-MTAB-12524): a sample can appear in the input yet contribute zero rows
# to per_cell_footprint if it is a BCR library (IGH/IGK/IGL only). E-MTAB-12524
# samples 13086148/13086149 are BCR (0 TCR contigs) and are correctly dropped by
# the chain filter — NOT by is_cell/high_confidence. Check chain composition first
# if a sample unexpectedly vanishes from the output.


def _to_bool(x):
    if isinstance(x, bool):
        return x
    if pd.isna(x):
        return np.nan
    s = str(x).strip().lower()
    if s in ("true", "1", "yes"):
        return True
    if s in ("false", "0", "no", "none"):
        return False
    return np.nan


def load_all_contigs(paths, sample_from_name=lambda p: os.path.basename(p).split("-all_contig")[0]):
    """Load one or more cellranger all_contig_annotations.csv into a single frame,
    tagging each with a `sample` column derived from the filename."""
    frames = []
    for p in paths:
        df = pd.read_csv(p)
        df["sample"] = sample_from_name(p)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def per_cell_footprint(contigs, chains=TCR_CHAINS,
                       require_cell=False, require_high_conf=False, min_umis=None):
    """Collapse a cellranger ALL-contig table to one row per (sample, barcode).

    Defaults deliberately apply NO is_cell/high_confidence filter (see module docstring).
    """
    df = contigs.copy()
    df["productive_b"] = df["productive"].map(_to_bool)
    df = df[df["chain"].isin(chains)].copy()
    if require_cell and "is_cell" in df:
        df = df[df["is_cell"] == True]  # noqa: E712
    if require_high_conf and "high_confidence" in df:
        df = df[df["high_confidence"] == True]  # noqa: E712
    if min_umis is not None and "umis" in df:
        df = df[df["umis"] >= min_umis]
    g = df.groupby(["sample", "barcode"])
    out = pd.DataFrame({
        "n_tcr_contigs":   g.size(),
        "n_productive":    g["productive_b"].apply(lambda s: (s == True).sum()),   # noqa: E712
        "n_nonproductive": g["productive_b"].apply(lambda s: (s == False).sum()),  # noqa: E712
        "max_umis":        g["umis"].max() if "umis" in df else np.nan,
        "chains":          g["chain"].apply(lambda s: ",".join(sorted(set(s)))),
    }).reset_index()
    out["has_productive"] = out["n_productive"] > 0
    out["has_nonproductive"] = out["n_nonproductive"] > 0
    out["nonprod_only"] = out["has_nonproductive"] & (~out["has_productive"])
    return out


def filter_sensitivity(contigs, chains=TCR_CHAINS):
    """Report nonprod_only counts across filter levels — the QC-destroys-signal check."""
    rows = []
    for rc, rh, label in [(False, False, "raw all-contig"),
                          (True, False, "is_cell"),
                          (True, True, "is_cell+high_confidence")]:
        fp = per_cell_footprint(contigs, chains=chains, require_cell=rc, require_high_conf=rh)
        rows.append({"filter": label, "cells": len(fp),
                     "nonprod_only": int(fp["nonprod_only"].sum()),
                     "rate": round(fp["nonprod_only"].mean(), 4) if len(fp) else 0.0})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    import sys
    paths = sorted(glob.glob(sys.argv[1])) if len(sys.argv) > 1 else []
    contigs = load_all_contigs(paths)
    print(filter_sensitivity(contigs).to_string(index=False))
