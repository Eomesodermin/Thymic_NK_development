"""
harvest_footprint.py — pool all cellranger vdj outputs, join to the Domínguez Conde
atlas NK annotations by barcode, and compute the TCR-footprint rate in
CD56-bright vs CD56-dim NK.  This is the Q1/Q2 test.

Inputs:
  - VDJ outputs: <work>/vdj_<lib>/outs/all_contig_annotations.csv (one per library)
  - harvest_map.json: [{vdj_lib, gex_prefix, atlas_prefix, join_tier, kind, sample, run}]
  - atlas H5AD: tildc.h5ad (T & ILC compartment, cell_type carries ontology NK labels)

Join key (validated on real data, 75.6% match):
  VDJ barcode  '<16bp>-1'  ->  strip '-1'  ->  '<atlas_prefix>_<16bp>'
  direct tier   : atlas_prefix = GEX library id (Pan_T####)
  czi_pooled    : atlas_prefix = paired prefix 'CZI-IAaaa+CZI-IAbbb' (first element = this GEX id)

Footprint definitions (from footprint_pipeline.py):
  nonprod_only = >=1 non-productive TCR contig AND no productive contig  <- aborted-relic target
"""
import pandas as pd
import numpy as np
import h5py
import json
import os
import glob
import sys

sys.path.insert(0, os.path.dirname(__file__))
from footprint_pipeline import per_cell_footprint, TCR_CHAINS  # noqa: E402

# Classify NK subsets by substring pattern, NOT exact-string allow-lists — the atlas
# ontology labels carry a ', human' suffix (e.g. 'CD16-negative, CD56-bright natural
# killer cell, human') and suffix variants differ; substring matching is robust to that.
def nk_class_of(label):
    """bright / dim / other from a cell_type ontology string."""
    if not isinstance(label, str):
        return "other"
    L = label.lower()
    if "natural killer" not in L:
        return "other"
    if "cd56-bright" in L or "cd56 bright" in L:
        return "bright"
    if "cd56-dim" in L or "cd56 dim" in L:
        return "dim"
    return "other"


def load_atlas_nk(h5ad_path):
    """Return DataFrame [key, cell_type, nk_class] for atlas cells, key='<prefix>_<16bp>'."""
    f = h5py.File(h5ad_path, "r")
    idx = np.array([x.decode() for x in f["obs"]["_index"][:]])
    ct = f["obs"]["cell_type"]
    cats = np.array([x.decode() for x in ct["categories"][:]])
    cell_type = cats[ct["codes"][:]]
    df = pd.DataFrame({"key": idx, "cell_type": cell_type})
    df["nk_class"] = df.cell_type.map(nk_class_of)
    return df


def harvest(work_dir, harvest_map_path, h5ad_path, out_prefix="q1"):
    hmap = json.load(open(harvest_map_path))
    atlas = load_atlas_nk(h5ad_path)
    atlas_keys = set(atlas.key)

    percell_all = []
    per_lib_rows = []
    for h in hmap:
        lib = h["vdj_lib"]
        # accept three layouts: cellranger run dir (vdj_<lib>/outs/...),
        # smoke-test dir (smoke_<lib>/outs/...), and the flattened harvest_export
        # tarball layout (<lib>__all_contig_annotations.csv in work_dir).
        candidates = [
            os.path.join(work_dir, f"vdj_{lib}", "outs", "all_contig_annotations.csv"),
            os.path.join(work_dir, f"smoke_{lib}", "outs", "all_contig_annotations.csv"),
            os.path.join(work_dir, f"{lib}__all_contig_annotations.csv"),
        ]
        fp_csv = next((c for c in candidates if os.path.exists(c)), None)
        if fp_csv is None:
            per_lib_rows.append({"vdj_lib": lib, "status": "MISSING"})
            continue
        contigs = pd.read_csv(fp_csv)
        contigs["sample"] = lib
        pc = per_cell_footprint(contigs)  # raw all-contig, no is_cell/high_conf filter
        # build join key
        pc["key"] = h["atlas_prefix"] + "_" + pc["barcode"].str.split("-").str[0]
        pc["vdj_lib"] = lib
        pc["join_tier"] = h["join_tier"]
        pc["matched"] = pc["key"].isin(atlas_keys)
        percell_all.append(pc)
        per_lib_rows.append({"vdj_lib": lib, "status": "OK",
                             "cells": len(pc), "matched": int(pc["matched"].sum()),
                             "tier": h["join_tier"]})

    percell = pd.concat(percell_all, ignore_index=True)
    # join to atlas NK class
    percell = percell.merge(atlas[["key", "cell_type", "nk_class"]], on="key", how="left")

    # ---- Q1/Q2: footprint rate in bright vs dim NK ----
    nk = percell[percell.nk_class.isin(["bright", "dim"])].copy()
    # NOTE: nk here counts only NK cells that carry >=1 TCR contig (the join is
    # keyed on VDJ barcodes). For a rate you must divide by the full aligned-NK
    # denominator (atlas NK from the aligned library prefixes), computed below.
    aligned_prefixes = set(h["atlas_prefix"] for h in hmap)
    atlas_pref = atlas.assign(prefix=atlas.key.str.rsplit("_", n=1).str[0])
    denom = (atlas_pref[atlas_pref.prefix.isin(aligned_prefixes)
             & atlas_pref.nk_class.isin(["bright", "dim"])]
             .nk_class.value_counts().to_dict())
    summ = nk.groupby("nk_class").agg(
        n_nk_with_tcr_contig=("key", "size"),
        n_nonprod_only=("nonprod_only", "sum"),
        n_has_nonprod=("has_nonproductive", "sum"),
        n_has_prod=("has_productive", "sum"),
    ).reset_index()
    summ["n_aligned_nk"] = summ.nk_class.map(denom)
    summ["footprint_rate"] = summ.n_nonprod_only / summ.n_aligned_nk

    pd.DataFrame(per_lib_rows).to_csv(f"{out_prefix}_per_library.csv", index=False)
    nk.to_csv(f"{out_prefix}_nk_percell.csv", index=False)
    summ.to_csv(f"{out_prefix}_bright_vs_dim.csv", index=False)
    return summ, nk, pd.DataFrame(per_lib_rows)


if __name__ == "__main__":
    work = sys.argv[1] if len(sys.argv) > 1 else "."
    hm = sys.argv[2] if len(sys.argv) > 2 else "harvest_map.json"
    h5 = sys.argv[3] if len(sys.argv) > 3 else "tildc.h5ad"
    summ, nk, perlib = harvest(work, hm, h5)
    print(summ.to_string(index=False))
