#!/usr/bin/env python
"""ingest_redeem.py — build schema-conforming AnnData for the Weng ReDeeM arm of Plan 2.

For each donor we join two sources the GEO deposit could not:
  1. `Consensus.final/` folder (Figshare 24418966) -> per-cell heteroplasmy + real depth,
     via mtclone.io.read_redeem_consensus. Barcodes here are in 10x-Multiome ATAC space.
  2. Annotated Multiome Seurat object (Figshare 23290004), exported by export_redeem_seurat.R
     to <prefix>.{cellmeta,nk_markers,wnn_umap}.tsv.gz. Barcodes here are in RNA space.

The ATAC<->RNA translation (redeemR ATACWhite/RNAWhite, 736,320 paired barcodes) is applied
inside read_redeem_consensus(barcode_map=...) so heteroplasmy cell_ids come out in RNA space,
matching the Seurat annotation. We then attach STD.CellType (NK identity), the bright/dim
marker panel, meanCov, and the authors' ClonalGroup reference; optionally restrict to NK.

Writes one .h5ad per donor (full marrow hierarchy, tagged) — NK subsetting is a downstream
classify step, but STD.CellType/markers are attached so it is one selection away.

Usage: ingest_redeem.py <consensus_dir> <seurat_prefix> <out.h5ad> --donor D --sample S [--nk-only]
"""
import sys, argparse
import numpy as np
import pandas as pd

MTCLONE = "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts"
sys.path.insert(0, MTCLONE)
from mtclone import io

WHITELIST_MAP = ("/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/raw/"
                 "GSE219014_redeem/whitelists/multiome_barcode_map.tsv.gz")
# Seurat sample-suffix -> tissue: -1 BMMC (where NK/mature live), -2 HSPC, -3 HSC
SUFFIX_TISSUE = {"1": "BMMC", "2": "HSPC", "3": "HSC"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("consensus_dir")
    ap.add_argument("seurat_prefix")
    ap.add_argument("out")
    ap.add_argument("--donor", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--thr", default="S")
    ap.add_argument("--nk-only", action="store_true")
    a = ap.parse_args()

    # ---- annotation (RNA-space barcodes)
    cm = pd.read_csv(f"{a.seurat_prefix}.cellmeta.tsv.gz", sep="\t").set_index("barcode")
    mk = pd.read_csv(f"{a.seurat_prefix}.nk_markers.tsv.gz", sep="\t").set_index("barcode")
    try:
        um = pd.read_csv(f"{a.seurat_prefix}.wnn_umap.tsv.gz", sep="\t").set_index("barcode")
    except FileNotFoundError:
        um = None

    # ---- heteroplasmy (ATAC-space -> translated to RNA space so it matches the Seurat annot)
    bmap = pd.read_csv(WHITELIST_MAP, sep="\t")
    adata = io.read_redeem_consensus(
        a.consensus_dir, thr=a.thr, dataset="GSE219014",
        sample=a.sample, donor=a.donor, tissue="BMMC", site="marrow",
        barcode_map=bmap, barcode_space="atac",
    )
    # obs_names carry the (already RNA-translated) barcode, optionally prefixed with "<sample>_".
    # Sample names contain underscores, so strip the exact prefix rather than split("_", 1).
    pref = a.sample + "_"
    raw = [c[len(pref):] if c.startswith(pref) else c for c in adata.obs_names]
    adata.obs["barcode_rna"] = raw

    # the Consensus barcodes have NO sample suffix; Seurat barcodes do (-1/-2/-3). The BMMC
    # consensus folder maps to the -1 suffix cells. Attach annotation by RNA barcode + suffix.
    cm_nosuf = cm.copy()
    cm_nosuf["suffix"] = cm_nosuf.index.str.rsplit("-", n=1).str[-1]
    cm_nosuf["barcode_rna"] = cm_nosuf.index.str.rsplit("-", n=1).str[0]
    # keep only the BMMC (-1) rows for a BMMC consensus folder
    cm_bmmc = cm_nosuf[cm_nosuf["suffix"] == "1"].set_index("barcode_rna")
    mk_bmmc = mk.copy()
    mk_bmmc["barcode_rna"] = mk_bmmc.index.str.rsplit("-", n=1).str[0]
    mk_bmmc = mk_bmmc[mk.index.str.endswith("-1")].set_index("barcode_rna")

    matched = adata.obs["barcode_rna"].isin(cm_bmmc.index)
    print(f"[redeem] {adata.n_obs} consensus cells; {matched.sum()} match Seurat BMMC annotation")
    adata = adata[matched.values].copy()

    idx = adata.obs["barcode_rna"].values
    for col in ["STD.CellType", "STD_Cat", "STD_Cat2", "meanCov", "ClonalGroup", "ClonalGroup.Prob"]:
        if col in cm_bmmc.columns:
            adata.obs[col] = cm_bmmc.loc[idx, col].values
    # bright/dim + NK marker panel into obsm
    mk_al = mk_bmmc.reindex(idx)
    adata.obsm["rna_markers"] = mk_al.to_numpy(dtype=np.float32)
    adata.uns["rna_marker_names"] = list(mk_bmmc.columns)
    if um is not None:
        um2 = um.copy(); um2["barcode_rna"] = um2.index.str.rsplit("-", n=1).str[0]
        um2 = um2[um.index.str.endswith("-1")].set_index("barcode_rna").reindex(idx)
        adata.obsm["wnn_umap"] = um2[["wnnUMAP_1", "wnnUMAP_2"]].to_numpy(dtype=np.float32)

    if a.nk_only:
        nk = adata.obs["STD.CellType"] == "NK"
        print(f"[redeem] restricting to NK: {nk.sum()} / {adata.n_obs}")
        adata = adata[nk.values].copy()

    print(f"[redeem] final {adata.n_obs} cells x {adata.n_vars} variants")
    print(f"[redeem] cell types: {dict(adata.obs['STD.CellType'].value_counts().head(8))}")
    print(f"[redeem] schema problems: {io.validate_schema(adata, strict=False)}")
    adata.write_h5ad(a.out)
    print(f"[redeem] -> {a.out}")


if __name__ == "__main__":
    main()
