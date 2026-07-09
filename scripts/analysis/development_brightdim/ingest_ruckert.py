#!/usr/bin/env python
"""ingest_ruckert.py — build schema-conforming AnnData for the Ruckert ASAP-seq arm.

For each ASAP-seq sample (mtASAP2/3/5) we have, after export_ruckert_mgatk.R:
  <prefix>.heteroplasmy.tsv.gz  (long: cell, variant, heteroplasmy, strand_cor)
  <prefix>.coverage.tsv.gz      (cell, coverage  — real mgatk per-cell mean depth)
plus the deposited demux + protein files:
  *_hashtags.csv.gz / *_donor.csv.gz  (barcode -> donor: CMVpos/neg individuals)
  *_ADT_{counts.txt,barcodes.tsv,tags.tsv}.gz  (surface protein, incl CD56 & CD16)

This driver ingests heteroplasmy+coverage through mtclone.io.read_mtscatac_heteroplasmy
(assay tag stays 'mtscATAC'; Ruckert is mtscATAC/ASAP), stamps per-cell DONOR from the
demux (clones are donor-private), and stashes the surface-protein matrix in
obsm['protein'] + uns['protein_names'] so classify can gate NK and bright/dim on real
CD56/CD16 protein. Writes one .h5ad per sample.
"""
import sys, os, gzip
import numpy as np
import pandas as pd
import scipy.io, scipy.sparse as sp

MTCLONE = "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts"
sys.path.insert(0, MTCLONE)
from mtclone import io

RAW = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/raw/GSE197008_ruckert"
PROC = "/Users/dilloncorvino/Documents/HPC_data/Thymic_NK_development/processed/GSE197008"

# sample -> (mgatk prefix, demux file, demux kind, ADT stem)
SAMPLES = {
    "mtASAP2": dict(prefix="mtASAP2", demux="GSM5906759_mtASAP2_hashtags.csv.gz",
                    adt="GSM5906759_mtASAP2"),
    "mtASAP3": dict(prefix="mtASAP3", demux="GSM5906760_mtASAP3_hashtags.csv.gz",
                    adt="GSM5906760_mtASAP3"),
    "mtASAP5": dict(prefix="mtASAP5", demux="GSM6413442_mtASAP5_donor.csv.gz",
                    adt="GSM6413442_mtASAP5"),
}


def load_demux(path):
    d = pd.read_csv(path)
    d.columns = ["barcode", "donor"]
    return dict(zip(d["barcode"], d["donor"]))


def load_protein(stem):
    """ADT: MatrixMarket (proteins x cells) + barcodes + tags -> (df cells x proteins)."""
    M = scipy.io.mmread(f"{RAW}/{stem}_ADT_counts.txt.gz").tocsr()      # proteins x cells
    bc = pd.read_csv(f"{RAW}/{stem}_ADT_barcodes.tsv.gz", header=None)[0].tolist()
    tags = pd.read_csv(f"{RAW}/{stem}_ADT_tags.tsv.gz", header=None)[0].tolist()
    df = pd.DataFrame(M.T.toarray(), index=bc, columns=tags)            # cells x proteins
    return df


def ingest_sample(name, cfg):
    het = pd.read_csv(f"{PROC}/{cfg['prefix']}.heteroplasmy.tsv.gz", sep="\t")
    cov = pd.read_csv(f"{PROC}/{cfg['prefix']}.coverage.tsv.gz", sep="\t").set_index("cell")["coverage"]
    demux = load_demux(f"{RAW}/{cfg['demux']}")

    a = io.read_mtscatac_heteroplasmy(
        het[["cell", "variant", "heteroplasmy"]],
        dataset="GSE197008", sample=name, donor="MULTIPLEXED",
        tissue="blood", site="blood",
        coverage=cov, coverage_source="mgatk_depth",
    )
    # per-cell donor from demux (strip the sample_ prefix mtclone added to cell_id)
    raw_bc = [c.split("_", 1)[1] if c.startswith(name + "_") else c for c in a.obs_names]
    a.obs["barcode"] = raw_bc
    a.obs["donor"] = [demux.get(b, "unassigned") for b in raw_bc]
    a = a[a.obs["donor"] != "unassigned"].copy()          # drop undemuxed

    # surface protein -> obsm['protein'] aligned to cells
    prot = load_protein(cfg["adt"])
    prot = prot.reindex([b for b in a.obs["barcode"]])
    a.obsm["protein"] = prot.to_numpy(dtype=np.float32)
    a.uns["protein_names"] = list(prot.columns)

    out = f"{PROC}/{name}.h5ad"
    a.write_h5ad(out)
    print(f"[{name}] cells={a.n_obs} variants={a.n_vars} donors={dict(a.obs['donor'].value_counts())}")
    print(f"[{name}] proteins={a.uns['protein_names']}")
    print(f"[{name}] schema problems: {io.validate_schema(a, strict=False)}")
    print(f"[{name}] -> {out}")
    return a


if __name__ == "__main__":
    which = sys.argv[1:] or list(SAMPLES)
    for name in which:
        ingest_sample(name, SAMPLES[name])
