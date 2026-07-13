"""Step 2 — ingest GSE302113 per-cell heteroplasmy into per-donor AnnData (Plan 3).

Reads the author-called per-cell heteroplasmy + variant_stats TSVs (downloaded to RAW) via
`mtclone.io.read_mtscatac_heteroplasmy`, concatenates per donor across sites (outer join, missing
variants -> 0), recomputes joint pseudobulk stats, restores the mtclone schema, and writes one
<donor>.h5ad. Fragments are NOT needed here (they are Step-3 NK-labeling input).

Inputs:
  RAW/<GSM>...cell_heteroplasmic_df.tsv.gz, RAW/<GSM>...variant_stats.tsv.gz
  download_manifest.json (gsm, donor, dx, site, tissue, sort, het, vs per GSM) — saved as artifact.

Usage:
    python ingest_gse302113.py --raw <RAW> --manifest download_manifest.json --out <processed>
"""
from __future__ import annotations
import argparse, os, sys, json, warnings, collections
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad

sys.path.insert(0, "/Users/dilloncorvino/Documents/Github/Eomesodermin/Thymic_NK_development/scripts")
import mtclone
from mtclone import io


def recompute_var_stats(a):
    het = a.X
    a.var["pseudobulk_heteroplasmy"] = np.asarray(het.mean(axis=0)).ravel().astype(np.float32)
    a.var["n_cells_detected"] = np.asarray((het > 0).sum(axis=0)).ravel().astype(int)
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    manifest = json.load(open(args.manifest))

    by_donor = collections.defaultdict(list)
    sample_log = []
    for m in manifest:
        het = os.path.join(args.raw, m["het"]); vs = os.path.join(args.raw, m["vs"])
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                a = io.read_mtscatac_heteroplasmy(het, vs, dataset="GSE302113", sample=m["gsm"],
                                                  donor=m["donor"], tissue=m["tissue"], site=m["site"])
            a.obs["sort"] = m["sort"]; a.obs["diagnosis"] = m["dx"]
            by_donor[m["donor"]].append(a)
            sample_log.append(dict(gsm=m["gsm"], donor=m["donor"], site=m["site"],
                                   sort=m["sort"], n_cells=a.n_obs, n_vars=a.n_vars))
        except Exception as e:
            sample_log.append(dict(gsm=m["gsm"], donor=m["donor"], error=str(e)))
            print("ERR", m["gsm"], e, flush=True)

    from mtclone.io import parse_variant_string
    for donor, objs in by_donor.items():
        a = ad.concat(objs, join="outer", label="_src", fill_value=0.0, index_unique=None)
        a.X = sp.csr_matrix(a.X, dtype=np.float32); a.X.data = np.clip(a.X.data, 0.0, 1.0)
        a.var["variant_id"] = a.var_names.astype(str)
        parsed = [parse_variant_string(v) for v in a.var_names]
        a.var["pos"] = [p[0] if p else -1 for p in parsed]
        a.var["ref"] = [p[1] if p else "N" for p in parsed]
        a.var["alt"] = [p[2] if p else "N" for p in parsed]
        a = recompute_var_stats(a)
        if "strand_correlation" not in a.var: a.var["strand_correlation"] = np.nan
        a.obs["cell_id"] = a.obs_names.astype(str)
        a.uns["mtclone_schema_version"] = mtclone.io.SCHEMA_VERSION
        mtclone.io.validate_schema(a, strict=False)
        for c in ["cell_id","donor","sample","tissue","site","assay","dataset",
                  "coverage_source","sort","diagnosis"]:
            if c in a.obs: a.obs[c] = a.obs[c].astype(str)
        a.write_h5ad(os.path.join(args.out, f"{donor}.h5ad"))
        print(f"{donor}: {a.n_obs} cells, {a.n_vars} vars, sites={sorted(a.obs['site'].unique())}", flush=True)

    pd.DataFrame(sample_log).to_csv(os.path.join(args.out, "sample_ingest_log.csv"), index=False)


if __name__ == "__main__":
    main()
