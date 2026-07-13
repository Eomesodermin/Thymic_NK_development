"""Gene-activity scoring from 10x ATAC fragments for a marker panel (GRCh38).

Plan 3 Step 3 companion to `mtclone.classify.label_nk(modality="atac")`, which scores a
cells x genes activity matrix but does not build one. This does exactly that, for the marker
panel only (NK identification does not need a genome-wide gene-activity matrix).

Model (ArchR-style, simplified): per gene, count fragment insertions (both ends) falling in a
window = gene body extended `upstream` bp past the TSS and `downstream` bp past the TES
(strand-aware). Counts are per called-cell barcode; cells are restricted to the barcodes in
that sample's mtclone heteroplasmy object. Output = one AnnData per sample (cells x genes,
layer 'counts'), plus a concatenated per-donor object. Values are raw insertion counts; the
classifier z-scores per gene internally, so no normalization is required here — but we also
store a depth-normalized log layer for QC.

Usage:
    python gene_activity_from_fragments.py --raw <dir> --het-dir <dir> --panel <json> \
        --manifest <download_manifest.json> --out <dir> [--jobs 8]
"""
from __future__ import annotations
import argparse, gzip, json, os, sys, time
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad


def build_windows(coords, upstream=5000, downstream=1000):
    """chrom -> sorted list of (win_start, win_end, gene) with strand-aware extension."""
    by_chrom = defaultdict(list)
    for g, c in coords.items():
        s, e, strand = c["start"], c["end"], c.get("strand", 1)
        if strand >= 0:                       # + strand: TSS = start
            ws, we = s - upstream, e + downstream
        else:                                  # - strand: TSS = end
            ws, we = s - downstream, e + upstream
        by_chrom[c["chrom"]].append((max(0, ws), we, g))
    for ch in by_chrom:
        by_chrom[ch].sort()
    return dict(by_chrom)


def score_sample(frag_path, called_barcodes, windows, genes):
    """Stream a fragments.tsv.gz, count insertions per (barcode, gene) for called cells."""
    gidx = {g: i for i, g in enumerate(genes)}
    bset = set(called_barcodes)
    counts = defaultdict(lambda: np.zeros(len(genes), dtype=np.float32))
    marker_chroms = set(windows)
    op = gzip.open if str(frag_path).endswith(".gz") else open
    n_lines = 0
    with op(frag_path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            n_lines += 1
            # split only what we need
            p = line.rstrip("\n").split("\t")
            ch = p[0]
            if ch not in marker_chroms:
                continue
            bc = p[3]
            if bc not in bset:
                continue
            start, end = int(p[1]), int(p[2])
            # both Tn5 insertion sites: start and end
            for site in (start, end):
                for ws, we, g in windows[ch]:
                    if site < ws:
                        break            # windows sorted; no later window can contain it... 
                    if ws <= site <= we:
                        counts[bc][gidx[g]] += 1.0
    return counts, n_lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--het-dir", required=True, help="dir with per-donor <donor>.h5ad (for called barcodes)")
    ap.add_argument("--panel", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()

    panel = json.load(open(args.panel))
    coords = panel["coords"]
    genes = sorted(coords)
    windows = build_windows(coords)
    manifest = json.load(open(args.manifest))
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    # called barcodes per SAMPLE (obs_names in donor objects are '<GSM>_<barcode>')
    hetdir = Path(args.het_dir)
    called = defaultdict(set)   # gsm -> {barcode}
    for h in hetdir.glob("*.h5ad"):
        a = ad.read_h5ad(h, backed="r")
        for cid, samp in zip(a.obs_names, a.obs["sample"]):
            bc = str(cid)[len(str(samp)) + 1:] if str(cid).startswith(str(samp)) else str(cid)
            called[str(samp)].add(bc)
    print({k: len(v) for k, v in called.items()}, flush=True)

    from concurrent.futures import ThreadPoolExecutor, as_completed

    def job(m):
        frag = os.path.join(args.raw, m["frag"])
        cts, nl = score_sample(frag, called[m["gsm"]], windows, genes)
        bcs = list(cts)
        if bcs:
            X = sp.csr_matrix(np.vstack([cts[b] for b in bcs]))
        else:
            X = sp.csr_matrix((0, len(genes)), dtype=np.float32)
        a = ad.AnnData(X=X, obs=pd.DataFrame(index=[f"{m['gsm']}_{b}" for b in bcs]),
                       var=pd.DataFrame(index=genes))
        a.obs["sample"] = m["gsm"]; a.obs["donor"] = m["donor"]
        a.obs["site"] = m["site"]; a.obs["sort"] = m["sort"]
        a.layers["counts"] = a.X.copy()
        a.write_h5ad(out / f"{m['gsm']}.geneactivity.h5ad")
        return m["gsm"], a.n_obs, int(nl)

    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(job, m): m["gsm"] for m in manifest}
        for i, f in enumerate(as_completed(futs)):
            gsm, n, nl = f.result()
            results.append(dict(gsm=gsm, n_cells=n, n_frag_lines=nl))
            print(f"[{i+1}/{len(manifest)}] {gsm}: {n} cells, {nl:,} frag lines", flush=True)
            json.dump(results, open(out / "geneactivity_progress.json", "w"), indent=1)

    print("DONE", len(results), "samples")


if __name__ == "__main__":
    main()
