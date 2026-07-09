"""io.py — ingest adapters + schema validation.

Every adapter returns an AnnData conforming to schema.md. Downstream code depends only on
that contract, never on the source format. See schema.md for the full field reference.
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.io
import anndata as ad

SCHEMA_VERSION = "1.0"

_REQUIRED_OBS = [
    "cell_id", "donor", "sample", "tissue", "site", "assay", "dataset",
    "coverage", "coverage_source",
]
_REQUIRED_VAR = [
    "variant_id", "pos", "ref", "alt", "pseudobulk_heteroplasmy", "n_cells_detected",
]


# ----------------------------------------------------------------------------- helpers
def canonical_variant_id(pos, ref, alt) -> str:
    """Normalize to `chrM:POS:REF>ALT` (1-based rCRS, uppercase)."""
    return f"chrM:{int(pos)}:{str(ref).upper()}>{str(alt).upper()}"


_VARIANT_PATTERNS = [
    # mgatk / common: "3243A>G", "3243_A_G", "chrM:3243:A:G", "chrM:3243:A>G"
    re.compile(r"^(?:chrM[:_])?(\d+)[:_]?([ACGTN])\s*[>:_]\s*([ACGTN])$", re.I),
    # ReDeeM style: "Variants" like "3243 A G" or "m.3243A>G"
    re.compile(r"^m?\.?(\d+)\s*([ACGTN])\s*>?\s*([ACGTN])$", re.I),
]


def parse_variant_string(s: str):
    """Return (pos, ref, alt) parsed from a heterogeneous variant string, or None."""
    s = str(s).strip()
    for pat in _VARIANT_PATTERNS:
        m = pat.match(s)
        if m:
            return int(m.group(1)), m.group(2).upper(), m.group(3).upper()
    return None


def _finalize(X, obs: pd.DataFrame, var: pd.DataFrame, *, assay: str, dataset: str,
              coverage=None, coverage_source: str = "derived_from_matrix") -> ad.AnnData:
    """Assemble a schema-conforming AnnData from a heteroplasmy matrix + frames."""
    X = sp.csr_matrix(X, dtype=np.float32)
    # clip to [0,1] defensively
    if X.max() > 1.0 + 1e-6:
        warnings.warn("heteroplasmy values >1 detected; clipping to [0,1]")
        X.data = np.clip(X.data, 0.0, 1.0)

    obs = obs.copy()
    var = var.copy()

    # ---- per-cell coverage
    if coverage is None:
        # derive: mean read support proxy unavailable -> use n informative variants detected
        # per cell as a *relative* coverage proxy (documented in schema as derived).
        n_det = np.asarray((X > 0).sum(axis=1)).ravel()
        coverage = n_det.astype(np.float32)
        coverage_source = "derived_from_matrix"
    obs["coverage"] = np.asarray(coverage, dtype=np.float32)
    obs["coverage_source"] = coverage_source

    obs["assay"] = assay
    obs["dataset"] = dataset
    for col in ("donor", "sample", "tissue", "site"):
        if col not in obs.columns:
            obs[col] = "unknown"
    if "cell_id" not in obs.columns:
        obs["cell_id"] = obs.index.astype(str)

    # ---- per-variant stats
    het = X
    var["pseudobulk_heteroplasmy"] = np.asarray(het.mean(axis=0)).ravel().astype(np.float32)
    var["n_cells_detected"] = np.asarray((het > 0).sum(axis=0)).ravel().astype(int)
    if "strand_correlation" not in var.columns:
        var["strand_correlation"] = np.nan

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names = var["variant_id"].astype(str).values
    adata.obs_names = obs["cell_id"].astype(str).values
    adata.uns["mtclone_schema_version"] = SCHEMA_VERSION
    return adata


# ----------------------------------------------------------------------------- validator
def validate_schema(adata: ad.AnnData, *, strict: bool = True) -> list[str]:
    """Return a list of schema violations. Raises on the first if strict."""
    problems: list[str] = []

    def _fail(msg):
        problems.append(msg)
        if strict:
            raise ValueError(f"mtclone schema violation: {msg}")

    if not sp.issparse(adata.X):
        _fail("X must be sparse")
    else:
        if adata.X.dtype != np.float32:
            _fail(f"X dtype must be float32, got {adata.X.dtype}")
        if adata.X.nnz and (adata.X.data.min() < -1e-6 or adata.X.data.max() > 1 + 1e-6):
            _fail("X values must lie in [0,1]")
    for c in _REQUIRED_OBS:
        if c not in adata.obs.columns:
            _fail(f"missing required obs column '{c}'")
        elif adata.obs[c].isna().all():
            _fail(f"required obs column '{c}' is all-NaN")
    for c in _REQUIRED_VAR:
        if c not in adata.var.columns:
            _fail(f"missing required var column '{c}'")
    if adata.obs_names.has_duplicates:
        _fail("cell_id (obs_names) not unique")
    if adata.var_names.has_duplicates:
        _fail("variant_id (var_names) not unique")
    if "variant_id" in adata.var.columns and not np.array_equal(
        adata.var_names.astype(str).values, adata.var["variant_id"].astype(str).values
    ):
        _fail("var_names must equal var['variant_id']")
    if "donor" in adata.obs and adata.obs["donor"].isna().any():
        _fail("donor is null for some cells")
    if adata.uns.get("mtclone_schema_version") != SCHEMA_VERSION:
        _fail("uns['mtclone_schema_version'] not set correctly")
    return problems


# ----------------------------------------------------------------------------- adapters
def read_mtscatac_heteroplasmy(cell_heteroplasmic_df, variant_stats=None, *,
                               dataset="GSE302113", sample="unknown", donor="unknown",
                               tissue="unknown", site="unknown",
                               coverage=None, coverage_source="derived_from_matrix") -> ad.AnnData:
    """Ingest a mtscATAC per-cell heteroplasmy table (Liu / GSE302113 style).

    `cell_heteroplasmic_df`: long or wide table of cell x variant heteroplasmy. This adapter
    accepts a long table with columns among {cell, barcode, variant, heteroplasmy, af} or a
    wide cell x variant matrix (index = cells, columns = variants).
    `variant_stats`: optional per-variant stats table (strand correlation, coverage).
    `coverage`: OPTIONAL per-cell mean mtDNA coverage. Pass a mapping/Series indexed by the
        matrix's cell ids, or an array aligned to matrix row order. Used for the ASAP-seq /
        mgatk path (GSE197008) where R exports true per-cell depth alongside heteroplasmy;
        leave None to derive a relative proxy from the matrix (Liu GSE302113 default).
    """
    df = _load_table(cell_heteroplasmic_df)
    Xdf = _long_or_wide_to_wide(df)  # cells (rows) x variants (cols), values in [0,1]

    var = pd.DataFrame(index=Xdf.columns)
    parsed = [parse_variant_string(v) for v in Xdf.columns]
    var["variant_id"] = [
        canonical_variant_id(*p) if p else str(v)
        for p, v in zip(parsed, Xdf.columns)
    ]
    var["pos"] = [p[0] if p else -1 for p in parsed]
    var["ref"] = [p[1] if p else "N" for p in parsed]
    var["alt"] = [p[2] if p else "N" for p in parsed]

    if variant_stats is not None:
        vs = _load_table(variant_stats)
        sc_col = next((c for c in vs.columns if "strand" in c.lower()), None)
        id_col = next((c for c in vs.columns if c.lower() in ("variant", "variant_id", "var")), None)
        if sc_col and id_col:
            m = {parse_variant_string(r[id_col]): r[sc_col] for _, r in vs.iterrows()}
            var["strand_correlation"] = [
                m.get((p[0], p[1], p[2]), np.nan) if p else np.nan for p in parsed
            ]

    obs = pd.DataFrame(index=Xdf.index.astype(str))
    obs["cell_id"] = [f"{sample}_{c}" if not str(c).startswith(sample) else str(c)
                      for c in Xdf.index]
    obs["donor"] = donor
    obs["sample"] = sample
    obs["tissue"] = tissue
    obs["site"] = site

    cov = coverage
    cov_src = coverage_source
    if coverage is not None:
        if isinstance(coverage, (pd.Series, dict)):
            cov = np.asarray([float(pd.Series(coverage).get(c, np.nan)) for c in Xdf.index],
                             dtype=np.float32)
        else:
            cov = np.asarray(coverage, dtype=np.float32)
            if cov.shape[0] != Xdf.shape[0]:
                raise ValueError(f"coverage length {cov.shape[0]} != n_cells {Xdf.shape[0]}")
        cov_src = coverage_source if coverage_source != "derived_from_matrix" else "depth_file"

    return _finalize(Xdf.values, obs, var, assay="mtscATAC", dataset=dataset,
                     coverage=cov, coverage_source=cov_src)


def read_redeem_mtx(mtx, barcodes, features, depth=None, *,
                    dataset="GSE219014", sample="unknown", donor="unknown",
                    tissue="BMMC", site="marrow") -> ad.AnnData:
    """Ingest a ReDeeM variant matrix (GSE219014 style).

    `mtx`     : MatrixMarket cell x variant (or variant x cell) file.
    `barcodes`: one cell barcode per line.
    `features`: one variant string per line (parsed to canonical id).
    `depth`   : OPTIONAL per-cell coverage vector/file. GSE219014 ships NO depth file, so this
                is None there and coverage is derived from the matrix (documented in schema).
    """
    M = scipy.io.mmread(str(mtx)).tocsr().astype(np.float32)
    bc = _read_lines(barcodes)
    feat = _read_lines(features)

    # orient to cells x variants
    if M.shape[0] == len(feat) and M.shape[1] == len(bc):
        M = M.T.tocsr()
    if not (M.shape[0] == len(bc) and M.shape[1] == len(feat)):
        raise ValueError(f"matrix {M.shape} incompatible with {len(bc)} barcodes / {len(feat)} features")

    # ReDeeM matrices may be counts; if values exceed 1, treat as alt-allele fraction via
    # row-normalization is NOT correct (no ref counts here) -> if integer-valued, binarize to
    # presence as a conservative heteroplasmy proxy and warn.
    if M.nnz and M.data.max() > 1.0 + 1e-6:
        warnings.warn("ReDeeM matrix has values >1 (looks like counts); using presence "
                      "(alt>0 -> heteroplasmy set to min(1, alt/(alt))) as a proxy. Provide "
                      "a true heteroplasmy/depth pair for exact fractions.")
        M = M.copy()
        M.data = np.minimum(M.data, 1.0)  # presence proxy in [0,1]

    parsed = [parse_variant_string(f) for f in feat]
    var = pd.DataFrame({
        "variant_id": [canonical_variant_id(*p) if p else str(f) for p, f in zip(parsed, feat)],
        "pos": [p[0] if p else -1 for p in parsed],
        "ref": [p[1] if p else "N" for p in parsed],
        "alt": [p[2] if p else "N" for p in parsed],
    })

    obs = pd.DataFrame(index=[str(b) for b in bc])
    obs["cell_id"] = [f"{sample}_{b}" if not str(b).startswith(sample) else str(b) for b in bc]
    obs["donor"] = donor
    obs["sample"] = sample
    obs["tissue"] = tissue
    obs["site"] = site

    cov, cov_src = None, "derived_from_matrix"
    if depth is not None:
        cov = _load_depth(depth, bc)
        cov_src = "depth_file"

    return _finalize(M, obs, var, assay="redeem", dataset=dataset,
                     coverage=cov, coverage_source=cov_src)


# 14-column layout of RawGenotypes.<thr>.StrandBalance (from redeemR::redeemR.read GiveName)
_REDEEM_RAWGT_COLS = ["UMI", "Cell", "Pos", "Variants", "Call", "Ref", "FamSize",
                      "GT_Cts", "CSS", "DB_Cts", "SG_Cts", "Plus", "Minus", "Depth"]
_REDEEM_THR = {"T": "Total", "LS": "VerySensitive", "S": "Sensitive", "VS": "Specific"}


def _redeem_variant_to_canonical(v: str):
    """ReDeeM `Variants` string -> (pos, ref, alt). Handles `Pos_Ref_Alt` (e.g. `3243_A_G`),
    `PosRefAlt` (`3243A>G`), and the `chrM` variants; falls back to the shared parser."""
    s = str(v).strip()
    parts = s.split("_")
    if len(parts) == 3 and parts[0].isdigit() and parts[1].isalpha() and parts[2].isalpha():
        return int(parts[0]), parts[1].upper(), parts[2].upper()
    return parse_variant_string(s)


def read_redeem_consensus(path, thr="S", *, dataset="GSE219014", sample="unknown",
                          donor="unknown", tissue="BMMC", site="marrow",
                          cell_whitelist=None, min_depth=1, barcode_map=None,
                          barcode_space="atac") -> ad.AnnData:
    """Ingest a ReDeeM-V `*.Consensus.final/` folder (Weng GSE219014 / Figshare 24418966).

    This is the adapter that resolves the GEO barcode gap: the GEO `.mtx` ships no
    barcodes/features, but the authors' `Consensus.final` folders carry the full per-molecule
    genotype table WITH cell barcodes, variant identities, and per-(cell,pos) depth — so true
    heteroplasmy (not the presence-proxy of `read_redeem_mtx`) is recoverable.

    Folder contents used (see redeemR::redeemR.read):
      - `RawGenotypes.<thr>.StrandBalance` — one row per (UMI, cell, variant); 14 cols
        (`_REDEEM_RAWGT_COLS`). Alt-allele molecule count per (Cell, Variant) = row count
        (== redeemR `Freq`); `Depth` = per-(cell,pos) coverage. hetero = Freq / Depth.
      - `QualifiedTotalCts` — per (cell, pos) total coverage; used for per-cell MEAN coverage
        (schema `obs['coverage']`, source `depth_file`).

    `thr` : consensus threshold — one of {'T','LS','S','VS'} (redeemR naming). Default 'S'
            (Sensitive), the setting used throughout the Weng reproducibility notebooks.
    `cell_whitelist` : optional iterable of barcodes to keep (e.g. cells present in the paired
            annotated Seurat object) — clones are called per real cell, so restrict early.
            Compared AFTER barcode translation, so pass whitelist barcodes in the SAME space
            the map translates INTO (RNA space if barcode_map is the ATAC->RNA multiome map).
    `min_depth` : drop (cell,variant) observations whose depth < this before computing hetero.
    `barcode_map` : optional 10x-Multiome ATAC<->RNA barcode translation. The Consensus.final
            cell barcodes are in **ATAC** whitelist space, but the annotated Seurat objects use
            **RNA** whitelist space — they do NOT overlap without translation. Pass a mapping
            (dict or 2-col DataFrame/path) so returned cell_ids match the Seurat annotation.
            Accepts: dict {atac: rna}; a DataFrame/tsv path with columns 'atac','rna'.
    `barcode_space` : which space the folder's barcodes are in ('atac', default). If a map is
            given, barcodes are translated to the OTHER space; unmapped barcodes are dropped.
    """
    d = Path(path)
    thr_name = _REDEEM_THR.get(thr, thr)
    # match ONLY the requested threshold (do not silently fall back to another thr's file)
    gt_candidates = ([d / f"RawGenotypes.{thr_name}.StrandBalance"]
                     + list(d.glob(f"RawGenotypes.{thr_name}.*"))
                     + list(d.glob(f"RawGenotypes.{thr_name}")))
    gt_path = next((p for p in gt_candidates if p.exists()), None)
    if gt_path is None:
        raise FileNotFoundError(
            f"no RawGenotypes.{thr_name}[.StrandBalance] in {d} (thr={thr!r}). "
            f"Present: {[p.name for p in d.iterdir()][:12]}")

    # ---- per-molecule genotype rows -> (Cell, Variant) alt-molecule counts + depth
    raw = pd.read_csv(gt_path, sep="\t", header=None, names=_REDEEM_RAWGT_COLS,
                      dtype={"Cell": str, "Variants": str, "Pos": int, "Depth": float})

    # ---- 10x Multiome barcode translation (ATAC folder space -> RNA Seurat space)
    xmap = None
    if barcode_map is not None:
        if isinstance(barcode_map, dict):
            bm = pd.DataFrame({"atac": list(barcode_map), "rna": list(barcode_map.values())})
        elif isinstance(barcode_map, (str, Path)):
            bm = pd.read_csv(barcode_map, sep="\t")
        else:
            bm = barcode_map.copy()
        src, dst = ("atac", "rna") if barcode_space == "atac" else ("rna", "atac")
        xmap = dict(zip(bm[src].astype(str), bm[dst].astype(str)))
    if xmap is not None:
        raw["Cell"] = raw["Cell"].map(xmap)
        raw = raw[raw["Cell"].notna()]                 # drop barcodes absent from the whitelist
        if raw.empty:
            raise ValueError(f"no barcodes translated via barcode_map for {gt_path} "
                             f"(barcode_space={barcode_space!r} — wrong direction?)")

    if cell_whitelist is not None:
        wl = set(str(x) for x in cell_whitelist)
        raw = raw[raw["Cell"].isin(wl)]
    if min_depth > 1:
        raw = raw[raw["Depth"] >= min_depth]
    if raw.empty:
        raise ValueError(f"no genotype rows left in {gt_path} after filtering "
                         f"(cell_whitelist / min_depth={min_depth})")

    # alt-molecule count per (Cell, Variant) == redeemR Freq
    freq = (raw.groupby(["Cell", "Variants"]).size().rename("freq").reset_index())
    # depth per (Cell, Variant): depth is per (cell,pos); take the max over the rows (all equal per pos)
    depth_cv = (raw.groupby(["Cell", "Variants"])["Depth"].max().rename("depth").reset_index())
    g = freq.merge(depth_cv, on=["Cell", "Variants"])
    g["hetero"] = np.minimum(1.0, g["freq"] / g["depth"].clip(lower=1.0))

    cells = pd.Index(sorted(g["Cell"].unique()), name="cell")
    variants = pd.Index(sorted(g["Variants"].unique()), name="variant")
    ci = {c: i for i, c in enumerate(cells)}
    vi = {v: i for i, v in enumerate(variants)}
    X = sp.csr_matrix(
        (g["hetero"].to_numpy(np.float32),
         (g["Cell"].map(ci).to_numpy(), g["Variants"].map(vi).to_numpy())),
        shape=(len(cells), len(variants)))

    # ---- per-cell mean coverage from QualifiedTotalCts (cell, pos, total_cts ...)
    cov = None
    cov_src = "derived_from_matrix"
    qtc_path = d / "QualifiedTotalCts"
    if qtc_path.exists():
        qtc = pd.read_csv(qtc_path, sep="\t", header=None)
        # cols: V1=Cell, V2=Pos, V3..=coverage flavours (Total/... per redeemR DepthSummary)
        qtc = qtc.rename(columns={0: "Cell", 1: "Pos", 2: "cov"})
        qtc["Cell"] = qtc["Cell"].astype(str)
        if xmap is not None:                            # same ATAC->RNA translation as genotypes
            qtc["Cell"] = qtc["Cell"].map(xmap)
            qtc = qtc[qtc["Cell"].notna()]
        cov_per_cell = qtc.groupby("Cell")["cov"].mean()
        cov = np.asarray([cov_per_cell.get(c, np.nan) for c in cells], dtype=np.float32)
        cov_src = "depth_file"

    parsed = [_redeem_variant_to_canonical(v) for v in variants]
    var = pd.DataFrame({
        "variant_id": [canonical_variant_id(*p) if p else str(v) for p, v in zip(parsed, variants)],
        "pos": [p[0] if p else -1 for p in parsed],
        "ref": [p[1] if p else "N" for p in parsed],
        "alt": [p[2] if p else "N" for p in parsed],
    })
    obs = pd.DataFrame(index=[str(c) for c in cells])
    obs["cell_id"] = [f"{sample}_{c}" if not str(c).startswith(sample) else str(c) for c in cells]
    obs["donor"] = donor
    obs["sample"] = sample
    obs["tissue"] = tissue
    obs["site"] = site

    return _finalize(X, obs, var, assay="redeem", dataset=dataset,
                     coverage=cov, coverage_source=cov_src)


def read_mgatk(mgatk_dir, *, dataset="unknown", sample="unknown", donor="unknown",
               tissue="unknown", site="unknown") -> ad.AnnData:
    """Ingest a generic mgatk output directory.

    Expects the standard mgatk final/ outputs: `*.A/C/G/T.txt.gz` allele-count matrices and
    `*.coverage.txt.gz`, or the summarized `*.heteroplasmy` / `*.variant` products. This is
    the target for the optional re-calling path and for third-party mtscATAC datasets.
    Minimal implementation: read a precomputed cell x variant heteroplasmy tsv if present.
    """
    d = Path(mgatk_dir)
    het = list(d.glob("*heteroplasm*")) + list(d.glob("*cell*variant*"))
    if not het:
        raise FileNotFoundError(
            f"no heteroplasmy product found in {d}. Full allele-count assembly "
            "(A/C/G/T + coverage -> heteroplasmy) is a TODO for the recall path; for now "
            "point read_mgatk at a precomputed cell x variant heteroplasmy table."
        )
    return read_mtscatac_heteroplasmy(het[0], dataset=dataset, sample=sample, donor=donor,
                                      tissue=tissue, site=site)


# ----------------------------------------------------------------------------- io utils
def _load_table(x) -> pd.DataFrame:
    if isinstance(x, pd.DataFrame):
        return x
    p = str(x)
    sep = "\t" if (p.endswith(".tsv") or p.endswith(".tsv.gz") or p.endswith(".txt")
                   or p.endswith(".txt.gz")) else ","
    return pd.read_csv(p, sep=sep, index_col=0)


def _long_or_wide_to_wide(df: pd.DataFrame) -> pd.DataFrame:
    """Return cells x variants wide matrix in [0,1]."""
    cols = {c.lower(): c for c in df.columns}
    het_col = next((cols[c] for c in ("heteroplasmy", "af", "vaf", "value") if c in cols), None)
    cell_col = next((cols[c] for c in ("cell", "barcode", "cell_id") if c in cols), None)
    var_col = next((cols[c] for c in ("variant", "variant_id", "var", "mutation") if c in cols), None)
    if het_col and cell_col and var_col:  # long
        wide = df.pivot_table(index=cell_col, columns=var_col, values=het_col,
                              aggfunc="max", fill_value=0.0)
        return wide.astype(np.float32)
    # already wide (index = cells, columns = variants)
    return df.astype(np.float32)


def _read_lines(path) -> list[str]:
    if isinstance(path, (list, tuple, np.ndarray, pd.Index)):
        return [str(x) for x in path]
    import gzip
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        return [ln.split("\t")[0].strip() for ln in fh if ln.strip()]


def _load_depth(depth, barcodes) -> np.ndarray:
    if isinstance(depth, (np.ndarray, list, pd.Series)):
        return np.asarray(depth, dtype=np.float32)
    df = _load_table(depth)
    # assume a single coverage column, index by barcode
    col = df.columns[-1]
    d = df[col]
    return np.asarray([d.get(b, np.nan) for b in barcodes], dtype=np.float32)
