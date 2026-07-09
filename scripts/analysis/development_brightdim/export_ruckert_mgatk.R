#!/usr/bin/env Rscript
# export_ruckert_mgatk.R
# ----------------------------------------------------------------------------
# Boundary step for the Ruckert (GSE197008/197037) ASAP-seq arm of Plan 2.
#
# The deposited GSM*_mgatk.rds is a Seurat `Assay` holding the RAW mgatk allele
# counts: 132,552 rows = {A,C,G,T} x 16,569 positions x {fwd,rev} strands,
# columns = cell barcodes. It is one step upstream of a heteroplasmy matrix, so
# this script performs the standard mgatk/Signac variant call:
#   1. collapse strands -> per-(base,pos,cell) counts (fwd+rev)
#   2. per (pos,cell) coverage = sum over 4 bases
#   3. ref allele per position from rCRS; alt = the 3 non-ref bases
#   4. heteroplasmy(cell,variant) = alt_count / coverage_at_pos
#   5. strand concordance per variant (fwd vs rev correlation) as a QC stat
#   6. keep variants seen in >= min_cells at >= min_het (informative set),
#      export a LONG cell x variant heteroplasmy table + per-cell mean coverage.
#
# Emits portable tsv the pure-Python mtclone.io.read_mtscatac_heteroplasmy reads
# (with coverage=), so mtclone stays R-free downstream.
#
# Usage:
#   Rscript export_ruckert_mgatk.R <mgatk.rds> <rCRS.fasta> <out_prefix> [min_cells] [min_het]
# ----------------------------------------------------------------------------
suppressMessages({library(SeuratObject); library(Matrix); library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
rds_path <- args[1]; fasta_path <- args[2]; out_prefix <- args[3]
min_cells <- ifelse(length(args) >= 4, as.integer(args[4]), 5L)
min_het   <- ifelse(length(args) >= 5, as.numeric(args[5]), 0.01)
strand_cor_min <- 0.65   # Lareau mtscATAC default strand-concordance filter

cat(sprintf("[mgatk] reading %s\n", rds_path))
# GEO's *_mgatk.rds.gz is externally gzipped in a form R's gzfile()/readRDS reject
# ("unknown input format"). Robust path: decompress to a tempfile, then readRDS.
read_rds_any <- function(p) {
  ok <- tryCatch({ obj <- readRDS(p); TRUE }, error = function(e) FALSE)
  if (ok) return(obj)
  if (!grepl("\\.gz$", p)) stop(sprintf("readRDS failed and %s is not .gz", p))
  tmp <- tempfile(fileext = ".rds")
  system2("gunzip", c("-c", shQuote(p)), stdout = tmp)
  on.exit(unlink(tmp), add = TRUE)
  readRDS(tmp)
}
a <- read_rds_any(rds_path)
M <- tryCatch(LayerData(a, layer = "counts"),
              error = function(e) GetAssayData(a, layer = "counts"))
if (is.null(M) || nrow(M) == 0) M <- GetAssayData(a, layer = "data")
rn <- rownames(M)
stopifnot(length(rn) > 0)

# parse rowname "{base}-{pos}-{strand}"
parts   <- tstrsplit(rn, "-", fixed = TRUE)
base_v  <- parts[[1]]; pos_v <- as.integer(parts[[2]]); strand_v <- parts[[3]]
stopifnot(all(base_v %in% c("A","C","G","T")))

# rCRS ref allele per position (1-based)
fa <- readLines(fasta_path)
ref_seq <- paste0(fa[!grepl("^>", fa)], collapse = "")
ref_seq <- toupper(ref_seq)
stopifnot(nchar(ref_seq) == 16569)
ref_at  <- function(p) substr(ref_seq, p, p)

cells <- colnames(M); n_cells <- length(cells)
cat(sprintf("[mgatk] %d rows x %d cells; collapsing strands\n", nrow(M), n_cells))

# For each (base,pos) we need fwd-row and rev-row indices into M. Build lookup keyed
# by base_pos so we can address fwd/rev directly (no 66k-way split).
key_bp   <- paste(base_v, pos_v, sep = "_")
row_fwd  <- integer(0); row_rev <- integer(0)
fwd_mask <- strand_v == "fwd"; rev_mask <- strand_v == "rev"
fwd_key  <- key_bp[fwd_mask]; rev_key <- key_bp[rev_mask]
fwd_idx  <- seq_along(rn)[fwd_mask]; rev_idx <- seq_along(rn)[rev_mask]
# align fwd and rev to a common (base,pos) universe
uniq_bp  <- sort(unique(key_bp))
ub_base  <- sub("_.*$", "", uniq_bp)
ub_pos   <- as.integer(sub("^.*_", "", uniq_bp))
fmap <- match(uniq_bp, fwd_key); rmap <- match(uniq_bp, rev_key)

Mfwd <- M[fwd_idx[fmap], , drop = FALSE]                 # (base,pos) x cell, fwd only
Mrev <- M[rev_idx[rmap], , drop = FALSE]                 # (base,pos) x cell, rev only
BP   <- Mfwd + Mrev                                      # strand-summed counts, sparse
rownames(BP) <- uniq_bp

# ---- per-position coverage = sum over 4 bases at that pos (sparse aggregate)
pos_fac    <- factor(ub_pos)
G <- sparseMatrix(i = as.integer(pos_fac), j = seq_len(nrow(BP)), x = 1,
                  dims = c(nlevels(pos_fac), nrow(BP)))
cov_by_pos <- as(G %*% BP, "CsparseMatrix")               # (pos) x cell, sparse
pos_levels <- as.integer(levels(pos_fac))
pos_row_of <- setNames(seq_along(pos_levels), pos_levels)  # pos -> row index in cov_by_pos

# ---- VECTORIZED variant call on sparse triplets (alt base != ref only)
cat("[mgatk] calling variants (alt != rCRS ref), vectorized\n")
ref_of_bp <- vapply(ub_pos, ref_at, character(1))
alt_rows  <- which(ub_base != ref_of_bp)
BPa <- BP[alt_rows, , drop = FALSE]                        # (n_alt) x cell, sparse alt counts
trip <- summary(as(BPa, "CsparseMatrix"))                  # data.frame i,j,x over nonzeros
trip <- as.data.table(trip)                                # i = alt-row idx, j = cell idx, x = alt count
# coverage at each (alt-row's position, cell)
trip[, pos := ub_pos[alt_rows][i]]
trip[, prow := pos_row_of[as.character(pos)]]
cov_at <- cov_by_pos[cbind(trip$prow, trip$j)]             # sparse elementwise lookup
trip[, cov := as.numeric(cov_at)]
trip[, het := ifelse(cov > 0, x / cov, 0)]
trip <- trip[het >= min_het & cov > 0]
# n cells per alt row; keep informative variants
nc <- trip[, .(ncell = uniqueN(j)), by = i]
keep_rows <- nc[ncell >= min_cells, i]
trip <- trip[i %in% keep_rows]
cat(sprintf("[mgatk] %d candidate alt rows -> %d pass min_cells=%d\n",
            length(alt_rows), length(keep_rows), min_cells))

# ---- strand concordance for surviving variants, VECTORIZED across cells.
# Pearson r between fwd and rev alt-count row-vectors, computed for all keep_rows at once:
#   r = (n*sum(xy) - sum(x)sum(y)) / sqrt((n*sum(x^2)-sum(x)^2)(n*sum(y^2)-sum(y)^2))
kr   <- keep_rows
Fk   <- Mfwd[alt_rows[kr], , drop = FALSE]                # (n_keep) x cell
Rk   <- Mrev[alt_rows[kr], , drop = FALSE]
n    <- ncol(Fk)
sx   <- Matrix::rowSums(Fk);           sy  <- Matrix::rowSums(Rk)
sxx  <- Matrix::rowSums(Fk * Fk);      syy <- Matrix::rowSums(Rk * Rk)
sxy  <- Matrix::rowSums(Fk * Rk)
num  <- n * sxy - sx * sy
den  <- sqrt(pmax(0, n * sxx - sx^2) * pmax(0, n * syy - sy^2))
scor <- ifelse(den > 0, num / den, NA_real_)
scor_of <- setNames(scor, kr)
trip[, strand_cor := scor_of[as.character(i)]]
# NA strand_cor = fwd or rev has zero variance = single-strand-only support = strand-bias
# artifact (the exact failure mode the concordance filter targets) -> DROP, do not keep.
n_na <- trip[is.na(strand_cor), uniqueN(i)]
trip <- trip[!is.na(strand_cor) & strand_cor >= strand_cor_min]
cat(sprintf("[mgatk] dropped %d single-strand-only variants (NA strand_cor)\n", n_na))

# ---- assemble long table with canonical Pos_Ref_Alt ids
trip[, variant := sprintf("%d_%s_%s", pos, ref_of_bp[alt_rows][i], ub_base[alt_rows][i])]
trip[, cell := cells[j]]
long <- trip[, .(cell, variant, heteroplasmy = het, strand_cor)]
cat(sprintf("[mgatk] %d informative variants, %d (cell,variant) obs\n",
            length(unique(long$variant)), nrow(long)))

# ---- per-cell mean coverage (mean per-position depth) -> real coverage
cell_cov <- data.table(cell = cells, coverage = as.numeric(Matrix::colMeans(cov_by_pos)))

fwrite(long,     paste0(out_prefix, ".heteroplasmy.tsv.gz"), sep = "\t")
fwrite(cell_cov, paste0(out_prefix, ".coverage.tsv.gz"),     sep = "\t")
# per-variant strand cor (unique)
vstats <- unique(long[, .(variant, strand_cor)])
fwrite(vstats,   paste0(out_prefix, ".variant_stats.tsv.gz"), sep = "\t")
cat(sprintf("[mgatk] wrote %s.{heteroplasmy,coverage,variant_stats}.tsv.gz\n", out_prefix))
