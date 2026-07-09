#!/usr/bin/env Rscript
# export_redeem_seurat.R
# ----------------------------------------------------------------------------
# Boundary step for the Weng ReDeeM arm of Plan 2. The annotated Multiome Seurat
# objects (Figshare 23290004) carry the per-cell cell-type call + RNA that the
# Consensus.final mtDNA folders lack. This script exports, per donor object:
#   <prefix>.cellmeta.tsv.gz  — barcode, STD.CellType, STD_Cat, STD_Cat2, Sample,
#                               meanCov, ClonalGroup, ClonalGroup.Prob (author ref)
#   <prefix>.nk_markers.tsv.gz — per-cell RNA expression of the bright/dim + NK
#                               marker genes (log-normalized), for classify.
#   <prefix>.wnn_umap.tsv.gz   — barcode, wnnUMAP_1/2 (paired-modality embedding)
#
# mtclone's ReDeeM ingest (read_redeem_consensus) then joins these onto the
# heteroplasmy matrix BY BARCODE to (a) restrict to NK, (b) label bright/dim.
#
# Usage: Rscript export_redeem_seurat.R <Seurat.RDS> <out_prefix>
# ----------------------------------------------------------------------------
suppressMessages({library(Seurat); library(SeuratObject); library(data.table)})

args <- commandArgs(trailingOnly = TRUE)
rds <- args[1]; out_prefix <- args[2]
# NK identity + bright/dim + progenitor markers (mtclone.classify vocabulary)
NK_GENES <- c("NCAM1","NCR1","GNLY","KLRD1","NKG7","KLRF1",          # NK identity
              "GZMK","SELL","XCL1","IL7R","GZMB","FGFBP2","PRF1","CX3CR1","FCGR3A")  # bright/dim
PROG_GENES <- c("CD34","SPINK2","MLLT3","AVP","HLF")                  # HSPC context

cat(sprintf("[seurat] reading %s\n", rds))
s <- readRDS(rds)
md <- s@meta.data
bc <- colnames(s)

meta_cols <- intersect(c("STD.CellType","STD_Cat","STD_Cat2","Sample","meanCov",
                         "ClonalGroup","ClonalGroup.Prob","seurat_clusters"),
                       colnames(md))
cm <- data.table(barcode = bc, as.data.table(md[, meta_cols, drop = FALSE]))
fwrite(cm, paste0(out_prefix, ".cellmeta.tsv.gz"), sep = "\t")
cat(sprintf("[seurat] %d cells; STD.CellType levels: %s\n", length(bc),
            paste(sort(unique(as.character(md$STD.CellType))), collapse = ",")))
cat(sprintf("[seurat] NK cells: %d\n", sum(md$STD.CellType == "NK", na.rm = TRUE)))

# RNA marker expression (log-normalized). SeuratObject 5 uses layers; get counts then
# normalize ourselves to avoid depending on assay-specific NormalizeData internals / Signac.
rna <- s[["RNA"]]
get_layer <- function(assay, layer) {
  tryCatch(SeuratObject::LayerData(assay, layer = layer), error = function(e) NULL)
}
dat <- get_layer(rna, "data")
if (is.null(dat)) {                                  # no normalized layer -> make it from counts
  cnt <- get_layer(rna, "counts")
  if (is.null(cnt)) cnt <- SeuratObject::GetAssayData(rna, layer = "counts")
  libsize <- Matrix::colSums(cnt)
  dat <- cnt
  dat@x <- log1p(dat@x / rep(libsize, diff(dat@p)) * 1e4)   # log-normalize (CP10k)
}
rna_genes <- rownames(dat)
genes <- intersect(c(NK_GENES, PROG_GENES), rna_genes)
expr <- as.matrix(dat[genes, , drop = FALSE])
em <- data.table(barcode = colnames(dat), as.data.table(t(expr)))
fwrite(em, paste0(out_prefix, ".nk_markers.tsv.gz"), sep = "\t")
cat(sprintf("[seurat] exported %d/%d marker genes present\n", length(genes),
            length(c(NK_GENES, PROG_GENES))))

# wnn.umap embedding if present
if ("wnn.umap" %in% names(s@reductions)) {
  emb <- Embeddings(s, "wnn.umap")
  ue <- data.table(barcode = rownames(emb), wnnUMAP_1 = emb[,1], wnnUMAP_2 = emb[,2])
  fwrite(ue, paste0(out_prefix, ".wnn_umap.tsv.gz"), sep = "\t")
}
cat(sprintf("[seurat] wrote %s.{cellmeta,nk_markers,wnn_umap}.tsv.gz\n", out_prefix))
