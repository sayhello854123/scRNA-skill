#!/usr/bin/env Rscript
# Summarize an unknown single-cell object: dims, assays, layers, metadata,
# and which analysis steps have already been run.
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.here, "_common.R"))

HELP <- "Usage: Rscript inspect_data.R <input> [--max-meta N]
Summarize a Seurat object (or any loadable input). Read-only."

opts <- parse_args(commandArgs(trailingOnly = TRUE),
                   spec = list(`max-meta` = "int"),
                   defaults = list(`max-meta` = 20),
                   help_text = HELP)
if (length(opts$`_`) < 1) die("need an input path. Try --help")

obj <- load_seurat(opts$`_`[[1]])
cat(sprintf("\n== Seurat object ==\n%d cells x %d features\n", ncol(obj), nrow(obj)))
cat("Active assay:", DefaultAssay(obj), "\n")
cat("Assays:", paste(Seurat::Assays(obj), collapse = ", "), "\n")
for (a in Seurat::Assays(obj)) {
  ly <- tryCatch(Layers(obj[[a]]), error = function(e) character(0))
  cat(sprintf("  [%s] layers: %s\n", a, paste(ly, collapse = ", ")))
}
cat("Reductions:", paste(Reductions(obj), collapse = ", "), "\n")
cat("Variable features:", length(VariableFeatures(obj)), "\n")
cat("Idents levels:", paste(head(levels(Idents(obj)), 30), collapse = ", "), "\n")

md <- obj@meta.data
cat(sprintf("\n== Metadata (%d columns) ==\n", ncol(md)))
for (col in head(colnames(md), opts$`max-meta`)) {
  x <- md[[col]]
  if (is.numeric(x)) {
    cat(sprintf("  %-22s numeric  range [%.3g, %.3g] median %.3g\n",
                col, min(x, na.rm = TRUE), max(x, na.rm = TRUE), median(x, na.rm = TRUE)))
  } else {
    lv <- length(unique(x))
    cat(sprintf("  %-22s factor   %d levels\n", col, lv))
  }
}
cat("\nComputed so far:",
    paste(c(if ("percent.mt" %in% colnames(md)) "QC",
            if (length(VariableFeatures(obj))) "HVG",
            if ("pca" %in% Reductions(obj)) "PCA",
            if ("umap" %in% Reductions(obj)) "UMAP",
            if (any(grepl("clusters|res", colnames(md)))) "clusters"),
          collapse = ", "), "\n")
