#!/usr/bin/env Rscript
# =============================================================================
# Single-cell RNA-seq analysis template (Seurat v5)
# Copy this file, edit the PARAMETERS block, and run:  Rscript my_analysis.R
# Prefer the scripts/ CLI tools for routine runs; use this when you want a
# single editable script for a bespoke analysis.
# =============================================================================

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})
set.seed(42)

## ---- PARAMETERS -------------------------------------------------------------
input_path    <- "filtered_feature_bc_matrix/"  # 10x dir, .h5, .rds, or .csv
output_path   <- "results/processed.rds"
mt_pattern    <- "^MT-"     # "^mt-" for mouse
min_features  <- 200
max_features  <- 6000
max_mt        <- 10
norm_method   <- "lognorm"  # "lognorm" or "sct"
n_features    <- 2000
dims          <- 1:30
resolution    <- 0.5
fig_dir       <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load ----------------------------------------------------------------
if (dir.exists(input_path)) {
  obj <- CreateSeuratObject(Read10X(input_path), project = "scRNA",
                            min.cells = 3, min.features = 200)
} else if (grepl("\\.rds$", input_path)) {
  obj <- readRDS(input_path)
} else if (grepl("\\.h5$", input_path)) {
  obj <- CreateSeuratObject(Read10X_h5(input_path), min.cells = 3, min.features = 200)
} else {
  stop("edit the loader for your input format")
}

## ---- 2. QC ------------------------------------------------------------------
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
ggsave(file.path(fig_dir, "qc_violin.png"),
       VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0),
       width = 10, height = 4)
obj <- subset(obj, subset = nFeature_RNA > min_features &
                            nFeature_RNA < max_features & percent.mt < max_mt)

## ---- 3. Normalize -----------------------------------------------------------
if (norm_method == "sct") {
  obj <- SCTransform(obj, variable.features.n = n_features, verbose = FALSE)
} else {
  obj <- NormalizeData(obj) |>
    FindVariableFeatures(nfeatures = n_features) |>
    ScaleData()
}

## ---- 4. Dimensionality reduction --------------------------------------------
obj <- RunPCA(obj, verbose = FALSE)
ggsave(file.path(fig_dir, "elbow.png"), ElbowPlot(obj, ndims = 50))
obj <- FindNeighbors(obj, dims = dims)
obj <- RunUMAP(obj, dims = dims)

## ---- 5. Cluster -------------------------------------------------------------
obj <- FindClusters(obj, resolution = resolution)
ggsave(file.path(fig_dir, "umap_clusters.png"), DimPlot(obj, label = TRUE))

## ---- 6. Markers -------------------------------------------------------------
if ("RNA" %in% Assays(obj)) obj <- JoinLayers(obj)
if (DefaultAssay(obj) == "SCT") obj <- PrepSCTFindMarkers(obj)
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dir.create("results/markers", showWarnings = FALSE, recursive = TRUE)
write.csv(markers, "results/markers/all_markers.csv", row.names = FALSE)
top <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5) %>% pull(gene)
ggsave(file.path(fig_dir, "marker_dotplot.png"),
       DotPlot(obj, features = unique(top)) + RotatedAxis(), width = 12)

## ---- 7. Annotate (edit after inspecting markers) ----------------------------
# new_ids <- c("0" = "CD4 T", "1" = "CD14+ Mono", "2" = "B", "3" = "CD8 T")
# obj <- RenameIdents(obj, new_ids)
# obj$cell_type <- Idents(obj)
# ggsave(file.path(fig_dir, "umap_celltypes.png"), DimPlot(obj, label = TRUE))

## ---- 8. Save ----------------------------------------------------------------
saveRDS(obj, output_path)
write.csv(obj@meta.data, "results/cell_metadata.csv")
cat("Done. Object ->", output_path, "\n")
