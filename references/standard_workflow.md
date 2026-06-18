# Standard Seurat v5 Workflow for Single-Cell Analysis

A complete, annotated pipeline from raw counts to annotated cell types. The
`scripts/` CLI tools automate each step; this document is the reference for the
underlying calls and for customizing beyond the script flags.

## 1. Loading data

```r
library(Seurat)

# 10x Cell Ranger output (a directory with barcodes/features/matrix)
counts <- Read10X(data.dir = "filtered_feature_bc_matrix/")
obj <- CreateSeuratObject(counts, project = "scRNA", min.cells = 3, min.features = 200)

# 10x HDF5
counts <- Read10X_h5("filtered_feature_bc_matrix.h5")

# Existing object
obj <- readRDS("data.rds")
```

`min.cells`/`min.features` apply a light filter at object creation. `.h5ad`
inputs must be converted first — see `integration.md`.

## 2. Quality control

```r
obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^MT-")   # mouse: "^mt-"
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
FeatureScatter(obj, "nCount_RNA", "percent.mt")
FeatureScatter(obj, "nCount_RNA", "nFeature_RNA")

obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
```

Set the upper `nFeature_RNA` bound to remove likely doublets. For dedicated
doublet detection use `DoubletFinder` or `scDblFinder` (separate packages).
Thresholds are dataset-specific — always inspect the violin/scatter plots.

## 3. Normalization

**LogNormalize** (classic three-step):

```r
obj <- NormalizeData(obj)                          # log(CPM/100 + 1)-style
obj <- FindVariableFeatures(obj, nfeatures = 2000) # HVGs by variance-stabilized dispersion
obj <- ScaleData(obj)                              # z-score; can regress covariates
```

**SCTransform** (recommended; regularized negative-binomial model, replaces all
three steps and stores results in a new `SCT` assay):

```r
obj <- SCTransform(obj, vars.to.regress = "percent.mt")
```

Install `glmGamPoi` to make SCTransform much faster. After SCTransform the
default assay becomes `SCT`; switch with `DefaultAssay(obj) <- "RNA"` when you
need raw/log counts.

## 4. Dimensionality reduction

```r
obj <- RunPCA(obj, npcs = 50)
ElbowPlot(obj, ndims = 50)          # pick where the curve flattens
DimHeatmap(obj, dims = 1:15, cells = 500, balanced = TRUE)

obj <- FindNeighbors(obj, dims = 1:30)
obj <- RunUMAP(obj, dims = 1:30)
DimPlot(obj, reduction = "umap")
# Optional: obj <- RunTSNE(obj, dims = 1:30)
```

The chosen `dims` propagate through neighbors, UMAP, and clustering. Too few
dims loses structure; too many adds noise.

## 5. Clustering

```r
obj <- FindNeighbors(obj, dims = 1:30)          # if not already done
obj <- FindClusters(obj, resolution = 0.5)      # algorithm 1 = Louvain (default)
# Leiden: FindClusters(obj, resolution = 0.5, algorithm = 4)  # needs leidenalg
DimPlot(obj, label = TRUE)

# Scan resolutions into separate columns
for (res in c(0.3, 0.5, 0.8, 1.0))
  obj <- FindClusters(obj, resolution = res, cluster.name = paste0("res", res))
```

Higher resolution → more, finer clusters. Cross-check with `clustree` if unsure.

## 6. Marker genes (exploratory)

Seurat v5 keeps per-sample layers; merge them before testing:

```r
obj <- JoinLayers(obj)                          # RNA assay: combine layers
# SCT objects: PrepSCTFindMarkers(obj) before marker tests

# One cluster vs the rest
FindMarkers(obj, ident.1 = "0", min.pct = 0.25, logfc.threshold = 0.25)

# All clusters at once (install 'presto' for speed)
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10)
```

`FindAllMarkers` p-values are **exploratory**: cells within a cluster are not
independent, so per-cell tests are anti-conservative. For condition/treatment
comparisons, pseudobulk (Section 9) and use a bulk DE tool.

## 7. Cell-type annotation

```r
panel <- c("CD3D", "CD14", "MS4A1", "NKG7", "FCGR3A", "LYZ", "PPBP")
DotPlot(obj, features = panel) + RotatedAxis()
FeaturePlot(obj, features = panel)
VlnPlot(obj, features = panel, pt.size = 0)

new_ids <- c("0" = "CD4 T", "1" = "CD14+ Mono", "2" = "B", "3" = "CD8 T",
             "4" = "NK", "5" = "FCGR3A+ Mono", "6" = "DC", "7" = "Platelet")
obj <- RenameIdents(obj, new_ids)
obj$cell_type <- Idents(obj)
DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE)
```

For automated reference-based annotation, use **Azimuth** (`RunAzimuth`) or
correlate cluster averages against a labeled reference. See `integration.md`.

## 8. Gene-set / module scoring

```r
features <- list(Tcell = c("CD3D", "CD3E", "CD3G"))
obj <- AddModuleScore(obj, features = features, name = "Tcell")
FeaturePlot(obj, "Tcell1")

# Cell cycle scoring with built-in gene lists
obj <- CellCycleScoring(obj, s.features = cc.genes$s.genes,
                        g2m.features = cc.genes$g2m.genes)
```

## 9. Pseudobulk + condition differential expression

```r
# Aggregate raw counts by sample x cell type, then run a bulk DE tool
pb <- AggregateExpression(obj, assays = "RNA", return.seurat = TRUE,
                          group.by = c("sample", "cell_type"))
# Export pb counts and analyze with DESeq2 / edgeR / pydeseq2 (see pydeseq2 skill)
```

Pseudobulk respects biological replication, giving valid p-values for
between-condition comparisons.

## 10. Saving

```r
saveRDS(obj, "results/processed.rds")
write.csv(obj@meta.data, "results/cell_metadata.csv")
```

## End-to-end via scripts

```bash
Rscript scripts/run_pipeline.R raw.rds -o processed.rds --resolution 0.5
# or step by step — see SKILL.md "Script Toolkit".
```
