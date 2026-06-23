---
name: seurat
description: Standard single-cell RNA-seq analysis in R with Seurat v5. Use for QC, normalization (LogNormalize/SCTransform), dimensionality reduction (PCA/UMAP/t-SNE), clustering, marker gene identification, multi-sample integration (Harmony/CCA/RPCA), cell-type annotation, visualization, and reading 10x / .rds / .h5ad inputs. Best for R-based scRNA-seq workflows. For Python use the scanpy skill; for deep-learning models use scvi-tools; for AnnData I/O use anndata.
license: MIT
metadata: {"version": "0.1", "skill-author": "xpq"}
---

# Seurat: Single-Cell Analysis in R

## Overview

Seurat is the most widely used R toolkit for single-cell RNA-seq analysis. Apply this skill for complete workflows including quality control, normalization, dimensionality reduction, clustering, marker gene identification, multi-sample integration, cell-type annotation, and publication-quality visualization. Target release: **Seurat v5** (v5.x), which introduces the `Assay5` layer-based data model and on-disk `BPCells` matrices for large datasets.

## Installation

Requires R **≥ 4.2**. Install Seurat v5 from CRAN plus the recommended performance packages (`presto` for fast marker detection, `BPCells` for on-disk matrices, `glmGamPoi` for fast SCTransform). Full instructions for every install scenario (v4, older versions, dev build, Docker) live in **`install.md`** — read it for environment setup.

```r
install.packages("Seurat")
# Performance + integration extras:
setRepositories(ind = 1:3, addURLs = c(
  "https://satijalab.r-universe.dev", "https://bnprks.r-universe.dev/"))
install.packages(c("BPCells", "presto", "glmGamPoi", "harmony"))
```

For R-native interop and converting between `.rds` and `.h5ad`, see `references/integration.md`. For Python-based analysis of the same data, use the **scanpy** skill; for AnnData format questions use **anndata**.

## When to Use This Skill

This skill should be used when:
- Analyzing single-cell RNA-seq data in **R** (10x Cell Ranger output, `.rds`, `.h5`, `.h5ad`, CSV/MTX)
- Performing quality control and doublet filtering on scRNA-seq datasets
- Normalizing with LogNormalize or SCTransform
- Creating UMAP, t-SNE, or PCA visualizations with Seurat plotting functions
- Identifying cell clusters and finding marker genes (`FindMarkers`/`FindAllMarkers`)
- Integrating multiple samples/batches (Harmony, CCA, RPCA, or Seurat v5 `IntegrateLayers`)
- Annotating cell types from known markers or reference mapping (Azimuth)
- Producing publication-quality single-cell figures in R

## Script Toolkit (prefer these over writing code from scratch)

This skill bundles ready-to-run `Rscript` CLI tools in `scripts/` for every common step. **Run these instead of hand-writing Seurat code** — they handle input loading by extension, sensible defaults, raw-count preservation, figure setup, and progress logging. Each reads and writes a Seurat object as `.rds`, so they chain together, and each supports `--help`.

All scripts source a shared `scripts/_common.R` helper (loading, saving, arg parsing, logging) — keep it alongside the others. Run from the skill directory or pass full paths; figures default to `./figures/`.

| Script | Purpose | Typical call |
|--------|---------|--------------|
| `run_pipeline.R` | **Full workflow in one command**: load → QC → normalize → PCA → (integrate) → UMAP → cluster → markers | `Rscript scripts/run_pipeline.R raw.rds -o processed.rds` |
| `inspect_data.R` | Summarize an unknown object (dims, assays, layers, metadata, what's computed) | `Rscript scripts/inspect_data.R data.rds` |
| `load_data.R` | Load any format (10x dir/.h5, csv, mtx, .h5ad) → write `.rds` | `Rscript scripts/load_data.R 10x_dir/ -o data.rds` |
| `qc_analysis.R` | QC metrics (`percent.mt`), before/after plots, filtering, optional doublet flag | `Rscript scripts/qc_analysis.R raw.rds -o qc.rds` |
| `normalize.R` | LogNormalize or SCTransform + variable features + scaling | `Rscript scripts/normalize.R qc.rds -o norm.rds --method sct` |
| `reduce_dimensions.R` | PCA + elbow plot, neighbors, UMAP, optional t-SNE | `Rscript scripts/reduce_dimensions.R norm.rds -o red.rds --dims 30` |
| `integrate.R` | Multi-sample integration: harmony / cca / rpca | `Rscript scripts/integrate.R red.rds -o int.rds --method harmony --batch sample` |
| `cluster.R` | `FindNeighbors` + `FindClusters` at one or many resolutions | `Rscript scripts/cluster.R red.rds -o clu.rds --resolution 0.3 0.6 1.0` |
| `find_markers.R` | `FindAllMarkers` (presto) + per-cluster CSVs + marker plots | `Rscript scripts/find_markers.R clu.rds -o clu.rds --groupby seurat_clusters` |
| `annotate.R` | Map clusters → cell types from JSON/CSV; optional marker dotplot | `Rscript scripts/annotate.R clu.rds -o ann.rds --mapping map.json` |
| `plot.R` | Generate umap/tsne/violin/dotplot/feature/heatmap plots from an object | `Rscript scripts/plot.R ann.rds --kind dotplot --genes CD3D CD14 --groupby cell_type` |

### One-shot end-to-end run

```bash
# Counts → clustered, marker-annotated object + figures + marker CSVs
Rscript scripts/run_pipeline.R raw.rds -o processed.rds --resolution 0.5 --method lognorm
# With multi-sample integration:
Rscript scripts/run_pipeline.R raw.rds -o processed.rds --batch sample --integrate harmony
# Reproducible parameters via JSON (keys mirror flag names):
Rscript scripts/run_pipeline.R raw.rds -o processed.rds --config assets/pipeline_config.json
```

### Step-by-step chain (when you need to inspect/iterate between stages)

```bash
Rscript scripts/qc_analysis.R       raw.rds  -o qc.rds
Rscript scripts/normalize.R         qc.rds   -o norm.rds --method lognorm
Rscript scripts/reduce_dimensions.R norm.rds -o red.rds  --dims 30
Rscript scripts/cluster.R           red.rds  -o clu.rds  --resolution 0.3 0.5 0.8
Rscript scripts/find_markers.R      clu.rds  -o clu.rds  --groupby seurat_clusters
# inspect results/markers/*.csv, decide labels, write a mapping JSON, then:
Rscript scripts/annotate.R          clu.rds  -o ann.rds  --mapping assets/celltype_mapping.json
```

The sections below document the underlying Seurat calls each script performs — read them when customizing beyond the script flags.

## Quick Start

```r
library(Seurat)

# Load 10x Cell Ranger output (filtered_feature_bc_matrix/)
counts <- Read10X(data.dir = "filtered_feature_bc_matrix/")
obj <- CreateSeuratObject(counts = counts, project = "scRNA",
                          min.cells = 3, min.features = 200)

# Other inputs:
# obj <- readRDS("data.rds")                         # existing Seurat object
# counts <- Read10X_h5("filtered_feature_bc_matrix.h5")
# .h5ad → convert first; see references/integration.md
```

### Seurat v5 object structure

```r
obj                       # prints assays, layers, dims, active assay
obj[["RNA"]]              # the RNA Assay5 object
Layers(obj[["RNA"]])      # e.g. "counts", "data", "scale.data"
obj@meta.data             # per-cell metadata (data.frame); also obj[[]]
obj$seurat_clusters       # a metadata column
Idents(obj)               # active cell identities
Reductions(obj)           # e.g. "pca", "umap"
VariableFeatures(obj)     # selected HVGs
```

In v5, an assay holds **layers** (multiple count/data matrices, one per sample before joining). After per-sample processing, call `JoinLayers(obj)` to merge them for differential expression.

## Standard Analysis Workflow

### 1. Quality Control

**Single sample** — inspect violin/scatter plots, then filter:

```r
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")  # "^mt-" for mouse
VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
```

**Multi-sample** — configure per-sample thresholds, then run DoubletFinder per sample before merging. The recommended pipeline (in `references/standard_workflow.md` Section 2b) is:

1. Per-sample `sample_config` with individual `nFeature_min/max` and `percent_mt_max`
2. Hard threshold filter → `NormalizeData` + `RunPCA` + `RunUMAP` (needed for DoubletFinder)
3. `paramSweep` + `find.pK` automated pK selection; doublet rate ≈ 0.8 % per 1 000 cells
4. `doubletFinder` → remove doublets → rebuild clean `CreateSeuratObject` from filtered counts
5. `merge(..., add.cell.ids = names(processed))` to combine samples

Use the script: `Rscript scripts/qc_analysis.R raw.rds -o qc.rds --min-features 200 --max-mt 10`

### 2. Normalization

**LogNormalize** (classic) or **SCTransform** (recommended; models technical noise):

```r
# LogNormalize path
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj)

# SCTransform path (replaces the three steps above)
obj <- SCTransform(obj, vars.to.regress = "percent.mt")
```

### 3. Dimensionality Reduction

```r
obj <- RunPCA(obj)
ElbowPlot(obj, ndims = 50)          # choose dims
obj <- FindNeighbors(obj, dims = 1:30)
obj <- RunUMAP(obj, dims = 1:30)
DimPlot(obj, reduction = "umap")
```

### 4. Clustering

```r
obj <- FindClusters(obj, resolution = 0.5)   # try 0.3–1.2
DimPlot(obj, label = TRUE)
```

### 5. Marker Genes

```r
obj <- JoinLayers(obj)               # v5: merge layers before DE
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Top markers per cluster:
library(dplyr); markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10)
```

`presto` accelerates the default Wilcoxon test substantially. Per-cell tests are for **exploratory** cluster markers; for rigorous condition comparisons, pseudobulk by sample × cell type and use a bulk DE tool (see `pydeseq2` / `references/integration.md`).

### 6. Cell-Type Annotation

```r
markers_panel <- c("CD3D", "CD14", "MS4A1", "NKG7", "FCGR3A")
DotPlot(obj, features = markers_panel) + RotatedAxis()
FeaturePlot(obj, features = markers_panel)

new_ids <- c("0" = "CD4 T", "1" = "CD14+ Mono", "2" = "B", "3" = "CD8 T")
obj <- RenameIdents(obj, new_ids)
obj$cell_type <- Idents(obj)
DimPlot(obj, group.by = "cell_type", label = TRUE)
```

Automated reference mapping is available via **Azimuth** (`RunAzimuth`). See `references/integration.md`.

### 7. Save

```r
saveRDS(obj, "results/processed.rds")
write.csv(obj@meta.data, "results/cell_metadata.csv")
```

## Multi-Sample Integration

When combining multiple samples/batches, normalize per sample (layers in v5), then integrate to remove batch effects while preserving biology:

```r
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA()
obj <- IntegrateLayers(obj, method = HarmonyIntegration,
                       orig.reduction = "pca", new.reduction = "harmony")
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30)
obj <- JoinLayers(obj)
```

Methods: `HarmonyIntegration` (fast, scalable), `CCAIntegration` (classic, good for divergent datasets), `RPCAIntegration` (faster CCA, conservative). Use the script: `Rscript scripts/integrate.R red.rds -o int.rds --method harmony --batch sample`. Full guidance in `references/integration.md`.

## Key Parameters to Adjust

- **QC**: `min.features` (200–500), `percent.mt` cutoff (5–20%), upper `nFeature_RNA` to drop doublets
- **HVGs**: `nfeatures` (2000–3000)
- **Dims**: `dims = 1:N` chosen from `ElbowPlot` (typically 20–50)
- **Clustering**: `resolution` (0.3–1.2; higher → more clusters)
- **Normalization**: prefer SCTransform for most datasets; LogNormalize for speed/compatibility

## Common Pitfalls and Best Practices

1. **v5 layers**: call `JoinLayers()` before `FindMarkers`/`FindAllMarkers` on multi-layer objects.
2. **Mito pattern**: `^MT-` for human, `^mt-` for mouse — check `rownames(obj)`.
3. **SCTransform DE**: run `PrepSCTFindMarkers()` before marker detection on SCT-integrated objects.
4. **Install `presto`**: dramatically speeds up `FindAllMarkers`.
5. **Integrate, don't merge**: a plain `merge()` does not remove batch effects.
6. **Pseudobulk for condition DE**: per-cell p-values are anti-conservative; aggregate by sample for rigorous comparisons.
7. **Choose dims from ElbowPlot**: too few loses structure, too many adds noise.
8. **Validate annotations** with multiple markers, not one.
9. **Set a seed**: clustering/UMAP are stochastic; scripts set `set.seed()` for reproducibility.
10. **Save checkpoints**: long workflows can fail partway — the scripts write `.rds` at each step.

## Bundled Resources

### scripts/ (CLI toolkit)
Composable `.rds`-in/`.rds`-out `Rscript` tools covering the whole workflow plus a one-command pipeline. See the **Script Toolkit** table above. Each has `--help`. Files:
- `_common.R` — shared loading/saving/arg-parsing/figure helpers sourced by the others (not a CLI)
- `run_pipeline.R` — full pipeline in one command (flags or `--config` JSON)
- `inspect_data.R`, `load_data.R` — explore and load/convert any input format
- `qc_analysis.R`, `normalize.R`, `reduce_dimensions.R`, `integrate.R`, `cluster.R` — pipeline steps
- `find_markers.R`, `annotate.R` — markers and annotation
- `plot.R` — generate any standard plot

**Default to these scripts before writing Seurat code from scratch.**

### references/
- `standard_workflow.md` — complete step-by-step workflow with explanations and code.
- `api_reference.md` — quick lookup of Seurat functions by category (I/O, preprocessing, tools, plotting, object accessors).
- `plotting_guide.md` — publication-quality visualization recipes (DimPlot, FeaturePlot, VlnPlot, DotPlot, DoHeatmap, styling, multi-panel).
- `integration.md` — multi-sample integration methods, SCTransform specifics, Azimuth reference mapping, and `.rds`↔`.h5ad` interop with the scanpy/anndata skills.

### install.md
Installation runbook: Seurat v5 from CRAN, performance packages, Signac/SeuratData/Azimuth/SeuratWrappers, Seurat v4, older versions, dev build, and Docker.

### assets/
- `analysis_template.R` — full end-to-end template to copy and customize.
- `pipeline_config.json` — parameter set for `run_pipeline.R --config`.
- `celltype_mapping.json` — cluster → cell-type map for `annotate.R --mapping`.
- `gene_signatures.json` — marker panels per cell type for annotation/scoring.

## Additional Resources

- **Seurat documentation**: https://satijalab.org/seurat/
- **Seurat v5 vignettes**: https://satijalab.org/seurat/articles/get_started_v5_new
- **Integration vignette**: https://satijalab.org/seurat/articles/seurat5_integration
- **Azimuth**: https://azimuth.hubmapconsortium.org/
- **Best practices**: Luecken & Theis (2019) "Current best practices in single-cell RNA-seq"
