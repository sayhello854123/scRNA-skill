# Integration, SCTransform & R↔Python Interop

## Multi-sample integration (Seurat v5)

A plain `merge()` concatenates objects but does **not** remove batch effects.
To combine samples for joint clustering/annotation, split the RNA assay into
per-sample layers, run the standard preprocessing, then `IntegrateLayers`.

```r
library(Seurat)

# Combine samples into one object (layers stay separate in v5)
obj <- merge(s1, y = c(s2, s3), add.cell.ids = c("s1", "s2", "s3"))
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)   # one layer per sample

obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA()

obj <- IntegrateLayers(
  obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)

obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- JoinLayers(obj)        # merge layers before DE/markers
```

### Choosing a method

| Method | `new.reduction` | Notes |
|--------|-----------------|-------|
| `HarmonyIntegration` | `harmony` | Fast, scalable, strong default for many samples |
| `CCAIntegration` | `integrated.cca` | Classic anchors; good when datasets are quite different; slower |
| `RPCAIntegration` | `integrated.rpca` | Faster, more conservative CCA variant |
| `JointPCAIntegration` | `integrated.dr` | Joint PCA anchors |
| `scVIIntegration` | `integrated.scvi` | Deep-learning; needs `scVI`/reticulate (see scvi-tools skill) |

The integrated reduction feeds `FindNeighbors`/`RunUMAP`; the original `pca`
stays available for comparison. Use the script:
`Rscript scripts/integrate.R red.rds -o int.rds --method harmony --batch sample`.

## SCTransform with integration

```r
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)
obj <- SCTransform(obj)
obj <- RunPCA(obj)
obj <- IntegrateLayers(obj, method = RPCAIntegration,
                       normalization.method = "SCT", new.reduction = "integrated.rpca")
# Before markers on SCT data:
obj <- PrepSCTFindMarkers(obj)
```

## Automated annotation (Azimuth)

```r
library(Azimuth)
obj <- RunAzimuth(obj, reference = "pbmcref")   # downloads a curated reference
DimPlot(obj, group.by = "predicted.celltype.l2", label = TRUE)
```

References include `pbmcref`, `lungref`, `kidneyref`, `bonemarrowref`, etc.
Predictions land in `predicted.*` metadata columns with confidence scores.

## R ↔ Python interop (.rds ↔ .h5ad)

Seurat is R-native; scanpy/anndata are Python-native. Convert at the file level.

### .h5ad → Seurat (`.rds`)

```r
# Option A: schard (lightweight, reads h5ad directly to Seurat)
# install.packages("schard")  # or remotes::install_github("cellgeni/schard")
obj <- schard::h5ad2seurat("data.h5ad")

# Option B: zellkonverter (Bioconductor) via SingleCellExperiment
# BiocManager::install("zellkonverter")
sce <- zellkonverter::readH5AD("data.h5ad")
obj <- Seurat::as.Seurat(sce, counts = "X", data = NULL)
```

`scripts/load_data.R` and `_common.R` use `schard` first, falling back to
`zellkonverter`.

### Seurat (`.rds`) → .h5ad

```r
library(Seurat); library(SeuratDisk)
SaveH5Seurat(obj, "data.h5Seurat")
Convert("data.h5Seurat", dest = "h5ad")        # writes data.h5ad

# Or via SingleCellExperiment + zellkonverter:
sce <- Seurat::as.SingleCellExperiment(obj)
zellkonverter::writeH5AD(sce, "data.h5ad")
```

Then analyze in Python with the **scanpy** skill (`sc.read_h5ad("data.h5ad")`).
Preserve raw counts, full metadata, and gene identifiers across the conversion;
verify cell/gene counts match afterward.

## When to hand off to other skills

- **scanpy** — Python scRNA-seq analysis on the converted `.h5ad`.
- **anndata** — AnnData structure/I/O questions.
- **scvi-tools** — deep-learning models (scVI/scANVI) for integration & label transfer.
- **pydeseq2** — rigorous pseudobulk differential expression between conditions.
- **cellxgene-census** — pulling public reference atlases for annotation.
