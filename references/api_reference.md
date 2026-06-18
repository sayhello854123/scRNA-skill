# Seurat v5 API Quick Reference

Function lookup by category. See https://satijalab.org/seurat/ for full docs.

## Reading / writing

| Function | Purpose |
|----------|---------|
| `Read10X(data.dir)` | Load 10x matrix directory → sparse counts |
| `Read10X_h5(file)` | Load 10x `.h5` → sparse counts |
| `CreateSeuratObject(counts, min.cells, min.features, project)` | Build object |
| `readRDS()` / `saveRDS()` | Load / save a Seurat object as `.rds` |
| `ReadMtx(mtx, cells, features)` | Load arbitrary MTX triplet |
| `as.Seurat(sce)` / `as.SingleCellExperiment(obj)` | Convert to/from SCE |

## Object accessors

| Expression | Returns |
|------------|---------|
| `obj[["RNA"]]` | the assay object |
| `Layers(obj[["RNA"]])` | layer names (`counts`, `data`, `scale.data`) |
| `LayerData(obj, "counts")` | a specific matrix |
| `obj@meta.data`, `obj[[]]` | per-cell metadata data.frame |
| `obj$col`, `obj[["col"]]` | one metadata column |
| `Idents(obj)` / `Idents(obj) <- "col"` | get/set active identities |
| `Assays(obj)`, `DefaultAssay(obj)` | assay names / active assay |
| `Reductions(obj)`, `Graphs(obj)` | reduction / graph names |
| `VariableFeatures(obj)` | selected HVGs |
| `Cells(obj)`, `rownames(obj)`, `ncol(obj)`, `nrow(obj)` | cells / genes / counts |
| `JoinLayers(obj)` | merge per-sample layers (v5) before DE |
| `split(obj[["RNA"]], f = ...)` | split layers by a factor (for integration) |

## Preprocessing

| Function | Purpose |
|----------|---------|
| `PercentageFeatureSet(obj, pattern)` | % counts matching a gene pattern (e.g. mito) |
| `subset(obj, subset = ...)` | filter cells/genes by expression on metadata |
| `NormalizeData(obj)` | log-normalize |
| `FindVariableFeatures(obj, nfeatures)` | select HVGs |
| `ScaleData(obj, vars.to.regress)` | z-score, optional covariate regression |
| `SCTransform(obj, vars.to.regress)` | regularized NB normalization (→ SCT assay) |

## Dimensionality reduction & clustering

| Function | Purpose |
|----------|---------|
| `RunPCA(obj, npcs)` | principal components |
| `RunUMAP(obj, dims, reduction)` | UMAP embedding |
| `RunTSNE(obj, dims)` | t-SNE embedding |
| `FindNeighbors(obj, dims, reduction)` | KNN/SNN graph |
| `FindClusters(obj, resolution, algorithm)` | graph clustering (1=Louvain, 4=Leiden) |
| `IntegrateLayers(obj, method, orig.reduction, new.reduction)` | batch integration |

Integration `method` values: `HarmonyIntegration`, `CCAIntegration`,
`RPCAIntegration`, `JointPCAIntegration`, `scVIIntegration`.

## Differential expression / markers

| Function | Purpose |
|----------|---------|
| `FindMarkers(obj, ident.1, ident.2)` | DE between two groups |
| `FindAllMarkers(obj, only.pos, min.pct, logfc.threshold)` | markers per cluster |
| `FindConservedMarkers(obj, ident.1, grouping.var)` | markers consistent across conditions |
| `PrepSCTFindMarkers(obj)` | required before marker tests on SCT objects |
| `AggregateExpression(obj, group.by, return.seurat)` | pseudobulk |

Test methods (`test.use`): `wilcox` (default, fast with `presto`), `MAST`,
`DESeq2`, `roc`, `LR`, `negbinom`, `poisson`.

## Annotation / scoring

| Function | Purpose |
|----------|---------|
| `RenameIdents(obj, ...)` | relabel identities |
| `AddModuleScore(obj, features, name)` | gene-set score per cell |
| `CellCycleScoring(obj, s.features, g2m.features)` | cell-cycle phase |
| `RunAzimuth(obj, reference)` | automated reference mapping (Azimuth package) |

## Plotting

| Function | Purpose |
|----------|---------|
| `DimPlot(obj, reduction, group.by, label)` | scatter of an embedding by group |
| `FeaturePlot(obj, features)` | expression overlaid on embedding |
| `VlnPlot(obj, features, group.by)` | violin per group |
| `DotPlot(obj, features, group.by)` | mean expr + % expressing |
| `DoHeatmap(obj, features, group.by)` | single-cell heatmap |
| `FeatureScatter(obj, feature1, feature2)` | scatter of two features |
| `ElbowPlot(obj)`, `VariableFeaturePlot(obj)`, `DimHeatmap(obj)` | diagnostics |
| `RotatedAxis()`, `NoLegend()`, `+ ggtitle()` | ggplot-style modifiers |

All plotting functions return `ggplot` objects; combine with `patchwork`
(`p1 | p2`, `p1 / p2`) and save with `ggplot2::ggsave()`.
