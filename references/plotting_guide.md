# Seurat Plotting Guide

All Seurat plotting functions return `ggplot` objects, so they can be themed,
combined with `patchwork`, and saved with `ggsave()`. Use `scripts/plot.R` for
one-off plots from a processed object.

## Saving figures

```r
library(ggplot2)
p <- DimPlot(obj, label = TRUE)
ggsave("figures/umap.png", p, width = 8, height = 6, dpi = 300)
ggsave("figures/umap.pdf", p, width = 8, height = 6)   # vector for publication
```

## Embeddings (UMAP / t-SNE / PCA)

```r
DimPlot(obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE)
DimPlot(obj, reduction = "umap", split.by = "condition")     # one panel per condition
DimPlot(obj, reduction = "umap", group.by = "sample")        # check batch mixing
```

## Expression on the embedding

```r
FeaturePlot(obj, features = c("CD3D", "CD14", "MS4A1", "NKG7"))
FeaturePlot(obj, "CD3D", split.by = "condition", order = TRUE)  # high cells on top
FeaturePlot(obj, "CD3D", min.cutoff = "q10", max.cutoff = "q90")  # clip outliers
```

## Marker visualization

```r
# Dot plot: dot size = % cells expressing, color = scaled mean expression
DotPlot(obj, features = markers, group.by = "cell_type") + RotatedAxis()

# Violin plots
VlnPlot(obj, features = c("CD3D", "CD14"), group.by = "cell_type", pt.size = 0)
VlnPlot(obj, "CD3D", group.by = "cell_type", split.by = "condition")

# Heatmap of top markers (uses scale.data)
library(dplyr)
top <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 10) %>% pull(gene)
DoHeatmap(obj, features = unique(top)) + NoLegend()
```

## Diagnostics

```r
VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
FeatureScatter(obj, "nCount_RNA", "nFeature_RNA")
ElbowPlot(obj, ndims = 50)               # choose PCs
VariableFeaturePlot(obj)                 # HVG mean-variance
DimHeatmap(obj, dims = 1:9, cells = 500, balanced = TRUE)
```

## Composing and styling

```r
library(patchwork)
p1 <- DimPlot(obj, group.by = "cell_type")
p2 <- FeaturePlot(obj, "CD3D")
(p1 | p2)                                 # side by side
(p1 / p2) + plot_annotation(tag_levels = "A")

# Theme tweaks (ggplot)
DimPlot(obj) +
  ggtitle("PBMC clusters") +
  theme(legend.position = "bottom") +
  NoAxes()

# Custom colors
DimPlot(obj, cols = c("CD4 T" = "#E41A1C", "B" = "#377EB8"))
FeaturePlot(obj, "CD3D", cols = c("lightgrey", "darkred"))
```

## Publication tips

- Export PDF/SVG (vector) for figures; PNG at `dpi = 300` for rasters.
- `label = TRUE, repel = TRUE` for readable on-plot cluster labels.
- Use consistent palettes across panels (define a named color vector once).
- `order = TRUE` in `FeaturePlot` plots expressing cells on top.
- `RotatedAxis()` keeps long gene names legible in `DotPlot`.
- `pt.size = 0` removes jitter dots from violins when there are many cells.
