#!/usr/bin/env Rscript
# Generate a standard single-cell plot from a processed object.
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.here, "_common.R"))

HELP <- "Usage: Rscript plot.R <input> --kind KIND [options]
Kinds: umap | tsne | violin | feature | dotplot | heatmap
Options:
  --kind KIND         plot type (required)
  --genes G [G ...]   genes (for violin/feature/dotplot/heatmap)
  --groupby COL       grouping/identity column (default 'seurat_clusters')
  --out FILE          figure filename (default '<kind>.png')
  --figdir DIR        figure output dir (default figures)
  --width N --height N inches"

opts <- parse_args(commandArgs(trailingOnly = TRUE),
  spec = list(kind = "string", genes = "multi", groupby = "string",
              out = "string", figdir = "string", width = "double", height = "double"),
  defaults = list(groupby = "seurat_clusters", figdir = "figures",
                  width = 8, height = 6),
  help_text = HELP)
if (length(opts$`_`) < 1) die("need an input path. Try --help")
if (is.null(opts$kind)) die("need --kind")

suppressMessages(library(Seurat))
obj <- load_seurat(opts$`_`[[1]])
g <- opts$genes
needs_genes <- opts$kind %in% c("violin", "feature", "dotplot", "heatmap")
if (needs_genes && is.null(g)) die("--kind ", opts$kind, " needs --genes")

p <- switch(opts$kind,
  umap    = DimPlot(obj, reduction = "umap", group.by = opts$groupby, label = TRUE),
  tsne    = DimPlot(obj, reduction = "tsne", group.by = opts$groupby, label = TRUE),
  violin  = VlnPlot(obj, features = g, group.by = opts$groupby, pt.size = 0),
  feature = FeaturePlot(obj, features = g),
  dotplot = DotPlot(obj, features = g, group.by = opts$groupby) + RotatedAxis(),
  heatmap = DoHeatmap(obj, features = g, group.by = opts$groupby),
  die("unknown --kind: ", opts$kind))

out <- if (is.null(opts$out)) sprintf("%s.png", opts$kind) else opts$out
save_plot(p, out, opts$figdir, width = opts$width, height = opts$height)
info("done")
