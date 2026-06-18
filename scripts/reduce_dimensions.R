#!/usr/bin/env Rscript
# PCA + elbow plot, neighborhood graph, UMAP, and optional t-SNE.
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.here, "_common.R"))

HELP <- "Usage: Rscript reduce_dimensions.R <input> -o out.rds [options]
Options:
  -o, --output FILE   output .rds (required)
  --dims N            number of PCs to use downstream (default 30)
  --npcs N            number of PCs to compute (default 50)
  --tsne              also run t-SNE
  --reduction NAME    base reduction for neighbors/UMAP (default 'pca')
  --figdir DIR        figure output dir (default figures)
  --seed N            random seed (default 42)"

opts <- parse_args(commandArgs(trailingOnly = TRUE),
  spec = list(output = "string", dims = "int", npcs = "int", tsne = "bool",
              reduction = "string", figdir = "string", seed = "int"),
  defaults = list(dims = 30, npcs = 50, reduction = "pca", figdir = "figures", seed = 42),
  help_text = HELP)
if (length(opts$`_`) < 1) die("need an input path. Try --help")
if (is.null(opts$output)) die("need -o/--output")

suppressMessages(library(Seurat))
set.seed(opts$seed)
obj <- load_seurat(opts$`_`[[1]])

if (!("pca" %in% Reductions(obj))) {
  info("running PCA (", opts$npcs, " PCs)")
  obj <- RunPCA(obj, npcs = opts$npcs, verbose = FALSE)
}
save_plot(ElbowPlot(obj, ndims = opts$npcs), "elbow.png", opts$figdir)

d <- 1:opts$dims
info("FindNeighbors + RunUMAP on '", opts$reduction, "' dims 1:", opts$dims)
obj <- FindNeighbors(obj, reduction = opts$reduction, dims = d, verbose = FALSE)
obj <- RunUMAP(obj, reduction = opts$reduction, dims = d, verbose = FALSE)
save_plot(DimPlot(obj, reduction = "umap"), "umap.png", opts$figdir)

if (isTRUE(opts$tsne)) {
  obj <- RunTSNE(obj, reduction = opts$reduction, dims = d)
  save_plot(DimPlot(obj, reduction = "tsne"), "tsne.png", opts$figdir)
}

save_seurat(obj, opts$output)
info("done")
