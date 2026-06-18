#!/usr/bin/env Rscript
# Graph-based clustering: FindNeighbors (if needed) + FindClusters at one or
# more resolutions. Each resolution is stored in its own metadata column.
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.here, "_common.R"))

HELP <- "Usage: Rscript cluster.R <input> -o out.rds [options]
Options:
  -o, --output FILE         output .rds (required)
  --resolution R [R ...]    one or more resolutions (default 0.5)
  --dims N                  dims if neighbors must be (re)built (default 30)
  --reduction NAME          reduction for neighbors (default 'pca')
  --algorithm N             1=Louvain, 4=Leiden (default 1)
  --figdir DIR              figure output dir (default figures)
  --seed N                  random seed (default 42)"

opts <- parse_args(commandArgs(trailingOnly = TRUE),
  spec = list(output = "string", resolution = "multi_double", dims = "int",
              reduction = "string", algorithm = "int", figdir = "string", seed = "int"),
  defaults = list(resolution = 0.5, dims = 30, reduction = "pca",
                  algorithm = 1, figdir = "figures", seed = 42),
  help_text = HELP)
if (length(opts$`_`) < 1) die("need an input path. Try --help")
if (is.null(opts$output)) die("need -o/--output")

suppressMessages(library(Seurat))
set.seed(opts$seed)
obj <- load_seurat(opts$`_`[[1]])

if (length(Seurat::Graphs(obj)) == 0) {
  info("no graph found; running FindNeighbors on '", opts$reduction, "'")
  obj <- FindNeighbors(obj, reduction = opts$reduction, dims = 1:opts$dims, verbose = FALSE)
}

for (res in opts$resolution) {
  key <- sprintf("clusters_res%s", res)
  obj <- FindClusters(obj, resolution = res, algorithm = opts$algorithm,
                      cluster.name = key, verbose = FALSE)
  n <- length(unique(obj@meta.data[[key]]))
  info(sprintf("resolution %.2f -> %d clusters (%s)", res, n, key))
  if ("umap" %in% Reductions(obj))
    save_plot(DimPlot(obj, group.by = key, label = TRUE) + ggplot2::ggtitle(key),
              sprintf("umap_%s.png", key), opts$figdir)
}
# Default identity = last (often highest) resolution requested
last <- sprintf("clusters_res%s", tail(opts$resolution, 1))
obj$seurat_clusters <- obj@meta.data[[last]]
Idents(obj) <- "seurat_clusters"

save_seurat(obj, opts$output)
info("done — active idents set to ", last)
