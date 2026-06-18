#!/usr/bin/env Rscript
# Map clusters -> cell-type labels from a JSON/CSV mapping, store as a
# 'cell_type' metadata column, and plot the annotated UMAP.
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.here, "_common.R"))

HELP <- "Usage: Rscript annotate.R <input> -o out.rds --mapping FILE [options]
Mapping file: JSON object {\"0\": \"CD4 T\", \"1\": \"B\", ...} or CSV with
columns cluster,cell_type.
Options:
  -o, --output FILE   output .rds (required)
  --mapping FILE      cluster -> cell-type map, JSON or CSV (required)
  --groupby COL       cluster column to map from (default 'seurat_clusters')
  --label COL         new metadata column name (default 'cell_type')
  --figdir DIR        figure output dir (default figures)"

opts <- parse_args(commandArgs(trailingOnly = TRUE),
  spec = list(output = "string", mapping = "string", groupby = "string",
              label = "string", figdir = "string"),
  defaults = list(groupby = "seurat_clusters", label = "cell_type", figdir = "figures"),
  help_text = HELP)
if (length(opts$`_`) < 1) die("need an input path. Try --help")
if (is.null(opts$output)) die("need -o/--output")
if (is.null(opts$mapping)) die("need --mapping")

suppressMessages(library(Seurat))
obj <- load_seurat(opts$`_`[[1]])
if (!(opts$groupby %in% colnames(obj@meta.data))) die("no column: ", opts$groupby)

ext <- tolower(tools::file_ext(opts$mapping))
mapping <- if (ext == "json") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) die("need 'jsonlite' for JSON mapping")
  unlist(jsonlite::read_json(opts$mapping, simplifyVector = TRUE))
} else {
  df <- read.csv(opts$mapping, colClasses = "character")
  setNames(df$cell_type, df$cluster)
}

clusters <- as.character(obj@meta.data[[opts$groupby]])
labels <- mapping[clusters]
missing <- setdiff(unique(clusters), names(mapping))
if (length(missing)) {
  info("WARNING: no mapping for clusters: ", paste(missing, collapse = ", "),
       " (left as-is)")
  labels[is.na(labels)] <- clusters[is.na(labels)]
}
obj@meta.data[[opts$label]] <- labels
Idents(obj) <- opts$label

save_plot(DimPlot(obj, group.by = opts$label, label = TRUE, repel = TRUE),
          "umap_celltypes.png", opts$figdir, width = 9)

save_seurat(obj, opts$output)
info("done — annotated ", length(unique(labels)), " cell types into '", opts$label, "'")
