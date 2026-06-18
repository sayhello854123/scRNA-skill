#!/usr/bin/env Rscript
# Load any supported input (10x dir/.h5, csv, .rds, .h5ad) and write a .rds
# Seurat object. Optionally tag a sample/batch name into metadata.
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.here, "_common.R"))

HELP <- "Usage: Rscript load_data.R <input> -o out.rds [options]
Inputs: 10x directory, .h5 (10x), .csv (genes x cells), .rds, .h5ad
Options:
  -o, --output FILE    output .rds (required)
  --project NAME       project name (default 'scRNA')
  --sample NAME        write a 'sample' metadata column with this value
  --min-cells N        min cells per gene (default 3)
  --min-features N     min features per cell (default 200)"

opts <- parse_args(commandArgs(trailingOnly = TRUE),
  spec = list(output = "string", project = "string", sample = "string",
              `min-cells` = "int", `min-features` = "int"),
  defaults = list(project = "scRNA", `min-cells` = 3, `min-features` = 200),
  help_text = HELP)
if (length(opts$`_`) < 1) die("need an input path. Try --help")
if (is.null(opts$output)) die("need -o/--output")

obj <- load_seurat(opts$`_`[[1]])
obj@project.name <- opts$project
if (!is.null(opts$sample)) obj$sample <- opts$sample
save_seurat(obj, opts$output)
info("done")
