#!/usr/bin/env Rscript
# Normalize and select variable features. Supports LogNormalize (default)
# and SCTransform. LogNormalize path also scales the data.
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.here, "_common.R"))

HELP <- "Usage: Rscript normalize.R <input> -o out.rds [options]
Options:
  -o, --output FILE     output .rds (required)
  --method M            'lognorm' (default) or 'sct'
  --nfeatures N         number of variable features (default 2000)
  --regress VARS        vars.to.regress, e.g. 'percent.mt nCount_RNA'
  --figdir DIR          figure output dir (default figures)"

opts <- parse_args(commandArgs(trailingOnly = TRUE),
  spec = list(output = "string", method = "string", nfeatures = "int",
              regress = "multi", figdir = "string"),
  defaults = list(method = "lognorm", nfeatures = 2000, figdir = "figures"),
  help_text = HELP)
if (length(opts$`_`) < 1) die("need an input path. Try --help")
if (is.null(opts$output)) die("need -o/--output")

suppressMessages(library(Seurat))
obj <- load_seurat(opts$`_`[[1]])
regress <- if (is.null(opts$regress)) NULL else opts$regress

if (tolower(opts$method) == "sct") {
  info("running SCTransform")
  obj <- SCTransform(obj, variable.features.n = opts$nfeatures,
                     vars.to.regress = regress, verbose = FALSE)
} else {
  info("running LogNormalize + FindVariableFeatures + ScaleData")
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = opts$nfeatures, verbose = FALSE)
  obj <- ScaleData(obj, vars.to.regress = regress, verbose = FALSE)
}

top <- head(VariableFeatures(obj), 15)
info("top variable features: ", paste(top, collapse = ", "))
p <- VariableFeaturePlot(obj)
save_plot(p, "variable_features.png", opts$figdir)

save_seurat(obj, opts$output)
info("done")
