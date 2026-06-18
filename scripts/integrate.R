#!/usr/bin/env Rscript
# Multi-sample integration via Seurat v5 IntegrateLayers (harmony / cca / rpca).
# Splits the RNA layers by a batch key, integrates, then rebuilds neighbors/UMAP
# on the integrated reduction.
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.here, "_common.R"))

HELP <- "Usage: Rscript integrate.R <input> -o out.rds --batch COL [options]
Options:
  -o, --output FILE   output .rds (required)
  --batch COL         metadata column identifying sample/batch (required)
  --method M          'harmony' (default), 'cca', or 'rpca'
  --dims N            dims for neighbors/UMAP on integrated space (default 30)
  --figdir DIR        figure output dir (default figures)
  --seed N            random seed (default 42)"

opts <- parse_args(commandArgs(trailingOnly = TRUE),
  spec = list(output = "string", batch = "string", method = "string",
              dims = "int", figdir = "string", seed = "int"),
  defaults = list(method = "harmony", dims = 30, figdir = "figures", seed = 42),
  help_text = HELP)
if (length(opts$`_`) < 1) die("need an input path. Try --help")
if (is.null(opts$output)) die("need -o/--output")
if (is.null(opts$batch)) die("need --batch (metadata column with sample IDs)")

suppressMessages(library(Seurat))
set.seed(opts$seed)
obj <- load_seurat(opts$`_`[[1]])
if (!(opts$batch %in% colnames(obj@meta.data))) die("no metadata column: ", opts$batch)

method_fun <- switch(tolower(opts$method),
  harmony = "HarmonyIntegration",
  cca     = "CCAIntegration",
  rpca    = "RPCAIntegration",
  die("unknown --method: ", opts$method))
new_red <- switch(tolower(opts$method),
  harmony = "harmony", cca = "integrated.cca", rpca = "integrated.rpca")

info("splitting RNA layers by '", opts$batch, "' and normalizing per layer")
obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[opts$batch]])
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)

info("integrating with ", method_fun)
obj <- IntegrateLayers(obj, method = get(method_fun),
                       orig.reduction = "pca", new.reduction = new_red,
                       verbose = FALSE)

d <- 1:opts$dims
obj <- FindNeighbors(obj, reduction = new_red, dims = d, verbose = FALSE)
obj <- RunUMAP(obj, reduction = new_red, dims = d, verbose = FALSE)
obj <- JoinLayers(obj)

save_plot(DimPlot(obj, reduction = "umap", group.by = opts$batch),
          "umap_integrated_by_batch.png", opts$figdir, width = 9)

save_seurat(obj, opts$output)
info("done — integrated reduction: ", new_red)
