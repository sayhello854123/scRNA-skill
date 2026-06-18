#!/usr/bin/env Rscript
# Quality control: compute percent.mt (+ ribo), plot before/after QC,
# and filter cells by feature counts and mito fraction.
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.here, "_common.R"))

HELP <- "Usage: Rscript qc_analysis.R <input> -o out.rds [options]
Options:
  -o, --output FILE      output .rds (required)
  --min-features N       min genes per cell (default 200)
  --max-features N       max genes per cell, drops likely doublets (default Inf)
  --max-mt PCT           max %% mitochondrial counts (default 10)
  --mt-pattern REGEX     mito gene pattern (default '^MT-'; mouse '^mt-')
  --figdir DIR           figure output dir (default figures)"

opts <- parse_args(commandArgs(trailingOnly = TRUE),
  spec = list(output = "string", `min-features` = "int", `max-features` = "double",
              `max-mt` = "double", `mt-pattern` = "string", figdir = "string"),
  defaults = list(`min-features` = 200, `max-features` = Inf, `max-mt` = 10,
                  `mt-pattern` = "^MT-", figdir = "figures"),
  help_text = HELP)
if (length(opts$`_`) < 1) die("need an input path. Try --help")
if (is.null(opts$output)) die("need -o/--output")

suppressMessages(library(Seurat))
obj <- load_seurat(opts$`_`[[1]])

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = opts$`mt-pattern`)
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

p_before <- VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                    ncol = 3, pt.size = 0)
save_plot(p_before, "qc_violin_before.png", opts$figdir, width = 10, height = 4)

n0 <- ncol(obj)
obj <- subset(obj, subset = nFeature_RNA >= opts$`min-features` &
                            nFeature_RNA <= opts$`max-features` &
                            percent.mt <= opts$`max-mt`)
info(sprintf("filtered %d -> %d cells (removed %d)", n0, ncol(obj), n0 - ncol(obj)))

p_after <- VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                   ncol = 3, pt.size = 0)
save_plot(p_after, "qc_violin_after.png", opts$figdir, width = 10, height = 4)

save_seurat(obj, opts$output)
info("done")
