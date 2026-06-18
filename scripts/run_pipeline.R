#!/usr/bin/env Rscript
# Full single-cell workflow in one command:
#   load -> QC -> normalize -> PCA -> (optional integrate) -> UMAP -> cluster -> markers
# Flags may also be supplied via a --config JSON (keys mirror flag names).
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.here, "_common.R"))

HELP <- "Usage: Rscript run_pipeline.R <input> -o out.rds [options]
Options:
  -o, --output FILE     output .rds (required)
  --config FILE         JSON of parameters (CLI flags override)
  --method M            normalization: 'lognorm' (default) or 'sct'
  --nfeatures N         variable features (default 2000)
  --max-mt PCT          max %% mito for QC (default 10)
  --min-features N      min genes/cell for QC (default 200)
  --mt-pattern REGEX    mito pattern (default '^MT-')
  --dims N              PCs for neighbors/UMAP (default 30)
  --resolution R [R..]  clustering resolution(s) (default 0.5)
  --batch COL           if set with --integrate, integrate by this column
  --integrate M         integration method: harmony / cca / rpca
  --groupby COL         marker grouping column (default 'seurat_clusters')
  --figdir DIR          figure output dir (default figures)
  --seed N              random seed (default 42)"

raw <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(raw,
  spec = list(output = "string", config = "string", method = "string",
              nfeatures = "int", `max-mt` = "double", `min-features` = "int",
              `mt-pattern` = "string", dims = "int", resolution = "multi_double",
              batch = "string", integrate = "string", groupby = "string",
              figdir = "string", seed = "int"),
  defaults = list(method = "lognorm", nfeatures = 2000, `max-mt` = 10,
                  `min-features` = 200, `mt-pattern` = "^MT-", dims = 30,
                  resolution = 0.5, groupby = "seurat_clusters",
                  figdir = "figures", seed = 42),
  help_text = HELP)
opts <- apply_config(opts, raw)
if (length(opts$`_`) < 1) die("need an input path. Try --help")
if (is.null(opts$output)) die("need -o/--output")

suppressMessages({ library(Seurat); library(dplyr) })
set.seed(opts$seed)
fig <- opts$figdir

obj <- load_seurat(opts$`_`[[1]])

## 1. QC ----------------------------------------------------------------------
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = opts$`mt-pattern`)
save_plot(VlnPlot(obj, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0),
          "qc_violin.png", fig, width = 10, height = 4)
n0 <- ncol(obj)
obj <- subset(obj, subset = nFeature_RNA >= opts$`min-features` & percent.mt <= opts$`max-mt`)
info(sprintf("QC: %d -> %d cells", n0, ncol(obj)))

do_integrate <- !is.null(opts$integrate) && !is.null(opts$batch)

## 2. Normalize ---------------------------------------------------------------
if (do_integrate) {
  info("splitting by '", opts$batch, "' for integration")
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[opts$batch]])
}
if (tolower(opts$method) == "sct") {
  obj <- SCTransform(obj, variable.features.n = opts$nfeatures, verbose = FALSE)
} else {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = opts$nfeatures, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
}

## 3. PCA + (integration) + UMAP ----------------------------------------------
obj <- RunPCA(obj, verbose = FALSE)
save_plot(ElbowPlot(obj, ndims = 50), "elbow.png", fig)
red <- "pca"
if (do_integrate) {
  mf <- switch(tolower(opts$integrate),
    harmony = "HarmonyIntegration", cca = "CCAIntegration", rpca = "RPCAIntegration",
    die("unknown --integrate: ", opts$integrate))
  red <- switch(tolower(opts$integrate),
    harmony = "harmony", cca = "integrated.cca", rpca = "integrated.rpca")
  info("integrating with ", mf)
  obj <- IntegrateLayers(obj, method = get(mf), orig.reduction = "pca",
                         new.reduction = red, verbose = FALSE)
}
d <- 1:opts$dims
obj <- FindNeighbors(obj, reduction = red, dims = d, verbose = FALSE)
obj <- RunUMAP(obj, reduction = red, dims = d, verbose = FALSE)

## 4. Cluster -----------------------------------------------------------------
for (res in opts$resolution) {
  key <- sprintf("clusters_res%s", res)
  obj <- FindClusters(obj, resolution = res, cluster.name = key, verbose = FALSE)
  info(sprintf("res %.2f -> %d clusters", res, length(unique(obj@meta.data[[key]]))))
}
obj$seurat_clusters <- obj@meta.data[[sprintf("clusters_res%s", tail(opts$resolution, 1))]]
Idents(obj) <- "seurat_clusters"
save_plot(DimPlot(obj, label = TRUE), "umap_clusters.png", fig)

## 5. Markers -----------------------------------------------------------------
if ("RNA" %in% Seurat::Assays(obj)) obj <- JoinLayers(obj)
if (DefaultAssay(obj) == "SCT") obj <- PrepSCTFindMarkers(obj)
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
dir.create("results/markers", showWarnings = FALSE, recursive = TRUE)
write.csv(markers, "results/markers/all_markers.csv", row.names = FALSE)
top <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 5) %>%
  pull(gene) %>% unique()
if (length(top))
  save_plot(DotPlot(obj, features = top) + RotatedAxis(), "marker_dotplot.png",
            fig, width = max(8, length(top) * 0.35))

save_seurat(obj, opts$output)
info("pipeline complete — object: ", opts$output,
     " | markers: results/markers/ | figures: ", fig, "/")
