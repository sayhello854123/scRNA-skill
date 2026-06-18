#!/usr/bin/env Rscript
# FindAllMarkers across groups (presto-accelerated), writing a combined CSV,
# per-group CSVs, and a top-marker dotplot. Use for exploratory cluster markers.
.here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(.here, "_common.R"))

HELP <- "Usage: Rscript find_markers.R <input> -o out.rds [options]
Options:
  -o, --output FILE     output .rds (markers also written to results/markers/)
  --groupby COL         identity column to test (default 'seurat_clusters')
  --min-pct X           min.pct (default 0.25)
  --logfc X             logfc.threshold (default 0.25)
  --top N               top genes per group for the dotplot (default 5)
  --only-pos            only positive markers (default TRUE; pass to keep)
  --outdir DIR          CSV output dir (default results/markers)
  --figdir DIR          figure output dir (default figures)"

opts <- parse_args(commandArgs(trailingOnly = TRUE),
  spec = list(output = "string", groupby = "string", `min-pct` = "double",
              logfc = "double", top = "int", `only-pos` = "bool",
              outdir = "string", figdir = "string"),
  defaults = list(groupby = "seurat_clusters", `min-pct` = 0.25, logfc = 0.25,
                  top = 5, `only-pos` = TRUE, outdir = "results/markers",
                  figdir = "figures"),
  help_text = HELP)
if (length(opts$`_`) < 1) die("need an input path. Try --help")
if (is.null(opts$output)) die("need -o/--output")

suppressMessages({ library(Seurat); library(dplyr) })
obj <- load_seurat(opts$`_`[[1]])
if (!(opts$groupby %in% colnames(obj@meta.data))) die("no column: ", opts$groupby)
Idents(obj) <- opts$groupby

# v5: merge layers before DE; SCT objects need PrepSCTFindMarkers.
if ("RNA" %in% Seurat::Assays(obj)) obj <- JoinLayers(obj)
if (DefaultAssay(obj) == "SCT") obj <- PrepSCTFindMarkers(obj)

if (!requireNamespace("presto", quietly = TRUE))
  info("note: 'presto' not installed — FindAllMarkers will be slower")

info("FindAllMarkers on '", opts$groupby, "'")
markers <- FindAllMarkers(obj, only.pos = opts$`only-pos`,
                          min.pct = opts$`min-pct`, logfc.threshold = opts$logfc,
                          verbose = FALSE)

dir.create(opts$outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(markers, file.path(opts$outdir, "all_markers.csv"), row.names = FALSE)
for (g in unique(markers$cluster)) {
  sub <- markers %>% filter(cluster == g) %>% arrange(desc(avg_log2FC))
  write.csv(sub, file.path(opts$outdir, sprintf("markers_%s.csv", g)), row.names = FALSE)
}
info("wrote marker CSVs to ", opts$outdir)

top <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = opts$top) %>%
  pull(gene) %>% unique()
if (length(top))
  save_plot(DotPlot(obj, features = top, group.by = opts$groupby) +
              Seurat::RotatedAxis(),
            "marker_dotplot.png", opts$figdir,
            width = max(8, length(top) * 0.35), height = 6)

save_seurat(obj, opts$output)
info("done")
