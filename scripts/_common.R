#!/usr/bin/env Rscript
# Shared helpers for the Seurat script toolkit.
#
# Every CLI script in this directory sources this file so that data loading,
# saving, argument parsing, figure setup, and logging behave consistently.
# This file is NOT a CLI itself; source it:
#
#     source(file.path(dirname(this_file()), "_common.R"))

suppressWarnings(suppressMessages({
  # Seurat is loaded lazily by callers that need it.
}))

# ---- logging -----------------------------------------------------------------

info <- function(...) cat(sprintf("[seurat] %s\n", paste0(...)), file = stderr())

die <- function(...) {
  cat(sprintf("Error: %s\n", paste0(...)), file = stderr())
  quit(status = 1, save = "no")
}

# ---- tiny argument parser ----------------------------------------------------
# Supports: positional args, --flag value, --flag value1 value2 (multi),
# and boolean --flag (no value). Returns a named list; positionals in $`_`.
#
# spec: named list where each name is a long flag (without --) and the value is
# one of "string", "int", "double", "bool", "multi" (space-separated values),
# "multi_int", "multi_double". Defaults are supplied separately.
parse_args <- function(args, spec, defaults = list(), help_text = NULL) {
  if (any(args %in% c("-h", "--help"))) {
    if (!is.null(help_text)) cat(help_text, "\n")
    quit(status = 0, save = "no")
  }
  out <- defaults
  pos <- c()
  i <- 1
  short <- c("-o" = "output")  # common alias
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% names(short)) a <- paste0("--", short[[a]])
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      type <- spec[[key]]
      if (is.null(type)) die("unknown flag: ", a)
      if (type == "bool") {
        out[[key]] <- TRUE; i <- i + 1; next
      }
      vals <- c()
      j <- i + 1
      while (j <= length(args) && !startsWith(args[[j]], "--") &&
             !(args[[j]] %in% names(short))) {
        vals <- c(vals, args[[j]]); j <- j + 1
      }
      if (length(vals) == 0) die("flag ", a, " expects a value")
      out[[key]] <- switch(type,
        string       = vals[[1]],
        int          = as.integer(vals[[1]]),
        double       = as.numeric(vals[[1]]),
        multi        = vals,
        multi_int    = as.integer(vals),
        multi_double = as.numeric(vals),
        die("bad spec type for ", key))
      i <- j
    } else {
      pos <- c(pos, a); i <- i + 1
    }
  }
  out[["_"]] <- pos
  out
}

# Merge values from a --config JSON file (flag names as keys) as lower-priority
# defaults than explicit CLI flags. Call after parse_args using the raw args to
# know which keys were set on the command line.
apply_config <- function(opts, args) {
  if (is.null(opts$config)) return(opts)
  if (!requireNamespace("jsonlite", quietly = TRUE))
    die("--config needs the 'jsonlite' package (install.packages('jsonlite'))")
  cfg <- jsonlite::read_json(opts$config, simplifyVector = TRUE)
  cli_keys <- sub("^--", "", grep("^--", args, value = TRUE))
  for (k in names(cfg)) {
    if (!(k %in% cli_keys)) opts[[k]] <- cfg[[k]]
  }
  opts
}

# ---- IO ----------------------------------------------------------------------

load_seurat <- function(path) {
  if (!file.exists(path) && !dir.exists(path)) die("input not found: ", path)
  suppressMessages(library(Seurat))
  ext <- tolower(tools::file_ext(path))
  if (dir.exists(path)) {
    info("loading 10x directory: ", path)
    counts <- Read10X(data.dir = path)
    return(CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200))
  }
  obj <- switch(ext,
    rds   = readRDS(path),
    h5    = { counts <- Read10X_h5(path)
              CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200) },
    h5ad  = load_h5ad(path),
    csv   = { m <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
              CreateSeuratObject(counts = m, min.cells = 3, min.features = 200) },
    mtx   = die("for MTX, pass the containing 10x directory instead"),
    die("unsupported extension: .", ext))
  if (!inherits(obj, "Seurat"))
    die("loaded object is not a Seurat object (got ", class(obj)[1], ")")
  info("loaded Seurat object: ", ncol(obj), " cells x ", nrow(obj), " features")
  obj
}

load_h5ad <- function(path) {
  if (!requireNamespace("schard", quietly = TRUE) &&
      !requireNamespace("zellkonverter", quietly = TRUE))
    die("reading .h5ad needs 'schard' or 'zellkonverter'. ",
        "See references/integration.md for conversion.")
  if (requireNamespace("schard", quietly = TRUE)) {
    return(schard::h5ad2seurat(path))
  }
  sce <- zellkonverter::readH5AD(path)
  Seurat::as.Seurat(sce, counts = "X", data = NULL)
}

save_seurat <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, path)
  info("wrote ", path)
}

# ---- figures -----------------------------------------------------------------

figdir <- function(d = NULL) {
  d <- if (is.null(d)) "figures" else d
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

save_plot <- function(plot, name, dir = "figures", width = 8, height = 6, dpi = 300) {
  d <- figdir(dir)
  fp <- file.path(d, name)
  suppressMessages(ggplot2::ggsave(fp, plot = plot, width = width,
                                   height = height, dpi = dpi))
  info("saved figure ", fp)
  invisible(fp)
}

# Path of the currently-running script, so sourcing _common.R works regardless
# of the caller's working directory.
this_file <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  "."
}
