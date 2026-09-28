#!/usr/bin/env Rscript

# Export cell-level metadata from the annotated Xenium Seurat RDS.
# The Python/Scanpy scripts use this table as the authority for cell IDs,
# sample/age labels, SnC calls, and cell type annotations, then subset the
# raw Xenium H5 matrices to exactly these cells.

suppressPackageStartupMessages({
  library(Seurat)
})

export_xenium_rds_metadata <- function(
    input_rds = file.path("data", "skin.rds"),
    outdir = file.path("results", "260727", "xenium_rds_metadata")) {
  outfile <- file.path(outdir, "xenium_rds_cell_metadata.csv.gz")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message("Reading Seurat RDS: ", input_rds)
  obj <- readRDS(input_rds)

  meta <- obj@meta.data
  meta$cell_id <- rownames(meta)
  meta <- meta[, c("cell_id", setdiff(colnames(meta), "cell_id")), drop = FALSE]

  if (!"sample" %in% colnames(meta)) {
    meta$sample <- sub("_.*$", "", meta$cell_id)
  }

  if (!"age" %in% colnames(meta)) {
    meta$age <- ifelse(grepl("-Y$", meta$sample, ignore.case = TRUE), "Young",
      ifelse(grepl("-O$", meta$sample, ignore.case = TRUE), "Old", NA_character_)
    )
  }

  message("Writing metadata: ", outfile)
  utils::write.csv(meta, gzfile(outfile), row.names = FALSE, quote = TRUE)

  summary <- data.frame(
    n_cells = nrow(meta),
    n_samples = length(unique(meta$sample)),
    samples = paste(sort(unique(meta$sample)), collapse = ", "),
    metadata_columns = paste(colnames(meta), collapse = ", ")
  )
  utils::write.csv(summary, file.path(outdir, "xenium_rds_metadata_export_summary.csv"), row.names = FALSE)

  message("Done. Cells: ", nrow(meta), "; samples: ", paste(sort(unique(meta$sample)), collapse = ", "))
  invisible(meta)
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  input_rds <- if (length(args) >= 1) args[[1]] else file.path("data", "skin.rds")
  outdir <- if (length(args) >= 2) args[[2]] else file.path("results", "260727", "xenium_rds_metadata")
  export_xenium_rds_metadata(input_rds, outdir)
}
