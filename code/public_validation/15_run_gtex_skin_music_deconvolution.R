#!/usr/bin/env Rscript

# MuSiC deconvolution of GTEx bulk skin using the same donor-balanced cells,
# Broad cell-type labels, and raw counts used for the custom CIBERSORTx input.

suppressPackageStartupMessages({
  library(MuSiC)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
reference_dir <- file.path(project_root, "results", "public_validation", "gtex_skin_cibersortx", "reference")
input_dir <- file.path(project_root, "results", "public_validation", "gtex_skin_cibersortx", "website_inputs")
output_dir <- file.path(project_root, "results", "public_validation", "gtex_skin_music")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

reference_path <- file.path(reference_dir, "gtex_skin_cibersortx_single_cell_reference.tsv")
cell_manifest_path <- file.path(reference_dir, "selected_cell_manifest.csv")
bulk_path <- file.path(input_dir, "gtex_bulk_skin_music_raw_counts.tsv")
bulk_manifest_path <- file.path(input_dir, "bulk_sample_manifest.csv")
stopifnot(file.exists(reference_path), file.exists(cell_manifest_path),
          file.exists(bulk_path), file.exists(bulk_manifest_path))

message("Reading single-cell raw-count reference")
sc_dt <- fread(reference_path, check.names = FALSE)
sc_genes <- sc_dt[[1]]
sc_columns <- names(sc_dt)[-1]
sc_counts <- Matrix(as.matrix(sc_dt[, -1]), sparse = TRUE)
rownames(sc_counts) <- sc_genes
colnames(sc_counts) <- sc_columns
rm(sc_dt); gc()

cell_manifest <- fread(cell_manifest_path)
cell_manifest <- cell_manifest[match(sc_columns, cibersortx_column)]
if (anyNA(cell_manifest$cibersortx_column)) stop("Reference columns do not match selected-cell manifest")
sc_sce <- SingleCellExperiment(
  assays = list(counts = sc_counts),
  colData = DataFrame(
    sampleID = as.character(cell_manifest$donor),
    cellType = as.character(cell_manifest$cell_type),
    sourceCellID = as.character(cell_manifest$source_cell_id),
    row.names = sc_columns
  )
)
rm(sc_counts); gc()

reference_audit <- as.data.table(as.data.frame(colData(sc_sce)))[,
  .(n_cells = .N), by = .(sampleID, cellType)]
fwrite(reference_audit, file.path(output_dir, "music_reference_donor_type_audit.csv"))

message("Reading aligned GTEx bulk raw counts")
bulk_dt <- fread(bulk_path, check.names = FALSE)
bulk_genes <- bulk_dt[[1]]
bulk_samples <- names(bulk_dt)[-1]
bulk_counts <- as.matrix(bulk_dt[, -1])
rownames(bulk_counts) <- bulk_genes
colnames(bulk_counts) <- bulk_samples
storage.mode(bulk_counts) <- "numeric"
rm(bulk_dt); gc()

common_genes <- intersect(rownames(sc_sce), rownames(bulk_counts))
if (length(common_genes) < 1000) stop("Too few shared genes: ", length(common_genes))
sc_sce <- sc_sce[common_genes, ]
bulk_counts <- bulk_counts[common_genes, , drop = FALSE]
message(sprintf("Running MuSiC with %d genes, %d cells, %d donors, %d cell types, and %d bulk samples",
                length(common_genes), ncol(sc_sce), length(unique(sc_sce$sampleID)),
                length(unique(sc_sce$cellType)), ncol(bulk_counts)))

set.seed(20260920)
music_result <- music_prop(
  bulk.mtx = bulk_counts,
  sc.sce = sc_sce,
  clusters = "cellType",
  samples = "sampleID",
  select.ct = NULL,
  verbose = FALSE
)

weighted <- as.data.frame(music_result$Est.prop.weighted, check.names = FALSE)
weighted <- cbind(sample_id = rownames(weighted), weighted)
allgene <- as.data.frame(music_result$Est.prop.allgene, check.names = FALSE)
allgene <- cbind(sample_id = rownames(allgene), allgene)
fwrite(weighted, file.path(output_dir, "gtex_skin_music_weighted_proportions.csv"))
fwrite(allgene, file.path(output_dir, "gtex_skin_music_allgene_nnls_proportions.csv"))

estimated_types <- colnames(music_result$Est.prop.weighted)
type_audit <- reference_audit[, .(n_donors = uniqueN(sampleID), n_cells = sum(n_cells)), by = cellType]
type_audit[, estimated_by_music := cellType %in% estimated_types]
type_audit[, exclusion_reason := fifelse(
  estimated_by_music, NA_character_,
  fifelse(n_donors < 2, "Only one reference donor; no cross-donor variance estimate",
          "Not returned by MuSiC")
)]
setorder(type_audit, cellType)
fwrite(type_audit, file.path(output_dir, "music_cell_type_inclusion_audit.csv"))

diagnostics <- data.table(
  sample_id = names(music_result$r.squared.full),
  r_squared_full = as.numeric(music_result$r.squared.full)
)
fwrite(diagnostics, file.path(output_dir, "gtex_skin_music_sample_diagnostics.csv"))

sample_manifest <- fread(bulk_manifest_path)
matched_manifest <- sample_manifest[match(weighted$sample_id, sample_id)]
if (anyNA(matched_manifest$sample_id)) stop("MuSiC output samples do not exactly match bulk manifest")
fwrite(cbind(matched_manifest, weighted[, -1, drop = FALSE]),
       file.path(output_dir, "gtex_skin_music_weighted_proportions_with_site.csv"))

saveRDS(
  list(
    music_result = music_result,
    reference_donor_type_audit = reference_audit,
    sample_manifest = sample_manifest,
    common_genes = common_genes,
    parameters = list(clusters = "cellType", samples = "sampleID", select.ct = NULL,
                      counts = "raw non-log counts", seed = 20260920)
  ),
  file.path(output_dir, "gtex_skin_music_result.rds"),
  compress = "xz"
)

session <- capture.output(sessionInfo())
writeLines(session, file.path(output_dir, "sessionInfo.txt"))
message("MuSiC outputs written to: ", output_dir)
