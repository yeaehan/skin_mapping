############################################################
# Xenium spatial scale diagnostics
#
# Quantifies why the 35 um hex size is a reasonable starting
# scale for Gi* analysis in this dataset.
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_file <- normalizePath(sub("^--file=", "", script_file_arg[1]), winslash = "/", mustWork = TRUE)
  project_root <- normalizePath(file.path(dirname(script_file), "../../.."), winslash = "/", mustWork = TRUE)
  setwd(project_root)
}

input_rds <- "data/skin.rds"
hex_bins_csv <- "results/260727/xenium_snc_hotspots/xenium_snc_hotspot_bins.csv"
output_dir <- "results/260727/xenium_spatial_scale_diagnostics"

hex_size_um <- 35
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

extract_fov_coords <- function(obj, fov_name) {
  coords <- GetTissueCoordinates(obj[[fov_name]])
  coords <- as.data.frame(coords)
  coords$coord_rowname <- rownames(coords)
  coords$fov <- fov_name

  if ("cell" %in% colnames(coords)) {
    coords$cell_id <- as.character(coords$cell)
  } else if ("cell_id" %in% colnames(coords)) {
    coords$cell_id <- as.character(coords$cell_id)
  } else {
    coords$cell_id <- as.character(coords$coord_rowname)
  }

  numeric_cols <- names(coords)[vapply(coords, is.numeric, logical(1))]
  if (length(numeric_cols) < 2) {
    stop("Could not find at least two numeric coordinate columns for ", fov_name)
  }

  x_col <- if ("x" %in% numeric_cols) "x" else numeric_cols[1]
  y_col <- if ("y" %in% numeric_cols) "y" else numeric_cols[2]

  coords %>%
    transmute(
      cell_id = cell_id,
      fov = fov,
      x = .data[[x_col]],
      y = .data[[y_col]]
    )
}

extract_fov_areas <- function(obj, fov_name) {
  segmentation <- slot(obj[[fov_name]], "boundaries")$segmentations
  polygons <- slot(segmentation, "polygons")

  data.frame(
    cell_id = vapply(polygons, function(p) slot(p, "ID"), character(1)),
    fov = fov_name,
    cell_area_um2 = vapply(polygons, function(p) slot(p, "area"), numeric(1))
  )
}

nearest_neighbor_distance <- function(df) {
  xy <- as.matrix(df[, c("x", "y")])
  n <- nrow(xy)
  if (n < 2) {
    return(rep(NA_real_, n))
  }

  vapply(seq_len(n), function(i) {
    dx <- xy[, 1] - xy[i, 1]
    dy <- xy[, 2] - xy[i, 2]
    d <- sqrt(dx^2 + dy^2)
    min(d[d > 0], na.rm = TRUE)
  }, numeric(1))
}

summary_stats <- function(x) {
  tibble(
    n = sum(!is.na(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    p05 = unname(quantile(x, 0.05, na.rm = TRUE)),
    p10 = unname(quantile(x, 0.10, na.rm = TRUE)),
    p25 = unname(quantile(x, 0.25, na.rm = TRUE)),
    median = median(x, na.rm = TRUE),
    p75 = unname(quantile(x, 0.75, na.rm = TRUE)),
    p90 = unname(quantile(x, 0.90, na.rm = TRUE)),
    p95 = unname(quantile(x, 0.95, na.rm = TRUE)),
    max = max(x, na.rm = TRUE)
  )
}

message("Reading ", input_rds)
obj <- readRDS(input_rds)

coords <- bind_rows(lapply(Images(obj), function(f) extract_fov_coords(obj, f)))
areas <- bind_rows(lapply(Images(obj), function(f) extract_fov_areas(obj, f)))
meta <- obj@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  select(cell_id, any_of(c("sample", "age", "cell_type_major", "cell_type_minor", "cell_type_final")))

cell_size <- coords %>%
  left_join(areas, by = c("cell_id", "fov")) %>%
  left_join(meta, by = "cell_id") %>%
  mutate(equivalent_diameter_um = 2 * sqrt(cell_area_um2 / pi))

cell_size$nearest_neighbor_um <- unlist(lapply(split(cell_size, cell_size$fov), nearest_neighbor_distance))

write_csv(cell_size, file.path(output_dir, "xenium_cell_size_and_spacing.csv"))

cell_area_summary <- bind_rows(
  summary_stats(cell_size$cell_area_um2) %>% mutate(group = "all_cells", metric = "cell_area_um2"),
  summary_stats(cell_size$equivalent_diameter_um) %>% mutate(group = "all_cells", metric = "equivalent_diameter_um"),
  summary_stats(cell_size$nearest_neighbor_um) %>% mutate(group = "all_cells", metric = "nearest_neighbor_um")
) %>%
  select(group, metric, everything())

cell_area_by_type <- cell_size %>%
  filter(!is.na(cell_type_major)) %>%
  group_by(cell_type_major) %>%
  summarise(
    summary_stats(cell_area_um2),
    median_equivalent_diameter_um = median(equivalent_diameter_um, na.rm = TRUE),
    median_nearest_neighbor_um = median(nearest_neighbor_um, na.rm = TRUE),
    .groups = "drop"
  )

hex_bins <- read_csv(hex_bins_csv, show_col_types = FALSE)
hex_area_um2 <- 3 * sqrt(3) / 2 * hex_size_um^2

hex_summary <- bind_rows(
  summary_stats(hex_bins$n_cells) %>% mutate(group = "retained_35um_hex_bins", metric = "cells_per_hex")
) %>%
  mutate(
    hex_size_um = hex_size_um,
    hex_area_um2 = hex_area_um2,
    median_cells_per_1000_um2 = median / hex_area_um2 * 1000,
    implied_area_per_cell_from_median_hex_um2 = hex_area_um2 / median,
    implied_diameter_from_median_hex_um = 2 * sqrt((hex_area_um2 / median) / pi)
  ) %>%
  select(group, metric, hex_size_um, hex_area_um2, everything())

write_csv(cell_area_summary, file.path(output_dir, "xenium_cell_area_spacing_summary.csv"))
write_csv(cell_area_by_type, file.path(output_dir, "xenium_cell_area_spacing_by_major_cell_type.csv"))
write_csv(hex_summary, file.path(output_dir, "xenium_35um_hex_cell_count_summary.csv"))

writeLines(
  c(
    paste0("35 um hex area = ", round(hex_area_um2, 1), " um2."),
    paste0("Retained 35 um hex bins = ", nrow(hex_bins), "."),
    paste0("Cells in retained bins = ", sum(hex_bins$n_cells), "."),
    "",
    "Cells per 35 um hex bin:",
    capture.output(print(hex_summary)),
    "",
    "Segmentation-derived cell size and spacing:",
    capture.output(print(cell_area_summary)),
    "",
    "By major cell type:",
    capture.output(print(cell_area_by_type))
  ),
  file.path(output_dir, "xenium_spatial_scale_diagnostics.txt")
)

message("Done. Results written to ", output_dir)
print(hex_summary)
print(cell_area_summary)
print(cell_area_by_type)
