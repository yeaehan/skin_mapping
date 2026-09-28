############################################################
# Cell types near Xenium SnC hot spots
#
# Uses SnC Getis-Ord Gi* hot spots from:
#   results/260727/xenium_snc_hotspots/xenium_snc_hotspot_results_long.csv
#
# Outputs:
#   results/260727/xenium_snc_hotspot_neighborhoods/
#     xenium_snc_hotspot_cell_neighborhoods.csv
#     xenium_snc_hotspot_celltype_composition.csv
#     xenium_snc_hotspot_celltype_composition_by_sample.csv
#     xenium_snc_hotspot_celltype_enrichment.csv
#     xenium_snc_hotspot_celltype_enrichment_by_sample.csv
#     xenium_snc_hotspot_celltype_<annotation>.pdf/.png
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_file <- normalizePath(sub("^--file=", "", script_file_arg[1]), winslash = "/", mustWork = TRUE)
  project_root <- normalizePath(file.path(dirname(script_file), "../../.."), winslash = "/", mustWork = TRUE)
  setwd(project_root)
}

input_rds <- "data/skin.rds"
hotspot_results_csv <- "results/260727/xenium_snc_hotspots/xenium_snc_hotspot_results_long.csv"
output_dir <- "results/260727/xenium_snc_hotspot_neighborhoods"

hex_size_um <- 35
hotspot_metric <- "snc_prop"
hotspot_class <- "Hot spot"
annotation_cols <- c("cell_type_minor", "cell_type_final", "immune_subclusters1")
zone_levels <- c("hotspot_core", "neighbor_1", "neighbor_2", "background")

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

axial_round <- function(q, r) {
  x <- q
  z <- r
  y <- -x - z

  rx <- round(x)
  ry <- round(y)
  rz <- round(z)

  x_diff <- abs(rx - x)
  y_diff <- abs(ry - y)
  z_diff <- abs(rz - z)

  use_x <- x_diff > y_diff & x_diff > z_diff
  use_y <- !use_x & y_diff > z_diff
  use_z <- !use_x & !use_y

  rx[use_x] <- -ry[use_x] - rz[use_x]
  ry[use_y] <- -rx[use_y] - rz[use_y]
  rz[use_z] <- -rx[use_z] - ry[use_z]

  data.frame(q = as.integer(rx), r = as.integer(rz))
}

assign_hexes <- function(df, size_um) {
  x0 <- min(df$x, na.rm = TRUE)
  y0 <- min(df$y, na.rm = TRUE)
  x_shift <- df$x - x0
  y_shift <- df$y - y0

  q_float <- (sqrt(3) / 3 * x_shift - 1 / 3 * y_shift) / size_um
  r_float <- (2 / 3 * y_shift) / size_um
  axial <- axial_round(q_float, r_float)

  df$hex_q <- axial$q
  df$hex_r <- axial$r
  df$hex_id <- paste(df$fov, axial$q, axial$r, sep = "__")
  df
}

hex_distance <- function(q1, r1, q2, r2) {
  pmax(abs(q1 - q2), abs(r1 - r2), abs((-q1 - r1) - (-q2 - r2)))
}

nearest_hotspot_distance <- function(df, hot_hexes) {
  fov_hot <- hot_hexes %>% filter(fov == df$fov[1])
  if (nrow(fov_hot) == 0) {
    return(rep(NA_integer_, nrow(df)))
  }

  vapply(seq_len(nrow(df)), function(i) {
    as.integer(min(hex_distance(df$hex_q[i], df$hex_r[i], fov_hot$hex_q, fov_hot$hex_r)))
  }, integer(1))
}

zone_from_distance <- function(d) {
  case_when(
    is.na(d) ~ "background",
    d == 0 ~ "hotspot_core",
    d == 1 ~ "neighbor_1",
    d == 2 ~ "neighbor_2",
    TRUE ~ "background"
  )
}

sem <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) {
    return(NA_real_)
  }
  stats::sd(x) / sqrt(length(x))
}

safe_one_sample_t <- function(x, mu = 0) {
  x <- x[!is.na(x)]
  if (length(x) < 2 || stats::sd(x) == 0) {
    return(NA_real_)
  }
  stats::t.test(x, mu = mu)$p.value
}

make_sample_composition <- function(cell_df, annotation_col) {
  ann <- annotation_col
  base <- cell_df %>%
    mutate(
      annotation_col = ann,
      cell_type = as.character(.data[[ann]]),
      cell_type = if_else(is.na(cell_type) | cell_type == "", "unknown", cell_type)
    )

  sample_zones <- base %>%
    distinct(annotation_col, sample, age, zone)

  cell_types <- base %>%
    distinct(annotation_col, cell_type)

  counts <- base %>%
    group_by(annotation_col, sample, age, zone, cell_type) %>%
    summarise(n_cells = n(), .groups = "drop")

  sample_zones %>%
    left_join(cell_types, by = "annotation_col", relationship = "many-to-many") %>%
    left_join(
      counts,
      by = c("annotation_col", "sample", "age", "zone", "cell_type")
    ) %>%
    mutate(n_cells = replace_na(n_cells, 0L)) %>%
    group_by(annotation_col, sample, age, zone) %>%
    mutate(
      sample_zone_total = sum(n_cells),
      prop = if_else(sample_zone_total > 0, n_cells / sample_zone_total, NA_real_)
    ) %>%
    ungroup()
}

summarise_composition <- function(sample_composition) {
  sample_composition %>%
    group_by(annotation_col, zone, cell_type) %>%
    summarise(
      n_samples = n_distinct(sample),
      total_cells = sum(n_cells),
      mean_prop = mean(prop, na.rm = TRUE),
      sem_prop = sem(prop),
      mean_n_cells = mean(n_cells, na.rm = TRUE),
      .groups = "drop"
    )
}

make_sample_enrichment <- function(cell_df, annotation_col, near_zones) {
  ann <- annotation_col
  near_label <- paste(near_zones, collapse = "+")
  base <- cell_df %>%
    mutate(
      near_definition = near_label,
      zone_set = if_else(zone %in% near_zones, "near", "reference"),
      annotation_col = ann,
      cell_type = as.character(.data[[ann]]),
      cell_type = if_else(is.na(cell_type) | cell_type == "", "unknown", cell_type)
    )

  sample_zone_sets <- base %>%
    distinct(annotation_col, near_definition, sample, age, zone_set)

  cell_types <- base %>%
    distinct(annotation_col, near_definition, cell_type)

  counts <- base %>%
    group_by(annotation_col, near_definition, sample, age, zone_set, cell_type) %>%
    summarise(n_cells = n(), .groups = "drop")

  props <- sample_zone_sets %>%
    left_join(cell_types, by = c("annotation_col", "near_definition"), relationship = "many-to-many") %>%
    left_join(
      counts,
      by = c("annotation_col", "near_definition", "sample", "age", "zone_set", "cell_type")
    ) %>%
    mutate(n_cells = replace_na(n_cells, 0L)) %>%
    group_by(annotation_col, near_definition, sample, age, zone_set) %>%
    mutate(
      zone_set_total = sum(n_cells),
      prop = if_else(zone_set_total > 0, n_cells / zone_set_total, NA_real_)
    ) %>%
    ungroup()

  props %>%
    select(annotation_col, near_definition, sample, age, zone_set, cell_type, n_cells, zone_set_total, prop) %>%
    pivot_wider(
      names_from = zone_set,
      values_from = c(n_cells, zone_set_total, prop),
      values_fill = list(n_cells = 0L, zone_set_total = 0L, prop = 0)
    ) %>%
    mutate(
      near_prop = prop_near,
      reference_prop = prop_reference,
      prop_diff = near_prop - reference_prop,
      near_prop_pseudo = near_prop + 1e-6,
      reference_prop_pseudo = reference_prop + 1e-6,
      log2_ratio = log2(near_prop_pseudo / reference_prop_pseudo)
    )
}

summarise_enrichment <- function(sample_enrichment) {
  sample_enrichment %>%
    group_by(annotation_col, near_definition, cell_type) %>%
    summarise(
      n_samples = n_distinct(sample),
      mean_near_prop = mean(near_prop, na.rm = TRUE),
      sem_near_prop = sem(near_prop),
      mean_reference_prop = mean(reference_prop, na.rm = TRUE),
      sem_reference_prop = sem(reference_prop),
      mean_prop_diff = mean(prop_diff, na.rm = TRUE),
      sem_prop_diff = sem(prop_diff),
      mean_log2_ratio = mean(log2_ratio, na.rm = TRUE),
      sem_log2_ratio = sem(log2_ratio),
      p_value = safe_one_sample_t(prop_diff, mu = 0),
      .groups = "drop"
    ) %>%
    group_by(annotation_col, near_definition) %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    arrange(annotation_col, near_definition, fdr, desc(mean_log2_ratio))
}

plot_annotation <- function(composition, enrichment, annotation_col, out_dir) {
  ann <- annotation_col
  comp <- composition %>%
    filter(.data$annotation_col == ann) %>%
    group_by(cell_type) %>%
    mutate(total = sum(total_cells)) %>%
    ungroup() %>%
    mutate(cell_type = fct_reorder(cell_type, total))

  keep_types <- comp %>%
    group_by(cell_type) %>%
    summarise(total = sum(total_cells), .groups = "drop") %>%
    slice_max(total, n = 16) %>%
    pull(cell_type)

  comp_plot <- comp %>%
    mutate(cell_type_plot = if_else(cell_type %in% keep_types, as.character(cell_type), "Other")) %>%
    group_by(zone, cell_type_plot) %>%
    summarise(mean_prop = sum(mean_prop), .groups = "drop") %>%
    ungroup() %>%
    mutate(zone = factor(zone, levels = zone_levels))

  enrich_plot <- enrichment %>%
    filter(
      .data$annotation_col == ann,
      near_definition == "hotspot_core+neighbor_1",
      is.finite(mean_log2_ratio)
    ) %>%
    arrange(desc(abs(mean_log2_ratio))) %>%
    slice_head(n = 18) %>%
    mutate(
      cell_type = fct_reorder(cell_type, mean_log2_ratio),
      significant = fdr < 0.05
    )

  p1 <- ggplot(comp_plot, aes(zone, mean_prop, fill = cell_type_plot)) +
    geom_col(width = 0.8, color = "white", linewidth = 0.1) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      title = paste0(ann, " composition by SnC hot-spot distance, sample mean"),
      x = NULL,
      y = "Mean sample cell fraction",
      fill = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  p2 <- ggplot(enrich_plot, aes(mean_log2_ratio, cell_type, fill = significant)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#D7191C"), drop = FALSE) +
    labs(
      title = "Sample-level enrichment: core + neighbor_1 vs reference",
      x = "Mean sample log2 ratio",
      y = NULL,
      fill = "FDR < 0.05"
    ) +
    theme_classic(base_size = 10)

  p <- p1 / p2
  safe_col <- gsub("[^A-Za-z0-9_]+", "_", ann)
  ggsave(file.path(out_dir, paste0("xenium_snc_hotspot_celltype_", safe_col, ".pdf")), p, width = 10, height = 11)
  ggsave(file.path(out_dir, paste0("xenium_snc_hotspot_celltype_", safe_col, ".png")), p, width = 10, height = 11, dpi = 300)
}

message("Reading ", input_rds)
hsskin.ST <- readRDS(input_rds)

message("Reading SnC hot spots from ", hotspot_results_csv)
hotspot_results <- read_csv(hotspot_results_csv, show_col_types = FALSE)
hot_hexes <- hotspot_results %>%
  filter(metric == hotspot_metric, hotspot_class == !!hotspot_class) %>%
  distinct(fov, hex_id, hex_q, hex_r)

if (nrow(hot_hexes) == 0) {
  stop("No hot-spot hexes found for metric = ", hotspot_metric)
}
message("SnC hot-spot hexes: ", nrow(hot_hexes))

missing_ann <- setdiff(annotation_cols, colnames(hsskin.ST@meta.data))
if (length(missing_ann) > 0) {
  warning("Missing annotation columns skipped: ", paste(missing_ann, collapse = ", "))
}
annotation_cols <- intersect(annotation_cols, colnames(hsskin.ST@meta.data))

coords_all <- bind_rows(lapply(Images(hsskin.ST), function(f) extract_fov_coords(hsskin.ST, f)))

meta_cols <- unique(c(
  "sample", "age", "snc", "sasp_num", "immune_cell",
  "cell_type_major", annotation_cols
))

meta <- hsskin.ST@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  select(cell_id, any_of(meta_cols))

cell_df <- coords_all %>%
  left_join(meta, by = "cell_id") %>%
  filter(!is.na(x), !is.na(y))

cell_hex <- bind_rows(lapply(split(cell_df, cell_df$fov), function(x) assign_hexes(x, hex_size_um)))
cell_hex$hex_distance_to_snc_hotspot <- unlist(lapply(split(cell_hex, cell_hex$fov), nearest_hotspot_distance, hot_hexes = hot_hexes))
cell_hex$zone <- zone_from_distance(cell_hex$hex_distance_to_snc_hotspot)
cell_hex$zone <- factor(cell_hex$zone, levels = zone_levels)

write_csv(cell_hex, file.path(output_dir, "xenium_snc_hotspot_cell_neighborhoods.csv"))

sample_composition <- bind_rows(lapply(annotation_cols, function(col) make_sample_composition(cell_hex, col)))
composition <- summarise_composition(sample_composition)
write_csv(sample_composition, file.path(output_dir, "xenium_snc_hotspot_celltype_composition_by_sample.csv"))
write_csv(composition, file.path(output_dir, "xenium_snc_hotspot_celltype_composition.csv"))

near_sets <- list(
  "core" = c("hotspot_core"),
  "core_neighbor1" = c("hotspot_core", "neighbor_1"),
  "core_neighbor1_neighbor2" = c("hotspot_core", "neighbor_1", "neighbor_2")
)

sample_enrichment <- bind_rows(lapply(annotation_cols, function(col) {
  bind_rows(lapply(near_sets, function(zones) make_sample_enrichment(cell_hex, col, zones)))
}))
enrichment <- summarise_enrichment(sample_enrichment)
write_csv(sample_enrichment, file.path(output_dir, "xenium_snc_hotspot_celltype_enrichment_by_sample.csv"))
write_csv(enrichment, file.path(output_dir, "xenium_snc_hotspot_celltype_enrichment.csv"))

for (col in annotation_cols) {
  plot_annotation(composition, enrichment, col, output_dir)
}

message("Cells by zone:")
print(table(cell_hex$zone, useNA = "ifany"))

message("Top enrichments for core + neighbor_1:")
print(
  enrichment %>%
    filter(near_definition == "hotspot_core+neighbor_1", fdr < 0.05) %>%
    arrange(desc(mean_log2_ratio)) %>%
    select(annotation_col, cell_type, mean_near_prop, mean_reference_prop, mean_prop_diff, mean_log2_ratio, fdr) %>%
    head(25)
)

message("Done. Results written to ", output_dir)
