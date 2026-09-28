############################################################
# Average SnC/SASP expression in Xenium hot and cold spots
#
# Uses Getis-Ord Gi* bin calls from:
#   results/260727/xenium_snc_hotspots/xenium_snc_hotspot_results_long.csv
#
# Outputs:
#   results/260727/xenium_snc_hotspot_expression/
#     xenium_snc_hotspot_definition.txt
#     xenium_snc_hotspot_cell_expression_long.csv
#     xenium_snc_hotspot_expression_summary.csv
#     xenium_snc_hotspot_expression_by_sample.csv
#     xenium_snc_hotspot_expression_summary.pdf/.png
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(Matrix)
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
output_dir <- "results/260727/xenium_snc_hotspot_expression"

hex_size_um <- 35
fdr_cutoff <- 0.05
metrics_to_summarize <- c("snc_prop", "mean_sasp_num")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sasp_genes <- c(
  "CCL11", "CCL13", "CCL14", "CCL15", "CCL16", "CCL19", "CCL2", "CCL20",
  "CCL21", "CCL26", "CCL28", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7",
  "CCL8", "CXCL1", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL14",
  "CXCL16", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL9", "IL12B",
  "IL13", "IL15", "IL18", "IL1A", "IL1B", "IL2", "IL23A", "IL4",
  "IL5", "IL6", "MMP12", "MMP9", "MMP3", "VEGFA", "TNF", "TGFB1",
  "CSF1", "EGF", "FGF2", "IFNB1", "CMKLR1"
)

snc_marker_genes <- c(
  "CDKN1A", "CDKN2A", "PCNA", "MKI67", "CXCL10", "CXCR3", "BCL2", "HMGB1", "LMNB1"
)

get_layer_data <- function(obj, assay, layer) {
  out <- tryCatch(
    LayerData(obj, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    return(out)
  }

  tryCatch(
    GetAssayData(obj, assay = assay, slot = layer),
    error = function(e) {
      stop("Could not get ", assay, " ", layer, " data: ", conditionMessage(e))
    }
  )
}

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

summarise_expression <- function(df, group_cols) {
  df %>%
    group_by(across(all_of(group_cols)), feature) %>%
    summarise(
      n_cells = n(),
      mean_count = mean(count, na.rm = TRUE),
      median_count = median(count, na.rm = TRUE),
      pct_expressing = mean(count > 0, na.rm = TRUE),
      .groups = "drop"
    )
}

message("Reading ", input_rds)
hsskin.ST <- readRDS(input_rds)

assay_pos <- if ("Xenium" %in% Assays(hsskin.ST)) "Xenium" else DefaultAssay(hsskin.ST)
counts <- get_layer_data(hsskin.ST, assay = assay_pos, layer = "counts")

sasp_genes_present <- intersect(sasp_genes, rownames(counts))
snc_marker_genes_present <- intersect(snc_marker_genes, rownames(counts))

message("SASP genes present: ", length(sasp_genes_present))
message("SnC marker genes present: ", paste(snc_marker_genes_present, collapse = ", "))

definition_lines <- c(
  "SnC hot/cold spot definition used in this analysis",
  "",
  paste0("Input object: ", input_rds),
  paste0("Hot-spot result file: ", hotspot_results_csv),
  paste0("Spatial binning: hexagonal bins with size = ", hex_size_um, " um."),
  "Each bin was summarized by n_cells, n_snc, snc_prop, and mean_sasp_num.",
  "",
  "SnC+ cell definition used before Gi*:",
  "  SnC+ = (CDKN1A>0 OR CDKN2A>0) AND NOT(PCNA>0 OR MKI67>0) AND sasp_num >= 5",
  "  sasp_num = number of SASP genes with raw Xenium count > 0 in that cell.",
  "",
  "Getis-Ord Gi* hot/cold spot definition:",
  "  Gi* was calculated on hex bins within each FOV using first-order hex neighbors, including the focal bin.",
  paste0("  Hot spot = Gi* z-score > 0 and Benjamini-Hochberg FDR < ", fdr_cutoff, "."),
  paste0("  Cold spot = Gi* z-score < 0 and Benjamini-Hochberg FDR < ", fdr_cutoff, "."),
  paste0("  Non-hot spot / Not significant = FDR >= ", fdr_cutoff, " among tested bins."),
  "  Untested cells are cells in hex bins below the original minimum-cell bin threshold or without a Gi* result.",
  "",
  "Expression summary:",
  "  Individual genes are raw Xenium counts per cell averaged within each hot/cold/non-significant class.",
  "  SASP_total_count = sum of raw counts across all SASP genes present in the panel for each cell.",
  "  SASP_mean_count_per_gene = SASP_total_count divided by the number of SASP genes present.",
  paste0("  SASP genes present in panel: ", length(sasp_genes_present), " / ", length(sasp_genes)),
  paste0("  Present SASP genes: ", paste(sasp_genes_present, collapse = ", "))
)
writeLines(definition_lines, file.path(output_dir, "xenium_snc_hotspot_definition.txt"))

coords_all <- bind_rows(lapply(Images(hsskin.ST), function(f) extract_fov_coords(hsskin.ST, f)))

meta <- hsskin.ST@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  select(cell_id, any_of(c("sample", "age", "snc", "sasp_num", "cell_type_minor", "immune_subclusters1")))

cell_df <- coords_all %>%
  left_join(meta, by = "cell_id") %>%
  filter(!is.na(x), !is.na(y))

cell_hex <- bind_rows(lapply(split(cell_df, cell_df$fov), function(x) assign_hexes(x, hex_size_um)))

feature_mat <- counts[unique(c(snc_marker_genes_present, sasp_genes_present)), colnames(hsskin.ST), drop = FALSE]
feature_df <- as.data.frame(t(as.matrix(feature_mat)))
feature_df$cell_id <- rownames(feature_df)

sasp_total <- Matrix::colSums(counts[sasp_genes_present, colnames(hsskin.ST), drop = FALSE])
sasp_detected <- Matrix::colSums(counts[sasp_genes_present, colnames(hsskin.ST), drop = FALSE] > 0)

expr_df <- feature_df %>%
  mutate(
    SASP_total_count = as.numeric(sasp_total[cell_id]),
    SASP_mean_count_per_gene = SASP_total_count / length(sasp_genes_present),
    SASP_num_detected = as.numeric(sasp_detected[cell_id])
  ) %>%
  pivot_longer(
    cols = -cell_id,
    names_to = "feature",
    values_to = "count"
  )

hotspot_results <- read_csv(hotspot_results_csv, show_col_types = FALSE) %>%
  filter(metric %in% metrics_to_summarize) %>%
  select(
    metric, fov, hex_id, hotspot_class, gi_star_z, gi_star_fdr,
    n_cells_bin = n_cells, n_snc_bin = n_snc, snc_prop_bin = snc_prop,
    mean_sasp_num_bin = mean_sasp_num
  )

cell_hotspot <- cell_hex %>%
  select(cell_id, fov, sample, age, hex_id, snc, sasp_num) %>%
  inner_join(hotspot_results, by = c("fov", "hex_id"), relationship = "many-to-many")

cell_expr_long <- cell_hotspot %>%
  inner_join(expr_df, by = "cell_id", relationship = "many-to-many") %>%
  mutate(
    hotspot_class = factor(hotspot_class, levels = c("Cold spot", "Not significant", "Hot spot")),
    metric = factor(metric, levels = metrics_to_summarize)
  )

write_csv(cell_expr_long, file.path(output_dir, "xenium_snc_hotspot_cell_expression_long.csv.gz"))

summary_global <- summarise_expression(
  cell_expr_long,
  c("metric", "hotspot_class")
) %>%
  arrange(metric, feature, hotspot_class)

write_csv(summary_global, file.path(output_dir, "xenium_snc_hotspot_expression_summary.csv"))

summary_by_sample <- summarise_expression(
  cell_expr_long,
  c("metric", "sample", "age", "hotspot_class")
) %>%
  arrange(metric, sample, feature, hotspot_class)

write_csv(summary_by_sample, file.path(output_dir, "xenium_snc_hotspot_expression_by_sample.csv"))

plot_features <- c(
  "SASP_total_count",
  "SASP_num_detected",
  intersect(c("CDKN1A", "CDKN2A", "PCNA", "MKI67", "CXCL10", "CXCR3"), snc_marker_genes_present)
)

plot_df <- summary_global %>%
  filter(feature %in% plot_features) %>%
  mutate(
    feature = factor(feature, levels = plot_features),
    hotspot_class = factor(hotspot_class, levels = c("Cold spot", "Not significant", "Hot spot"))
  )

p <- ggplot(plot_df, aes(hotspot_class, mean_count, fill = hotspot_class)) +
  geom_col(width = 0.75) +
  facet_grid(feature ~ metric, scales = "free_y") +
  scale_fill_manual(
    values = c("Cold spot" = "#2C7BB6", "Not significant" = "grey75", "Hot spot" = "#D7191C"),
    drop = FALSE
  ) +
  labs(
    title = "Average raw Xenium SnC/SASP signal by Gi* class",
    x = NULL,
    y = "Mean raw count per cell",
    fill = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text.y = element_text(angle = 0),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "xenium_snc_hotspot_expression_summary.pdf"), p, width = 9, height = 11)
ggsave(file.path(output_dir, "xenium_snc_hotspot_expression_summary.png"), p, width = 9, height = 11, dpi = 300)

message("Global expression summary:")
print(
  summary_global %>%
    filter(feature %in% plot_features) %>%
    select(metric, hotspot_class, feature, n_cells, mean_count, median_count, pct_expressing) %>%
    arrange(metric, feature, hotspot_class)
)

message("Done. Results written to ", output_dir)
