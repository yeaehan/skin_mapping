############################################################
# Xenium SnC/SASP hot spot and cold spot analysis
#
# Input:
#   data/skin.rds  Seurat Xenium object
#
# Outputs:
#   results/260727/xenium_snc_hotspots/
#     xenium_snc_hotspot_bins.csv
#     xenium_snc_hotspot_results_long.csv
#     xenium_snc_hotspot_summary.csv
#     xenium_snc_hotspots_<metric>.pdf/.png
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(Matrix)
  library(patchwork)
})

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_file <- normalizePath(sub("^--file=", "", script_file_arg[1]), winslash = "/", mustWork = TRUE)
  project_root <- normalizePath(file.path(dirname(script_file), "../../.."), winslash = "/", mustWork = TRUE)
  setwd(project_root)
}

# -------------------------
# User-facing parameters
# -------------------------
input_rds <- "data/skin.rds"
output_dir <- "results/260727/xenium_snc_hotspots"

hex_size_um <- 35
hex_neighbor_order <- 1
min_cells_per_hex <- 5
n_permutations <- 199
fdr_cutoff <- 0.05
set.seed(260726)

metrics_to_test <- c("snc_prop", "mean_sasp_num")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Helpers
# -------------------------
get_layer_data <- function(obj, assay, layer) {
  out <- tryCatch(
    LayerData(obj, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    return(out)
  }

  slot_name <- layer
  tryCatch(
    GetAssayData(obj, assay = assay, slot = slot_name),
    error = function(e) {
      stop("Could not get ", assay, " ", layer, " data: ", conditionMessage(e))
    }
  )
}

get_pos <- function(obj, feature, assay = "Xenium", layer = "counts") {
  if (!assay %in% Assays(obj)) {
    stop("Assay not found: ", assay)
  }
  if (!feature %in% rownames(obj[[assay]])) {
    return(setNames(rep(FALSE, ncol(obj)), colnames(obj)))
  }
  m <- get_layer_data(obj, assay = assay, layer = layer)[feature, ]
  setNames(as.numeric(m) > 0, names(m))
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
  x_col <- if ("x" %in% numeric_cols) "x" else numeric_cols[1]
  y_col <- if ("y" %in% numeric_cols) "y" else numeric_cols[2]

  if (length(numeric_cols) < 2) {
    stop("Could not find at least two numeric coordinate columns for ", fov_name)
  }

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

  center_x <- size_um * sqrt(3) * (axial$q + axial$r / 2) + x0
  center_y <- size_um * 3 / 2 * axial$r + y0

  df$hex_q <- axial$q
  df$hex_r <- axial$r
  df$hex_x <- center_x
  df$hex_y <- center_y
  df$hex_id <- paste(df$fov, axial$q, axial$r, sep = "__")
  df
}

axial_offsets <- function(order = 1) {
  offsets <- expand.grid(dq = -order:order, dr = -order:order)
  offsets$ds <- -offsets$dq - offsets$dr
  offsets %>%
    filter(pmax(abs(dq), abs(dr), abs(ds)) <= order) %>%
    select(dq, dr)
}

calc_gistar_for_group <- function(df, metric, order = 1, n_perm = 0) {
  df <- df %>% arrange(hex_q, hex_r)
  x <- df[[metric]]
  n <- length(x)

  out <- df
  out$gi_star_z <- NA_real_
  out$gi_star_p_analytic <- NA_real_
  out$gi_star_p_perm <- NA_real_
  out$n_neighbors <- NA_integer_

  if (n < 3 || stats::sd(x, na.rm = TRUE) == 0 || all(is.na(x))) {
    return(out)
  }

  offsets <- axial_offsets(order)
  key_to_index <- seq_len(n)
  names(key_to_index) <- paste(df$hex_q, df$hex_r, sep = "__")
  neighbor_index <- vector("list", n)

  for (i in seq_len(n)) {
    keys <- paste(df$hex_q[i] + offsets$dq, df$hex_r[i] + offsets$dr, sep = "__")
    neighbor_index[[i]] <- unname(key_to_index[keys[!is.na(key_to_index[keys])]])
  }

  xbar <- mean(x, na.rm = TRUE)
  s <- sqrt(mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))
  z <- rep(NA_real_, n)

  for (i in seq_len(n)) {
    idx <- neighbor_index[[i]]
    wij <- length(idx)
    denom <- s * sqrt((n * wij - wij^2) / (n - 1))
    if (wij > 0 && denom > 0) {
      z[i] <- (sum(x[idx], na.rm = TRUE) - xbar * wij) / denom
    }
  }

  p_analytic <- 2 * stats::pnorm(-abs(z))
  p_perm <- rep(NA_real_, n)

  if (n_perm > 0) {
    perm_more_extreme <- rep(0L, n)
    
    # Pre-calculate constants once OUTSIDE the permutation loop
    xbar <- mean(x, na.rm = TRUE)
    s_pop <- sqrt(mean((x - xbar)^2, na.rm = TRUE))

    for (b in seq_len(n_perm)) {
      xp <- sample(x, size = n, replace = FALSE)
      zp <- rep(NA_real_, n)

      for (i in seq_len(n)) {
        idx <- neighbor_index[[i]]
        wij <- length(idx)
        denom <- s_pop * sqrt((n * wij - wij^2) / (n - 1))
        if (wij > 0 && !is.na(denom) && denom > 0) {
          zp[i] <- (sum(xp[idx], na.rm = TRUE) - xbar * wij) / denom
        }
      }
      
      # Safely evaluate comparisons, ignoring NAs
      valid_cmp <- !is.na(zp) & !is.na(z)
      perm_more_extreme[valid_cmp] <- perm_more_extreme[valid_cmp] + 
        as.integer(abs(zp[valid_cmp]) >= abs(z[valid_cmp]))
    }
    
    p_perm <- (perm_more_extreme + 1) / (n_perm + 1)
  }

  out$gi_star_z <- z
  out$gi_star_p_analytic <- p_analytic
  out$gi_star_p_perm <- p_perm
  out$n_neighbors <- vapply(neighbor_index, length, integer(1))
  out
}

plot_hotspots <- function(results_long, metric, out_dir, fdr = 0.05) {
  df <- results_long %>%
    filter(.data$metric == metric) %>%
    mutate(
      hotspot_class = case_when(
        gi_star_fdr < fdr & gi_star_z > 0 ~ "Hot spot",
        gi_star_fdr < fdr & gi_star_z < 0 ~ "Cold spot",
        TRUE ~ "Not significant"
      ),
      hotspot_class = factor(
        hotspot_class,
        levels = c("Cold spot", "Not significant", "Hot spot")
      )
    )

  p_z <- ggplot(df, aes(hex_x, hex_y, color = gi_star_z)) +
    geom_point(size = 0.75, alpha = 0.9) +
    facet_wrap(~ sample) +
    coord_fixed() +
    scale_y_reverse() +
    scale_color_gradient2(
      low = "#2C7BB6",
      mid = "grey92",
      high = "#D7191C",
      midpoint = 0,
      na.value = "grey75"
    ) +
    labs(
      title = paste0("Getis-Ord Gi* z-score: ", metric),
      x = NULL,
      y = NULL,
      color = "Gi* z"
    ) +
    theme_void(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )

  p_class <- ggplot(df, aes(hex_x, hex_y, color = hotspot_class)) +
    geom_point(size = 0.75, alpha = 0.9) +
    facet_wrap(~ sample) +
    coord_fixed() +
    scale_y_reverse() +
    scale_color_manual(
      values = c(
        "Cold spot" = "#2C7BB6",
        "Not significant" = "grey84",
        "Hot spot" = "#D7191C"
      ),
      drop = FALSE
    ) +
    labs(
      title = paste0("Significant hot/cold spots: ", metric),
      x = NULL,
      y = NULL,
      color = NULL
    ) +
    theme_void(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )

  p <- p_z / p_class
  safe_metric <- gsub("[^A-Za-z0-9_]+", "_", metric)
  ggsave(
    file.path(out_dir, paste0("xenium_snc_hotspots_", safe_metric, ".pdf")),
    p,
    width = 12,
    height = 10
  )
  ggsave(
    file.path(out_dir, paste0("xenium_snc_hotspots_", safe_metric, ".png")),
    p,
    width = 12,
    height = 10,
    dpi = 300
  )
}

# -------------------------
# Load object and annotate SnC
# -------------------------
message("Reading ", input_rds)
hsskin.ST <- readRDS(input_rds)

assay_pos <- if ("Xenium" %in% Assays(hsskin.ST)) "Xenium" else DefaultAssay(hsskin.ST)
layer_pos <- "counts"

sasp_genes <- c(
  "CCL11", "CCL13", "CCL14", "CCL15", "CCL16", "CCL19", "CCL2", "CCL20",
  "CCL21", "CCL26", "CCL28", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7",
  "CCL8", "CXCL1", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL14",
  "CXCL16", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL9", "IL12B",
  "IL13", "IL15", "IL18", "IL1A", "IL1B", "IL2", "IL23A", "IL4",
  "IL5", "IL6", "MMP12", "MMP9", "MMP3", "VEGFA", "TNF", "TGFB1",
  "CSF1", "EGF", "FGF2", "IFNB1", "CMKLR1"
)

sasp_genes_present <- intersect(sasp_genes, rownames(hsskin.ST[[assay_pos]]))
message("SASP genes detected in panel: ", length(sasp_genes_present))

counts <- get_layer_data(hsskin.ST, assay = assay_pos, layer = layer_pos)
sasp_matrix <- counts[sasp_genes_present, , drop = FALSE]
hsskin.ST$sasp_num <- Matrix::colSums(sasp_matrix > 0)

p21 <- get_pos(hsskin.ST, "CDKN1A", assay_pos, layer_pos)
p16 <- get_pos(hsskin.ST, "CDKN2A", assay_pos, layer_pos)
PCNA <- get_pos(hsskin.ST, "PCNA", assay_pos, layer_pos)
MKI67 <- get_pos(hsskin.ST, "MKI67", assay_pos, layer_pos)

hsskin.ST$snc <- ifelse(
  ((p21 | p16) & !(PCNA | MKI67) & hsskin.ST$sasp_num >= 5),
  "SnC+",
  "SnC-"
)
hsskin.ST$snc <- factor(hsskin.ST$snc, levels = c("SnC-", "SnC+"))

message("SnC table:")
print(table(hsskin.ST$snc, useNA = "ifany"))

# -------------------------
# Extract coordinates
# -------------------------
fov_names <- Images(hsskin.ST)
message("Extracting coordinates from ", length(fov_names), " FOVs")

coords_all <- bind_rows(lapply(fov_names, function(f) extract_fov_coords(hsskin.ST, f)))

meta <- hsskin.ST@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  select(
    cell_id,
    any_of(c("sample", "age", "cell_type_major", "cell_type_minor", "cell_type_final")),
    sasp_num,
    snc
  )

cell_df <- coords_all %>%
  left_join(meta, by = "cell_id") %>%
  filter(!is.na(x), !is.na(y), !is.na(snc), !is.na(sasp_num))

if (!"sample" %in% colnames(cell_df)) {
  cell_df$sample <- cell_df$fov
}

if (!"age" %in% colnames(cell_df)) {
  cell_df$age <- NA_character_
}

message("Cells with coordinates and metadata: ", nrow(cell_df))

cell_hex <- bind_rows(lapply(split(cell_df, cell_df$fov), function(x) {
  assign_hexes(x, hex_size_um)
}))

hex_bins <- cell_hex %>%
  group_by(fov, sample, age, hex_id, hex_q, hex_r, hex_x, hex_y) %>%
  summarise(
    n_cells = n(),
    n_snc = sum(snc == "SnC+", na.rm = TRUE),
    snc_prop = n_snc / n_cells,
    mean_sasp_num = mean(sasp_num, na.rm = TRUE),
    sum_sasp_num = sum(sasp_num, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_hex)

write_csv(hex_bins, file.path(output_dir, "xenium_snc_hotspot_bins.csv"))
message("Hex bins retained: ", nrow(hex_bins))

# -------------------------
# Run Getis-Ord Gi*
# -------------------------
results_long <- bind_rows(lapply(metrics_to_test, function(metric) {
  message("Running Gi* for ", metric)
  hex_bins %>%
    group_by(fov) %>%
    group_modify(~ calc_gistar_for_group(
      .x,
      metric = metric,
      order = hex_neighbor_order,
      n_perm = n_permutations
    )) %>%
    ungroup() %>%
    mutate(metric = metric)
}))

# Grouping by fov and metric ensures FDR is calculated per sample
results_long <- results_long %>%
  group_by(fov, metric) %>%
  mutate(
    gi_star_p = gi_star_p_analytic,
    gi_star_fdr = p.adjust(gi_star_p, method = "BH"),
    hotspot_class = case_when(
      gi_star_fdr < fdr_cutoff & gi_star_z > 0 ~ "Hot spot",
      gi_star_fdr < fdr_cutoff & gi_star_z < 0 ~ "Cold spot",
      TRUE ~ "Not significant"
    )
  ) %>%
  ungroup()

write_csv(results_long, file.path(output_dir, "xenium_snc_hotspot_results_long.csv"))

summary_tbl <- results_long %>%
  group_by(metric, fov, sample, age, hotspot_class) %>%
  summarise(
    n_hex = n(),
    n_cells = sum(n_cells, na.rm = TRUE),
    n_snc = sum(n_snc, na.rm = TRUE),
    median_gi_star_z = median(gi_star_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(metric, fov, hotspot_class)

write_csv(summary_tbl, file.path(output_dir, "xenium_snc_hotspot_summary.csv"))

for (metric in metrics_to_test) {
  plot_hotspots(results_long, metric, output_dir, fdr = fdr_cutoff)
}

message("Done. Results written to ", output_dir)
