############################################################
# Young vs old comparison of Xenium SnC hot spots
#
# Compares:
#   1) SnC hot-spot burden
#   2) Cell-type composition inside hot spots
#   3) SASP signal inside/near hot spots
#   4) Neighboring-cell composition around hot spots
############################################################

suppressPackageStartupMessages({
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

hotspot_results_csv <- "results/260727/xenium_snc_hotspots/xenium_snc_hotspot_results_long.csv"
cell_neighborhoods_csv <- "results/260727/xenium_snc_hotspot_neighborhoods/xenium_snc_hotspot_cell_neighborhoods.csv"
expression_by_sample_csv <- "results/260727/xenium_snc_hotspot_expression/xenium_snc_hotspot_expression_by_sample.csv"
output_dir <- "results/260727/xenium_snc_hotspot_young_old_comparison"

hotspot_metric <- "snc_prop"
annotation_cols <- c("cell_type_minor", "cell_type_final", "immune_subclusters1")
zone_levels <- c("hotspot_core", "neighbor_1", "neighbor_2", "background")
near_zones <- c("hotspot_core", "neighbor_1")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mean_sem <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(tibble(mean = NA_real_, sem = NA_real_, n_samples = 0L))
  }
  tibble(
    mean = mean(x),
    sem = ifelse(length(x) > 1, sd(x) / sqrt(length(x)), NA_real_),
    n_samples = length(x)
  )
}

sem <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) {
    return(NA_real_)
  }
  stats::sd(x) / sqrt(length(x))
}

make_sample_composition <- function(df, annotation_col, zones_to_keep = zone_levels) {
  ann <- annotation_col
  base <- df %>%
    filter(zone %in% zones_to_keep) %>%
    mutate(
      annotation_col = ann,
      cell_type = as.character(.data[[ann]]),
      cell_type = if_else(is.na(cell_type) | cell_type == "", "unknown", cell_type)
    )

  sample_zones <- base %>%
    distinct(annotation_col, age, sample, zone)

  cell_types <- base %>%
    distinct(annotation_col, cell_type)

  counts <- base %>%
    group_by(annotation_col, age, sample, zone, cell_type) %>%
    summarise(n_cells = n(), n_snc = sum(snc == "SnC+", na.rm = TRUE), .groups = "drop")

  sample_zones %>%
    left_join(cell_types, by = "annotation_col", relationship = "many-to-many") %>%
    left_join(
      counts,
      by = c("annotation_col", "age", "sample", "zone", "cell_type")
    ) %>%
    mutate(
      n_cells = replace_na(n_cells, 0L),
      n_snc = replace_na(n_snc, 0L)
    ) %>%
    group_by(annotation_col, age, sample, zone) %>%
    mutate(
      sample_zone_total = sum(n_cells),
      prop = if_else(sample_zone_total > 0, n_cells / sample_zone_total, NA_real_)
    ) %>%
    ungroup()
}

make_age_composition <- function(sample_composition) {
  sample_composition %>%
    group_by(annotation_col, age, zone, cell_type) %>%
    summarise(
      n_samples = n_distinct(sample),
      total_cells = sum(n_cells),
      mean_prop = mean(prop, na.rm = TRUE),
      sem_prop = sem(prop),
      mean_n_cells = mean(n_cells, na.rm = TRUE),
      .groups = "drop"
    )
}

safe_one_sample_t <- function(x, mu = 0) {
  x <- x[!is.na(x)]
  if (length(x) < 2 || stats::sd(x) == 0) {
    return(NA_real_)
  }
  stats::t.test(x, mu = mu)$p.value
}

safe_two_group_t <- function(value, group) {
  keep <- !is.na(value) & !is.na(group)
  value <- value[keep]
  group <- droplevels(factor(group[keep]))
  if (nlevels(group) != 2 || any(table(group) < 2)) {
    return(NA_real_)
  }
  stats::t.test(value ~ group)$p.value
}

composition_age_tests <- function(sample_composition) {
  sample_composition %>%
    group_by(annotation_col, zone, cell_type) %>%
    summarise(
      n_young = n_distinct(sample[age == "Young"]),
      n_old = n_distinct(sample[age == "Old"]),
      mean_prop_young = mean(prop[age == "Young"], na.rm = TRUE),
      mean_prop_old = mean(prop[age == "Old"], na.rm = TRUE),
      old_minus_young = mean_prop_old - mean_prop_young,
      p_value = safe_two_group_t(prop, age),
      .groups = "drop"
    ) %>%
    group_by(annotation_col, zone) %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    arrange(annotation_col, zone, fdr, desc(abs(old_minus_young)))
}

make_sample_enrichment <- function(df, annotation_col) {
  ann <- annotation_col
  base <- df %>%
    mutate(
      zone_set = if_else(zone %in% near_zones, "near", "reference"),
      annotation_col = ann,
      cell_type = as.character(.data[[ann]]),
      cell_type = if_else(is.na(cell_type) | cell_type == "", "unknown", cell_type)
    )

  sample_zone_sets <- base %>%
    distinct(annotation_col, age, sample, zone_set)

  cell_types <- base %>%
    distinct(annotation_col, cell_type)

  counts <- base %>%
    group_by(annotation_col, age, sample, zone_set, cell_type) %>%
    summarise(n_cells = n(), .groups = "drop")

  props <- sample_zone_sets %>%
    left_join(cell_types, by = "annotation_col", relationship = "many-to-many") %>%
    left_join(
      counts,
      by = c("annotation_col", "age", "sample", "zone_set", "cell_type")
    ) %>%
    mutate(n_cells = replace_na(n_cells, 0L)) %>%
    group_by(annotation_col, age, sample, zone_set) %>%
    mutate(
      zone_set_total = sum(n_cells),
      prop = if_else(zone_set_total > 0, n_cells / zone_set_total, NA_real_)
    ) %>%
    ungroup()

  props %>%
    select(annotation_col, age, sample, zone_set, cell_type, n_cells, zone_set_total, prop) %>%
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

summarise_sample_enrichment <- function(sample_enrichment) {
  sample_enrichment %>%
    group_by(annotation_col, age, cell_type) %>%
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
    group_by(annotation_col, age) %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    arrange(annotation_col, age, fdr, desc(mean_log2_ratio))
}

sample_enrichment_age_tests <- function(sample_enrichment) {
  sample_enrichment %>%
    group_by(annotation_col, cell_type) %>%
    summarise(
      n_young = n_distinct(sample[age == "Young"]),
      n_old = n_distinct(sample[age == "Old"]),
      mean_log2_ratio_young = mean(log2_ratio[age == "Young"], na.rm = TRUE),
      mean_log2_ratio_old = mean(log2_ratio[age == "Old"], na.rm = TRUE),
      old_minus_young_log2_ratio = mean_log2_ratio_old - mean_log2_ratio_young,
      p_value = safe_two_group_t(log2_ratio, age),
      .groups = "drop"
    ) %>%
    group_by(annotation_col) %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    arrange(annotation_col, fdr, desc(abs(old_minus_young_log2_ratio)))
}

plot_composition <- function(composition, annotation_col, zones_to_plot, file_stem) {
  ann <- annotation_col
  comp <- composition %>%
    filter(annotation_col == ann, zone %in% zones_to_plot)

  keep_types <- comp %>%
    group_by(cell_type) %>%
    summarise(total = sum(total_cells), .groups = "drop") %>%
    slice_max(total, n = 12) %>%
    pull(cell_type)

  plot_df <- comp %>%
    mutate(cell_type_plot = if_else(cell_type %in% keep_types, cell_type, "Other")) %>%
    group_by(age, zone, cell_type_plot) %>%
    summarise(mean_prop = sum(mean_prop), .groups = "drop") %>%
    ungroup() %>%
    mutate(zone = factor(zone, levels = zone_levels))

  p <- ggplot(plot_df, aes(age, mean_prop, fill = cell_type_plot)) +
    geom_col(width = 0.78, color = "white", linewidth = 0.1) +
    facet_wrap(~ zone, nrow = 1) +
    scale_y_continuous(labels = percent_format()) +
    labs(
      title = paste0(ann, " composition near SnC hot spots, sample mean"),
      x = NULL,
      y = "Mean sample cell fraction",
      fill = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(legend.position = "right")

  ggsave(file.path(output_dir, paste0(file_stem, ".pdf")), p, width = 11, height = 6)
  ggsave(file.path(output_dir, paste0(file_stem, ".png")), p, width = 11, height = 6, dpi = 300)
}

message("Reading inputs")
hotspot_results <- read_csv(hotspot_results_csv, show_col_types = FALSE)
cell_df <- read_csv(cell_neighborhoods_csv, show_col_types = FALSE) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    zone = factor(zone, levels = zone_levels),
    snc_binary = snc == "SnC+"
  )
expr_sample <- read_csv(expression_by_sample_csv, show_col_types = FALSE) %>%
  mutate(age = factor(age, levels = c("Young", "Old")))

# 1) Hot-spot burden by sample/age.
hotspot_burden_by_sample <- hotspot_results %>%
  filter(metric == hotspot_metric) %>%
  group_by(age, sample) %>%
  summarise(
    n_hex = n(),
    n_hot_hex = sum(hotspot_class == "Hot spot", na.rm = TRUE),
    n_not_sig_hex = sum(hotspot_class == "Not significant", na.rm = TRUE),
    n_cells_in_hot_hexes = sum(n_cells[hotspot_class == "Hot spot"], na.rm = TRUE),
    n_snc_in_hot_hexes = sum(n_snc[hotspot_class == "Hot spot"], na.rm = TRUE),
    mean_snc_prop_hot_hexes = mean(snc_prop[hotspot_class == "Hot spot"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hot_hex_prop = n_hot_hex / n_hex,
    snc_prop_in_hot_hexes = if_else(
      n_cells_in_hot_hexes > 0,
      n_snc_in_hot_hexes / n_cells_in_hot_hexes,
      NA_real_
    )
  )

hotspot_fraction_intensity <- hotspot_burden_by_sample %>%
  transmute(
    age,
    sample,
    n_hot_hex,
    n_hex,
    hot_hex_prop,
    n_cells_in_hot_hexes,
    n_snc_in_hot_hexes,
    snc_prop_in_hot_hexes,
    mean_snc_prop_hot_hexes
  )

hotspot_burden_by_age <- hotspot_burden_by_sample %>%
  group_by(age) %>%
  summarise(
    n_samples = n(),
    total_hot_hex = sum(n_hot_hex),
    total_hex = sum(n_hex),
    mean_hot_hex_prop = mean(hot_hex_prop),
    sem_hot_hex_prop = sd(hot_hex_prop) / sqrt(n()),
    total_cells_in_hot_hexes = sum(n_cells_in_hot_hexes),
    total_snc_in_hot_hexes = sum(n_snc_in_hot_hexes),
    pooled_snc_prop_in_hot_hexes = total_snc_in_hot_hexes / total_cells_in_hot_hexes,
    .groups = "drop"
  )

write_csv(hotspot_burden_by_sample, file.path(output_dir, "young_old_hotspot_burden_by_sample.csv"))
write_csv(hotspot_burden_by_age, file.path(output_dir, "young_old_hotspot_burden_by_age.csv"))
write_csv(hotspot_fraction_intensity, file.path(output_dir, "young_old_hotspot_fraction_intensity_by_sample.csv"))

# 2) Cell-type composition inside hot spots and neighboring zones.
sample_composition <- bind_rows(lapply(annotation_cols, function(col) make_sample_composition(cell_df, col)))
composition <- make_age_composition(sample_composition)
composition_age_test <- composition_age_tests(sample_composition)

write_csv(composition, file.path(output_dir, "young_old_celltype_composition_by_zone.csv"))
write_csv(sample_composition, file.path(output_dir, "young_old_celltype_composition_by_sample_zone.csv"))
write_csv(composition_age_test, file.path(output_dir, "young_old_celltype_composition_age_tests.csv"))

# 3) SASP by age/zone, using the same cell-level sasp_num used in the SnC gate.
sasp_by_sample_zone <- cell_df %>%
  group_by(age, sample, zone) %>%
  summarise(
    n_cells = n(),
    n_snc = sum(snc_binary, na.rm = TRUE),
    snc_prop = mean(snc_binary, na.rm = TRUE),
    mean_sasp_num = mean(sasp_num, na.rm = TRUE),
    median_sasp_num = median(sasp_num, na.rm = TRUE),
    pct_sasp_num_ge5 = mean(sasp_num >= 5, na.rm = TRUE),
    .groups = "drop"
  )

sasp_by_age_zone <- sasp_by_sample_zone %>%
  group_by(age, zone) %>%
  summarise(
    n_samples = n_distinct(sample),
    pooled_cells = sum(n_cells),
    pooled_snc = sum(n_snc),
    pooled_snc_prop = pooled_snc / pooled_cells,
    mean_sasp_num_sample_mean = mean(mean_sasp_num, na.rm = TRUE),
    sem_sasp_num_sample_mean = sd(mean_sasp_num, na.rm = TRUE) / sqrt(n()),
    mean_pct_sasp_ge5 = mean(pct_sasp_num_ge5, na.rm = TRUE),
    sem_pct_sasp_ge5 = sd(pct_sasp_num_ge5, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

write_csv(sasp_by_sample_zone, file.path(output_dir, "young_old_sasp_by_sample_zone.csv"))
write_csv(sasp_by_age_zone, file.path(output_dir, "young_old_sasp_by_age_zone.csv"))

# Raw summed SASP counts from the expression summary.
sasp_features <- c("SASP_total_count", "SASP_num_detected", "SASP_mean_count_per_gene")
sasp_raw_by_sample <- expr_sample %>%
  filter(metric == hotspot_metric, feature %in% sasp_features) %>%
  select(metric, age, sample, hotspot_class, feature, n_cells, mean_count, median_count, pct_expressing)

sasp_raw_by_age <- sasp_raw_by_sample %>%
  group_by(age, hotspot_class, feature) %>%
  summarise(
    n_samples = n_distinct(sample),
    pooled_cells_in_summary = sum(n_cells),
    mean_of_sample_means = mean(mean_count, na.rm = TRUE),
    sem_of_sample_means = sd(mean_count, na.rm = TRUE) / sqrt(n()),
    mean_pct_expressing = mean(pct_expressing, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(sasp_raw_by_sample, file.path(output_dir, "young_old_raw_sasp_by_sample_hotspot_class.csv"))
write_csv(sasp_raw_by_age, file.path(output_dir, "young_old_raw_sasp_by_age_hotspot_class.csv"))

# 4) Sample-level enrichment of neighboring cells within each age.
sample_enrichment <- bind_rows(lapply(annotation_cols, function(col) make_sample_enrichment(cell_df, col)))
enrichment_by_age <- summarise_sample_enrichment(sample_enrichment)
enrichment_age_test <- sample_enrichment_age_tests(sample_enrichment)

write_csv(sample_enrichment, file.path(output_dir, "young_old_neighbor_enrichment_by_sample.csv"))
write_csv(enrichment_by_age, file.path(output_dir, "young_old_neighbor_enrichment_by_age.csv"))
write_csv(enrichment_age_test, file.path(output_dir, "young_old_neighbor_enrichment_age_tests.csv"))

top_enrichment <- enrichment_by_age %>%
  filter(is.finite(mean_log2_ratio), fdr < 0.05, mean_log2_ratio > 0) %>%
  group_by(annotation_col, age) %>%
  arrange(desc(mean_log2_ratio), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

significant_neighbor_effects <- enrichment_by_age %>%
  filter(is.finite(mean_log2_ratio), fdr < 0.05) %>%
  group_by(annotation_col, age) %>%
  arrange(desc(abs(mean_log2_ratio)), .by_group = TRUE) %>%
  ungroup()

write_csv(top_enrichment, file.path(output_dir, "young_old_top_neighbor_enrichments_by_age.csv"))
write_csv(significant_neighbor_effects, file.path(output_dir, "young_old_significant_neighbor_effects_by_age.csv"))

# Plots.
for (ann in annotation_cols) {
  plot_composition(
    composition,
    annotation_col = ann,
    zones_to_plot = c("hotspot_core", "neighbor_1"),
    file_stem = paste0("young_old_", ann, "_composition_hotspot_neighbor")
  )
}

p_sasp <- sasp_by_age_zone %>%
  filter(zone %in% c("hotspot_core", "neighbor_1", "background")) %>%
  mutate(zone = factor(zone, levels = zone_levels)) %>%
  ggplot(aes(zone, mean_sasp_num_sample_mean, fill = age)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(
    aes(
      ymin = mean_sasp_num_sample_mean - sem_sasp_num_sample_mean,
      ymax = mean_sasp_num_sample_mean + sem_sasp_num_sample_mean
    ),
    position = position_dodge(width = 0.75),
    width = 0.2
  ) +
  scale_fill_manual(values = c("Young" = "#5DA5DA", "Old" = "#D7191C")) +
  labs(
    title = "SASP genes detected per cell by SnC hot-spot distance",
    x = NULL,
    y = "Mean sasp_num per cell",
    fill = NULL
  ) +
  theme_classic(base_size = 10)

ggsave(file.path(output_dir, "young_old_sasp_by_hotspot_distance.pdf"), p_sasp, width = 8, height = 5)
ggsave(file.path(output_dir, "young_old_sasp_by_hotspot_distance.png"), p_sasp, width = 8, height = 5, dpi = 300)

p_burden <- hotspot_burden_by_sample %>%
  ggplot(aes(age, hot_hex_prop, color = age)) +
  geom_point(size = 2.5, position = position_jitter(width = 0.08, height = 0)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black") +
  scale_y_continuous(labels = percent_format()) +
  scale_color_manual(values = c("Young" = "#5DA5DA", "Old" = "#D7191C")) +
  labs(
    title = "Fraction of 35 um hex bins called SnC hot spots",
    x = NULL,
    y = "Hot-spot hex fraction",
    color = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "young_old_hotspot_burden_by_sample.pdf"), p_burden, width = 5, height = 5)
ggsave(file.path(output_dir, "young_old_hotspot_burden_by_sample.png"), p_burden, width = 5, height = 5, dpi = 300)

p_fraction_intensity <- hotspot_fraction_intensity %>%
  ggplot(aes(hot_hex_prop, snc_prop_in_hot_hexes, color = age)) +
  geom_point(aes(size = n_snc_in_hot_hexes), alpha = 0.85) +
  geom_text(
    aes(label = sample),
    nudge_y = 0.006,
    size = 3,
    show.legend = FALSE,
    check_overlap = TRUE
  ) +
  scale_x_continuous(labels = percent_format()) +
  scale_y_continuous(labels = percent_format()) +
  scale_size_continuous(range = c(2.5, 7), breaks = pretty_breaks(n = 4)) +
  scale_color_manual(values = c("Young" = "#5DA5DA", "Old" = "#D7191C")) +
  labs(
    title = "SnC hot-spot extent and intensity by sample",
    x = "Hot-spot hex fraction",
    y = "SnC fraction inside hot-spot hexes",
    color = NULL,
    size = "SnC cells\nin hot spots"
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "young_old_hotspot_fraction_vs_intensity.pdf"), p_fraction_intensity, width = 6.5, height = 5.2)
ggsave(file.path(output_dir, "young_old_hotspot_fraction_vs_intensity.png"), p_fraction_intensity, width = 6.5, height = 5.2, dpi = 300)

message("Hot-spot burden by age:")
print(hotspot_burden_by_age)

message("SASP by age/zone:")
print(sasp_by_age_zone)

message("Significant neighbor effects:")
print(significant_neighbor_effects %>% select(annotation_col, age, cell_type, mean_near_prop, mean_reference_prop, mean_prop_diff, mean_log2_ratio, fdr))

message("Done. Results written to ", output_dir)
