##########################################################
# Xenium SnC CDKN subtype burden by sample
#
# Denominator = all cells in each sample
# Bars = p21+, p16+, and p16+p21+ SnC+ cells
# Fill = Xenium minor cell type
##########################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(grid)
  library(patchwork)
})

project_dir <- "//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2"
setwd(project_dir)

input_rds <- "data/skin.rds"
outdir <- "results/260830/xenium_snc_cdkn_celltype_burden_by_sample"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

hsskin.ST <- readRDS(input_rds)
required_metadata <- c("sample", "snc", "cell_type_minor", "cell_type_major")
missing_metadata <- setdiff(required_metadata, colnames(hsskin.ST@meta.data))
if (length(missing_metadata) > 0) {
  stop("Missing required Xenium metadata: ", paste(missing_metadata, collapse = ", "))
}

marker_genes <- c("CDKN1A", "CDKN2A")
if (!all(marker_genes %in% rownames(hsskin.ST))) {
  stop("CDKN1A and/or CDKN2A are absent from the Xenium object.")
}

# Use detected transcripts, not normalized expression, to define p21/p16 status.
marker_counts <- FetchData(
  hsskin.ST,
  vars = marker_genes,
  layer = "counts"
)

md <- hsskin.ST@meta.data %>%
  mutate(
    cell_id = rownames(hsskin.ST@meta.data),
    sample = as.character(sample)
  ) %>%
  left_join(
    marker_counts %>%
      mutate(cell_id = rownames(marker_counts)),
    by = "cell_id"
  ) %>%
  filter(!is.na(sample)) %>%
  mutate(
    # Xenium sample suffix is the authoritative age mapping.
    age = case_when(
      grepl("-Y$", sample, ignore.case = TRUE) ~ "Young",
      grepl("-O$", sample, ignore.case = TRUE) ~ "Old",
      TRUE ~ NA_character_
    ),
    age = factor(age, levels = c("Young", "Old")),
    p21_pos = CDKN1A > 0,
    p16_pos = CDKN2A > 0,
    cdkn_state = case_when(
      p21_pos & p16_pos ~ "p16+p21+",
      p21_pos ~ "p21+",
      p16_pos ~ "p16+",
      TRUE ~ NA_character_
    ),
    cell_type_minor = as.character(cell_type_minor),
    cell_type_major = as.character(cell_type_major)
  )

if (any(is.na(md$age))) {
  stop(
    "Samples without a valid -Y/-O suffix: ",
    paste(sort(unique(md$sample[is.na(md$age)])), collapse = ", ")
  )
}

age_audit <- md %>% distinct(sample, age) %>% arrange(age, sample)
if (nrow(age_audit) != 6 || sum(age_audit$age == "Young") != 3 || sum(age_audit$age == "Old") != 3) {
  stop("Expected exactly three Young (-Y) and three Old (-O) Xenium samples.")
}
write.csv(age_audit, file.path(outdir, "xenium_sample_age_audit.csv"), row.names = FALSE)
print(age_audit)

cdkn_order <- c("p21+", "p16+", "p16+p21+")

minor_colors <- c(
  # Epidermal lineages: anchored to CODEX/4i basal red and suprabasal blue.
  "Basal Keratinocytes" = "#b2182b",
  "Bulge Keratinocytes" = "#d6604d",
  "Spinous Keratinocytes" = "#2166ac",
  "Granular Keratinocytes" = "#67a9cf",
  "Hair Follicle Cells" = "#ef8a62",

  # Fibroblast lineages: shades around the CODEX/4i fibroblast navy.
  "Papillary Fibroblasts" = "#4575b4",
  "Reticular Fibroblasts" = "#053061",
  "Dermal Fibroblasts" = "#0b3c6f",

  # Structural/vascular lineages.
  "Endothelial Cells" = "#33a02c",
  "Muscle Cells" = "#4393c3",

  # Immune lineages use distinct complementary colors.
  "CD4+ T Cells" = "#984ea3",
  "CD8+ T Cells/NK Cells" = "#6a3d9a",
  "B Cells" = "#cab2d6",
  "Dendritic Cells" = "#ffbf00",
  "Langerhans Cells" = "#fdbf6f",
  "Macrophages" = "#01665e",
  "Mast Cells" = "#ff7f00",

  # Exact CODEX/4i melanocyte color.
  "Melanocytes" = "#d1e5f0"
)

observed_types <- sort(unique(md$cell_type_minor[!is.na(md$cell_type_minor)]))
missing_colors <- setdiff(observed_types, names(minor_colors))
if (length(missing_colors) > 0) {
  stop("Cell types missing palette colors: ", paste(missing_colors, collapse = ", "))
}
celltype_order <- names(minor_colors)

sample_totals <- md %>% count(sample, age, name = "total_cells")

snc_counts <- md %>%
  filter(
    snc == "SnC+",
    !is.na(cdkn_state),
    !is.na(cell_type_minor)
  ) %>%
  count(sample, age, cdkn_state, cell_type_minor, name = "snc_count")

snc_cdkn_sample <- snc_counts %>%
  complete(
    nesting(sample, age),
    cdkn_state = cdkn_order,
    cell_type_minor = celltype_order,
    fill = list(snc_count = 0)
  ) %>%
  left_join(sample_totals, by = c("sample", "age")) %>%
  mutate(
    proportion_total = snc_count / total_cells,
    cdkn_state = factor(cdkn_state, levels = cdkn_order),
    cell_type_minor = factor(cell_type_minor, levels = celltype_order)
  )

sample_order <- sample_totals %>%
  arrange(age, sample) %>%
  pull(sample) %>%
  unique()

snc_cdkn_sample <- snc_cdkn_sample %>%
  mutate(sample = factor(sample, levels = sample_order))

sample_totals <- sample_totals %>%
  mutate(sample = factor(sample, levels = sample_order))

bar_totals <- snc_cdkn_sample %>%
  group_by(sample, age, cdkn_state) %>%
  summarise(
    n_snc = sum(snc_count),
    total_cells = first(total_cells),
    proportion_total = sum(proportion_total),
    .groups = "drop"
  )

p_main <- ggplot(
  snc_cdkn_sample,
  aes(x = cdkn_state, y = proportion_total, fill = cell_type_minor)
) +
  geom_col(width = 0.75, color = "white", linewidth = 0.35) +
  facet_wrap(~sample, nrow = 1, strip.position = "bottom") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = minor_colors, breaks = celltype_order, drop = FALSE) +
  labs(x = NULL, y = "% of Total Cells", fill = "Cell Type") +
  theme_classic(base_size = 14) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 11, color = "black", margin = margin(t = 6)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 10, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(face = "bold"),
    panel.spacing.x = unit(0.8, "lines"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(t = 0, r = 15, b = 15, l = 15)
  )

age_strip_data <- sample_totals %>%
  arrange(sample) %>%
  mutate(x = seq_along(sample))

age_group_labels <- age_strip_data %>%
  group_by(age) %>%
  summarise(x = mean(x), .groups = "drop")

p_age_strip <- ggplot(age_strip_data, aes(x = x, y = 1, fill = age)) +
  geom_tile(width = 1, height = 0.7) +
  scale_fill_manual(values = c("Young" = "#BDBDBD", "Old" = "#525252")) +
  geom_text(
    data = age_group_labels,
    aes(x = x, y = 1, label = age),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4.5,
    color = "white"
  ) +
  scale_x_continuous(limits = c(0.5, length(sample_order) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 1.5), expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(t = 0, r = 15, b = 0, l = 15))

# Put the legend in its own column so the age strip aligns with sample facets.
main_legend <- cowplot::get_legend(p_main)
p_plot_column <-
  p_age_strip /
  (p_main + theme(legend.position = "none")) +
  plot_layout(heights = c(0.35, 6))

p_final <-
  p_plot_column |
  wrap_elements(main_legend) +
  plot_layout(widths = c(1, 0.20))

print(p_final)

ggsave(
  file.path(outdir, "xenium_snc_cdkn_celltype_burden_by_sample.png"),
  p_final,
  width = 20,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(outdir, "xenium_snc_cdkn_celltype_burden_by_sample.pdf"),
  p_final,
  width = 20,
  height = 7,
  units = "in",
  bg = "white"
)

write.csv(
  snc_cdkn_sample,
  file.path(outdir, "xenium_snc_cdkn_celltype_burden_by_sample.csv"),
  row.names = FALSE
)
write.csv(
  bar_totals,
  file.path(outdir, "xenium_snc_cdkn_subtype_burden_by_sample.csv"),
  row.names = FALSE
)


##########################################################
# Major-cell-type version of the same plot
##########################################################

major_colors <- c(
  # Shared colors from the current 4i/CODEX CDKN burden figures.
  # Keratinocytes use the Basal Keratinocyte red as the single major-lineage color.
  "Keratinocytes" = "#b2182b",
  "Fibroblasts" = "#053061",
  "Melanocytes" = "#d1e5f0",
  "Muscle Cells" = "#4393c3",
  # Xenium-only major lineages: kept distinct from the shared colors.
  "Blood Vessel Cells" = "#33a02c",
  "Immune Cells" = "#6a3d9a"
)
major_order <- names(major_colors)

observed_major <- sort(unique(md$cell_type_major[!is.na(md$cell_type_major)]))
missing_major_colors <- setdiff(observed_major, major_order)
if (length(missing_major_colors) > 0) {
  stop(
    "Major cell types missing palette colors: ",
    paste(missing_major_colors, collapse = ", ")
  )
}

snc_major_counts <- md %>%
  filter(
    snc == "SnC+",
    !is.na(cdkn_state),
    !is.na(cell_type_major)
  ) %>%
  count(sample, age, cdkn_state, cell_type_major, name = "snc_count")

snc_cdkn_major_sample <- sample_totals %>%
  select(sample, age) %>%
  crossing(
    cdkn_state = cdkn_order,
    cell_type_major = major_order
  ) %>%
  left_join(
    snc_major_counts,
    by = c("sample", "age", "cdkn_state", "cell_type_major")
  ) %>%
  mutate(snc_count = replace_na(snc_count, 0L)) %>%
  left_join(sample_totals, by = c("sample", "age")) %>%
  mutate(
    proportion_total = snc_count / total_cells,
    sample = factor(sample, levels = sample_order),
    cdkn_state = factor(cdkn_state, levels = cdkn_order),
    cell_type_major = factor(cell_type_major, levels = major_order)
  )

major_bar_totals <- snc_cdkn_major_sample %>%
  group_by(sample, age, cdkn_state) %>%
  summarise(
    n_snc = sum(snc_count),
    total_cells = first(total_cells),
    proportion_total = sum(proportion_total),
    .groups = "drop"
  )

p_major_main <- ggplot(
  snc_cdkn_major_sample,
  aes(x = cdkn_state, y = proportion_total, fill = cell_type_major)
) +
  geom_col(width = 0.75, color = "white", linewidth = 0.35) +
  facet_wrap(~sample, nrow = 1, strip.position = "bottom") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = major_colors, breaks = major_order, drop = FALSE) +
  labs(x = NULL, y = "% of Total Cells", fill = "Major Cell Type") +
  theme_classic(base_size = 14) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 11, color = "black", margin = margin(t = 6)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 10, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(face = "bold"),
    panel.spacing.x = unit(0.8, "lines"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(t = 0, r = 15, b = 15, l = 15)
  )

major_legend <- cowplot::get_legend(p_major_main)
p_major_column <-
  p_age_strip /
  (p_major_main + theme(legend.position = "none")) +
  plot_layout(heights = c(0.35, 6))
p_major_final <-
  p_major_column |
  wrap_elements(major_legend) +
  plot_layout(widths = c(1, 0.20))

print(p_major_final)

ggsave(
  file.path(outdir, "xenium_snc_cdkn_celltype_major_burden_by_sample.png"),
  p_major_final,
  width = 20,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(outdir, "xenium_snc_cdkn_celltype_major_burden_by_sample.pdf"),
  p_major_final,
  width = 20,
  height = 7,
  units = "in",
  bg = "white"
)
write.csv(
  snc_cdkn_major_sample,
  file.path(outdir, "xenium_snc_cdkn_celltype_major_burden_by_sample.csv"),
  row.names = FALSE
)
write.csv(
  major_bar_totals,
  file.path(outdir, "xenium_snc_cdkn_celltype_major_subtype_totals_by_sample.csv"),
  row.names = FALSE
)

cat("\nSaved Xenium figure and tables to:\n", outdir, "\n")
