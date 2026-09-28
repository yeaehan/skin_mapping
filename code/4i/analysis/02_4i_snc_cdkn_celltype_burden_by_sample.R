##########################################################
# 4i SnC CDKN subtype burden by sample
#
# Denominator = all cells in each sample
# Bars = p21+, p16+, and p16+p21+ SnC+ cells
# Fill = 4i minor skin cell type
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

input_rdata <- "data/4i/seurat_list_4i_annotated.RData"
outdir <- "results/260830/4i_snc_cdkn_celltype_burden_by_sample"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

loaded <- new.env(parent = emptyenv())
load(input_rdata, envir = loaded)
if (!exists("seurat_list", envir = loaded, inherits = FALSE)) {
  stop("seurat_list was not found in ", input_rdata)
}
seurat_list <- loaded$seurat_list

sample_id_from_name <- function(x) {
  case_when(
    grepl("CAMS18a", x, ignore.case = TRUE) ~ "CAMS18a",
    grepl("CAMS16a", x, ignore.case = TRUE) ~ "CAMS16a",
    grepl("CAMS5a", x, ignore.case = TRUE) ~ "CAMS5a",
    grepl("Subject1", x, ignore.case = TRUE) ~ "Subject1",
    grepl("CAMS38a", x, ignore.case = TRUE) ~ "CAMS38a",
    grepl("0201UA", x, ignore.case = TRUE) ~ "0201UA",
    grepl("202UA", x, ignore.case = TRUE) ~ "202UA",
    grepl("0203UA", x, ignore.case = TRUE) ~ "0203UA",
    grepl("0204UA", x, ignore.case = TRUE) ~ "0204UA",
    grepl("207UA", x, ignore.case = TRUE) ~ "207UA",
    TRUE ~ NA_character_
  )
}

young_samples <- c("0201UA", "202UA", "0203UA", "0204UA", "207UA")
old_samples <- c("CAMS18a", "CAMS16a", "CAMS5a", "Subject1", "CAMS38a")

md <- bind_rows(lapply(names(seurat_list), function(object_name) {
  obj <- seurat_list[[object_name]]
  required <- c("snc", "cdkn_state", "minorskincell")
  missing <- setdiff(required, colnames(obj@meta.data))
  if (length(missing) > 0) {
    stop(object_name, " is missing metadata: ", paste(missing, collapse = ", "))
  }
  obj@meta.data %>%
    transmute(
      cell_id = rownames(obj@meta.data),
      object_name = object_name,
      sample = sample_id_from_name(object_name),
      snc = as.character(snc),
      cdkn_state = as.character(cdkn_state),
      minorskincell = as.character(minorskincell)
    )
})) %>%
  mutate(
    age = case_when(
      sample %in% young_samples ~ "Young",
      sample %in% old_samples ~ "Old",
      TRUE ~ NA_character_
    ),
    age = factor(age, levels = c("Young", "Old"))
  )

if (any(is.na(md$sample)) || any(is.na(md$age))) {
  stop("At least one 4i object could not be assigned to the explicit sample/age map.")
}

age_audit <- md %>% distinct(object_name, sample, age) %>% arrange(age, sample)
if (
  nrow(age_audit) != 10 ||
  !setequal(age_audit$sample[age_audit$age == "Young"], young_samples) ||
  !setequal(age_audit$sample[age_audit$age == "Old"], old_samples)
) {
  stop("4i sample-age validation failed.")
}
write.csv(age_audit, file.path(outdir, "4i_sample_age_audit.csv"), row.names = FALSE)
print(age_audit %>% select(sample, age))

cdkn_order <- c("p21+", "p16+", "p16+p21+")
celltype_order <- c(
  "Fibroblast",
  "Melanocyte",
  "Basal Keratinocyte",
  "Suprabasal Keratinocyte",
  "Corneocyte",
  "none"
)
celltype_colors <- c(
  "Fibroblast" = "#053061",
  "Melanocyte" = "#d1e5f0",
  "Basal Keratinocyte" = "#b2182b",
  "Suprabasal Keratinocyte" = "#2166ac",
  "Corneocyte" = "#fddbc7",
  "none" = "#d9d9d9"
)

unknown_types <- setdiff(unique(md$minorskincell), celltype_order)
unknown_types <- unknown_types[!is.na(unknown_types)]
if (length(unknown_types) > 0) {
  stop("Unexpected 4i minor cell types: ", paste(unknown_types, collapse = ", "))
}

sample_totals <- md %>% count(sample, age, name = "total_cells")

snc_counts <- md %>%
  filter(
    snc == "SnC+",
    cdkn_state %in% cdkn_order,
    !is.na(minorskincell)
  ) %>%
  count(sample, age, cdkn_state, minorskincell, name = "snc_count")

# Start from every sample x CDKN state x cell type combination. This retains
# samples with zero SnC+ cells as explicit zero-height bars.
snc_cdkn_sample <- sample_totals %>%
  select(sample, age) %>%
  crossing(
    cdkn_state = cdkn_order,
    minorskincell = celltype_order
  ) %>%
  left_join(
    snc_counts,
    by = c("sample", "age", "cdkn_state", "minorskincell")
  ) %>%
  mutate(snc_count = replace_na(snc_count, 0L)) %>%
  left_join(sample_totals, by = c("sample", "age")) %>%
  mutate(
    proportion_total = snc_count / total_cells,
    cdkn_state = factor(cdkn_state, levels = cdkn_order),
    minorskincell = factor(minorskincell, levels = celltype_order)
  )

sample_order <- sample_totals %>%
  arrange(age, sample) %>%
  pull(sample)

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
  aes(x = cdkn_state, y = proportion_total, fill = minorskincell)
) +
  geom_col(width = 0.75, color = "white", linewidth = 0.35) +
  facet_wrap(~sample, nrow = 1, strip.position = "bottom") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = celltype_colors, breaks = celltype_order, drop = FALSE) +
  labs(x = NULL, y = "% of Total Cells", fill = "Cell Type") +
  theme_classic(base_size = 14) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 10, color = "black", margin = margin(t = 6)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 9, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(face = "bold"),
    panel.spacing.x = unit(0.65, "lines"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
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

# Place the legend outside the strip/main column to keep exact strip alignment.
main_legend <- cowplot::get_legend(p_main)
p_plot_column <-
  p_age_strip /
  (p_main + theme(legend.position = "none")) +
  plot_layout(heights = c(0.35, 6))
p_final <-
  p_plot_column |
  wrap_elements(main_legend) +
  plot_layout(widths = c(1, 0.14))

print(p_final)

ggsave(
  file.path(outdir, "4i_snc_cdkn_celltype_burden_by_sample.png"),
  p_final,
  width = 24,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(outdir, "4i_snc_cdkn_celltype_burden_by_sample.pdf"),
  p_final,
  width = 24,
  height = 7,
  units = "in",
  bg = "white"
)

write.csv(
  snc_cdkn_sample,
  file.path(outdir, "4i_snc_cdkn_celltype_burden_by_sample.csv"),
  row.names = FALSE
)
write.csv(
  bar_totals,
  file.path(outdir, "4i_snc_cdkn_subtype_burden_by_sample.csv"),
  row.names = FALSE
)

cat("\nSaved 4i figure and tables to:\n", outdir, "\n")
