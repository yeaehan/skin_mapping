##########################################################
# CODEX SnC CDKN subtype burden by sample
#
# Denominator = ALL CELLS IN EACH SAMPLE
#
# Three bars per sample:
#   p21+
#   p16+
#   p16+p21+
#
# Bar height = proportion of ALL cells in sample
# Fill = minor skin cell type
#
# Samples ordered:
#   Young -> Old
#
# Top annotation:
#   continuous Young / Old strip
##########################################################


##########################################################
# 0. Libraries + data
##########################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(grid)
library(patchwork)

setwd(
  "//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2"
)

load(
  "data/codex/seurat_list_codex_annotated.RData"
)

obj_codex <- merge(
  seurat_list_codex[[1]],
  y = seurat_list_codex[-1],
  merge.data = TRUE
)

cat("CODEX data loaded.\n")


##########################################################
# 1. Plot settings
##########################################################

cdkn_order <- c(
  "p21+",
  "p16+",
  "p16+p21+"
)

celltype_order <- c(
  "Fibroblast",
  "Melanocyte",
  "Basal Keratinocyte",
  "Suprabasal Keratinocyte",
  "Corneocyte",
  "Smooth Muscle",
  "none"
)

celltype_colors <- c(
  "Fibroblast"              = "#053061",
  "Melanocyte"              = "#d1e5f0",
  "Basal Keratinocyte"      = "#b2182b",
  "Suprabasal Keratinocyte" = "#2166ac",
  "Corneocyte"              = "#fddbc7",
  "Smooth Muscle"           = "#4393c3",
  "none"                    = "#d9d9d9"
)


##########################################################
# 2. Prepare metadata
##########################################################

md <- obj_codex@meta.data %>%
  filter(
    !is.na(sample)
  ) %>%
  mutate(
    sample = as.character(sample),
    # Assign age from sample identity. Do not trust inherited metadata:
    # only these two CODEX samples are Young; every other sample is Old.
    age = if_else(
      tolower(trimws(sample)) %in% c("207ua", "7503"),
      "Young",
      "Old"
    ),
    age = factor(
      age,
      levels = c(
        "Young",
        "Old"
      )
    )
  )

young_samples_detected <- md %>%
  filter(age == "Young") %>%
  distinct(sample) %>%
  pull(sample)

if (
  length(young_samples_detected) != 2 ||
  !setequal(tolower(young_samples_detected), c("207ua", "7503"))
) {
  stop(
    "Age-label validation failed. Young samples must be exactly 207UA and 7503."
  )
}


##########################################################
# Check samples
##########################################################

cat("\nSamples detected:\n")

print(
  md %>%
    distinct(
      sample,
      age
    ) %>%
    arrange(
      age,
      sample
    )
)


##########################################################
# 3. Total cells per sample
#
# This is the denominator.
##########################################################

sample_totals <- md %>%
  count(
    sample,
    age,
    name = "total_cells"
  )

cat("\nTotal cells per sample:\n")

print(sample_totals)


##########################################################
# 4. Count SnC+ cells
#
# Group by:
#   sample
#   age
#   CDKN subtype
#   minor skin cell type
##########################################################

snc_counts <- md %>%
  filter(
    snc == "SnC+",
    !is.na(cdkn_state),
    !is.na(minorskincell1),
    cdkn_state %in% cdkn_order
  ) %>%
  count(
    sample,
    age,
    cdkn_state,
    minorskincell1,
    name = "snc_count"
  )


##########################################################
# 5. Make sure every sample has all 3 CDKN bars
#
# Even if there are zero:
#
#   p21+
#   p16+
#   p16+p21+
##########################################################

snc_cdkn_sample <- snc_counts %>%

  complete(
    nesting(
      sample,
      age
    ),

    cdkn_state = cdkn_order,

    minorskincell1 = celltype_order,

    fill = list(
      snc_count = 0
    )
  ) %>%

  left_join(
    sample_totals,
    by = c(
      "sample",
      "age"
    )
  ) %>%

  mutate(

    # IMPORTANT:
    # denominator = ALL cells in sample
    proportion_total =
      snc_count / total_cells,

    cdkn_state = factor(
      cdkn_state,
      levels = cdkn_order
    ),

    minorskincell1 = factor(
      minorskincell1,
      levels = celltype_order
    )
  )


##########################################################
# 6. Order samples
#
# Young first
# Old second
##########################################################

sample_order <- sample_totals %>%
  arrange(
    age,
    sample
  ) %>%
  pull(sample) %>%
  unique()


snc_cdkn_sample <- snc_cdkn_sample %>%
  mutate(
    sample = factor(
      sample,
      levels = sample_order
    )
  )


sample_totals <- sample_totals %>%
  mutate(
    sample = factor(
      sample,
      levels = sample_order
    )
  )


##########################################################
# 7. Calculate total burden for each bar
##########################################################

bar_totals <- snc_cdkn_sample %>%

  group_by(
    sample,
    age,
    cdkn_state
  ) %>%

  summarise(
    n_snc = sum(snc_count),

    total_cells =
      first(total_cells),

    proportion_total =
      sum(proportion_total),

    .groups = "drop"
  )


cat("\nSnC subtype burden per sample:\n")

print(
  bar_totals %>%
    mutate(
      percent_total =
        100 * proportion_total
    )
)


##########################################################
# 8. Main plot
##########################################################

cat("\nGenerating main plot...\n")

p_main <- ggplot(
  snc_cdkn_sample,
  aes(
    x = cdkn_state,
    y = proportion_total,
    fill = minorskincell1
  )
) +

  geom_col(
    width = 0.75,
    color = "white",
    linewidth = 0.4
  ) +

  # One sample per facet
  # All samples in ONE row
  # Sample names underneath
  facet_wrap(
    ~ sample,
    nrow = 1,
    strip.position = "bottom"
  ) +

  # Keep all three CDKN categories
  scale_x_discrete(
    drop = FALSE
  ) +

  # Y axis = % of ALL cells
  scale_y_continuous(
    labels = percent_format(
      accuracy = 0.1
    ),

    expand = expansion(
      mult = c(
        0,
        0.05
      )
    )
  ) +

  scale_fill_manual(
    values = celltype_colors,
    breaks = celltype_order,
    drop = FALSE
  ) +

  labs(
    x = NULL,
    y = "% of Total Cells",
    fill = "Cell Type"
  ) +

  theme_classic(
    base_size = 14
  ) +

  theme(

    # Sample names below plot
    strip.placement = "outside",

    strip.background =
      element_blank(),

    strip.text.x = element_text(
      face = "bold",
      size = 11,
      color = "black",
      margin = margin(
        t = 6
      )
    ),

    # p21 / p16 / double+
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "bold",
      size = 10,
      color = "black"
    ),

    axis.text.y = element_text(
      color = "black"
    ),

    axis.title.y = element_text(
      face = "bold"
    ),

    # Space between samples
    panel.spacing.x = unit(
      0.8,
      "lines"
    ),

    legend.position = "right",

    legend.title = element_text(
      face = "bold"
    ),

    legend.text = element_text(
      size = 11
    ),

    plot.margin = margin(
      t = 0,
      r = 15,
      b = 15,
      l = 15
    )
  )


##########################################################
# 9. Prepare Young / Old annotation strip
##########################################################

age_strip_data <- sample_totals %>%

  arrange(sample) %>%

  mutate(
    x = seq_along(sample)
  )


##########################################################
# Calculate center of Young and Old sections
##########################################################

age_group_labels <- age_strip_data %>%

  group_by(age) %>%

  summarise(
    xmin = min(x),
    xmax = max(x),
    x = mean(x),
    .groups = "drop"
  )


##########################################################
# 10. Make Young / Old top strip
##########################################################

p_age_strip <- ggplot(
  age_strip_data,
  aes(
    x = x,
    y = 1,
    fill = age
  )
) +

  geom_tile(
    width = 1,
    height = 0.7
  ) +

  scale_fill_manual(
    values = c(
      "Young" = "#BDBDBD",
      "Old"   = "#525252"
    )
  ) +

  # Put labels directly on strip
  geom_text(
    data = age_group_labels,
    aes(
      x = x,
      y = 1,
      label = age
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4.5,
    color = "white"
  ) +

  scale_x_continuous(
    limits = c(
      0.5,
      length(sample_order) + 0.5
    ),
    expand = c(
      0,
      0
    )
  ) +

  scale_y_continuous(
    limits = c(
      0.5,
      1.5
    ),
    expand = c(
      0,
      0
    )
  ) +

  theme_void() +

  theme(
    legend.position = "none",

    plot.margin = margin(
      t = 0,
      r = 15,
      b = 0,
      l = 15
    )
  )


##########################################################
# 11. Combine age strip + main plot
##########################################################

# Keep the legend outside the age-strip/main-plot column. Previously the
# strip used the width occupied by the legend, making its two Young tiles
# visually span additional Old sample facets.
main_legend <- cowplot::get_legend(p_main)

p_plot_column <-
  p_age_strip /
  (p_main + theme(legend.position = "none")) +
  plot_layout(
    heights = c(
      0.35,
      6
    )
  )

p_final <-
  p_plot_column |
  wrap_elements(main_legend) +
  plot_layout(
    widths = c(
      1,
      0.16
    )
  )


##########################################################
# 12. Display
##########################################################

print(p_final)


##########################################################
# 13. Save PNG
##########################################################

ggsave(
  filename =
    "results/codex_snc_cdkn_celltype_burden_by_sample.png",

  plot = p_final,

  width = 20,
  height = 6.5,

  units = "in",

  dpi = 300,

  bg = "white"
)


##########################################################
# 14. Save detailed values
##########################################################

write.csv(
  snc_cdkn_sample,

  "results/codex_snc_cdkn_celltype_burden_by_sample.csv",

  row.names = FALSE
)


##########################################################
# 15. Save total subtype burden
##########################################################

write.csv(
  bar_totals,

  "results/codex_snc_cdkn_subtype_burden_by_sample.csv",

  row.names = FALSE
)


##########################################################
# Done
##########################################################

cat(
  "\nSaved figure:\n",
  "results/codex_snc_cdkn_celltype_burden_by_sample.png\n"
)

cat("\nFinished.\n")
