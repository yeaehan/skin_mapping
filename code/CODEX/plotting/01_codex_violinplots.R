library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)

# Plot single-cell violins with sample-level mean dots and sample-level statistics.
plot_seurat_violin_with_sampledots <- function(obj,
                                               feature,
                                               group_col = "snc",
                                               sample_col = "sample",
                                               assay = "RNA",
                                               layer = "data",
                                               cols = c("SnC-" = "#c8def9", "SnC+" = "#d7b0b0"),
                                               dot_cols = c("SnC-" = "#000080", "SnC+" = "#a00000")) {
  
  # 1) Base Seurat violin
  p <- VlnPlot(
    obj,
    features = feature,
    group.by = group_col,
    assay = assay,
    layer = layer,
    pt.size = 0,        # remove cell-level dots
    cols = cols
  )
  
  # 2) Fetch data for sample-level dots
  df <- FetchData(
    obj,
    vars = c(feature, group_col, sample_col),
    assay = assay,
    layer = layer
  )
  colnames(df) <- c("value", "group", "sample")
  
  df <- df %>%
    filter(!is.na(group), !is.na(sample)) %>%
    mutate(group = factor(group, levels = c("SnC-", "SnC+")))
  
  # 3) Compute sample means (Pseudobulk)
  df_sample <- df %>%
    group_by(sample, group) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  ymax <- max(df$value, na.rm = TRUE)
  
  # 4) Add sample dots and map explicit aesthetics to stats
  p <- p +
    geom_point(
      data = df_sample,
      aes(x = group, y = mean_value, color = group),
      position = position_jitter(width = 0.08, height = 0, seed = 42),
      size = 3,
      alpha = 0.9,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = dot_cols) +
    
    stat_summary(fun = mean, geom = "crossbar",
                 width = 0.4, color = "black", linewidth = 0.8) +
    
    # FIX: Add explicit mapping inside stat_compare_means so it prints on the plot
    stat_compare_means(
      data = df_sample,
      aes(x = group, y = mean_value),
      method = "t.test", 
      label = "p.format",     # Shows clean 'p = 0.045' format instead of wild asterisks
      label.x = 1.3,
      label.y = ymax * 1,  # Slightly above ymax to keep it safe from cutoffs
      size = 5,
      inherit.aes = FALSE     # Tells ggpubr to ignore raw single-cell dimensions
    ) +
    
    theme_classic(base_size = 18) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    labs(title = feature, y = "Expression")
  
  return(p)
}

obj_codex$age <- factor(obj_codex$age, levels = c("Young", "Old"))
obj_codex$snc <- factor(obj_codex$snc, levels = c("SnC-", "SnC+"))


features <- c("BCL2", "HMGB1", "LaminB1")

p_list <- lapply(features, function(f) {
  plot_seurat_violin_with_sampledots(
    obj = obj_codex,
    feature = f,
    group_col = "snc",
    sample_col = "sample",
    assay = "RNA",
    layer = "data"
  )
})

final_plot <- p_list[[1]] + p_list[[2]] + p_list[[3]]

final_plot

ggsave(
  "//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2/results/260523/codex_violins_subsenmarkers_sampledots_snc.png",
  plot = final_plot,
  width = 3,
  height = 12,
  dpi = 300
)

ggsave(
  "//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2/results/260523/codex_violins_subsenmarkers_sampledots_snc.pdf",
  plot = final_plot,
  width = 3,
  height = 12,
  dpi = 300
)
