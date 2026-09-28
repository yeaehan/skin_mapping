##########################################################
# CODEX cell analysis (copied from 4i logic)
##########################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(ggpubr)

setwd("//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2")
load("data/codex/seurat_list_codex.RData")  # loads seurat_list_codex
annotation_only <- "--annotate-only" %in% commandArgs(trailingOnly = TRUE)

##########################################################
# 1) Marker panel
##########################################################
marker_panel_codex <- c(
  "BCL2", "SMA", "p21", "Keratin 14", "Vimentin", "Keratin 5", "p16", "CD34",
  "IP10", "CD8", "hIL-1alpha", "CD14", "LaminB1", "Ki67", "Filaggrin", "CXCR3",
  "HMGB1", "Sirt1", "IDO1", "FOXP3", "IL6", "MelanA", "IL8", "CD4",
  "CD68", "HLA-DR", "CD45", "PCNA", "ICOS", "Pan-Cytokeratin", "CD3e", "CD11c"
)

##########################################################
# 2) Build marker_list_codex (counts > 0)
##########################################################
marker_list_codex <- list()

for (sample_name in names(seurat_list_codex)) {
  cat("Processing:", sample_name, "\n")
  obj <- seurat_list_codex[[sample_name]]
  
  sample_markers <- list()
  for (marker in marker_panel_codex) {
    if (marker %in% rownames(obj[["RNA"]]$counts)) {
      sample_markers[[marker]] <- obj[["RNA"]]$counts[marker, ] > 0
    } else {
      warning(paste("Marker", marker, "not found in", sample_name))
      sample_markers[[marker]] <- NULL
    }
  }
  marker_list_codex[[sample_name]] <- sample_markers
}

get_marker <- function(marker_vecs, marker, cells) {
  v <- marker_vecs[[marker]]
  if (is.null(v)) return(setNames(rep(FALSE, length(cells)), cells))
  v[cells]
}

##########################################################
# 3) Annotate cdkn_state + snc + major/minor + immune
##########################################################
for (sample_name in names(seurat_list_codex)) {
  cat("Annotating:", sample_name, "\n")
  
  obj <- seurat_list_codex[[sample_name]]
  cells <- colnames(obj)
  m <- marker_list_codex[[sample_name]]
  
  # --- markers used in rules ---
  p21   <- get_marker(m, "p21", cells)
  p16   <- get_marker(m, "p16", cells)
  PCNA  <- get_marker(m, "PCNA", cells)
  IP10  <- get_marker(m, "IP10", cells)
  IL6  <- get_marker(m, "IL6", cells)
  IL8  <- get_marker(m, "IL8", cells)
  CXCR3 <- get_marker(m, "CXCR3", cells)
  bcl2  <- get_marker(m, "BCL2", cells)
  
  MelanA    <- get_marker(m, "MelanA", cells)
  Vimentin  <- get_marker(m, "Vimentin", cells)
  SMA       <- get_marker(m, "SMA", cells)
  Keratin5  <- get_marker(m, "Keratin 5", cells)
  Keratin14 <- get_marker(m, "Keratin 14", cells)
  Filaggrin <- get_marker(m, "Filaggrin", cells)
  PanCyto <- get_marker(m, "Pan-Cytokeratin", cells)
  
  CD4   <- get_marker(m, "CD4", cells)
  CD68  <- get_marker(m, "CD68", cells)
  CD8   <- get_marker(m, "CD8", cells)
  FOXP3 <- get_marker(m, "FOXP3", cells)
  CD14  <- get_marker(m, "CD14", cells)
  CD45  <- get_marker(m, "CD45", cells)
  
  # --- NEW: cdkn_state (same as 4i) ---
  obj$cdkn_state <- dplyr::case_when(
    p16 & !p21 ~ "p16+",
    p21 & !p16 ~ "p21+",
    p16 & p21  ~ "p16+p21+",
    TRUE       ~ "p16-p21-"
  )
  obj$cdkn_state <- factor(obj$cdkn_state,
                           levels = c("p16-p21-", "p21+", "p16+", "p16+p21+"))
  
  # --- SnC label (same logic as your 4i w/ bcl2 + IP10|CXCR3) ---
  # obj$snc <- ifelse(((p21 | p16) & (!PCNA) & bcl2 & (IP10 |IL6 | IL8 )), "SnC+", "SnC-")
  # obj$snc <- factor(obj$snc, levels = c("SnC-", "SnC+"))
  
  obj$snc <- ifelse(((p21 | p16) & (!PCNA) & (IP10 |IL6 | IL8 )), "SnC+", "SnC-")
  obj$snc <- factor(obj$snc, levels = c("SnC-", "SnC+"))
  
  # --- Major/minor skin cell types ---
  obj$majorskincell <- dplyr::case_when(
    MelanA ~ "Melanocyte",
    SMA ~ "Smooth Muscle",
    Vimentin ~ "Fibroblast",
    (Keratin5 | Keratin14 | Filaggrin) ~ "Keratinocyte",
    TRUE ~ "none"
  )
  
  obj$minorskincell1 <- dplyr::case_when(
    MelanA ~ "Melanocyte",
    SMA ~ "Smooth Muscle",
    Vimentin ~ "Fibroblast",
    (Keratin5 | Keratin14) ~ "Basal Keratinocyte",
    Filaggrin ~ "Corneocyte",
    !Filaggrin & !Keratin5 & !Keratin14 & PanCyto ~ "Suprabasal Keratinocyte",
    TRUE ~ "none"
  )
  
  obj$majorskincell <- factor(obj$majorskincell,
                              levels = c("none","Smooth Muscle","Fibroblast","Melanocyte","Keratinocyte"))
  obj$minorskincell1 <- factor(obj$minorskincell1,
                              levels = c("none","Smooth Muscle","Fibroblast","Melanocyte",
                                         "Basal Keratinocyte", "Suprabasal Keratinocyte", "Corneocyte"))
  
  # --- Immune cell types (same as before) ---
  obj$immcelltype <- dplyr::case_when(
    
    CD4 & FOXP3 ~ "CD4+ FOXP3+",
    CD4 & CXCR3 ~ "CD4+ CXCR3+",
    CD14 ~ "CD14+",
    CD68 ~ "CD68+",
    CD8 ~ "CD8+",
    CD4 ~ "CD4+",
    CD45 ~ "CD45+",
    TRUE ~ "none"
  )
  
  # Grab the spatial matrix for this specific object
  spatial_matrix <- obj@reductions$spatial@cell.embeddings
  
  # Assign the first column (X) and second column (Y) to the metadata
  # (Adjust the [,1] and [,2] if your matrix has different dimensions)
  obj$Spatial_X <- spatial_matrix[, 1]
  obj$Spatial_Y <- spatial_matrix[, 2]
  
  # Ensure age order
  obj$age <- factor(obj$age, levels = c("Young", "Old"))
  
  # Create “replicate” id for plotting (if missing)
  if (!"source_file" %in% colnames(obj@meta.data)) obj$source_file <- sample_name
  
  obj$age_major <- paste(obj$age, obj$majorskincell, sep = " ")
  obj$age_minor <- paste(obj$age, obj$minorskincell1, sep = " ")
  
  seurat_list_codex[[sample_name]] <- obj
}

# Save immediately after the reproducible annotation stage. This avoids
# requiring all downstream exploratory plotting blocks to run successfully.
save(
  seurat_list_codex,
  file = "data/codex/seurat_list_codex_annotated.RData"
)

if (annotation_only) {
  cat("CODEX annotation-only run complete.\n")
  quit(save = "no", status = 0)
}


##########################################################
# 4) SnC proportion vs age (bar plot with sample dots)
##########################################################
library(dplyr)
library(ggplot2)
library(scales)
library(dplyr)
library(ggplot2)
library(scales)

blue_col <- "#296d9e"
snc_columns <- c("snc") # The two metadata columns you want to plot

# Loop through both "snc" (with BCL2) and "snc2" (without BCL2)
for (col_name in snc_columns) {
  
  # 1. Extract raw counts for the specific column
  df_codex <- bind_rows(lapply(names(seurat_list_codex), function(sample_name) {
    md <- seurat_list_codex[[sample_name]]@meta.data
    
    # Filter out NAs for the specific target column
    md_clean <- md[!is.na(md[[col_name]]), ]
    
    data.frame(
      sample = sample_name,
      age = unique(md_clean$age)[1],
      n_snc = sum(md_clean[[col_name]] == "SnC+", na.rm = TRUE),
      n_total = nrow(md_clean)
    )
  })) %>% 
    mutate(
      age = factor(age, levels = c("Young","Old")),
      snc_prop = n_snc / n_total  # Individual sample fraction
    )
  
  # 2. Calculate the pooled proportions and Standard Error
  df_pooled <- df_codex %>%
    group_by(age) %>%
    summarise(
      pooled_prop = sum(n_snc) / sum(n_total),            # True pooled fraction
      se = sd(snc_prop) / sqrt(n()),                      # SEM of the biological replicates
      .groups = "drop"
    )
  
  # 3. Dynamically set ymax and plot
  ymax <- max(c(df_codex$snc_prop, df_pooled$pooled_prop + df_pooled$se), na.rm = TRUE) * 1.15
  
  p_age <- ggplot() +
    geom_col(
      data = df_pooled, aes(x = age, y = pooled_prop),
      fill = blue_col, color = "black", width = 0.6
    ) +
    geom_errorbar(
      data = df_pooled, aes(x = age, ymin = pooled_prop - se, ymax = pooled_prop + se),
      width = 0.15, linewidth = 0.8, color = "black"
    ) +
    geom_point(
      data = df_codex, aes(x = age, y = snc_prop),
      position = position_jitter(width = 0.07, height = 0, seed = 42),
      size = 3, color = "black"
    ) +
    scale_y_continuous(
      limits = c(0, ymax),
      breaks = pretty(c(0, ymax), n = 4),
      labels = percent_format(accuracy = 0.1), # <--- Change this from 1 to 0.1
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 18)
  
  # Print to console/viewer
  print(p_age)
  
  # 4. Dynamically name and save the file based on the column used
  file_name <- paste0("results/260303/codex_", col_name, "_age.pdf")
  ggsave(filename = file_name, plot = p_age, width = 3.6, height = 4.2, dpi = 500)
  
  cat("Successfully saved:", file_name, "\n")
}



##########################################################
# 5) Sub-senescence marker violins by SnC (regular violins like 4i)
##########################################################
plot_violin_marker <- function(obj,
                               feature,
                               group_col = "snc2",
                               assay = "RNA",
                               layer = "data",
                               cols = c("SnC-" = "#296d9e", "SnC+" = "#cc0b01")) {
  
  df <- FetchData(obj, vars = c(feature, group_col), assay = assay, layer = layer)
  colnames(df) <- c("value", "group")
  
  df <- df %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::mutate(group = factor(group, levels = c("SnC-", "SnC+")))
  
  ymax <- max(df$value, na.rm = TRUE)
  
  ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.9) +
    stat_summary(fun = mean, geom = "crossbar",
                 width = 0.4, color = "black", linewidth = 0.8) +
    stat_compare_means(method = "wilcox.test", label = "p.signif",
                       label.x = 1.5, label.y = ymax * 0.98, size = 6) +
    scale_fill_manual(values = cols) +
    theme_classic(base_size = 18) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    labs(y = "Expression", title = toupper(feature))
}

obj_codex <- merge(seurat_list_codex[[1]], y = seurat_list_codex[-1], merge.data = TRUE)

features <- c("HMGB1","LaminB1", "BCL2")
p_list <- lapply(features, function(f) plot_violin_marker(obj_codex, f, group_col="snc2", assay="RNA", layer="data"))
final_plot <- p_list[[1]] + p_list[[2]] + p_list[[3]]
final_plot

ggsave("results/260506/codex_violins_subsenmarkers(nobcl2).pdf", plot = final_plot, width = 6, height = 3, dpi = 300)

##########################################################
# SnC Marker Distribution
##########################################################
VlnPlot(
  obj_codex, 
  features = c("p21", "p16", "PCNA", "BCL2", "IP10", "IL6", "IL8"), 
  group.by = "snc", 
  pt.size = 0 # Turn off points if you have thousands of cells
)


##########################################################
# 6) SnC frequency by cell type (sample-level) + bar facets
##########################################################
major_anatomical_order <- c("none","Smooth Muscle","Fibroblast","Melanocyte","Keratinocyte")
minor_anatomical_order <- c("none","Smooth Muscle","Fibroblast","Melanocyte","Basal Keratinocyte","Suprabasal Keratinocyte","Corneocyte")

freq_major <- obj_codex@meta.data %>%
  filter(!is.na(age) & !is.na(majorskincell) & !is.na(snc) & !is.na(source_file)) %>%
  group_by(source_file, majorskincell, age) %>%
  summarise(
    Total_Cells = n(),
    SnC_Plus_Cells = sum(snc == "SnC+"),
    Freq_SnC_Plus = 100 * SnC_Plus_Cells / Total_Cells,
    .groups = "drop"
  ) %>%
  mutate(majorskincell = factor(majorskincell, levels = major_anatomical_order),
         age = factor(age, levels = c("Young","Old")))

freq_minor <- obj_codex@meta.data %>%
  filter(!is.na(age) & !is.na(minorskincell1) & !is.na(snc) & !is.na(source_file)) %>%
  group_by(source_file, minorskincell1, age) %>%
  summarise(
    Total_Cells = n(),
    SnC_Plus_Cells = sum(snc == "SnC+"),
    Freq_SnC_Plus = 100 * SnC_Plus_Cells / Total_Cells,
    .groups = "drop"
  ) %>%
  mutate(minorskincell1 = factor(minorskincell1, levels = minor_anatomical_order),
         age = factor(age, levels = c("Young","Old"))) %>%
  filter(!is.na(minorskincell1))

p_major_bar <- ggplot(freq_major, aes(x = age, y = Freq_SnC_Plus)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6,
               fill = "#296d9e", color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15,
               linewidth = 0.8, color = "black") +
  geom_point(position = position_jitter(width = 0.07, height = 0),
             size = 3, color = "black") +
  facet_wrap(~ majorskincell, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "% SnC+ Cells") +
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  )

p_minor_bar <- ggplot(freq_minor, aes(x = age, y = Freq_SnC_Plus)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6,
               fill = "#296d9e", color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15,
               linewidth = 0.8, color = "black") +
  geom_point(position = position_jitter(width = 0.07, height = 0),
             size = 3, color = "black") +
  facet_wrap(~ minorskincell1, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "% SnC+ Cells") +
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  )

ggsave("results/260303/codex_barplot_snc_bcl2_major.pdf", plot = p_major_bar, width = 3, height = 11, dpi = 300)
ggsave("results/260303/codex_barplot_snc_bcl2_minor.pdf", plot = p_minor_bar, width = 3, height = 12, dpi = 300)

##########################################################
# 7) Donut: SnC composition by (minor) cell type (Young vs Old)
##########################################################

minor_colors <- c(
  "Basal Keratinocyte" = "#1f78b4",
  "Fibroblast" = "#33a02c",
  "Melanocyte" = "#6a3d9a",
  "Corneocyte" = "#fb9a99",
  "Smooth Muscle" = "#ff7f00",
  "none" = "grey80"
)

snc_comp <- obj_codex@meta.data %>%
  filter(snc == "SnC+" & !is.na(age) & !is.na(minorskincell1)) %>%
  group_by(age, minorskincell1) %>%
  summarise(
    cells = n(),
    .groups = "drop"
  ) %>%
  group_by(age) %>%
  mutate(
    pct = 100 * cells / sum(cells)
  )

snc_comp <- snc_comp %>%
  mutate(age = factor(age, levels = c("Young", "Old")))

donut_one <- function(df, age_value) {
  ggplot(df %>% filter(age == age_value),
         aes(x = 2, y = pct, fill = minorskincell1)) +
    geom_col(width = 1, color = "white", linewidth = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    scale_fill_manual(values = minor_colors, drop = FALSE) +
    theme_void() +
    ggtitle(as.character(age_value)) +
    theme(
      text = element_text(family = "Helvetica", size = 14),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}



p_young <- donut_one(snc_comp, "Young")
p_old   <- donut_one(snc_comp, "Old")

p_donut_minor <- (p_young | p_old) + plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("results/260617/codex_snc_donut_minor.pdf", p_donut_minor, width = 7.5, height = 3.5)
ggsave("results/260303/codex_snc_donut_major.pdf", p_donut_major, width = 7.5, height = 3.5)

##########################################################
# 8) Donut: p16/p21 distribution within SnC (Young vs Old) with ONE legend
##########################################################
cdkn_colors <- c(
  "p16-p21-" = "#d9d9d9",
  "p21+" = "#296d9e",
  "p16+" = "#cc0b01",
  "p16+p21+" = "#7f7f7f"
)

snc_marker_comp <- obj_codex@meta.data %>%
  filter(snc == "SnC+") %>%
  count(age, cdkn_state) %>%
  group_by(age) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    age = factor(age, levels = c("Young","Old")),
    cdkn_state = factor(cdkn_state, levels = c("p16-p21-", "p21+", "p16+", "p16+p21+"))
  )

donut_cdkn_noscale <- function(df, age_value) {
  ggplot(df %>% dplyr::filter(age == age_value),
         aes(x = 2, y = pct, fill = cdkn_state)) +
    geom_col(width = 1, color = "white", linewidth = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    theme_void() +
    ggtitle(as.character(age_value)) +
    theme(
      text = element_text(family = "Helvetica", size = 14),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.margin = margin(10, 10, 10, 10)
    )
}

p_young <- donut_cdkn_noscale(snc_marker_comp, "Young") + theme(legend.position = "none")
p_old   <- donut_cdkn_noscale(snc_marker_comp, "Old")

p_cdkn <- (p_young | p_old) +
  plot_layout(guides = "collect") &
  scale_fill_manual(values = cdkn_colors, drop = TRUE) &   # <- key change
  guides(fill = guide_legend(title = "cdkn_state")) &
  theme(legend.position = "right")

p_cdkn

ggsave("results/260303/codex_snc_donut_cdkn.pdf", p_cdkn, width = 7.5, height = 3.5)

#########################
# SnC+ cell type distributions in young and old
#########################

library(ggplot2)
library(dplyr)
library(scales)

# Define your colors
major_colors <- c(
  "none" = "grey80",
  "Fibroblast" = "#33a02c",
  "Melanocyte" = "#6a3d9a",
  "Keratinocyte" = "#1f78b4"
)

minor_colors <- c(
  "Suprabasal Keratinocyte" = "yellow",
  "Basal Keratinocyte" = "#1f78b4",
  "Fibroblast" = "#33a02c",
  "Melanocyte" = "#6a3d9a",
  "Corneocyte" = "#fb9a99",
  "Smooth Muscle" = "#ff7f00",
  "none" = "grey80"
)
# --- MAJOR CELL TYPES (Out of Total Tissue) ---
df_major_total <- obj_codex@meta.data %>%
  filter(!is.na(age) & !is.na(majorskincell) & !is.na(snc)) %>%
  # 1. Group by age to get the true denominator (all cells)
  group_by(age) %>%
  mutate(total_cells_in_age = n()) %>%
  # 2. Filter down to only the SnC+ cells
  filter(snc2 == "SnC+") %>%
  # 3. Count how many SnC+ cells belong to each cell type
  group_by(age, majorskincell, total_cells_in_age) %>%
  summarise(snc_count = n(), .groups = "drop") %>%
  # 4. Divide by the total tissue cell count
  mutate(pct_of_total = snc_count / total_cells_in_age) %>%
  mutate(age = factor(age, levels = c("Young", "Old")))

p_major_total <- ggplot(df_major_total, aes(x = age, y = pct_of_total, fill = majorskincell)) +
  geom_col(width = 0.6, color = "black") +
  # Use percent_format to show 1%, 2%, 3% etc. on the y-axis
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_fill_manual(values = major_colors) +
  theme_classic(base_size = 16) +
  labs(
    title = "SnC+ Cells as % of Total Tissue (Major)", 
    x = "Age", 
    y = "% of Total Cells", 
    fill = "Cell Type"
  )

# --- MINOR CELL TYPES (Out of Total Tissue) ---
df_minor_total <- obj_codex@meta.data %>%
  filter(!is.na(age) & !is.na(minorskincell1) & !is.na(snc)) %>%
  group_by(age) %>%
  mutate(total_cells_in_age = n()) %>%
  filter(snc2 == "SnC+") %>%
  group_by(age, minorskincell1, total_cells_in_age) %>%
  summarise(snc_count = n(), .groups = "drop") %>%
  mutate(pct_of_total = snc_count / total_cells_in_age) %>%
  mutate(age = factor(age, levels = c("Young", "Old")))

p_minor_total <- ggplot(df_minor_total, aes(x = age, y = pct_of_total, fill = minorskincell1)) +
  geom_col(width = 0.6, color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_fill_manual(values = minor_colors) +
  theme_classic(base_size = 16) +
  labs(
    title = "SnC+ Cells as % of Total Tissue (Minor)", 
    x = "Age", 
    y = "% of Total Cells", 
    fill = "Cell State"
  )

# Save the plots
ggsave("results/260303/codex_snc_total_tissue_major.pdf", p_major_total, width = 5, height = 5)
ggsave("results/260303/codex_snc_total_tissue_minor(nobcl).pdf", p_minor_total, width = 6, height = 5)

#########################
# SnC+ cell type distributions in young and old (normalized by total SnC)
#########################

major_colors <- c(
  "none" = "grey80",
  "Fibroblast" = "#33a02c",
  "Melanocyte" = "#6a3d9a",
  "Keratinocyte" = "#1f78b4",
  "Smooth Muscle" = "#ff7f00" 
)

minor_colors <- c(
  "Basal Keratinocyte" = "#1f78b4",
  "Fibroblast" = "#33a02c",
  "Melanocyte" = "#6a3d9a",
  "Corneocyte" = "#fb9a99",
  "Smooth Muscle" = "#ff7f00",
  "none" = "grey80"
)

RdBu color scheme

# --- MAJOR CELL TYPES (Composition out of 100% SnC+ cells) ---
df_major_comp <- obj_codex@meta.data %>%
  filter(snc2 == "SnC+" & !is.na(age) & !is.na(majorskincell)) %>%
  count(age, majorskincell) %>%
  group_by(age) %>%
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  mutate(age = factor(age, levels = c("Young", "Old")))

p_major_comp <- ggplot(df_major_comp, aes(x = age, y = pct, fill = majorskincell)) +
  geom_col(width = 0.6, color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = major_colors) +
  theme_classic(base_size = 16) +
  labs(
    # title = "SnC+ Cell Type Composition (Major)", 
    x = "Age", 
    y = "Proportion of SnC+ Cells", 
    fill = "Cell Type"
  )
library(dplyr)
library(ggplot2)
library(scales)

# ------------------------------------------------------------------
# 1. MASTER COLOR PALETTE & STACK ORDERS
# ------------------------------------------------------------------
# This elegant palette covers both major and minor cell categories
elegant_colors <- c(
  # Structural/Keratinocyte family (Consistent Reds/Oranges)
  "Keratinocyte"            = "#67001f",
  "Basal Keratinocyte"      = "#b2182b",
  "Suprabasal Keratinocyte" = "#2166ac",
  "Corneocyte"              = "#fddbc7",
  
  # Structural/Stromal (Consistent Blues)
  "Fibroblast"              = "#053061",
  "Smooth Muscle"           = "#4393c3",
  "Melanocyte"              = "#d1e5f0",
  
  # Neutral/None
  "none"                    = "#f7f7f7" # Light grey for none
)



# Lock in the vertical stack orders so legends and bars match
major_order <- c("Keratinocyte", "Melanocyte", "Fibroblast", "Smooth Muscle", "none")
minor_order <- c(
  "Suprabasal Keratinocyte", "Basal Keratinocyte", "Corneocyte", 
  "Melanocyte", "Fibroblast", "Smooth Muscle", "none"
)

# ------------------------------------------------------------------
# 2. MAJOR CELL TYPES PLOT
# ------------------------------------------------------------------
df_major_total <- obj_codex@meta.data %>%
  filter(!is.na(age) & !is.na(majorskincell) & !is.na(snc)) %>%
  group_by(age) %>%
  mutate(total_cells_in_age = n()) %>%
  filter(snc2 == "SnC+") %>%
  group_by(age, majorskincell, total_cells_in_age) %>%
  summarise(snc_count = n(), .groups = "drop") %>%
  mutate(pct_of_total = snc_count / total_cells_in_age) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    majorskincell = factor(majorskincell, levels = major_order)
  )

p_major <- ggplot(df_major_total, aes(x = age, y = pct_of_total, fill = majorskincell)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = elegant_colors) +
  theme_minimal(base_size = 15, base_family = "sans") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", margin = margin(t = 5)),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    legend.key = element_blank(),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 13),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  labs(x = NULL, y = NULL, fill = "Major Cell Type")

print(p_major)

# ------------------------------------------------------------------
# 3. MINOR CELL TYPES PLOT
# ------------------------------------------------------------------
df_minor_total <- obj_codex@meta.data %>%
  filter(!is.na(age) & !is.na(minorskincell1) & !is.na(snc)) %>%
  group_by(age) %>%
  mutate(total_cells_in_age = n()) %>%
  filter(snc2 == "SnC+") %>%
  group_by(age, minorskincell1, total_cells_in_age) %>%
  summarise(snc_count = n(), .groups = "drop") %>%
  mutate(pct_of_total = snc_count / total_cells_in_age) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    minorskincell1 = factor(minorskincell1, levels = minor_order)
  )

p_minor <- ggplot(df_minor_total, aes(x = age, y = pct_of_total, fill = minorskincell1)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.7) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = elegant_colors) +
  theme_minimal(base_size = 15, base_family = "sans") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", margin = margin(t = 5)),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    legend.key = element_blank(),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 13),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  labs(x = NULL, y = NULL, fill = "Minor Cell State")

print(p_minor)



# Save the plots
ggsave("results/260506/codex_snc_composition_stacked_major.pdf", plot = p_major, width = 4, height = 5)
ggsave("results/260506/codex_snc_composition_stacked_minor.pdf", plot = p_minor, width = 6, height = 10)




major_order <- c("Fibroblast", "Melanocyte", "Keratinocyte", "Smooth Muscle", "none")
minor_order <- c(
  "Fibroblast",
  "Melanocyte",
  "Basal Keratinocyte",
  "Suprabasal Keratinocyte",
  "Corneocyte",
  "Smooth Muscle",
  "none"
)

age_colors <- c("Young" = "#9ecae1", "Old" = "#08519c")

# ------------------------------------------------------------------
# 2. MAJOR CELL TYPES PLOT (CODEX)
#    SnC+ burden out of total cells in each age group
# ------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(scales)
library(grid)

# define which major labels count as skin cells
major_skin_levels <- c("Fibroblast", "Melanocyte", "Keratinocyte", "Smooth Muscle")

# denominator = all skin cells per age group
denom_major_skin <- obj_codex@meta.data %>%
  filter(!is.na(age), !is.na(majorskincell)) %>%
  filter(majorskincell %in% major_skin_levels) %>%
  group_by(age) %>%
  summarise(total_skin_cells_in_age = n(), .groups = "drop")

# numerator = SnC+ cells per major cell type per age group
df_major_total_codex <- obj_codex@meta.data %>%
  filter(!is.na(age), !is.na(majorskincell), !is.na(snc2)) %>%
  filter(majorskincell %in% major_skin_levels) %>%
  filter(snc2 == "SnC+") %>%
  group_by(age, majorskincell) %>%
  summarise(snc_count = n(), .groups = "drop") %>%
  left_join(denom_major_skin, by = "age") %>%
  mutate(
    pct_of_total = snc_count / total_skin_cells_in_age,
    age = factor(age, levels = c("Young", "Old")),
    majorskincell = factor(majorskincell, levels = major_order)
  )

p_major_total_codex <- ggplot(df_major_total_codex, aes(x = age, y = pct_of_total, fill = majorskincell)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.7) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = elegant_colors, drop = FALSE) +
  theme_minimal(base_size = 15, base_family = "sans") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", margin = margin(t = 5)),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    legend.key = element_blank(),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 13),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  labs(x = NULL, y = NULL, fill = "Major Cell Type")

print(p_major_total_codex)

# ------------------------------------------------------------------
# 3. MINOR CELL TYPES PLOT (CODEX)
#    SnC+ burden out of total cells in each age group
# ------------------------------------------------------------------
# define which minor labels count as skin cells
minor_skin_levels <- c(
  "Fibroblast",
  "Melanocyte",
  "Basal Keratinocyte",
  "Suprabasal Keratinocyte",
  "Corneocyte",
  "Smooth Muscle"
)

# denominator = all skin cells per age group
denom_minor_skin <- obj_codex@meta.data %>%
  filter(!is.na(age), !is.na(minorskincell1)) %>%
  filter(minorskincell1 %in% minor_skin_levels) %>%
  group_by(age) %>%
  summarise(total_skin_cells_in_age = n(), .groups = "drop")

# numerator = SnC+ cells per minor cell type per age group
df_minor_total_codex <- obj_codex@meta.data %>%
  filter(!is.na(age), !is.na(minorskincell1), !is.na(snc2)) %>%
  filter(minorskincell1 %in% minor_skin_levels) %>%
  filter(snc2 == "SnC+") %>%
  group_by(age, minorskincell1) %>%
  summarise(snc_count = n(), .groups = "drop") %>%
  left_join(denom_minor_skin, by = "age") %>%
  mutate(
    pct_of_total = snc_count / total_skin_cells_in_age,
    age = factor(age, levels = c("Young", "Old")),
    minorskincell1 = factor(minorskincell1, levels = minor_order)
  )

p_minor_total_codex <- ggplot(df_minor_total_codex, aes(x = age, y = pct_of_total, fill = minorskincell1)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.7) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = elegant_colors, drop = FALSE) +
  theme_minimal(base_size = 15, base_family = "sans") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", margin = margin(t = 5)),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    legend.key = element_blank(),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 13),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  labs(x = NULL, y = NULL, fill = "Minor Cell State")

print(p_minor_total_codex)

ggsave("results/260303/codex_snc_skin_stacked_major.pdf", plot = p_major_total_codex, width = 5, height = 5, dpi = 300)
ggsave("results/260303/codex_snc_skin_stacked_minor.pdf", plot = p_minor_total_codex, width = 6, height = 5, dpi = 300)

# ------------------------------------------------------------------
# 4. MAJOR CELL TYPES: Lineage-Specific SnC+ Burden (CODEX)
#    % SnC+ within each lineage
# ------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(scales)
library(grid)

# define which major labels count as skin cells
major_skin_levels <- c("Fibroblast", "Melanocyte", "Keratinocyte", "Smooth Muscle")

# denominator = all skin cells per age group
denom_major_skin <- obj_codex@meta.data %>%
  filter(!is.na(age), !is.na(majorskincell)) %>%
  filter(majorskincell %in% major_skin_levels) %>%
  group_by(age) %>%
  summarise(total_skin_cells_in_age = n(), .groups = "drop")

# numerator = SnC+ cells per major cell type per age group
df_major_total_codex <- obj_codex@meta.data %>%
  filter(!is.na(age), !is.na(majorskincell), !is.na(snc2)) %>%
  filter(majorskincell %in% major_skin_levels) %>%
  filter(snc2 == "SnC+") %>%
  group_by(age, majorskincell) %>%
  summarise(snc_count = n(), .groups = "drop") %>%
  left_join(denom_major_skin, by = "age") %>%
  mutate(
    pct_of_total = snc_count / total_skin_cells_in_age,
    age = factor(age, levels = c("Young", "Old")),
    majorskincell = factor(majorskincell, levels = major_order)
  )

p_major_total_codex <- ggplot(df_major_total_codex, aes(x = age, y = pct_of_total, fill = majorskincell)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.7) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(values = elegant_colors, drop = FALSE) +
  theme_minimal(base_size = 15, base_family = "sans") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", margin = margin(t = 5)),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    legend.key = element_blank(),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 13),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  labs(x = NULL, y = NULL, fill = "Major Cell Type")

print(p_major_total_codex)

print(p_major_intrinsic_codex)

# ------------------------------------------------------------------
# 5. MINOR CELL TYPES: Lineage-Specific SnC+ Burden (CODEX)
#    % SnC+ within each lineage
# ------------------------------------------------------------------
df_minor_intrinsic_codex <- obj_codex@meta.data %>%
  filter(!is.na(age) & !is.na(minorskincell1) & !is.na(snc2)) %>%
  filter(minorskincell1 != "none") %>%
  group_by(minorskincell1, age) %>%
  summarise(
    total_lineage_cells = n(),
    snc_cells = sum(snc2 == "SnC+"),
    pct_snc = snc_cells / total_lineage_cells,
    .groups = "drop"
  ) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    minorskincell1 = factor(
      minorskincell1,
      levels = c("Fibroblast", "Melanocyte", "Basal Keratinocyte", "Suprabasal Keratinocyte", "Corneocyte", "Smooth Muscle")
    )
  )

p_minor_intrinsic_codex <- ggplot(df_minor_intrinsic_codex, aes(x = minorskincell1, y = pct_snc, fill = age)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_manual(values = age_colors) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(
    axis.text.x = element_text(color = "black", face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = NULL
  )

print(p_minor_intrinsic_codex)

ggsave("results/260303/codex_snc_intrinsic_major.pdf", plot = p_major_intrinsic_codex, width = 6, height = 5, dpi = 300)
ggsave("results/260303/codex_snc_intrinsic_minor.pdf", plot = p_minor_intrinsic_codex, width = 7, height = 5, dpi = 300)



####################################################################
# 10) Spatial Tissue Mapping by Sample (`cdkn_state`)
####################################################################
library(ggplot2)
library(dplyr)
library(patchwork)

# Ensure you use your exact matching color palette configuration
cdkn_colors <- c(
  "p16-p21-" = "#d3d3d3",
  "p21+"     = "#010087",
  "p16+"     = "#8d0000",
  "p16+p21+" = "#006900"
)

# Container list to store each separate specimen plot
for (sample_name in names(seurat_list_codex)) {
  cat("Exporting individual spatial map for:", sample_name, "\n")
  
  # Extract the metadata for just this specific sample slice
  meta_df <- seurat_list_codex[[sample_name]]@meta.data
  
  p_individual <- ggplot(meta_df, aes(x = Spatial_X, y = Spatial_Y, color = cdkn_state)) +
    # shape = 16 (solid circle) with a slight alpha looks very crisp for tissue dots
    geom_point(size = 0.5, alpha = 0.8, shape = 16) +
    scale_color_manual(values = cdkn_colors, drop = FALSE, name = "CDKN State") +
    
    # CRITICAL: This freezes the true micro-meter physical aspect ratio of the skin biopsy
    coord_fixed() + 
    
    # Clean up standard plotting elements to make it look like a slide scan
    theme_void(base_size = 14) + 
    labs(title = sample_name) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 10)),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 10)
    )
  
  # Generate a clean, standardized filename for each specimen
  # (replaces spaces/special characters in sample name with underscores just in case)
  clean_filename <- gsub("[^A-Za-z0-9_-]", "_", sample_name)
  output_path <- paste0("results/260506/spatial_", clean_filename, ".pdf")
  
  # Save the individual plot
  ggsave(
    filename = output_path,
    plot = p_individual,
    width = 6,   # Perfect square-ish aspect ratio that lets coord_fixed work beautifully
    height = 5,
    dpi = 300
  )
}

cat("Spatial tissue graphics exported successfully to target folders.\n")


####################################################################
# 11) Polished Dot Plot: Size = Prevalence, Color = Expression Strength
####################################################################
library(Seurat)
library(dplyr)
library(ggplot2)

# ----------------------------------------------------------------------
# 1) Subset Object to Keep "Old" and Screen Out Background
# ----------------------------------------------------------------------
# Isolate your target subset to avoid global calculating shifts
obj_old_pos <- subset(
  x = obj_joined,
  subset = age == "Old" & 
    !is.na(minorskincell1) & 
    !is.na(cdkn_state)
)

# ----------------------------------------------------------------------
# 2) Enforce Your Anatomical Tissue Stack Ordering on Y-Axis
# ----------------------------------------------------------------------
custom_axis_order <- c(
  "none",                    # Absolute bottom row
  "Smooth Muscle",           
  "Fibroblast",              
  "Melanocyte",              
  "Basal Keratinocyte",      
  "Suprabasal Keratinocyte", 
  "Corneocyte"               # Absolute top row
)

# Filter identity factors down and enforce custom level sequencing
obj_old_pos$minorskincell1 <- factor(
  obj_old_pos$minorskincell1, 
  levels = custom_axis_order
)

# Set the active identity column for Seurat's DotPlot plotting engine
Idents(obj_old_pos) <- "minorskincell1"


# ----------------------------------------------------------------------
# 3) Generate Independent Seurat DotPlots to Decompress Color Ranges
# ----------------------------------------------------------------------
p_base <- DotPlot(
  object = obj_old_pos,
  features = c("p16", "p21"),
  assay = "RNA",
  col.min = NULL, 
  col.max = NULL  
) + scale_size_continuous(
  range = c(4, 14),          # Increases minimum and maximum dot circle diameters
  labels = function(x) paste0(round(x, 1), "%"),
  name = "Percent\nExpressing"
) + 
  scale_color_gradientn(
    colors = c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b")) 


print(p_base)
# ----------------------------------------------------------------------
# 4) Save Clean Clean Layout Out to Your Direct Path
# ----------------------------------------------------------------------
ggsave(
  filename = "results/260523/codex_old_skin_cdkn_direct_seurat_dotplot.pdf",
  plot = p_base,
  width = 6.2,
  height = 4.5,
  dpi = 300
)


####################################################################
# Immune Cell Proportions
####################################################################


# 1. FetchData safely queries across all split Seurat v5 layers automatically
marker_data <- FetchData(
  object = obj_joined, 
  vars = c("CD45", "CD8", "CD4", "CD68", "CD11c", "CD14", "FOXP3", "CXCR3", "Ki67"),
  layer = "counts"
)

# 2. Extract logical boolean vectors from the cleanly compiled dataframe
CD45  <- marker_data[["CD45"]] > 0
CD8   <- marker_data[["CD8"]] > 0
CD4   <- marker_data[["CD4"]] > 0
CD68  <- marker_data[["CD68"]] > 0
CD14  <- marker_data[["CD14"]] > 0
FOXP3 <- marker_data[["FOXP3"]] > 0
CXCR3 <- marker_data[["CXCR3"]] > 0
ki67  <- marker_data[["Ki67"]] > 0
CD11c  <- marker_data[["CD11c"]] > 0

# 2. Build metadata column with strict priority-ordered logic
obj_joined@meta.data$immcelltype4 <- case_when(
  # Only evaluate if cell is CD45+ (optional, remove condition if you want loose gating)
  
  # Combinatorial/Specific States First
  CD4 & ki67   ~ "CD4+ Ki67+",
  CD4 & FOXP3  ~ "CD4+ FOXP3+",
  CD4 & !FOXP3 ~ "CD4+ FOXP3-",
  
  # Single Marker/Lineage Rest-of-Pool
  CD8          ~ "CD8+",
  CD4          ~ "CD4+",
  CD14         ~ "CD14+",
  CD68         ~ "CD68+",
  
  TRUE         ~ "none"
)

target_order <- c(
  "CD4+ FOXP3+",
  "CD4+ CXCR3+",
  "CD4+",
  "CD8+",
  "CD68+",
  "CD14+",
  "CD11c+"
)

# 2. Blue and Red Scheme Custom Palette mapping your exact string names
# Maps T cells to vibrant/deep reds and Myeloid/MHC groups to steel/ice blues
elegant_colors <- c(
  # --- GROUP 1: CD4+ T Cells (Anchored by your Deep Velvet Red) ---
  "CD4+ FOXP3+"  = "#8d0000",  # Exact requested anchor Deep Velvet Red
  "CD4+ CXCR3+"  = "#c96464",  # Harmonized mid-tone rose terracotta
  "CD4+"         = "#e8b0b0",  # Harmonized pale pastel rose
  
  # --- GROUP 2: CD8+, Myeloid & Core Lineages (Anchored by your Ultramarine Blue) ---
  "CD8+"         = "#4263EB",  # Exact requested anchor Ultramarine Blue
  "CD68+"        = "#437cb4",  # Distinct light slate/steel blue
  "CD14+"        = "#9bc3e6",  # Distinct crisp light sky blue
  "CD11c+"        = "#768594"   # Muted mineral indigo-grey (clean transitional boundary)
)



obj_joined@meta.data$immcelltype5 <- case_when(
  # Only evaluate if cell is CD45+ (optional, remove condition if you want loose gating)
  
  CD8          ~ "CD8+",
  CD4          ~ "CD4+",
  CD14         ~ "CD14+",
  CD68         ~ "CD68+",
  
  TRUE         ~ "none"
)

target_order1 <- c(
  "CD4+",
  "CD8+",
  "CD68+",
  "CD14+"
)

# 2. Blue and Red Scheme Custom Palette mapping your exact string names
# Maps T cells to vibrant/deep reds and Myeloid/MHC groups to steel/ice blues
elegant_colors1 <- c(
  "CD4+"         = "#c96464",
  
  # --- GROUP 2: CD8+, Myeloid & Core Lineages (Anchored by your Ultramarine Blue) ---
  "CD8+"         = "#e8b0b0",  # Exact requested anchor Ultramarine Blue
  "CD68+"        = "#437cb4",  # Distinct light slate/steel blue
  "CD14+"        = "#9bc3e6"  # Distinct crisp light sky blue
)

obj_joined@meta.data$immcelltype2 <- case_when(
  
  CD4  & FOXP3 ~ "CD4+ FOXP3+",
  CD4  & CXCR3 ~ "CD4+ CXCR3+",
  CD8 ~ "CD8+",
  CD14 ~ "CD14+",
  CD11c ~ "CD11c+",
  #CD45 ~ "CD45+",
  CD68 ~ "CD68+",
  CD4 ~ "CD4+",
  TRUE ~ 'none'
)

# 3. Process data out of ALL tissue cells (denominator captures background "none" cells)
plot_df <- obj_joined@meta.data %>%
  filter(!is.na(age)) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    immcelltype2 = as.character(immcelltype2)
  ) %>%
  count(age, immcelltype2) %>%           # Generates a column named 'n' automatically!
  group_by(age) %>%
  mutate(percentage = n / sum(n)) %>%   # Denominator includes background structural cells!
  ungroup() %>%
  filter(immcelltype2 != "none") %>%     # Drop background structural cells before tracking proportions
  mutate(
    immcelltype2 = factor(immcelltype2, levels = target_order),
    
    # Generate labels. Suppresses micro-fractions under 0.1% to prevent text overlaps
    bar_label = ifelse(percentage >= 0.001, percent(percentage, accuracy = 0.1), "")
  )

# 4. Generate the exact PDF canvas matching your beautiful custom style template
pdf("results/codex_snc_composition_stacked_immune2.pdf", width = 6.5, height = 7)

ggplot(plot_df, aes(x = age, y = percentage, fill = immcelltype2)) +
  # Mirrored your custom bar formatting geometry exactly
  geom_col(width = 0.6, color = "white", linewidth = 0.7) +
  
  # Inject percentage text strings cleanly centered inside each stacked block layer
  geom_text(
    aes(label = bar_label),
    position = position_stack(vjust = 0.5),
    color = "black",
    size = 3
  ) +
  
  # Mirrored your y-scale expansion parameter rules exactly
  scale_y_continuous(labels = percent_format(accuracy = 0.1), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = elegant_colors, breaks = target_order) + 
  
  # Mirrored your theme architecture and background dashed grid formatting exactly
  theme_minimal(base_size = 15, base_family = "sans") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.line.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", margin = margin(t = 5)),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    legend.key = element_blank(),
    legend.key.size = unit(0.6, "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 13),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +
  labs(x = NULL, y = NULL, fill = "Immune Subset")

dev.off()

####################################################################
# Immune Cell Proportions
####################################################################

# obj_joined@meta.data$immcelltype3 <- case_when(
# 
#   # Combinatorial/Specific States First
#   CD4 & FOXP3  ~ "CD4+ FOXP3+",
#   CD4 & CXCR3 ~ "CD4+ CXCR3+",
#   
#   # Single Marker/Lineage Rest-of-Pool
#   CD8          ~ "CD8+",
#   CD4          ~ "CD4+",
#   CD4          ~ "CD11c+",
#   CD14         ~ "CD14+",
#   CD68         ~ "CD68+",
#   CD11c          ~ "CD11c+",
#   
#   TRUE         ~ "none"
# )


obj_joined@meta.data$immcelltype2 <- case_when(
  
  CD4  & FOXP3 ~ "CD4+ FOXP3+",
  CD4  & CXCR3 ~ "CD4+ CXCR3+",
  CD8 ~ "CD8+",
  CD14 ~ "CD14+",
  CD11c ~ "CD11c+",
  CD45 ~ "CD45+",
  CD68 ~ "CD68+",
  CD4 ~ "CD4+",
  TRUE ~ 'none'
)


mycols <- c(
  "SnC-" = "#bebebe",  # Exact Slate Charcoal Gray
  "SnC+" = "#3156a6"   # Exact Classic Steel Blue
)

age_groups <- c('Old', 'Young')

for (current_age in age_groups) {
  
  cat(paste("Processing data for age group:", current_age, "\n"))
  
  # 1. Subset the Seurat object for the current age group
  seurat_subset <- subset(obj_joined, subset = age == current_age)
  
  # 2. Calculate the cell type vs. native snc contingency table directly
  immsen <- as.data.frame(table(seurat_subset$immcelltype2, seurat_subset$snc))
  print(immsen)
  # 3. Calculate relative frequency (proportion)
  your_data <- immsen %>%
    group_by(Var1) %>%
    # Calculate frequency as proportion of the total for that cell type
    mutate(Freq = Freq / sum(Freq)) %>%
    ungroup()
  
  # 4. Clean up columns and filter 'none'
  colnames(your_data) <- c('cell_types', 'SnC', 'Freq')
  your_data <- your_data[your_data$cell_types != "none", ]
  
  # 5. Enforce specific immune cell subset sequencing layout
  your_data$cell_types <- factor(your_data$cell_types, levels = target_order)
  
  # 6. Calculate percentages and text positioning
  your_data$Freq <- round(your_data$Freq * 100, 2)
  
  your_data <- your_data %>%
    arrange(cell_types, desc(SnC)) %>% # Arrange correctly for stacking
    group_by(cell_types) %>%
    # Calculate y-position for text labels dead-center inside the segments
    mutate(lab.ypos = cumsum(Freq) - 0.5 * Freq) %>%
    ungroup()
  
  # 7. Define the vector PDF output file name path
  output_filename <- paste0("//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2/results/codex_imm22_", tolower(current_age), "_immsnc.pdf")
  
  # 8. Start vector graphics device (PDF) with your aspect ratio dimension scaling
  pdf(output_filename, width = 30, height = 8)
  
  # 9. Generate the faceted layout matching your template architecture
  p <- ggplot(your_data, aes(x = "", y = Freq, group = SnC, color = SnC, fill = SnC)) +
    geom_bar(width = 1, stat = "identity", color = "white", linewidth = 0.3) + 
    coord_polar("y", start = 0) + # Creates the pie effect
    geom_text(aes(y = lab.ypos, label = paste0(Freq, "%")), color = "black", size = 5, fontface = "bold") +
    scale_fill_manual(values = mycols) +
    scale_color_manual(values = mycols) +
    facet_wrap(~ cell_types, ncol = 10) + # Dynamic side-by-side positioning grid
    labs(
      fill = "Senescent Status"
    ) +
    theme_void() +
    theme(
      strip.text = element_text(size = 15, face = "bold"),
      plot.title = element_text(size = 15, hjust = 0.5), # Standard centered placement alignment
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15)
    )
  
  print(p)
  
  # 10. Close vector graphics pipeline device driver and flush files
  dev.off()
  
  cat(paste("Successfully saved vector plot panel to:", output_filename, "\n"))
}
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Enforce your custom target facet layout positioning order using exact labels
target_order <- c(
  "CD45+",
  "CD4+",
  "CD4+ FOXP3+",
  "CD4+ CXCR3+",
  "CD8+",
  "CD68+",
  "CD14+"
)

# Your exact requested base color definitions
mycols <- c(
  "SnC-" = "#bebebe",  # Slate Charcoal Gray
  "SnC+" = "#3156a6"   # Classic Steel Blue
)

age_groups <- c('Old', 'Young')

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

# 1. Enforce your exact facet order layout
target_order <- c(
  "CD45+",
  "CD4+",
  "CD4+ FOXP3+",
  "CD4+ CXCR3+",
  "CD8+",
  "CD68+",
  "CD14+"
)

# Your exact two-color definitions: Blue for SnC+, Grey for the rest (SnC-)
mycols <- c(
  "SnC+" = "#3156a6",
  "SnC-" = "#bebebe"
)

age_groups <- c('Old', 'Young')

for (current_age in age_groups) {
  
  cat(paste("Processing data for age group:", current_age, "\n"))
  
  # 1. Subset the Seurat object for the current age group
  seurat_subset <- subset(obj_joined, subset = age == current_age)
  
  # 2. Extract underlying cell logical marker profiles cleanly
  cell_matrix <- data.frame(
    snc   = as.character(seurat_subset$snc),
    CD45  = seurat_subset@assays$RNA$counts["CD45", ] > 0,
    CD4   = seurat_subset@assays$RNA$counts["CD4", ] > 0,
    CD8   = seurat_subset@assays$RNA$counts["CD8", ] > 0,
    CD14  = seurat_subset@assays$RNA$counts["CD14", ] > 0,
    CD68  = seurat_subset@assays$RNA$counts["CD68", ] > 0,
    FOXP3 = seurat_subset@assays$RNA$counts["FOXP3", ] > 0,
    CXCR3 = seurat_subset@assays$RNA$counts["CXCR3", ] > 0
  ) %>%
    filter(!is.na(snc))
  
  # 3. Create clean, independent summaries calculating the true SnC+ percentage directly
  calc_snc_pct <- function(df, label) {
    if(nrow(df) == 0) return(data.frame(cell_types = label, snc_pct = 0))
    pct <- mean(df$snc == "SnC+") * 100
    data.frame(cell_types = label, snc_pct = round(pct, 2))
  }
  
  summary_df <- bind_rows(
    calc_snc_pct(cell_matrix %>% filter(CD45 == TRUE), "CD45+"),
    calc_snc_pct(cell_matrix %>% filter(CD4 == TRUE), "CD4+"),
    calc_snc_pct(cell_matrix %>% filter(CD8 == TRUE), "CD8+"),
    calc_snc_pct(cell_matrix %>% filter(CD14 == TRUE), "CD14+"),
    calc_snc_pct(cell_matrix %>% filter(CD68 == TRUE), "CD68+"),
    calc_snc_pct(cell_matrix %>% filter(CD4 == TRUE & FOXP3 == TRUE), "CD4+ FOXP3+"),
    calc_snc_pct(cell_matrix %>% filter(CD4 == TRUE & CXCR3 == TRUE), "CD4+ CXCR3+")
  )
  
  # 4. Pivot to generate explicit SnC+ and SnC- rows so ggplot can build the pie slices
  your_data <- summary_df %>%
    mutate(SnC_plus = snc_pct, SnC_minus = 100 - snc_pct) %>%
    select(-snc_pct) %>%
    pivot_longer(cols = c(SnC_plus, SnC_minus), names_to = "SnC", values_to = "Freq") %>%
    mutate(SnC = ifelse(SnC == "SnC_plus", "SnC+", "SnC-"))
  
  # Enforce factor levels and ordering layers cleanly
  your_data$cell_types <- factor(your_data$cell_types, levels = target_order)
  your_data$SnC        <- factor(your_data$SnC, levels = c("SnC-", "SnC+"))
  
  # Calculate text label center positions
  your_data <- your_data %>%
    arrange(cell_types, desc(SnC)) %>%
    group_by(cell_types) %>%
    mutate(lab.ypos = cumsum(Freq) - 0.5 * Freq) %>%
    ungroup()
  
  # Print the clean data frame to the console to verify the math is flawless
  print(your_data)
  
  # 5. Save vector PDF output panel
  output_filename <- paste0("//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2/results/codex_im_", tolower(current_age), "_immsnc.pdf")
  num_cols <- length(target_order)
  pdf_width <- max(6, num_cols * 4.0)
  
  pdf(output_filename, width = pdf_width, height = 8)
  
  p <- ggplot(your_data, aes(x = "", y = Freq, group = SnC, color = SnC, fill = SnC)) +
    geom_bar(width = 1, stat = "identity", color = "white", linewidth = 0.3) + 
    coord_polar("y", start = 0) + 
    # Places text label metrics only on segments that have real values (> 0%)
    geom_text(aes(y = lab.ypos, label = ifelse(Freq > 0, paste0(Freq, "%"), "")), color = "black", size = 5, fontface = "bold") +
    scale_fill_manual(values = mycols) +
    scale_color_manual(values = mycols) +
    facet_wrap(~ cell_types, ncol = num_cols) + 
    labs(fill = "Senescent Status") +
    theme_void() +
    theme(
      strip.text = element_text(size = 15, face = "bold"),
      plot.title = element_text(size = 15, hjust = 0.5), 
      legend.text = element_text(size = 15),
      legend.title = element_text(size = 15)
    )
  
  print(p)
  dev.off()
  
  cat(paste("Successfully saved vector plot panel to:", output_filename, "\n"))
}


#######################################
#Pie charts for snc+
#####################################

# 1. Define the colors (using your existing palette structure)
minor_colors <- c(
  "Basal Keratinocyte" = "#1f78b4",
  "Fibroblast" = "#33a02c",
  "Melanocyte" = "#6a3d9a",
  "Corneocyte" = "#fb9a99",
  "Smooth Muscle" = "#ff7f00",
  "none" = "grey80"
)

# 2. Prepare the data
# We filter by 'minorskincell1' to stay consistent with your CODEX script logic
snc_comp <- obj_codex@meta.data %>%
  filter(snc == "SnC+" & !is.na(age) & !is.na(minorskincell1)) %>%
  group_by(age, minorskincell1) %>%
  summarise(
    cells = n(),
    .groups = "drop"
  ) %>%
  group_by(age) %>%
  mutate(
    pct = 100 * cells / sum(cells)
  )

library(ggplot2)
library(dplyr)
library(patchwork)

# snc_comp must have: age, majorskincell, pct (or cells)
# Ensure consistent order of ages
snc_comp <- snc_comp %>%
  mutate(age = factor(age, levels = c("Young", "Old")))

# 3. Use your existing donut function
# (This remains unchanged from your 4i script)
donut_one <- function(df, age_value) {
  ggplot(df %>% filter(age == age_value),
         aes(x = 2, y = pct, fill = minorskincell1)) +
    geom_col(width = 1, color = "white", linewidth = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    scale_fill_manual(values = minor_colors, drop = FALSE) +
    theme_void() +
    ggtitle(as.character(age_value)) +
    theme(
      text = element_text(family = "Helvetica", size = 14),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# 4. Generate and save
p_young <- donut_one(snc_comp, "Young")
p_old   <- donut_one(snc_comp, "Old")

p_codex_donut <- (p_young | p_old) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# Save as PDF
ggsave("results/260617/codex_snc_donut_minor.pdf", p_codex_donut, width = 7.5, height = 3.5)

snc_comp1 <- obj_codex@meta.data %>%
  filter(!is.na(age) & !is.na(minorskincell1)) %>%
  group_by(age, minorskincell1) %>%
  summarise(
    cells = n(),
    .groups = "drop"
  ) %>%
  group_by(age) %>%
  mutate(
    pct = 100 * cells / sum(cells)
  )

library(ggplot2)
library(dplyr)
library(patchwork)

# snc_comp must have: age, majorskincell, pct (or cells)
# Ensure consistent order of ages
snc_comp1 <- snc_comp1 %>%
  mutate(age = factor(age, levels = c("Young", "Old")))
# 3. Use your existing donut function
# (This remains unchanged from your 4i script)
donut_one1 <- function(df, age_value) {
  ggplot(df %>% filter(age == age_value),
         aes(x = 2, y = pct, fill = minorskincell1)) +
    geom_col(width = 1, color = "white", linewidth = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    scale_fill_manual(values = minor_colors, drop = FALSE) +
    theme_void() +
    ggtitle(as.character(age_value)) +
    theme(
      text = element_text(family = "Helvetica", size = 14),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# 4. Generate and save
p_young <- donut_one(snc_comp1, "Young")
p_old   <- donut_one(snc_comp1, "Old")

p_codex_donut <- (p_young | p_old) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# Save as PDF
ggsave("results/260617/codex_total_donut_minor.pdf", p_codex_donut, width = 7.5, height = 3.5)

##########################################################
# SnC marker dot plot by immune cell type
##########################################################

features_snc <- c("p21", "p16", "PCNA", "BCL2", "IP10", "IL6", "IL8", "CXCR3")

obj_codex$immcelltype <- factor(
  obj_codex$immcelltype,
  levels = c(
    "CD4+ FOXP3+",
    "CD4+ CXCR3+",
    "CD4+",
    "CD8+",
    "CD14+",
    "CD68+",
    "CD45+",
    "none"
  )
)

p_snc_marker_dot <- DotPlot(
  object = obj_codex,
  features = features_snc,
  group.by = "immcelltype",
  assay = "RNA",
  col.min = NULL,
  col.max = NULL
) +
  scale_size_continuous(
    range = c(4, 14),
    labels = function(x) paste0(round(x, 1), "%"),
    name = "Percent\nExpressing"
  ) +
  scale_color_gradientn(
    colors = c("#2166ac", "#67a9cf", "#f7f7f7", "#ef8a62", "#b2182b"),
    name = "Average\nExpression"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  labs(
    title = "Senescence Marker Expression Across CODEX Immune Cell Types"
  )

p_snc_marker_dot

ggsave(
  filename = "results/codex_snc_marker_dotplot_by_immcelltype.pdf",
  plot = p_snc_marker_dot,
  width = 9,
  height = 5.5
)

##########################################################
# 9) Save annotated list
##########################################################
save.image("data/codex/seurat_list_codex_annotated.RData") #260617 260708
