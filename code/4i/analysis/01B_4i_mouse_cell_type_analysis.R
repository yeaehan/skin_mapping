##########################################################
# 4i MOUSE ANALYSIS + PLOTTING (FINAL CLEAN VERSION)
##########################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(ggpubr)

setwd("//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2")
load("data/4i_mouse/260416.RData")
load("data/4i_mouse/seurat_list_4i_mouse.RData")

##########################################################
# Marker panel (corrected)
##########################################################
marker_panel <- c(
  "p21","p16","PCNA","HMGB1","LaminB1","BCL2",
  "IL1a","IL1b","RANTES",
  "Keratin 14","Keratin 10","Involucrin",
  "MelanA","DCT",
  "Vimentin","CD3","F4_80"
)

marker_list <- list()

for (sample_name in names(seurat_list_mouse)) {
  
  obj <- seurat_list_mouse[[sample_name]]
  
  sample_markers <- list()
  
  for (marker in marker_panel) {
    if (marker %in% rownames(obj[["RNA"]]$counts)) {
      sample_markers[[marker]] <- obj[["RNA"]]$counts[marker, ] > 0
    } else {
      sample_markers[[marker]] <- NULL
    }
  }
  
  marker_list[[sample_name]] <- sample_markers
}

get_marker <- function(marker_vecs, marker, cells) {
  v <- marker_vecs[[marker]]
  if (is.null(v)) return(setNames(rep(FALSE, length(cells)), cells))
  v[cells]
}

##########################################################
# ANNOTATION
##########################################################
for (sample_name in names(seurat_list_mouse)) {
  
  obj <- seurat_list_mouse[[sample_name]]
  cells <- colnames(obj)
  m <- marker_list[[sample_name]]
  
  p21 <- get_marker(m,"p21",cells)
  p16 <- get_marker(m,"p16",cells)
  PCNA <- get_marker(m,"PCNA",cells)
  BCL2 <- get_marker(m,"BCL2",cells)
  
  IL1a <- get_marker(m,"IL1a",cells)
  IL1b <- get_marker(m,"IL1b",cells)
  RANTES <- get_marker(m,"RANTES",cells)
  
  K14 <- get_marker(m,"Keratin 14",cells)
  K10 <- get_marker(m,"Keratin 10",cells)
  INV <- get_marker(m,"Involucrin",cells)
  
  MelanA <- get_marker(m,"MelanA",cells)
  DCT <- get_marker(m,"DCT",cells)
  
  Vimentin <- get_marker(m,"Vimentin",cells)
  
  # cdkn state
  obj$cdkn_state <- case_when(
    p16 & !p21 ~ "p16+",
    p21 & !p16 ~ "p21+",
    p16 & p21 ~ "p16+p21+",
    TRUE ~ "p16-p21-"
  )
  
  # SnC definition (FIXED)
  obj$snc <- ifelse(
    (p21 | p16) & (!PCNA) & (IL1a | IL1b | RANTES),
    "SnC+","SnC-"
  )
  
  obj$snc2 <- ifelse(
    (p21 | p16) & (!PCNA) & BCL2 & (IL1a | IL1b | RANTES),
    "SnC+","SnC-"
  )
  
  obj$snc <- factor(obj$snc, levels=c("SnC-","SnC+"))
  obj$snc2 <- factor(obj$snc2, levels=c("SnC-","SnC+"))
  
  # cell types
  obj$majorskincell <- case_when(
    MelanA | DCT ~ "Melanocyte",
    Vimentin ~ "Fibroblast",
    (K14 | K10 | INV) ~ "Keratinocyte",
    TRUE ~ "none"
  )
  
  obj$minorskincell <- case_when(
    MelanA | DCT ~ "Melanocyte",
    Vimentin ~ "Fibroblast",
    K14 ~ "Basal Keratinocyte",
    K10 ~ "Suprabasal Keratinocyte",
    INV ~ "Corneocyte",
    TRUE ~ "none"
  )
  
  obj$age <- factor(obj$age, levels=c("Young","Old"))
  
  seurat_list_mouse[[sample_name]] <- obj
}

##########################################################
# MERGE
##########################################################
obj_4i_mouse <- merge(seurat_list_mouse[[1]], y = seurat_list_mouse[-1])

##########################################################
# SnC proportion vs age
##########################################################
blue_col <- "#296d9e"

df_mouse <- bind_rows(lapply(names(seurat_list_mouse), function(sample_name) {
  
  md <- seurat_list_mouse[[sample_name]]@meta.data
  
  data.frame(
    sample = sample_name,
    age = unique(md$age)[1],
    snc_prop = mean(md$snc == "SnC+")
  )
})) %>%
  mutate(age=factor(age,levels=c("Young","Old")))

p_age <- ggplot(df_mouse, aes(x=age,y=snc_prop)) +
  stat_summary(fun=mean,geom="bar",fill=blue_col,color="black") +
  stat_summary(fun.data=mean_se,geom="errorbar", width=0.1) +
  geom_point(position=position_jitter(width=0.07),size=3) +
  scale_y_continuous(labels=percent_format(accuracy=0.1)) +
  theme_classic(base_size=18)

print(p_age)

##########################################################
# VIOLIN PLOTS (correct markers)
##########################################################

plot_seurat_violin_with_sampledots <- function(obj,
                                               feature,
                                               group_col = "age",
                                               sample_col = "sample",
                                               assay = "RNA",
                                               layer = "data",
                                               cols = c("Young" = "#c8def9", "Old" = "#d7b0b0"),
                                               dot_cols = c("Young" = "#000080", "Old" = "#a00000")) {
  
  # ---------------------------------------
  # 1) Base Seurat violin
  # ---------------------------------------
  p <- VlnPlot(
    obj,
    features = feature,
    group.by = group_col,
    assay = assay,
    layer = layer,
    pt.size = 0,
    cols = cols
  )
  
  # ---------------------------------------
  # 2) Fetch data for sample-level dots
  # ---------------------------------------
  df <- FetchData(
    obj,
    vars = c(feature, group_col, sample_col),
    assay = assay,
    layer = layer
  )
  colnames(df) <- c("value", "group", "sample")
  
  df <- df %>%
    filter(!is.na(group), !is.na(sample)) %>%
    mutate(group = factor(group, levels = c("Young", "Old")))
  
  # ---------------------------------------
  # 3) Compute sample means
  # ---------------------------------------
  df_sample <- df %>%
    group_by(sample, group) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  ymax <- max(df$value, na.rm = TRUE)
  
  # ---------------------------------------
  # 4) Add sample dots (colored by age)
  # ---------------------------------------
  p <- p +
    geom_point(
      data = df_sample,
      aes(
        x = as.numeric(factor(group, levels = c("Young", "Old"))),
        y = mean_value,
        color = group
      ),
      position = position_jitter(width = 0.05, height = 0, seed = 42),
      size = 3,
      alpha = 0.9,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = dot_cols) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      label.x = 1.5,
      label.y = ymax * 0.98,
      size = 6
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

# If needed, merge the 4i list first
# obj_4i <- merge(seurat_list[[1]], y = seurat_list[-1], merge.data = TRUE)

# Use source_file if sample column is missing
if (!"sample" %in% colnames(obj_4i_mouse@meta.data)) {
  if ("source_file" %in% colnames(obj_4i_mouse@meta.data)) {
    obj_4i_mouse$sample <- obj_4i_mouse$source_file
  } else {
    obj_4i_mouse$sample <- obj_4i_mouse$orig.ident
  }
}

# Make sure age is formatted consistently
obj_4i_mouse$age <- factor(obj_4i_mouse$age, levels = c("Young", "Old"))

features <- c("BCL2", "HMGB1", "LaminB1")

p_list <- lapply(features, function(f) {
  plot_seurat_violin_with_sampledots(
    obj = obj_4i_mouse,
    feature = f,
    group_col = "age",
    sample_col = "sample",
    assay = "RNA",
    layer = "data"
  )
})

final_plot <- p_list[[1]] + p_list[[2]] + p_list[[3]]

final_plot

ggsave(
  "//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2/results/260303/4i_mouse_subsenmarkers_sampledots_by_age.png",
  plot = final_plot,
  width = 10,
  height = 3,
  dpi = 300
)

ggsave(
  "//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2/results/260303/4i_mouse_subsenmarkers_sampledots_by_age.pdf",
  plot = final_plot,
  width = 10,
  height = 3,
  dpi = 300
)

##########################################################
# SnC distribution by cell type
##########################################################
freq_major <- obj_4i_mouse@meta.data %>%
  group_by(majorskincell, age) %>%
  summarise(
    snc_pct = mean(snc=="SnC+")*100,
    .groups="drop"
  )

ggplot(freq_major,aes(x=snc_pct,y=majorskincell,color=age))+
  geom_point(size=4)+
  theme_bw()

##########################################################
# UMAP (low marker panel fix)
##########################################################
obj_4i_mouse <- NormalizeData(obj_4i_mouse)
obj_4i_mouse <- ScaleData(obj_4i_mouse)

obj_4i_mouse <- RunPCA(obj_4i_mouse, npcs=5)
obj_4i_mouse <- FindNeighbors(obj_4i_mouse,dims=1:5)
obj_4i_mouse <- FindClusters(obj_4i_mouse,resolution=0.5)
obj_4i_mouse <- RunUMAP(obj_4i_mouse,dims=1:5)

p1 <- DimPlot(obj_4i_mouse,label=TRUE)
p2 <- DimPlot(obj_4i_mouse,group.by="majorskincell")
p3 <- DimPlot(obj_4i_mouse,group.by="snc",cols=c("grey80","red"))

(p1|p2)/(p3)



##########################################################
# Cell Type Distribution
##########################################################

# ------------------------------------------------------------------
# 1. MASTER COLOR PALETTE & STACK ORDERS FOR obj_4i_mouse
# ------------------------------------------------------------------
elegant_colors <- c(
  "Keratinocyte"            = "#67001f",
  "Basal Keratinocyte"      = "#b2182b",
  "Corneocyte"              = "#fddbc7",
  "none"                    = "#f7f7f7",
  "Melanocyte"              = "#d1e5f0",
  "Suprabasal Keratinocyte" = "#2166ac",
  "Fibroblast"              = "#053061"
)

major_order <- c("Keratinocyte", "Melanocyte", "Fibroblast", "none")
minor_order <- c(
  "Suprabasal Keratinocyte",
  "Basal Keratinocyte",
  "Corneocyte",
  "Melanocyte",
  "Fibroblast",
  "none"
)

# ------------------------------------------------------------------
# 2. SnC+ COMPOSITION PLOTS (among all cells in each age group)
# ------------------------------------------------------------------
df_major_snc <- obj_4i_mouse@meta.data %>%
  filter(!is.na(age), !is.na(majorskincell), !is.na(snc)) %>%
  group_by(age) %>%
  mutate(total_cells_in_age = n()) %>%
  filter(snc == "SnC+") %>%
  group_by(age, majorskincell, total_cells_in_age) %>%
  summarise(snc_count = n(), .groups = "drop") %>%
  mutate(
    pct_of_total = snc_count / total_cells_in_age,
    age = factor(age, levels = c("Young", "Old")),
    majorskincell = factor(majorskincell, levels = major_order)
  )

p_major_snc <- ggplot(df_major_snc, aes(x = age, y = pct_of_total, fill = majorskincell)) +
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

df_minor_snc <- obj_4i_mouse@meta.data %>%
  filter(!is.na(age), !is.na(minorskincell), !is.na(snc)) %>%
  group_by(age) %>%
  mutate(total_cells_in_age = n()) %>%
  filter(snc == "SnC+") %>%
  group_by(age, minorskincell, total_cells_in_age) %>%
  summarise(snc_count = n(), .groups = "drop") %>%
  mutate(
    pct_of_total = snc_count / total_cells_in_age,
    age = factor(age, levels = c("Young", "Old")),
    minorskincell = factor(minorskincell, levels = minor_order)
  )

p_minor_snc <- ggplot(df_minor_snc, aes(x = age, y = pct_of_total, fill = minorskincell)) +
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

print(p_major_snc)
print(p_minor_snc)

# ------------------------------------------------------------------
# 3. OVERALL CELL TYPE COMPOSITION PLOTS (all cells, not just SnC+)
# ------------------------------------------------------------------
df_major_cells <- obj_4i_mouse@meta.data %>%
  filter(!is.na(age), !is.na(majorskincell)) %>%
  count(age, majorskincell, name = "n_cells") %>%
  group_by(age) %>%
  mutate(
    total_cells_in_age = sum(n_cells),
    pct_of_total = n_cells / total_cells_in_age
  ) %>%
  ungroup() %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    majorskincell = factor(majorskincell, levels = major_order)
  )

p_major_cells <- ggplot(df_major_cells, aes(x = age, y = pct_of_total, fill = majorskincell)) +
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

df_minor_cells <- obj_4i_mouse@meta.data %>%
  filter(!is.na(age), !is.na(minorskincell)) %>%
  count(age, minorskincell, name = "n_cells") %>%
  group_by(age) %>%
  mutate(
    total_cells_in_age = sum(n_cells),
    pct_of_total = n_cells / total_cells_in_age
  ) %>%
  ungroup() %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    minorskincell = factor(minorskincell, levels = minor_order)
  )

p_minor_cells <- ggplot(df_minor_cells, aes(x = age, y = pct_of_total, fill = minorskincell)) +
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

print(p_major_cells)
print(p_minor_cells)

# ------------------------------------------------------------------
# 4. SAVE ALL PLOTS
# ------------------------------------------------------------------
ggsave("results/260303/4i_mouse_snc_composition_stacked_major.pdf", plot = p_major_snc, width = 5, height = 8)
ggsave("results/260303/4i_mouse_snc_composition_stacked_minor.pdf", plot = p_minor_snc, width = 6, height = 8)

ggsave("results/260303/4i_mouse_celltype_composition_stacked_major.pdf", plot = p_major_cells, width = 5, height = 8)
ggsave("results/260303/4i_mouse_celltype_composition_stacked_minor.pdf", plot = p_minor_cells, width = 6, height = 8)



##########################################################
# p16/p21 by cell types
##########################################################

# --------------------------------------------------
# helper: compute per-sample % marker+ cells
# marker positivity = expression > 0
# --------------------------------------------------
make_marker_prop_df <- function(obj,
                                marker_feature,
                                celltype_col = "minorskincell",
                                sample_col = "sample_id",
                                age_col = "age",
                                assay = "RNA",
                                layer = "counts") {
  
  # fetch expression + metadata together
  df <- FetchData(
    obj,
    vars = c(marker_feature, celltype_col, sample_col, age_col),
    assay = assay,
    layer = layer
  )
  
  colnames(df) <- c("marker_value", "minor_celltype", "sample_use", "age")
  
  df <- df %>%
    mutate(
      age = factor(age, levels = c("Young", "Old")),
      minor_celltype = factor(minor_celltype, levels = minor_order),
      marker_pos = marker_value > 0
    ) %>%
    filter(
      !is.na(age),
      !is.na(minor_celltype),
      !is.na(sample_use),
      minor_celltype != "none"
    ) %>%
    group_by(sample_use, age, minor_celltype) %>%
    summarise(
      total_cells = n(),
      marker_pos_cells = sum(marker_pos, na.rm = TRUE),
      prop_marker_pos = marker_pos_cells / total_cells,
      pct_marker_pos = 100 * prop_marker_pos,
      .groups = "drop"
    )
  
  df
}

# --------------------------------------------------
# build data frames
# --------------------------------------------------
df_p16 <- make_marker_prop_df(obj_4i_mouse, marker_feature = "p16")
df_p21 <- make_marker_prop_df(obj_4i_mouse, marker_feature = "p21")

# --------------------------------------------------
# plotting function
# --------------------------------------------------
plot_marker_minor <- function(df, y_lab = "% positive cells", title_text = NULL) {
  
  ggplot(df, aes(x = age, y = pct_marker_pos)) +
    stat_summary(
      fun = mean,
      geom = "bar",
      width = 0.6,
      fill = "#296d9e",
      color = "black"
    ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.05,
      linewidth = 0.4,
      alpha = 0.7,
      color = "black"
    ) +
    geom_point(
      position = position_jitter(width = 0.07, height = 0, seed = 42),
      size = 2.8,
      color = "black"
    ) +
    facet_wrap(~ minor_celltype, ncol = 1, scales = "free_y") +
    labs(
      x = NULL,
      y = y_lab,
      title = title_text
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 12)
    )
}

# --------------------------------------------------
# make plots
# --------------------------------------------------
p_p16_minor <- plot_marker_minor(
  df_p16,
  y_lab = "% p16+ cells",
  title_text = "p16+ cells by minor skin cell type"
)

p_p21_minor <- plot_marker_minor(
  df_p21,
  y_lab = "% p21+ cells",
  title_text = "p21+ cells by minor skin cell type"
)

print(p_p16_minor)
print(p_p21_minor)

# --------------------------------------------------
# save
# --------------------------------------------------
ggsave(
  "results/260303/4i_mouse_p16_positive_by_minor_celltype_age.pdf",
  plot = p_p16_minor,
  width = 4,
  height = 10,
  dpi = 300
)

ggsave(
  "results/260303/4i_mouse_p21_positive_by_minor_celltype_age.pdf",
  plot = p_p21_minor,
  width = 4,
  height = 10,
  dpi = 300
)

##########################################################
# Spatial Plots
##########################################################

library(colorspace)
library(scales)
library(patchwork)

# --------------------------------------------------
# Inputs
# --------------------------------------------------
# Assumes you already have:
#   seurat_list_mouse
# and each object already contains:
#   minorskincell
#   cdkn_state
#   sample_id / source_file
#   spatial reduction OR Spatial_X / Spatial_Y metadata

# --------------------------------------------------
# Color palettes
# --------------------------------------------------
minor_colors <- c(
  "Keratinocyte"            = "#67001f",
  "Basal Keratinocyte"      = "#b2182b",
  "Corneocyte"              = "#fddbc7",
  "none"                    = "#f7f7f7",
  "Melanocyte"              = "orange",
  "Suprabasal Keratinocyte" = "#2166ac",
  "Fibroblast"              = "#053061"
)

cdkn_colors <- c(
  "p16-p21-" = "#d9d9d9",
  "p21+" = "#296d9e",
  "p16+" = "#cc0b01",
  "p16+p21+" = "#7f7f7f"
)

minor_order <- c(
  "Fibroblast",
  "Melanocyte",
  "Basal Keratinocyte",
  "Suprabasal Keratinocyte",
  "Corneocyte",
  "none"
)

cdkn_order <- c("p16-p21-", "p21+", "p16+", "p16+p21+")

# --------------------------------------------------
# Helper: extract spatial coordinates
# --------------------------------------------------
get_spatial_df <- function(obj) {
  md <- obj@meta.data
  md$cell_id <- rownames(md)
  
  # Case 1: already stored in metadata
  if (all(c("Spatial_X", "Spatial_Y") %in% colnames(md))) {
    df <- md %>%
      mutate(
        x = Spatial_X,
        y = Spatial_Y
      )
    return(df)
  }
  
  # Case 2: stored in spatial reduction
  if ("spatial" %in% names(obj@reductions)) {
    emb <- as.data.frame(obj@reductions$spatial@cell.embeddings)
    emb$cell_id <- rownames(emb)
    
    # take first 2 columns regardless of names
    colnames(emb)[1:2] <- c("x", "y")
    
    df <- md %>%
      left_join(emb[, c("cell_id", "x", "y")], by = "cell_id")
    return(df)
  }
  
  stop("No spatial coordinates found: expected Spatial_X/Spatial_Y metadata or a spatial reduction.")
}

# --------------------------------------------------
# Helper: one spatial plot
# --------------------------------------------------
plot_spatial_categorical <- function(df, color_col, palette, levels_use, title_text,
                                     point_size = 0.45, alpha_val = 0.95) {
  
  df[[color_col]] <- factor(df[[color_col]], levels = levels_use)
  
  ggplot(df, aes(x = x, y = y, color = .data[[color_col]])) +
    geom_point(size = point_size, alpha = alpha_val, stroke = 0) +
    coord_fixed() +
    scale_y_reverse() +
    scale_color_manual(values = palette, drop = FALSE, na.value = "grey85") +
    labs(title = title_text, color = NULL) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right",
      legend.text = element_text(size = 10)
    )
}

# --------------------------------------------------
# Loop through samples and make two plots each
# --------------------------------------------------
dir.create("results/260303/4i_mouse_spatial", recursive = TRUE, showWarnings = FALSE)

for (sample_name in names(seurat_list_mouse)) {
  
  obj <- seurat_list_mouse[[sample_name]]
  
  # skip sample if required metadata are missing
  needed_cols <- c("centroid_x_pixels", "centroid_y_pixels", "minorskincell", "cdkn_state")
  missing_cols <- setdiff(needed_cols, colnames(obj@meta.data))
  
  if (length(missing_cols) > 0) {
    message("Skipping ", sample_name, " because missing: ", paste(missing_cols, collapse = ", "))
    next
  }
  
  obj$minorskincell <- factor(obj$minorskincell, levels = minor_order)
  obj$cdkn_state <- factor(obj$cdkn_state, levels = cdkn_order)
  
  df_plot <- obj@meta.data %>%
    mutate(
      x = centroid_x_pixels,
      y = centroid_y_pixels
    ) %>%
    filter(!is.na(x), !is.na(y))
  
  # minor: none at bottom
  df_plot_minor <- df_plot %>%
    mutate(
      minorskincell = factor(minorskincell, levels = minor_order),
      draw_order = ifelse(minorskincell == "none", 0, 1)
    ) %>%
    arrange(draw_order)
  
  # cdkn: double negatives at bottom
  df_plot_cdkn <- df_plot %>%
    mutate(
      cdkn_state = factor(cdkn_state, levels = cdkn_order),
      draw_order = ifelse(cdkn_state == "p16-p21-", 0, 1)
    ) %>%
    arrange(draw_order)
  p_minor <- plot_spatial_categorical(
    df = df_plot_minor,
    color_col = "minorskincell",
    palette = minor_colors,
    levels_use = minor_order,
    title_text = paste0(sample_name, " - Minor skin cell types"),
    point_size = 0.5
  )
  
  p_cdkn <- plot_spatial_categorical(
    df = df_plot_cdkn,
    color_col = "cdkn_state",
    palette = cdkn_colors,
    levels_use = cdkn_order,
    title_text = paste0(sample_name, " - p16/p21 state"),
    point_size = 0.5
  )
  
  print(p_minor)
  print(p_cdkn)
  
  ggsave(
    filename = paste0("results/260303/4i_mouse_spatial/", sample_name, "_minor_skin_celltypes.pdf"),
    plot = p_minor,
    width = 7,
    height = 6,
    dpi = 600
  )
  
  ggsave(
    filename = paste0("results/260303/4i_mouse_spatial/", sample_name, "_cdkn_state.pdf"),
    plot = p_cdkn,
    width = 7,
    height = 6,
    dpi = 600
  )
  
  p_combined <- p_minor + p_cdkn + plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "right")
  
  ggsave(
    filename = paste0("results/260303/4i_mouse_spatial/", sample_name, "_minor_and_cdkn_combined.pdf"),
    plot = p_combined,
    width = 12,
    height = 6,
    dpi = 600
  )
}

save.image("data/4i_mouse/260416.RData")





snc_markers <- c(
  "p21", "p16", "PCNA", "BCL2",
  "IL1a", "IL1b", "RANTES",
  "HMGB1", "LaminB1"
)

snc_markers <- snc_markers[snc_markers %in% rownames(obj_4i_mouse)]

# --------------------------------------------------
# make sure sample column exists
# --------------------------------------------------
if (!"sample" %in% colnames(obj_4i_mouse@meta.data)) {
  if ("sample_id" %in% colnames(obj_4i_mouse@meta.data)) {
    obj_4i_mouse$sample <- obj_4i_mouse$sample_id
  } else if ("source_file" %in% colnames(obj_4i_mouse@meta.data)) {
    obj_4i_mouse$sample <- obj_4i_mouse$source_file
  } else {
    obj_4i_mouse$sample <- obj_4i_mouse$orig.ident
  }
}

obj_4i_mouse$age <- factor(obj_4i_mouse$age, levels = c("Young", "Old"))

# --------------------------------------------------
# fetch expression + metadata
# --------------------------------------------------
df_expr <- FetchData(
  obj_4i_mouse,
  vars = c(snc_markers, "sample", "age"),
  assay = "RNA",
  layer = "data"
)

# --------------------------------------------------
# per-sample mean intensity for each marker
# --------------------------------------------------
df_sample_marker <- df_expr %>%
  tidyr::pivot_longer(
    cols = all_of(snc_markers),
    names_to = "marker",
    values_to = "intensity"
  ) %>%
  group_by(sample, age, marker) %>%
  summarise(
    mean_intensity = mean(intensity, na.rm = TRUE),
    median_intensity = median(intensity, na.rm = TRUE),
    pct_positive = mean(intensity > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# inspect table
print(df_sample_marker)

# --------------------------------------------------
# 1) sample-level mean intensity plot
# bars = age-group mean, dots = samples
# --------------------------------------------------
p_mean <- ggplot(df_sample_marker, aes(x = age, y = mean_intensity)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    width = 0.6,
    fill = "#296d9e",
    color = "black"
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.05,
    linewidth = 0.4,
    alpha = 0.7,
    color = "black"
  ) +
  geom_point(
    aes(color = age),
    position = position_jitter(width = 0.07, height = 0, seed = 42),
    size = 2.8,
    alpha = 0.9
  ) +
  scale_color_manual(values = c("Young" = "#000080", "Old" = "#a00000")) +
  facet_wrap(~ marker, ncol = 3, scales = "free_y") +
  labs(
    x = NULL,
    y = "Mean intensity",
    title = "SnC marker intensities by sample"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11)
  )

print(p_mean)

# --------------------------------------------------
# 2) optional heatmap-style sample x marker plot
# --------------------------------------------------
sample_order <- df_sample_marker %>%
  distinct(sample, age) %>%
  arrange(age, sample) %>%
  pull(sample)

p_heat <- df_sample_marker %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    marker = factor(marker, levels = snc_markers)
  ) %>%
  ggplot(aes(x = marker, y = sample, fill = mean_intensity)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "grey95", high = "#08519c") +
  labs(
    x = NULL,
    y = NULL,
    fill = "Mean\nintensity",
    title = "Sample-level SnC marker intensity heatmap"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

print(p_heat)

















library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)

# ---------------------------------------
# 1) Define markers (exactly your panel)
# ---------------------------------------
sen_markers <- c(
  "p21", "p16", "PCNA", "HMGB1", "LaminB1", "BCL2"
) #removed "IL1a", "IL1b"

skin_markers <- c(
  "Keratin 14", "Keratin 10", "MelanA", "Vimentin", "Involucrin")

# keep only those present
sen_markers <- sen_markers[sen_markers %in% rownames(seu)]
skin_markers <- skin_markers[skin_markers %in% rownames(seu)]

# ---------------------------------------
# 2) Fetch expression
# ---------------------------------------
expr_df <- FetchData(
  seu,
  vars = c(sen_markers, skin_markers),
  assay = "RNA",
  layer = "data"
)

# ---------------------------------------
# 3–4) Build grid with TRUE co-expression
# ---------------------------------------
plot_df <- expand.grid(
  sen_marker = sen_markers,
  skin_marker = skin_markers,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    # TRUE co-expression (% cells expressing BOTH)
    pct_plot = mean(
      expr_df[[sen_marker]] > 0 &
        expr_df[[skin_marker]] > 0,
      na.rm = TRUE
    ) * 100,
    
    # average expression (for color)
    avg_plot = mean(
      (expr_df[[sen_marker]] + expr_df[[skin_marker]]) / 2,
      na.rm = TRUE
    )
  ) %>%
  ungroup()
# ---------------------------------------
# 5) Plot
# ---------------------------------------
p_grid <- ggplot(plot_df, aes(x = skin_marker, y = sen_marker)) +
  geom_point(aes(size = pct_plot, color = avg_plot)) +
  
  scale_size(range = c(2, 10), name = "% expressed") +
  scale_color_gradientn(
    colours = c("lightgrey", "#9ecae1", "#3182bd", "#08519c"),
    name = "Avg expression"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = "Skin / lineage markers",
    y = "Senescence markers"
  )

print(p_grid)

ggsave(
  "results/260506/4i_mouse_marker_grid_allcells_fixed.pdf",
  p_grid,
  width = 6,
  height = 5,
  dpi = 300
)

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)






library(Seurat)
library(dplyr)
library(ggplot2)

# ---------------------------------------
# colors (same style as your example)
# ---------------------------------------
mycols_p16 <- c("#3156a6", "#cd8729")
mycols_p21 <- c("#3156a6", "#cd8729")

# ---------------------------------------
# skin markers to include
# ---------------------------------------
skin_markers <- c(
  "Keratin 14", "Keratin 10", "MelanA",
  "DCT", "Vimentin", "Involucrin"
)

skin_markers <- skin_markers[skin_markers %in% rownames(seu)]

# ---------------------------------------
# helper to build pie-chart dataframe
# ---------------------------------------
make_skin_marker_pie_df <- function(seu, sen_marker = "p16", skin_markers) {
  
  df <- FetchData(
    seu,
    vars = c(skin_markers, sen_marker),
    assay = "RNA",
    layer = "counts"
  )
  
  out <- lapply(skin_markers, function(sm) {
    
    df_sub <- df[df[[sm]] > 0, , drop = FALSE]
    
    if (nrow(df_sub) == 0) {
      return(data.frame(
        cell_types = sm,
        state = c("neg", "pos"),
        Freq = c(0, 0)
      ))
    }
    
    tab <- table(ifelse(df_sub[[sen_marker]] > 0, "pos", "neg"))
    
    # ensure both levels exist
    n_neg <- if ("neg" %in% names(tab)) as.numeric(tab["neg"]) else 0
    n_pos <- if ("pos" %in% names(tab)) as.numeric(tab["pos"]) else 0
    
    data.frame(
      cell_types = sm,
      state = c("neg", "pos"),
      Freq = c(n_neg, n_pos)
    )
  })
  
  your_data <- bind_rows(out) %>%
    group_by(cell_types) %>%
    mutate(Freq = Freq / sum(Freq)) %>%
    ungroup()
  
  your_data$Freq <- round(your_data$Freq * 100, 2)
  
  your_data <- your_data %>%
    arrange(desc(cell_types)) %>%
    group_by(cell_types) %>%
    mutate(lab.ypos = cumsum(Freq) - 0.5 * Freq) %>%
    ungroup()
  
  your_data
}

# ---------------------------------------
# make dataframes
# ---------------------------------------
your_data_p16 <- make_skin_marker_pie_df(
  seu = seu,
  sen_marker = "p16",
  skin_markers = skin_markers
)

your_data_p21 <- make_skin_marker_pie_df(
  seu = seu,
  sen_marker = "p21",
  skin_markers = skin_markers
)

# optional: order facets biologically
facet_order <- c("Vimentin", "MelanA", "DCT", "Keratin 14", "Keratin 10", "Involucrin")
facet_order <- facet_order[facet_order %in% unique(your_data_p16$cell_types)]

your_data_p16$cell_types <- factor(your_data_p16$cell_types, levels = facet_order)
your_data_p21$cell_types <- factor(your_data_p21$cell_types, levels = facet_order)

# ---------------------------------------
# p16 pies
# ---------------------------------------
p_p16 <- ggplot(your_data_p16, aes(x = "", y = Freq, group = state, color = state, fill = state)) +
  geom_bar(width = 2, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(y = lab.ypos, label = Freq), color = "black", size = 4) +
  scale_fill_manual(values = mycols_p16) +
  scale_color_manual(values = mycols_p16) +
  facet_grid(. ~ cell_types) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# ---------------------------------------
# p21 pies
# ---------------------------------------
p_p21 <- ggplot(your_data_p21, aes(x = "", y = Freq, group = state, color = state, fill = state)) +
  geom_bar(width = 2, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(y = lab.ypos, label = Freq), color = "black", size = 4) +
  scale_fill_manual(values = mycols_p21) +
  scale_color_manual(values = mycols_p21) +
  facet_grid(. ~ cell_types) +
  theme_void() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

print(p_p16)
print(p_p21)

ggsave("results/260303/4i_mouse_skinmarker_p16_pies_stylematch.pdf",
       plot = p_p16, width = 10, height = 4, dpi = 300)

ggsave("results/260303/4i_mouse_skinmarker_p21_pies_stylematch.pdf",
       plot = p_p21, width = 10, height = 4, dpi = 300)