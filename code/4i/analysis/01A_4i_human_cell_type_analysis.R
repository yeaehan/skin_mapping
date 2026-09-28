
##########################################################
#I haven't integrated 4i samples yet!!!
##########################################################



library(Seurat)
library(dplyr)
library(ggplot2)

#source("code/plotfunctions.R")
setwd("//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2")

save.image('data/4i/260225.RData') #260326 260331 260403 260617
load("data/4i/260225.RData")

load("data/4i/seurat_list_4i_annotated.RData")

load("data/4i/seurat_list_4i.RData")

marker_panel <- c(
  "p21", "p16", "PCNA", "HMGB1", "LaminB1", "BCL2", "IP10", "CXCR3",
  "Keratin 5", "Keratin 10", "MelanA", "Filaggrin", "Vimentin"
)

# This will store results for all samples
marker_list <- list()

for (sample_name in names(seurat_list)) {
  
  cat("Processing:", sample_name, "\n")
  
  obj <- seurat_list[[sample_name]]
  
  # Store markers for this object
  sample_markers <- list()
  
  for (marker in marker_panel) {
    
    if (marker %in% rownames(obj@assays$RNA$counts)) {
      
      sample_markers[[marker]] <- 
        obj@assays$RNA$counts[marker, ] > 0
      
    } else {
      
      warning(paste("Marker", marker, "not found in", sample_name))
      sample_markers[[marker]] <- NULL
      
    }
  }
  
  marker_list[[sample_name]] <- sample_markers
}

library(dplyr)

# Helper: safely fetch a marker logical vector for an object (returns all FALSE if missing)
get_marker <- function(marker_vecs, marker, cells) {
  v <- marker_vecs[[marker]]
  if (is.null(v)) return(setNames(rep(FALSE, length(cells)), cells))
  # ensure same order as cells
  v[cells]
}

for (sample_name in names(seurat_list)) {
  cat("Annotating:", sample_name, "\n")
  
  obj <- seurat_list[[sample_name]]
  cells <- colnames(obj)
  m <- marker_list[[sample_name]]  # list of logical vectors per marker
  
  # markers used in rules
  p21       <- get_marker(m, "p21", cells)
  p16       <- get_marker(m, "p16", cells)
  PCNA      <- get_marker(m, "PCNA", cells)
  IP10      <- get_marker(m, "IP10", cells)
  bcl2      <- get_marker(m, "BCL2", cells)
  MelanA    <- get_marker(m, "MelanA", cells)
  Vimentin  <- get_marker(m, "Vimentin", cells)
  Keratin5  <- get_marker(m, "Keratin 5", cells)
  Keratin10 <- get_marker(m, "Keratin 10", cells)
  Filaggrin <- get_marker(m, "Filaggrin", cells)
  Cxcr3 <- get_marker(m, "CXCR3", cells)
  
  obj$cdkn_state <- dplyr::case_when(
    p16 & !p21 ~ "p16+",
    p21 & !p16 ~ "p21+",
    p16 & p21  ~ "p16+p21+",
    TRUE       ~ "p16-p21-"
  )
  obj$cdkn_state <- factor(
    obj$cdkn_state,
    levels = c("p16-p21-", "p21+", "p16+", "p16+p21+")
  )
  
  # ---- Add metadata ----
  obj$snc <- ifelse(((p21 | p16) & (!PCNA) & (IP10)), "SnC+", "SnC-")
  #obj$snc2 <- ifelse(((p21 | p16) & (!PCNA) & (IP10)), "SnC+", "SnC-")
  
  obj$majorskincell <- dplyr::case_when(
    MelanA ~ "Melanocyte",
    Vimentin ~ "Fibroblast",
    (Keratin5 | Keratin10 | Filaggrin) ~ "Keratinocyte",
    TRUE ~ "none"
  )
  
  obj$minorskincell <- dplyr::case_when(
    MelanA ~ "Melanocyte",
    Vimentin ~ "Fibroblast",
    Keratin5 ~ "Basal Keratinocyte",
    Keratin10 ~ "Suprabasal Keratinocyte",
    Filaggrin ~ "Corneocyte",
    TRUE ~ "none"
  )
  
  # (optional) set factors for consistent plotting order
  obj$snc <- factor(obj$snc, levels = c("SnC-", "SnC+"))
  obj$snc2 <- factor(obj$snc2, levels = c("SnC-", "SnC+"))
  
  obj$majorskincell <- factor(obj$majorskincell,
                              levels = c("none", "Fibroblast", "Melanocyte", "Keratinocyte"))
  
  # Minor skin cells: Ordered Deep (Dermis) to Superficial (Stratum Corneum)
  obj$minorskincell <- factor(obj$minorskincell,
                              levels = c("none", 
                                         "Fibroblast", 
                                         "Melanocyte", 
                                         "Basal Keratinocyte", 
                                         "Suprabasal Keratinocyte", 
                                         "Corneocyte"))
  
  # Grab the spatial matrix for this specific object
  spatial_matrix <- obj@reductions$spatial@cell.embeddings
  
  # Assign the first column (X) and second column (Y) to the metadata
  # (Adjust the [,1] and [,2] if your matrix has different dimensions)
  obj$Spatial_X <- spatial_matrix[, 1]
  obj$Spatial_Y <- spatial_matrix[, 2]
  
  seurat_list[[sample_name]] <- obj
}

seurat_list <- lapply(seurat_list, function(obj) {
  
  # 1. Rename the Age values
  obj$age[obj$age == "old"] <- "Old"
  obj$age[obj$age == "young"] <- "Young"

  
  return(obj)
})

# Quick sanity check on one sample
table(seurat_list[[1]]$snc)
table(seurat_list[[1]]$majorskincell)
table(seurat_list[[1]]$minorskincell)

library(dplyr)
library(ggplot2)
library(scales)

blue_col <- "#296d9e"
snc_columns <- c("snc", "snc2") # The two metadata columns you want to plot

# Loop through both "snc" (with BCL2) and "snc2" (without BCL2)
for (col_name in snc_columns) {
  
  # 1. Extract raw counts for the specific column
  df_4i <- bind_rows(lapply(names(seurat_list), function(sample_name) {
    md <- seurat_list[[sample_name]]@meta.data
    
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
  df_pooled <- df_4i %>%
    group_by(age) %>%
    summarise(
      pooled_prop = sum(n_snc) / sum(n_total),            # True pooled fraction
      se = sd(snc_prop) / sqrt(n()),                      # SEM of the biological replicates
      .groups = "drop"
    )
  
  # 3. Dynamically set ymax and plot
  ymax <- max(c(df_4i$snc_prop, df_pooled$pooled_prop + df_pooled$se), na.rm = TRUE) * 1.15
  
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
      data = df_4i, aes(x = age, y = snc_prop),
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
  file_name <- paste0("results/260303/4i_", col_name, "_age.pdf")
  ggsave(filename = file_name, plot = p_age, width = 3.6, height = 4.2, dpi = 500)
  
  cat("Successfully saved:", file_name, "\n")
}

#######################################################################
# Fixing the snc assignemnt
#######################################################################
# 1. Transfer snc2 data to snc
obj_4i$snc <- obj_4i$snc2

# 2. Delete the redundant snc2 column
obj_4i$snc2 <- NULL

# 3. Update the seurat_list objects as well, so everything stays in sync
seurat_list <- lapply(seurat_list, function(x) {
  x$snc <- x$snc2
  x$snc2 <- NULL
  return(x)
})

# 4. If you have any dataframes in your global environment 
# created from obj_4i@meta.data (like df_4i), re-create them:
# This ensures they no longer contain the snc2 column.
df_4i <- obj_4i@meta.data
#######################################################################

##########################################################
#Subsenescence marker expression by SnC
##########################################################
install.packages("patchwork")
install.packages("devtools")

# Make sure you have ggpubr loaded for stat_compare_means
library(ggpubr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

plot_violin_marker <- function(obj,
                               feature,
                               group_col = "snc",
                               assay = "RNA",
                               layer = "data",
                               cols = c("SnC-" = "#296d9e", "SnC+" = "#cc0b01")) {
  
  df <- FetchData(obj, vars = c(feature, group_col), assay = assay, layer = layer)
  colnames(df) <- c("value", "group")
  
  df <- df %>%
    dplyr::filter(!is.na(group)) %>%
    dplyr::mutate(group = factor(group, levels = c("SnC-", "SnC+")))
  
  ymax <- max(df$value, na.rm = TRUE)
  
  p <- ggplot(df, aes(x = group, y = value, fill = group)) +
    
    geom_violin(trim = TRUE, scale = "width", alpha = 0.9) +
    
    stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.4,
      color = "black",
      linewidth = 0.8
    ) +
    
    stat_compare_means(
      method = "wilcox.test",
      label = "p.signif",
      label.x = 1.5,
      label.y = ymax * 0.98,
      size = 6
    ) +
    
    scale_fill_manual(values = cols) +
    
    theme_classic(base_size = 18) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    ) +
    
    labs(
      y = "Expression",
      title = toupper(feature)
    )
  
  return(p)
}

# Run the list exactly as you had it
features <- c("HMGB1", "LaminB1")

obj_4i <- merge(seurat_list[[1]], y = seurat_list[-1], merge.data = TRUE)

p_list <- lapply(features, function(f) plot_violin_marker(
  obj = obj_4i,
  feature = f,
  group_col = "snc",
  assay = "RNA",
  layer = "data"   
))

# Combine and print!
final_plot = p_list[[1]] + p_list[[2]]

final_plot
#ggsave("results/260225/4i_split_violins_subsenmarkers(cxcr3).png", plot = final_plot, width = 6, height = 3, dpi = 300)
ggsave("results/260303/4i_split_violins_subsenmarkers(cxcr3)bcl2.pdf", plot = final_plot, width = 6, height = 3, dpi = 300)

##########################################################
# SnC Marker Distribution
##########################################################
VlnPlot(
  obj_4i, 
  features = c("p21", "p16", "PCNA", "BCL2", "IP10"), 
  group.by = "snc", 
  pt.size = 0 # Turn off points if you have thousands of cells
)


##########################################################
# SnC Expression by Cell types
##########################################################
major_anatomical_order <- c(
  "none",
  "Fibroblast",    # Dermis (Bottom of plot)
  "Melanocyte",    # Basal layer
  "Keratinocyte"
)

minor_anatomical_order <- c(
  "none",
  "Fibroblast",    # Dermis (Bottom of plot)
  "Melanocyte",    # Basal layer
  "Basal Keratinocyte",         # Epidermis (Deep)
  "Suprabasal Keratinocyte",    # Epidermis (Mid)
  "Corneocyte"     # Epidermis (Surface / Top of plot)
)

#### calculate the frequency of SnC in each cell type by age

freq_major <- obj_4i@meta.data %>%
  filter(!is.na(age) & !is.na(majorskincell) & !is.na(snc)) %>%
  group_by(source_file, majorskincell, age) %>%
  summarise(
    Total_Cells = n(),
    SnC_Plus_Cells = sum(snc == "SnC+"),
    Freq_SnC_Plus = (SnC_Plus_Cells / Total_Cells) * 100,
    .groups = 'drop'
  ) %>%
  # Lock in your custom biological factor order
  mutate(majorskincell = factor(majorskincell, levels = major_anatomical_order))
  


freq_minor <- obj_4i@meta.data %>%
  # Filter out NAs to keep the data clean
  filter(!is.na(age) & !is.na(minorskincell) & !is.na(snc)) %>%
  group_by(source_file, minorskincell, age) %>%
  summarise(
    Total_Cells = n(),
    SnC_Plus_Cells = sum(snc == "SnC+"),
    Freq_SnC_Plus = (SnC_Plus_Cells / Total_Cells) * 100,
    .groups = 'drop'
  ) %>%
  # Lock in your custom biological factor order
  mutate(minorskincell = factor(minorskincell, levels = minor_anatomical_order)) %>%
  # Drop any minor cell types that aren't in your list (optional, but keeps plot clean)
  filter(!is.na(minorskincell))


major = ggplot(freq_major, aes(x = Freq_SnC_Plus, y = majorskincell)) +
  
  # Draw the connecting line between Young and Old
  geom_line(aes(group = majorskincell), color = "gray70", linewidth = 1.2) +
  
  # Add the dots on the ends
  geom_point(aes(color = age), size = 4.5) +
  
  # Colors (Assuming your age column has "Young" and "Old")
  scale_color_manual(values = c("Young" = "#296d9e", "Old" = "#cc0b01")) +
  theme_bw(base_size = 15) +
  labs(
    #title = "Senescence Shift by Major Skin Cell Type",
    x = "% SnC+ Cells",
    y = NULL, # Removes the redundant Y-axis label
    color = "Age"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
    axis.text.y = element_text(face = "bold", color = "black")
  )


minor = ggplot(freq_minor, aes(x = Freq_SnC_Plus, y = minorskincell)) +
  
  # The connecting line
  geom_line(aes(group = minorskincell), color = "gray70", linewidth = 1.2) +
  
  # The Young vs Old dots
  geom_point(aes(color = age), size = 4.5) +
  
  scale_color_manual(values = c("Young" = "#296d9e", "Old" = "#cc0b01")) +
  theme_bw(base_size = 15) +
  labs(
    #title = "Senescence Shift by Minor Cell State",
    #subtitle = "Ordered anatomically from Dermis to Skin Surface",
    x = "% SnC+ Cells",
    y = NULL, 
    color = "Age"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dashed", color = "gray90"),
    axis.text.y = element_text(face = "bold", color = "black")
  )

ggsave("results/260303/4i_barplot_snc_major.pdf", plot = major, width = 10, height = 3, dpi = 300)
ggsave("results/260303/4i_barplot_snc_minor.pdf", plot = minor, width = 10, height = 3, dpi = 300)


# Ensure Young first
freq_major <- freq_major %>%
  mutate(age = factor(age, levels = c("Young", "Old")))

p_major_bar <- ggplot(freq_major, aes(x = age, y = Freq_SnC_Plus)) +
  
  # mean bars
  stat_summary(
    fun = mean,
    geom = "bar",
    width = 0.6,
    fill = "#296d9e",
    color = "black"
  ) +
  
  # SEM error bars
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.15,
    linewidth = 0.8,
    color = "black"
  ) +
  
  # sample dots
  geom_point(
    position = position_jitter(width = 0.07, height = 0),
    size = 3,
    color = "black"
  ) +
  
  facet_wrap(~ majorskincell, ncol = 1, scales = "free_y") +
  
  labs(x = NULL, y = "% SnC+ Cells") +
  
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  )

print(p_major_bar)

freq_minor <- freq_minor %>%
  mutate(age = factor(age, levels = c("Young", "Old")))

p_minor_bar <- ggplot(freq_minor, aes(x = age, y = Freq_SnC_Plus)) +
  
  # mean bars
  stat_summary(
    fun = mean,
    geom = "bar",
    width = 0.6,
    fill = "#296d9e",
    color = "black"
  ) +
  
  # SEM error bars
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.15,
    linewidth = 0.8,
    color = "black"
  ) +
  
  # sample dots
  geom_point(
    position = position_jitter(width = 0.07, height = 0),
    size = 3,
    color = "black"
  ) +
  
  facet_wrap(~ minorskincell, ncol = 1, scales = "free_y") +
  
  labs(x = NULL, y = "% SnC+ Cells") +
  
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.8),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  )

print(p_minor_bar)

ggsave("results/260225/4i_barplot_snc_bcl2_major.pdf", plot = p_major_bar, width = 3, height = 11, dpi = 300)
ggsave("results/260225/4i_barplot_snc_bcl2_minor.pdf", plot = p_minor_bar, width = 3, height = 15, dpi = 300)




####################
# pie chart
#####################

minor_colors <- c(
  "none" = "grey80",
  "Fibroblast" = "#33a02c",
  "Melanocyte" = "#6a3d9a",
  "Basal Keratinocyte" = "#1f78b4",
  "Suprabasal Keratinocyte" = "#74add1",
  "Corneocyte" = "#fb9a99"
)

snc_comp <- obj_4i@meta.data %>%
  filter(snc == "SnC+" & !is.na(age) & !is.na(minorskincell)) %>%
  group_by(age, minorskincell) %>%
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

donut_one <- function(df, age_value) {
  ggplot(df %>% filter(age == age_value),
         aes(x = 2, y = pct, fill = minorskincell)) +
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

p <- (p_young | p_old) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("results/260617/4i_snc_donut_young_old_major.pdf", p, width = 7.5, height = 3.5)

snc_comp1 <- obj_4i@meta.data %>%
  filter(!is.na(age) & !is.na(minorskincell)) %>%
  group_by(age, minorskincell) %>%
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

donut_one <- function(df, age_value) {
  ggplot(df %>% filter(age == age_value),
         aes(x = 2, y = pct, fill = minorskincell)) +
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

p_young <- donut_one(snc_comp1, "Young")
p_old   <- donut_one(snc_comp1, "Old")

p <- (p_young | p_old) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave("results/260617/4i_total_donut_young_old_minor.pdf", p, width = 7.5, height = 3.5)

#########################
# p16/p21 distribution by SnC in young and old
#########################

cdkn_colors <- c(
  "p16-p21-" = "#d9d9d9",
  "p21+" = "#296d9e",
  "p16+" = "#cc0b01",
  "p16+p21+" = "#7f7f7f"
)

snc_marker_comp <- obj_4i@meta.data %>%
  dplyr::filter(snc == "SnC+") %>%
  dplyr::count(age, cdkn_state) %>%
  dplyr::group_by(age) %>%
  dplyr::mutate(pct = 100 * n / sum(n))

donut_marker <- function(df, age_value) {
  
  ggplot(df %>% dplyr::filter(age == age_value),
         aes(x = 2, y = pct, fill = cdkn_state)) +
    
    geom_col(width = 1, color = "white", linewidth = 1) +
    
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    
    scale_fill_manual(
      values = cdkn_colors,
      drop = FALSE
    ) +
    
    theme_void() +
    
    ggtitle(age_value) +
    
    theme(
      text = element_text(family = "Helvetica", size = 14),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10,10,10,10)
    )
}

library(patchwork)

p_young <- donut_marker(snc_marker_comp, "Young")
p_young <- p_young + theme(legend.position = "none")
p_old   <- donut_marker(snc_marker_comp, "Old")

p <- (p_young | p_old) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right") &
  scale_fill_manual(values = cdkn_colors, drop = FALSE)

ggsave("results/260303/4i_snc_donut_young_old_cdkn.pdf", p, width = 7.5, height = 3.5)




#########################
# SnC+ cell type distributions in young and old
#########################

library(ggplot2)
library(dplyr)
library(scales)
# ------------------------------------------------------------------
# 1. MASTER COLOR PALETTE & STACK ORDERS (Adapted for 4i)
# ------------------------------------------------------------------
elegant_colors <- c(
  "Keratinocyte"            = "#67001f",  
  "Basal Keratinocyte"      = "#b2182b",  
  "Corneocyte"              = "#fddbc7",  
  "Fibroblast"              = "#f7f7f7",  
  "Melanocyte"              = "#d1e5f0",  
  "Smooth Muscle"           = "#4393c3",  
  "Suprabasal Keratinocyte" = "#2166ac",  
  "none"                    = "#053061"   
)

# Lock in the vertical stack orders so legends and bars match perfectly
major_order <- c("Keratinocyte", "Melanocyte", "Fibroblast", "none")
minor_order <- c(
  "Suprabasal Keratinocyte", "Basal Keratinocyte", "Corneocyte", 
  "Melanocyte", "Fibroblast", "none"
)

# ------------------------------------------------------------------
# 2. MAJOR CELL TYPES PLOT (Out of Total Tissue)
# ------------------------------------------------------------------
df_major_total <- obj_4i@meta.data %>%
  filter(!is.na(age) & !is.na(majorskincell) & !is.na(snc)) %>%
  group_by(age) %>%
  mutate(total_cells_in_age = n()) %>%
  filter(snc == "SnC+") %>%
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
# 3. MINOR CELL TYPES PLOT (Out of Total Tissue)
# ------------------------------------------------------------------
df_minor_total <- obj_4i@meta.data %>%
  filter(!is.na(age) & !is.na(minorskincell) & !is.na(snc2)) %>%
  group_by(age) %>%
  mutate(total_cells_in_age = n()) %>%
  filter(snc2 == "SnC+") %>%
  group_by(age, minorskincell, total_cells_in_age) %>%
  summarise(snc_count = n(), .groups = "drop") %>%
  mutate(pct_of_total = snc_count / total_cells_in_age) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    minorskincell = factor(minorskincell, levels = minor_order)
  )

p_minor <- ggplot(df_minor_total, aes(x = age, y = pct_of_total, fill = minorskincell)) +
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

# ------------------------------------------------------------------
# 4. SAVE PLOTS (Tall aspect ratio applied)
# ------------------------------------------------------------------
ggsave("results/260303/4i_snc_total_tissue_major.pdf", plot = p_major, width = 5, height = 5)
ggsave("results/260303/4i_snc_total_tissue_minor(nobcl).pdf", plot = p_minor, width = 6, height = 5)

#########################
# SnC+ cell type distributions in young and old (normalized by total SnC)
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
  "none" = "grey80",
  "Fibroblast" = "#33a02c",
  "Melanocyte" = "#6a3d9a",
  "Basal Keratinocyte" = "#1f78b4",
  "Suprabasal Keratinocyte" = "#74add1",
  "Corneocyte" = "#fb9a99"
)

# --- MAJOR CELL TYPES (Composition out of 100% SnC+ cells) ---
df_major_comp <- obj_4i@meta.data %>%
  filter(snc2 == "SnC+" & !is.na(age) & !is.na(majorskincell)) %>%
  group_by(age, majorskincell) %>%
  summarise(cells = n(), .groups = "drop") %>%
  # Group by age to calculate the proportion out of total SnC+ cells in that age group
  group_by(age) %>%
  mutate(pct = cells / sum(cells)) %>%
  mutate(age = factor(age, levels = c("Young", "Old")))

p_major_comp <- ggplot(df_major_comp, aes(x = age, y = pct, fill = majorskincell)) +
  geom_col(width = 0.6, color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = major_colors) +
  theme_classic(base_size = 16) +
  labs(
    #title = "SnC+ Cell Type Composition (Major)", 
    x = "Age", 
    y = "Proportion of SnC+ Cells", 
    fill = "Cell Type"
  )

print(p_major_comp)


# --- MINOR CELL TYPES (Composition out of 100% SnC+ cells) ---
df_minor_comp <- obj_4i@meta.data %>%
  filter(snc2 == "SnC+" & !is.na(age) & !is.na(minorskincell)) %>%
  group_by(age, minorskincell) %>%
  summarise(cells = n(), .groups = "drop") %>%
  group_by(age) %>%
  mutate(pct = cells / sum(cells)) %>%
  mutate(age = factor(age, levels = c("Young", "Old")))

p_minor_comp <- ggplot(df_minor_comp, aes(x = age, y = pct, fill = minorskincell)) +
  geom_col(width = 0.6, color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = minor_colors) +
  theme_classic(base_size = 16) +
  labs(
    #title = "SnC+ Cell Type Composition (Minor)", 
    x = "Age", 
    y = "Proportion of SnC+ Cells", 
    fill = "Cell State"
  )

print(p_minor_comp)

# Save the plots
ggsave("results/260303/4i_snc_composition_stacked_major.pdf", plot = p_major_comp, width = 5, height = 5, dpi = 300)
ggsave("results/260303/4i_snc_composition_stacked_minor.pdf", plot = p_minor_comp, width = 6, height = 5, dpi = 300)



#########################
# SnC+ cell type distributions in young and old (normalized by total skin cells)
#########################


library(ggplot2)
library(dplyr)
library(scales)

elegant_colors <- c(
  "Keratinocyte"            = "#67001f",  
  "Basal Keratinocyte"      = "#b2182b",  
  "Corneocyte"              = "#fddbc7",  
  "Fibroblast"              = "#f7f7f7",  
  "Melanocyte"              = "#d1e5f0",  
  "Suprabasal Keratinocyte" = "#2166ac",  
  "none"                    = "lightgray"   
)

major_colors <- c(
  "none" = "grey80",
  "Fibroblast" = "#33a02c",
  "Melanocyte" = "#6a3d9a",
  "Keratinocyte" = "#1f78b4"
)

minor_colors <- c(
  "none" = "grey80",
  "Fibroblast" = "#33a02c",
  "Melanocyte" = "#6a3d9a",
  "Basal Keratinocyte" = "#1f78b4",
  "Suprabasal Keratinocyte" = "#74add1",
  "Corneocyte" = "#fb9a99"
)

major_order <- c("Fibroblast", "Melanocyte", "Keratinocyte", "none")
minor_order <- c(
  "Fibroblast",
  "Melanocyte",
  "Basal Keratinocyte",
  "Suprabasal Keratinocyte",
  "Corneocyte",
  "none"
)

# ------------------------------------------------------------------
# 2. MAJOR CELL TYPES PLOT (4i)
# ------------------------------------------------------------------
df_major_total_4i <- obj_4i@meta.data %>%
  filter(!is.na(age) & !is.na(majorskincell) & !is.na(snc2)) %>%
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

p_major_total <- ggplot(df_major_total_4i, aes(x = age, y = pct_of_total, fill = majorskincell)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.7) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
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

print(p_major_total)

# ------------------------------------------------------------------
# 3. MINOR CELL TYPES PLOT (4i)
# ------------------------------------------------------------------
df_minor_total_4i <- obj_4i@meta.data %>%
  filter(!is.na(age) & !is.na(minorskincell) & !is.na(snc2)) %>%
  group_by(age) %>%
  mutate(total_cells_in_age = n()) %>%
  filter(snc2 == "SnC+") %>%
  group_by(age, minorskincell, total_cells_in_age) %>%
  summarise(snc_count = n(), .groups = "drop") %>%
  mutate(pct_of_total = snc_count / total_cells_in_age) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    minorskincell = factor(minorskincell, levels = minor_order)
  )

p_minor_total <- ggplot(df_minor_total_4i, aes(x = age, y = pct_of_total, fill = minorskincell)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.7) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.05))
  ) +
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

print(p_minor_total)

ggsave("results/260303/4i_snc_skin_stacked_major.pdf", plot = p_major_total, width = 5, height = 5, dpi = 300)
ggsave("results/260303/4i_snc_skin_stacked_minor.pdf", plot = p_minor_total, width = 6, height = 5, dpi = 300)











library(ggplot2)
library(dplyr)
library(scales)

# 1. Define specific colors for Age
age_colors <- c("Young" = "#9ecae1", "Old" = "#08519c")

# ------------------------------------------------------------------
# 2. MAJOR CELL TYPES: Lineage-Specific SnC+ Burden
# ------------------------------------------------------------------
df_major_intrinsic <- obj_4i@meta.data %>%
  filter(!is.na(age) & !is.na(majorskincell) & !is.na(snc2)) %>%
  filter(majorskincell != "none") %>% # Drop unclassified cells
  # Group by cell type AND age
  group_by(majorskincell, age) %>%
  summarise(
    total_lineage_cells = n(),
    snc_cells = sum(snc2 == "SnC+"),
    pct_snc = snc_cells / total_lineage_cells, # The crucial math shift
    .groups = "drop"
  ) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    majorskincell = factor(majorskincell, levels = c("Fibroblast", "Melanocyte", "Keratinocyte"))
  )

p_major_intrinsic <- ggplot(df_major_intrinsic, aes(x = majorskincell, y = pct_snc, fill = age)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.5) +
  scale_y_continuous(
    labels = percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_manual(values = age_colors) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(
    axis.text.x = element_text(color = "black", face = "bold", margin = margin(t = 5)),
    axis.text.y = element_text(color = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    #title = "Intrinsic SnC+ Burden (Major Lineages)",
    x = NULL, 
    y = NULL
  )

print(p_major_intrinsic)

# ------------------------------------------------------------------
# 3. MINOR CELL TYPES: Lineage-Specific SnC+ Burden
# ------------------------------------------------------------------
df_minor_intrinsic <- obj_4i@meta.data %>%
  filter(!is.na(age) & !is.na(minorskincell) & !is.na(snc2)) %>%
  filter(minorskincell != "none") %>%
  group_by(minorskincell, age) %>%
  summarise(
    total_lineage_cells = n(),
    snc_cells = sum(snc2 == "SnC+"),
    pct_snc = snc_cells / total_lineage_cells,
    .groups = "drop"
  ) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    minorskincell = factor(minorskincell, levels = c("Fibroblast", "Melanocyte", "Basal Keratinocyte", "Suprabasal Keratinocyte", "Corneocyte"))
  )

p_minor_intrinsic <- ggplot(df_minor_intrinsic, aes(x = minorskincell, y = pct_snc, fill = age)) +
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
    #title = "Intrinsic SnC+ Burden (Minor States)",
    x = NULL, 
    y = NULL
  )

print(p_minor_intrinsic)

ggsave("results/260303/4i_snc_intrinsic_major.pdf", plot = p_major_intrinsic, width = 6, height = 5, dpi = 300)
ggsave("results/260303/4i_snc_intrinsic_minor.pdf", plot = p_minor_intrinsic, width = 7, height = 5, dpi = 300)












save.image("data/4i/seurat_list_4i_annotated.RData")

#########################
# Clustering
#########################
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 1. Normalize and identify variable features (if not already done)
obj_4i <- NormalizeData(obj_4i, assay = "RNA")
obj_4i <- FindVariableFeatures(obj_4i, selection.method = "vst", nfeatures = 2000)

# 2. Scale the data (Important: Scale all genes if you plan to plot heatmaps later)
all.genes <- rownames(obj_4i)
obj_4i <- ScaleData(obj_4i, features = all.genes)

# 3. Run PCA
obj_4i <- RunPCA(obj_4i, features = VariableFeatures(object = obj_4i), verbose = FALSE)

# (Optional Check: Look at the elbow plot to decide how many PCs to use. Usually 15-30 is good)
ElbowPlot(obj_4i)

# 4. Find Neighbors and Cluster the cells
# We will use 20 PCs here, adjust if your ElbowPlot suggests otherwise.
obj_4i <- FindNeighbors(obj_4i, dims = 1:5)

# The 'resolution' parameter dictates how granular the clusters are. 
# 0.4 - 0.8 is standard. Higher = more clusters.
obj_4i <- FindClusters(obj_4i, resolution = 0.5)

# 5. Run UMAP for 2D visualization
obj_4i <- RunUMAP(obj_4i, dims = 1:5)

# Plot 1: Standard Seurat Clusters
p1 <- DimPlot(obj_4i, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("Seurat Clusters")

# Plot 2: Major Cell Types overlaid on UMAP
p2 <- DimPlot(obj_4i, reduction = "umap", group.by = "majorskincell", pt.size = 0.5) + 
  ggtitle("Major Cell Types")

# Plot 3: Highlight the SnC+ cells
# We use cols to make SnC- gray and SnC+ stand out in red
p3 <- DimPlot(obj_4i, reduction = "umap", group.by = "snc", 
              cols = c("SnC-" = "gray80", "SnC+" = "red"), pt.size = 0.5) + 
  ggtitle("Senescent Cells (SnC+)")

# Plot 4: Feature plot of a continuous module score (e.g., SASP score or p21)
# Assuming you still have 'sasp_score' or 'p21' in the metadata/assay
p4 <- FeaturePlot(obj_4i, features = c("p21", "p16"), ncol = 2) 

# Display them together
(p1 | p2) / (p3 | p4)

