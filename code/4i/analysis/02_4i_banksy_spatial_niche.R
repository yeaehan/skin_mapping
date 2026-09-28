install.packages("R.utils")


remotes::install_github("prabhakarlab/Banksy")
remotes::install_github("satijalab/seurat-wrappers")

##########################################################
# BANKSY domain analysis for 4i Seurat objects
# - Runs BANKSY per sample first
# - Uses spatial coordinates from metadata
# - Uses a curated 4i marker panel
# - Stores BANKSY domains back into each Seurat object
# - Produces spatial domain plots and domain summaries
##########################################################

##########################################################
# Load packages
##########################################################
library(Seurat)
library(SeuratWrappers)
library(Banksy)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(purrr)
library(Matrix)

##########################################################
# Set working directory and load data
##########################################################
setwd("//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2")

load("data/4i/260225.RData")
load("data/4i/seurat_list_4i.RData")
load("data/4i/seurat_list_4i_annotated.RData")
load("data/4i/banksy_merged_4i_annotated.RData")



##########################################################
# Define 4i marker panel for BANKSY
# Keep this fairly compact and biologically interpretable
##########################################################
marker_panel <- c(
  "p21", "p16", "PCNA", "HMGB1", "LaminB1", "BCL2", "IP10", "CXCR3",
  "Keratin 5", "Keratin 10", "MelanA", "Filaggrin", "Vimentin"
)

##########################################################
# Helper: pick usable BANKSY features from the marker panel
# Filters out very sparse markers within each sample
##########################################################
# pick_banksy_features_4i <- function(
#     obj,
#     assay = "RNA",
#     marker_panel = marker_panel,
#     min_pct = 0.02,
#     min_var = 0
# ) {
#   DefaultAssay(obj) <- assay
#   
#   if (ncol(GetAssayData(obj, assay = assay, layer = "data")) == 0) {
#     obj <- NormalizeData(obj, assay = assay, verbose = FALSE)
#   }
#   
#   mat <- GetAssayData(obj, assay = assay, layer = "data")
#   feats <- intersect(marker_panel, rownames(mat))
#   
#   if (length(feats) == 0) {
#     warning("No BANKSY markers found in this object.")
#     return(character(0))
#   }
#   
#   pct_expr <- Matrix::rowMeans(mat[feats, , drop = FALSE] > 0)
#   vars <- apply(as.matrix(mat[feats, , drop = FALSE]), 1, var)
#   
#   keep <- feats[pct_expr >= min_pct & vars > min_var]
#   keep
# }

##########################################################
# Helper: run BANKSY on one sample
# IMPORTANT: do not run ScaleData after RunBanksy
##########################################################
run_banksy_one <- function(
    obj,
    sample_name = NULL,
    assay = "RNA",
    features = NULL,
    lambda = 0.8,
    k_geom = 12,
    resolution = 0.35,
    npcs = 10,
    dimx = "Spatial_X",
    dimy = "Spatial_Y",
    seed = 1234
) {
  DefaultAssay(obj) <- assay
  
  if (ncol(GetAssayData(obj, assay = assay, layer = "data")) == 0) {
    obj <- NormalizeData(obj, assay = assay, verbose = FALSE)
  }
  
  if (is.null(features) || length(features) == 0) {
    stop("No features supplied to run_banksy_one().")
  }
  
  set.seed(seed)
  
  obj <- RunBanksy(
    object = obj,
    lambda = lambda,
    assay = assay,
    slot = "data",
    dimx = dimx,
    dimy = dimy,
    features = features,
    k_geom = k_geom,
    verbose = TRUE
  )
  
  DefaultAssay(obj) <- "BANKSY"
  
  obj <- RunPCA(
    object = obj,
    assay = "BANKSY",
    reduction.name = "pca.banksy",
    features = rownames(obj[["BANKSY"]]),
    npcs = npcs,
    verbose = FALSE
  )
  
  obj <- FindNeighbors(
    object = obj,
    reduction = "pca.banksy",
    dims = 1:npcs,
    verbose = FALSE
  )
  
  domain_col <- paste0("banksy_domain_lam", gsub("\\.", "_", lambda), "_k", k_geom)
  
  obj <- FindClusters(
    object = obj,
    resolution = resolution,
    cluster.name = domain_col,
    verbose = FALSE
  )
  
  obj <- RunUMAP(
    object = obj,
    reduction = "pca.banksy",
    dims = 1:npcs,
    reduction.name = "umap.banksy",
    verbose = FALSE
  )
  
  if (!is.null(sample_name)) {
    obj$banksy_sample <- sample_name
  }
  
  obj
}

##########################################################
# Helper: spatial plot for BANKSY domains
##########################################################
plot_banksy_domains <- function(
    obj,
    domain_col = NULL,
    dimx = "Spatial_X",
    dimy = "Spatial_Y",
    pt_size = 0.4,
    title = NULL
) {
  if (is.null(domain_col)) {
    domain_col <- grep("^banksy_domain", colnames(obj@meta.data), value = TRUE)[1]
  }
  
  df <- obj@meta.data %>%
    mutate(
      x_plot = .data[[dimx]],
      y_plot = .data[[dimy]],
      domain_plot = .data[[domain_col]]
    )
  
  ggplot(df, aes(x = x_plot, y = y_plot, color = domain_plot)) +
    geom_point(size = pt_size) +
    scale_y_reverse() +
    coord_fixed() +
    theme_void() +
    labs(title = title, color = "BANKSY domain") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

##########################################################
# Helper: summarize BANKSY domains by biology
##########################################################
summarize_banksy_domains <- function(obj, domain_col = NULL) {
  if (is.null(domain_col)) {
    domain_col <- grep("^banksy_domain", colnames(obj@meta.data), value = TRUE)[1]
  }
  
  obj@meta.data %>%
    mutate(domain = .data[[domain_col]]) %>%
    group_by(domain) %>%
    summarise(
      n_cells = n(),
      pct_snc = mean(snc == "SnC+", na.rm = TRUE) * 100,
      pct_fib = mean(minorskincell == "Fibroblast", na.rm = TRUE) * 100,
      pct_mel = mean(minorskincell == "Melanocyte", na.rm = TRUE) * 100,
      pct_basal = mean(minorskincell == "Basal Keratinocyte", na.rm = TRUE) * 100,
      pct_supra = mean(minorskincell == "Suprabasal Keratinocyte", na.rm = TRUE) * 100,
      pct_corn = mean(minorskincell == "Corneocyte", na.rm = TRUE) * 100,
      pct_none = mean(minorskincell == "none", na.rm = TRUE) * 100,
      pct_young = mean(age == "Young", na.rm = TRUE) * 100,
      pct_old = mean(age == "Old", na.rm = TRUE) * 100,
      .groups = "drop"
    )
}

##########################################################
# Helper: plot domain composition heatmap
##########################################################
plot_domain_celltype_heatmap <- function(obj, domain_col = NULL) {
  if (is.null(domain_col)) {
    domain_col <- grep("^banksy_domain", colnames(obj@meta.data), value = TRUE)[1]
  }
  
  df <- obj@meta.data %>%
    mutate(domain = .data[[domain_col]]) %>%
    count(domain, minorskincell) %>%
    group_by(domain) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ungroup()
  
  ggplot(df, aes(x = minorskincell, y = domain, fill = pct)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(name = "% of domain") +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = "BANKSY domain") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

# ##########################################################
# # Pick BANKSY features per sample
# ##########################################################
# banksy_features_by_sample <- lapply(seurat_list, function(obj) {
#   pick_banksy_features_4i(
#     obj = obj,
#     assay = "RNA",
#     marker_panel = marker_panel,
#     min_pct = 0.02,
#     min_var = 0
#   )
# })
# 
# ##########################################################
# # Review feature availability per sample
# ##########################################################
# banksy_features_by_sample

##########################################################
# Run BANKSY per sample
# Recommended starting values for domain segmentation:
# lambda = 0.8
# k_geom = 12
##########################################################
banksy_results_4i <- list()

for (sample_name in names(seurat_list)) {
  message("Running BANKSY for sample: ", sample_name)
  
  obj <- seurat_list[[sample_name]]
  
  # Force the use of your entire marker panel.
  # intersect() simply ensures the gene name actually exists in the object 
  # to prevent Seurat from throwing a "feature not found" error.
  feats <- intersect(marker_panel, rownames(obj))
  
  # Quick check to make sure your crucial markers didn't get lost
  if (!all(c("p16", "p21") %in% feats)) {
    warning("Sample ", sample_name, " is missing p16 or p21 in its count matrix!")
  }
  
  if (length(feats) < 4) {
    warning("Skipping ", sample_name, " because it has fewer than 4 markers.")
    next
  }
  
  obj_b <- run_banksy_one(
    obj = obj,
    sample_name = sample_name,
    assay = "RNA",
    features = feats,   # <-- Now strictly using your forced panel!
    lambda = 0.8,
    k_geom = 12,
    resolution = 0.35,
    npcs = min(10, length(feats)),
    dimx = "Spatial_X",
    dimy = "Spatial_Y",
    seed = 1234
  )
  
  banksy_results_4i[[sample_name]] <- obj_b
}

##########################################################
# Save BANKSY results list
##########################################################
save(banksy_results_4i, file = "data/4i/banksy_results_4i_per_sample.RData")

##########################################################
# Plot BANKSY domains for each sample and save
##########################################################
dir.create("results/260303/banksy_4i", recursive = TRUE, showWarnings = FALSE)

for (sample_name in names(banksy_results_4i)) {
  obj <- banksy_results_4i[[sample_name]]
  domain_col <- grep("^banksy_domain", colnames(obj@meta.data), value = TRUE)[1]
  
  p_dom <- plot_banksy_domains(
    obj = obj,
    domain_col = domain_col,
    dimx = "Spatial_X",
    dimy = "Spatial_Y",
    pt_size = 0.4,
    title = paste0(sample_name, " | ", domain_col)
  )
  
  ggsave(
    filename = file.path("results/260303/banksy_4i", paste0(sample_name, "_", domain_col, "_spatial.pdf")),
    plot = p_dom,
    width = 4.5,
    height = 4.5,
    dpi = 300
  )
  
  ggsave(
    filename = file.path("results/260303/banksy_4i", paste0(sample_name, "_", domain_col, "_spatial.png")),
    plot = p_dom,
    width = 4.5,
    height = 4.5,
    dpi = 300
  )
}

##########################################################
# Plot BANKSY domain composition heatmap for each sample
##########################################################
for (sample_name in names(banksy_results_4i)) {
  obj <- banksy_results_4i[[sample_name]]
  domain_col <- grep("^banksy_domain", colnames(obj@meta.data), value = TRUE)[1]
  
  p_heat <- plot_domain_celltype_heatmap(obj, domain_col = domain_col)
  
  ggsave(
    filename = file.path("results/260303/banksy_4i", paste0(sample_name, "_", domain_col, "_celltype_heatmap.pdf")),
    plot = p_heat,
    width = 5.5,
    height = 3.5,
    dpi = 300
  )
}

##########################################################
# Create per-sample BANKSY domain summaries
##########################################################
banksy_domain_summaries <- lapply(names(banksy_results_4i), function(sample_name) {
  obj <- banksy_results_4i[[sample_name]]
  domain_col <- grep("^banksy_domain", colnames(obj@meta.data), value = TRUE)[1]
  
  summarize_banksy_domains(obj, domain_col = domain_col) %>%
    mutate(sample = sample_name, domain_col = domain_col)
})

banksy_domain_summaries <- bind_rows(banksy_domain_summaries)

##########################################################
# Save domain summary table
##########################################################
write.csv(
  banksy_domain_summaries,
  file = "results/260303/banksy_4i/banksy_domain_summaries_4i.csv",
  row.names = FALSE
)

##########################################################
# Quick view of domain summaries
##########################################################
banksy_domain_summaries

##########################################################
# Optional: parameter sweep on one representative sample
# Useful before deciding final lambda/k_geom
##########################################################
sample_for_tuning <- names(seurat_list)[1]
obj_tune <- seurat_list[[sample_for_tuning]]
feats_tune <- banksy_features_by_sample[[sample_for_tuning]]

param_grid <- expand.grid(
  lambda = c(0.2, 0.5, 0.8),
  k_geom = c(8, 12, 16),
  stringsAsFactors = FALSE
)

banksy_tuning_results <- list()

for (i in seq_len(nrow(param_grid))) {
  lam_i <- param_grid$lambda[i]
  kg_i  <- param_grid$k_geom[i]
  nm_i  <- paste0("lam", lam_i, "_k", kg_i)
  
  message("Tuning BANKSY: ", nm_i)
  
  obj_i <- run_banksy_one(
    obj = obj_tune,
    sample_name = sample_for_tuning,
    assay = "RNA",
    features = feats_tune,
    lambda = lam_i,
    k_geom = kg_i,
    resolution = 0.35,
    npcs = min(10, length(feats_tune)),
    dimx = "Spatial_X",
    dimy = "Spatial_Y",
    seed = 1234
  )
  
  banksy_tuning_results[[nm_i]] <- obj_i
}

##########################################################
# Save parameter sweep results
##########################################################
save(banksy_tuning_results, file = "data/4i/banksy_tuning_results_4i.RData")

##########################################################
# Plot parameter sweep results for visual comparison
##########################################################
dir.create("results/260303/banksy_4i/tuning", recursive = TRUE, showWarnings = FALSE)

for (nm in names(banksy_tuning_results)) {
  obj <- banksy_tuning_results[[nm]]
  domain_col <- grep("^banksy_domain", colnames(obj@meta.data), value = TRUE)[1]
  
  p_tune <- plot_banksy_domains(
    obj = obj,
    domain_col = domain_col,
    dimx = "Spatial_X",
    dimy = "Spatial_Y",
    pt_size = 0.4,
    title = paste0(sample_for_tuning, " | ", nm)
  )
  
  ggsave(
    filename = file.path("results/260303/banksy_4i/tuning", paste0(sample_for_tuning, "_", nm, ".pdf")),
    plot = p_tune,
    width = 4.5,
    height = 4.5,
    dpi = 300
  )
}

##########################################################
# Optional: merge all BANKSY-annotated samples after per-sample run
# This is useful only for downstream summary/plotting, not for initial discovery
##########################################################
banksy_merged_4i <- merge(
  x = banksy_results_4i[[1]],
  y = banksy_results_4i[-1],
  merge.data = TRUE
)

banksy_merged_4i$age <- factor(banksy_merged_4i$age, levels = c("Young", "Old"))

banksy_merged_4i <- ScaleData(
  object = banksy_merged_4i, 
  assay = "BANKSY",
  verbose = FALSE
)

DefaultAssay(banksy_merged_4i) <- "BANKSY"

# 4. Run a single, global PCA across all cells from all samples
banksy_merged_4i <- RunPCA(
  object = banksy_merged_4i,
  assay = "BANKSY",
  features = rownames(banksy_merged_4i[["BANKSY"]]),
  npcs = 10,
  reduction.name = "pca_joint",
  verbose = FALSE
)

# 5. Build the joint neighborhood graph
banksy_merged_4i <- FindNeighbors(
  object = banksy_merged_4i,
  reduction = "pca_joint",
  dims = 1:10,
  verbose = FALSE
)

# 6. Identify the universal domains 
# (We name the column explicitly so it doesn't get confused with the single-sample runs)
banksy_merged_4i <- FindClusters(
  object = banksy_merged_4i,
  resolution = 0.03, #0.01 would be four niches
  cluster.name = "joint_banksy_domain",
  verbose = FALSE
)

# 7. Run a global UMAP to visualize how all samples overlap
banksy_merged_4i <- RunUMAP(
  object = banksy_merged_4i,
  reduction = "pca_joint",
  dims = 1:10,
  reduction.name = "umap_joint",
  verbose = FALSE
)

##########################################################
# Save merged annotated object
##########################################################
save(banksy_merged_4i, file = "data/4i/banksy_merged_4i_annotated.RData")


##########################################################
# Composite Analysis
##########################################################


library(dplyr)
library(ggplot2)
library(dplyr)
library(ggplot2)

# Make sure the results folder exists
dir.create("results/260303/banksy_4i", recursive = TRUE, showWarnings = FALSE)

### ------------------------------------------------------------------------
### 1. SAVE THE COMPOSITION PLOTS
### ------------------------------------------------------------------------

# Plot 1: CDKN State Composition
joint_comp_cdkn <- banksy_merged_4i@meta.data %>%
  count(joint_banksy_domain, age, cdkn_state) %>%
  group_by(joint_banksy_domain, age) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p_cdkn <- ggplot(joint_comp_cdkn, aes(x = joint_banksy_domain, y = pct, fill = cdkn_state)) +
  geom_col(position = "stack") +
  facet_wrap(~ age, ncol = 1) + 
  theme_bw(base_size = 14) +
  labs(
    #title = "CDKN State Composition of Joint Spatial Domains",
    x = "Joint BANKSY Domain",
    y = "% of Cells in Domain",
    fill = "CDKN State"
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )

ggsave("results/260303/banksy_4i/joint_composition_CDKN.pdf", plot = p_cdkn, width = 6, height = 8)


# Plot 2: Minor Skin Cell Composition
joint_comp_minor <- banksy_merged_4i@meta.data %>%
  count(joint_banksy_domain, age, minorskincell) %>%
  group_by(joint_banksy_domain, age) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p_minor <- ggplot(joint_comp_minor, aes(x = joint_banksy_domain, y = pct, fill = minorskincell)) +
  geom_col(position = "stack") +
  facet_wrap(~ age, ncol = 1) + 
  theme_bw(base_size = 14) +
  labs(
    #title = "Cell Type Composition of Joint Spatial Domains",
    x = "Joint BANKSY Domain",
    y = "% of Cells in Domain",
    fill = "Cell Type"   # <-- Fixed your legend title here!
  ) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    panel.grid = element_blank()
  )

ggsave("results/260303/banksy_4i/joint_composition_minorskincell.pdf", plot = p_minor, width = 6, height = 8)


### ------------------------------------------------------------------------
### 2. PLOT THE SPATIAL NICHES IN EVERY SAMPLE
### ------------------------------------------------------------------------

niche_colors <- c(
  "0" = "#E69F00", # Orange
  "1" = "#56B4E9", # Light Blue
  "2" = "#009E73", # Green
  "3" = "#CC79A7", # Pink/Magenta
  "4" = "#0072B2", # Dark Blue
  "5" = "#D55E00", # Vermilion/Red
  "6" = "#F0E442"  # Yellow
)

# 2. Grab the list of unique samples
unique_samples <- unique(banksy_merged_4i$banksy_sample)

banksy_merged_4i$joint_banksy_domain <- factor(banksy_merged_4i$joint_banksy_domain, levels = c("0", "1", "2", "3", "4", "5", "6"))

for (sample_name in unique_samples) {
  
  df_plot <- banksy_merged_4i@meta.data %>%
    filter(banksy_sample == sample_name) %>%
    arrange(joint_banksy_domain) # <-- THIS SORTS THE ROWS 0 TO 6!
  
  sample_age <- df_plot$age[1]
  
  p_spatial <- ggplot(df_plot, aes(x = Spatial_X, y = Spatial_Y, color = joint_banksy_domain)) +
    geom_point(size = 0.6) +
    scale_y_reverse() +    
    coord_fixed() +        
    scale_color_manual(
      values = niche_colors, 
      drop = FALSE,
      breaks = c("0", "1", "2", "3", "4", "5", "6") # Forces legend order
    ) + 
    theme_void() +         
    labs(
      color = "Joint Domain"
    ) +
    theme(
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path("results/260303/banksy_4i", paste0(sample_name, "_realigned_spatial.pdf")),
    plot = p_spatial,
    width = 5, 
    height = 5, 
    dpi = 300
  )
}

# 
# # 3. Run the plotting loop
# for (sample_name in unique_samples) {
#   
#   df_plot <- banksy_merged_4i@meta.data %>%
#     filter(banksy_sample == sample_name)
#   
#   sample_age <- df_plot$age[1]
#   
#   p_spatial <- ggplot(df_plot, aes(x = Spatial_X, y = Spatial_Y, color = joint_banksy_domain)) +
#     geom_point(size = 0.6) +
#     scale_y_reverse() +    
#     coord_fixed() +        
#     # Apply your hardcoded colors and force the legend to show all 7
#     scale_color_manual(values = niche_colors, drop = FALSE, breaks = c("0", "1", "2", "3", "4", "5", "6")) + 
#     theme_void() +         
#     labs(
#       title = paste0(sample_name, " (", sample_age, ")"),
#       color = "Joint Domain"
#     ) +
#     theme(
#       plot.title = element_text(hjust = 0.5, face = "bold"),
#       legend.position = "right"
#     )
#   
#   ggsave(
#     filename = file.path("results/260303/banksy_4i", paste0(sample_name, "_joint_spatial_map.pdf")),
#     plot = p_spatial,
#     width = 5, 
#     height = 5, 
#     dpi = 300
#   )
# }
















##########################################################
# Session save
##########################################################
save.image("data/4i/260303_banksy_4i_workspace.RData")