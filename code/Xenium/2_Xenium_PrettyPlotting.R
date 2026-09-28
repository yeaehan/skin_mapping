library(Seurat)
library(ggplot2)
library(scales)
library(colorspace)
library(patchwork)

# ---------------------------------------
# Helpers
# ---------------------------------------

# Usage in your ggplot code:
# + scale_color_manual(values = minor_colors)

make_pastel_discrete <- function(vals) {
  vals <- vals[!is.na(vals)]
  vals <- unique(as.character(vals))
  base_cols <- hue_pal()(length(vals))   # ggplot-like default hues
  pastel_cols <- colorspace::lighten(base_cols, amount = 0.35)
  names(pastel_cols) <- vals
  pastel_cols
}

# make sure UMAP exists
if (!"umap" %in% names(hsskin.ST@reductions)) {
  stop("UMAP reduction not found in hsskin.ST")
}

# enforce age order
hsskin.ST$age <- factor(hsskin.ST$age, levels = c("Young", "Old"))

# palettes
major_cols <- make_pastel_discrete(hsskin.ST$cell_type_major)
minor_cols <- make_pastel_discrete(hsskin.ST$cell_type_minor)
age_cols   <- c("Young" = "#9ecae1", "Old" = "#f4a3a8")

# ---------------------------------------
# UMAP: major cell types
# ---------------------------------------
p_major <- DimPlot(
  hsskin.ST,
  reduction = "umap",
  group.by = "cell_type_major",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  scale_color_manual(values = major_colors, na.value = "grey80") +
  labs(title = "Cell type major") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  )

# ---------------------------------------
# UMAP: minor cell types
# ---------------------------------------
p_minor <- DimPlot(
  hsskin.ST,
  reduction = "umap",
  group.by = "cell_type_minor",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  scale_color_manual(values = minor_cols, na.value = "grey80") +
  labs(title = "Cell type minor") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  )

# ---------------------------------------
# UMAP: age
# ---------------------------------------
p_age <- DimPlot(
  hsskin.ST,
  reduction = "umap",
  group.by = "age",
  raster = FALSE
) +
  scale_color_manual(values = age_cols, na.value = "grey80") +
  labs(title = "Age") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  )

# print individually
p_major
p_minor
p_age

# combined
p_umap_all <- (p_major | p_minor | p_age)
p_umap_all

# save
ggsave("skin_mapping/results/260506/xenium_umap_celltype_major.pdf", p_major, width = 8, height = 5, dpi = 300)
ggsave("skin_mapping/results/260225/xenium_umap_celltype_minor_pastel.pdf", p_minor, width = 8, height = 5, dpi = 300)
ggsave("skin_mapping/results/260225/xenium_umap_age.pdf", p_age, width = 6, height = 5, dpi = 300)
ggsave("skin_mapping/results/260225/xenium_umap_major_minor_age_combined.pdf", p_umap_all, width = 22, height = 5, dpi = 300)


# ---------------------------------------
# 2) extract coordinates from all FOVs
# ---------------------------------------
fov_names <- Images(hsskin)

extract_fov_coords <- function(obj, fov_name) {
  coords <- GetTissueCoordinates(obj[[fov_name]])
  coords <- as.data.frame(coords)
  coords$cell_id <- rownames(coords)
  coords$fov <- fov_name
  
  # robustly find x/y columns
  numeric_cols <- names(coords)[sapply(coords, is.numeric)]
  if (length(numeric_cols) < 2) {
    stop(paste("Could not find at least two numeric coordinate columns for", fov_name))
  }
  
  coords <- coords %>%
    rename(
      x = all_of(numeric_cols[1]),
      y = all_of(numeric_cols[2])
    )
  
  coords
}

coords_all <- bind_rows(lapply(fov_names, function(f) extract_fov_coords(hsskin, f)))

coords_all <- coords_all %>%
  select(-cell_id) %>%
  rename(cell_id = cell)    # keep real cell names
# ---------------------------------------
# 3) join coordinates with metadata
# ---------------------------------------
md <- hsskin@meta.data %>%
  mutate(cell_id = rownames(.))

plot_df <- coords_all %>%
  left_join(md, by = "cell_id")

# sanity check
# head(plot_df)
# table(plot_df$sample, useNA = "ifany")

# ---------------------------------------
# 4) custom ggplot function
# ---------------------------------------
library(ggplot2)
library(dplyr)
library(patchwork)

minor_colors <- c(
  # --- Keratinocyte / Epidermal / Epithelial Family (Browns, Tans, Reds) ---
  "Basal Keratinocytes"           = "#dfc299", # Golden Tan
  "Bulge Keratinocytes"           = "#e7dbc9", # Light Bone Tan
  "Spinous Keratinocytes"         = "#a0816c", # Wood Brown
  "Granular Keratinocytes"        = "#6e5646", # Dark Bark Brown
  "Hair Follicle Cells"           = "#d14949", # Coral Red
  
  # --- Dermal Stromal / Structural Family (Greens, Blues, Purples) ---
  "Papillary Fibroblasts"         = "#a6dbb1", # Soft Mint Green (Light Dermal)
  "Reticular Fibroblasts"         = "#136432", # Forest Green (Deep Dermal)
  "Dermal Fibroblasts"            = "#54a868", # Intermediate Green (General Dermal)
  "Endothelial Cells"             = "#4393c3", # Steel Blue
  "Muscle Cells"                  = "#6a3d9a", # Deep structural Violet
  
  # --- Immune Lymphoid Lineages (High-Contrast Pinks, Magentas, Plums) ---
  "CD4+ T Cells"                  = "#ff69b4", # Bright Hot Pink (Unified map target)
  "CD8+ T Cells/NK Cells"         = "#e3d3e8", # Blackberry (Matches cluster 21)
  "B Cells"                       = "#984ea3", # Vibrant Orchid Purple
  
  # --- Immune Myeloid Lineages (Teals, Oranges, Golds) ---
  "Dendritic Cells"               = "#ffd700", # Bright Gold
  "Langerhans Cells"              = "#fcae91", # Warm Peach/Apricot (Fixed spelling typo)
  "Macrophages"                   = "#01665e", # Deep Dermal Teal
  "Mast Cells"                    = "#ff7f00", # Vivid Safety Orange (Assigned to cluster 23)
  
  # --- Outliers ---
  "Melanocytes"                   = "grey"     # Grey baseline profile
)

# 1. Define a "Professional" Spatial Theme
# Removes unnecessary borders, cleans up legends, and optimizes for tissue visualization
spatial_theme <- theme_void(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white", color = NA), # Black background makes colors pop
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
  )

plot_spatial_pretty <- function(df, sample_name, pt_size = 0.8) {
  
  df_sub <- df %>%
    filter(sample == sample_name, !is.na(cell_type_minor))
  
  # Ensure the factor order matches your intended palette order
  # This prevents random color swapping
  df_sub$cell_type_minor <- factor(df_sub$cell_type_minor, levels = names(minor_colors))
  
  ggplot(df_sub, aes(x = x, y = y, color = cell_type_minor)) +
    geom_point(size = pt_size, alpha = 0.8, stroke = 0) + # Stroke 0 for cleaner edges
    coord_fixed() +
    scale_y_reverse() +
    scale_color_manual(values = minor_colors, drop = FALSE) + # drop=FALSE keeps unused colors in legend
    labs(title = sample_name) +
    spatial_theme
}

# 2. Iterate and Save
sample_names <- unique(as.character(plot_df$sample[!is.na(plot_df$sample)]))

for (s in sample_names) {
  p <- plot_spatial_pretty(plot_df, s)
  
  ggsave(
    filename = paste0("results/260506/spatial_major_", gsub("[^A-Za-z0-9_-]", "_", s), ".pdf"),
    plot = p,
    width = 7, # Slightly larger for better clarity
    height = 6,
    device = "pdf",
    dpi = 600
  )
}



library(Seurat)
library(dplyr)
library(ggplot2)

# Ensure output path directory structure is initialized
dir.create("results/260506", recursive = TRUE, showWarnings = FALSE)

# Define your exact custom discrete senescence color palette configuration
cdkn_colors <- c(
  "p16-p21-" = "#d3d3d3",
  "p21+"     = "#010087",
  "p16+"     = "#8d0000",
  "p16+p21+" = "#006900"
)


# ----------------------------------------------------------------------
# 1) EXTRACT COORDINATES (Using your exact working logic)
# ----------------------------------------------------------------------
fov_names <- Images(hsskin.ST)

extract_fov_coords <- function(obj, fov_name) {
  coords <- GetTissueCoordinates(obj[[fov_name]])
  coords <- as.data.frame(coords)
  coords$cell_id <- rownames(coords)
  coords$fov <- fov_name
  
  # robustly find x/y columns
  numeric_cols <- names(coords)[sapply(coords, is.numeric)]
  if (length(numeric_cols) < 2) {
    stop(paste("Could not find at least two numeric coordinate columns for", fov_name))
  }
  
  coords <- coords %>%
    rename(
      x = all_of(numeric_cols[1]),
      y = all_of(numeric_cols[2])
    )
  
  coords
}

coords_all <- bind_rows(lapply(fov_names, function(f) extract_fov_coords(hsskin.ST, f)))

coords_all <- coords_all %>%
  select(-cell_id) %>%
  rename(cell_id = cell)    # keep real cell names

# ----------------------------------------------------------------------
# 2) FETCH DATA & DYNAMICALLY BUILD HISTOLOGY METADATA
# ----------------------------------------------------------------------
# Pull raw matrices from the Xenium assay layers to determine CDKN positivity thresholds
df_expression <- FetchData(
  object = hsskin.ST,
  vars = c("sample", "cell_type_minor", "CDKN2A", "CDKN1A"),
  assay = "Xenium",
  layer = "data"
)

# Map discrete boolean flags (>0 means transcript is present) into categorical state tiers
md <- df_expression %>%
  mutate(
    cell_id = rownames(.),
    p16_pos = CDKN2A > 0,
    p21_pos = CDKN1A > 0,
    cdkn_state = case_when(
      p16_pos & p21_pos  ~ "p16+p21+",
      p16_pos & !p21_pos ~ "p16+",
      !p16_pos & p21_pos ~ "p21+",
      TRUE               ~ "p16-p21-"
    ),
    cdkn_state = factor(cdkn_state, levels = names(cdkn_colors))
  )
  # Restrict data rows strictly to skin structural categories

# ----------------------------------------------------------------------
# 3) JOIN COORDINATES WITH INFERRED SENESCENCE METADATA
# ----------------------------------------------------------------------
plot_df <- coords_all %>%
  left_join(md, by = "cell_id") %>%
  filter(!is.na(sample)) # Clean out cells outside matching filtered bounds

# ----------------------------------------------------------------------
# 4) ITERATE THROUGH INDIVIDUAL SPECIMENS TO EXPORT PDFs
# ----------------------------------------------------------------------
all_samples <- unique(plot_df$sample)

for (sample_name in all_samples) {
  cat("Exporting skin structural spatial map for:", sample_name, "\n")
  
  # Subset down to the current active looping specimen slice
  df_sub <- plot_df %>% 
    filter(sample == sample_name, !is.na(cdkn_state))
  
  p_individual <- ggplot(df_sub, aes(x = x, y = y, color = cdkn_state)) +
    # solid circle (shape=16) optimized for dense spatial spot grids
    geom_point(size = 0.5, alpha = 0.8, shape = 16, stroke = 0) +
    
    # Preserves true micron aspect ratio and handles potential orientation flips
    coord_fixed() + 
    scale_y_reverse() +   
    
    scale_color_manual(values = cdkn_colors, drop = FALSE, name = "CDKN State") +
    labs(title = sample_name) +
    
    # Standard microscope slide visualization style overrides
    theme_void(base_size = 14) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 10)),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 10)
    )
  
  # Format specific name targets safely
  clean_filename <- gsub("[^A-Za-z0-9_-]", "_", sample_name)
  output_path <- paste0("skin_mapping/results/260523/spatial_xenium_structural_", clean_filename, ".pdf")
  
  ggsave(
    filename = output_path,
    plot = p_individual,
    width = 6,   
    height = 5,
    dpi = 300
  )
}

cat("All structural tissue spatial maps successfully exported to target folders.\n")

library(dplyr)
library(ggplot2)
library(patchwork)

# 1. Define custom colors for Xenium cell_type_major
xenium_colors <- c(
  "Keratinocytes" = "#1f78b4",
  "Fibroblasts"   = "#33a02c",
  "Melanocytes"   = "#6a3d9a",
  "Muscle Cells"  = "#ff7f00",
  "none"          = "grey80"
)

# 2. Prepare the Xenium data (Total cell proportions)
skin_comp <- hsskin@meta.data %>%
  # Ensure we have data for age and your major lineage categories
  filter(!is.na(age) & !is.na(cell_type_major)) %>%
  # Keep primary skin cell types, collapse others (Immune, Blood Vessel, etc.) into "none"
  mutate(
    skin_lineage = case_when(
      cell_type_major %in% c("Keratinocytes", "Fibroblasts", "Melanocytes", "Muscle Cells") ~ as.character(cell_type_major),
      TRUE ~ "none"
    )
  ) %>%
  # Count occurrences per age group
  count(age, skin_lineage) %>%
  group_by(age) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  # Convert variables into ordered factors matching your target layout
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    skin_lineage = factor(skin_lineage, levels = names(xenium_colors))
  )

# 3. Use your updated Xenium donut function
donut_one1 <- function(df, age_value) {
  ggplot(df %>% filter(age == age_value),
         aes(x = 2, y = pct, fill = skin_lineage)) +
    geom_col(width = 1, color = "white", linewidth = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    scale_fill_manual(values = xenium_colors, drop = FALSE) +
    theme_void() +
    ggtitle(as.character(age_value)) +
    theme(
      text = element_text(family = "Helvetica", size = 14),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# 4. Generate and save the composite layout
p_young <- donut_one1(skin_comp, "Young")
p_old   <- donut_one1(skin_comp, "Old")

p_xenium_donut <- (p_young | p_old) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# Save as PDF using your specific date path structure
ggsave("results/260617/xenium_total_donut_major.pdf", p_xenium_donut, width = 7.5, height = 3.5)


library(Seurat)
library(dplyr)
library(ggplot2)

# 1. Subset object to focus strictly on Immune Cells
immune_subset <- subset(hsskin, subset = cell_type_major == "Immune Cells")

# Ensure age is formatted correctly as a factor
immune_subset$age <- factor(immune_subset$age, levels = c("Young", "Old"))

# 2. Define your target feature list matching exact human Xenium nomenclature
target_features <- c(
  "GZMB", "PRF1",       # Cytotoxicity
  "LAG3", "PDCD1", "HAVCR2", # Exhaustion Markers (PD-1, TIM-3)
  "CD8A", "GZMK"        # CD8 T cell / Tissue-resident memory state
)

Idents(immune_subset) <- 'cell_type_minor'
# 3. Option A: Generate a clean split-violin plot comparing Young vs Old
p_violins <- VlnPlot(
  object = immune_subset,
  features = target_features,
  split.by = "age",          # Highlights Young vs Old comparison within each feature
  #group.by = "cell_type_major", # Aggregated across the filtered immune matrix
  assay = "Xenium",          # Target your spatial transcriptomics matrix assay
  layer = "data",
  pt.size = 0,             # Low point size to prevent overlapping data clouding
  combine = TRUE
) & 
  scale_fill_manual(values = c("Young" = "#296d9e", "Old" = "#cc0b01")) &
  theme(text = element_text(family = "Helvetica"))

# Display plot
print(p_violins)

# Save the violin comparison
ggsave("results/260617/immune_age_expression_violins.pdf", p_violins, width = 10, height = 7)


# 4. Option B: Generate a DotPlot (Alternative if you want to look at % cell abundance shifts)
p_immune_dotplot <- DotPlot(
  object = immune_subset,
  features = target_features,
  group.by = "age",            # Compares Young rows vs Old rows side-by-side
  assay = "Xenium",
  #layer = "data",
  dot.scale = 6,               # Your custom dot size scale floor
  scale = FALSE                # Your custom setting: Unscaled raw-normalized values
) +
  # Your exact light blue -> dark blue gradient
  scale_color_gradientn(
    colours = c("lightgray", "#9ecae1", "#3182bd", "#08519c")
  ) +
  # Your exact theme overrides
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), # Rotated marker names
    axis.text.y = element_text(size = 9, face = "bold"),
    legend.position = "right",
    legend.justification = c(0, 0.9),
    legend.key.height = unit(0.3, "cm"),
    legend.box.spacing = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(5.5, 30, 5.5, 5.5)   # Keeps right-hand clipping safe
  )

# Save the dotplot comparison
ggsave("results/260617/immune_age_expression_dotplot.pdf", p_immune_dotplot, width = 4, height = 2)


library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Filter down strictly to the Immune Cells compartment
immune_subset <- subset(hsskin, subset = cell_type_major == "Immune Cells")

# 2. Define your target feature markers
genes_to_plot <- c("GZMB", "PRF1", "LAG3", "PDCD1", "HAVCR2", "CD8A", "GZMK")

# 3. Extract the expression data layer matrix
expression_matrix <- GetAssayData(immune_subset, assay = "Xenium", layer = "data")[genes_to_plot, ]
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Filter down strictly to the Immune Cells compartment
immune_subset <- subset(hsskin, subset = cell_type_major == "Immune Cells")

immune_subset$CD4_pos <- GetAssayData(
  immune_subset,
  assay = "Xenium",
  layer = "counts"
)["CD4", ] > 0

immune_subset$CD8_pos <- GetAssayData(
  immune_subset,
  assay = "Xenium",
  layer = "counts"
)["CD8A", ] > 0

cd4_subset <- subset(immune_subset, subset = CD4_pos)
cd8_subset <- subset(immune_subset, subset = CD8_pos)

# 2. Define your target feature markers
genes_to_plot <- c("GZMB", "PRF1", "LAG3", "PDCD1", "HAVCR2", "GZMK")

# 3. Extract the normalized data layer matrix
expression_matrix <- GetAssayData(immune_subset, assay = "Xenium", layer = "data")[genes_to_plot, ]

# 4. Calculate per-sample positivity percentages (for dots and statistics)
sample_positivity <- as.data.frame(t(as.matrix(expression_matrix))) %>%
  mutate(
    age = immune_subset$age,
    sample = immune_subset$sample
  ) %>%
  pivot_longer(cols = all_of(genes_to_plot), names_to = "Gene", values_to = "Expression") %>%
  mutate(is_positive = Expression > 0) %>%
  group_by(age, sample, Gene) %>%
  summarize(
    total_cells = n(),
    positive_cells = sum(is_positive),
    pct_positive = (positive_cells / total_cells) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    Gene = factor(Gene, levels = genes_to_plot)
  )

# 5. Calculate global pooled positivity percentages (for the bar heights)
global_positivity <- as.data.frame(t(as.matrix(expression_matrix))) %>%
  mutate(age = immune_subset$age) %>%
  pivot_longer(cols = all_of(genes_to_plot), names_to = "Gene", values_to = "Expression") %>%
  mutate(is_positive = Expression > 0) %>%
  group_by(age, Gene) %>%
  summarize(
    total_cells = n(),
    positive_cells = sum(is_positive),
    pct_positive = (positive_cells / total_cells) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    Gene = factor(Gene, levels = genes_to_plot)
  )

# 6. Compute p-values across samples for each gene
p_values <- sample_positivity %>%
  group_by(Gene) %>%
  summarize(
    p_val = if(length(unique(age)) == 2) {
      t.test(pct_positive ~ age)$p.value
    } else {
      NA
    },
    .groups = "drop"
  ) %>%
  mutate(p_label = case_when(
    is.na(p_val) ~ "",
    p_val < 0.001 ~ "p < 0.001",
    p_val < 0.05  ~ sprintf("p = %.3f *", p_val), # Asterisk marks significance
    TRUE          ~ sprintf("p = %.3f", p_val)
  ))

# 7. Coordinate text positions smoothly above the highest elements
label_positions <- sample_positivity %>%
  group_by(Gene) %>%
  summarize(max_dot = max(pct_positive), .groups = "drop") %>%
  left_join(global_positivity %>% group_by(Gene) %>% summarize(max_bar = max(pct_positive), .groups = "drop"), by = "Gene") %>%
  rowwise() %>%
  mutate(upper_bound = max(max_dot, max_bar)) %>%
  ungroup() %>%
  left_join(p_values, by = "Gene") %>%
  mutate(y_pos = upper_bound + 4) # 4% margin headroom above tallest asset

# 8. Plot with global bars, sample dots, and p-value labels
p_positivity_bars <- ggplot() +
  # Global proportion bars
  geom_col(
    data = global_positivity,
    aes(x = Gene, y = pct_positive, fill = age),
    position = position_dodge(width = 0.75), 
    width = 0.65, 
    color = "black", 
    linewidth = 0.3
  ) +
  # Per-sample replicate dots (properly aligned/dodged over respective bars)
  geom_point(
    data = sample_positivity,
    aes(x = Gene, y = pct_positive, group = age),
    position = position_dodge(width = 0.75),
    color = "black",
    fill = "white", # High-contrast white fill with crisp black outline
    shape = 21,
    size = 1.4,
    stroke = 0.5,
    show.legend = FALSE
  ) +
  # P-value text blocks centered over each comparison bracket
  geom_text(
    data = label_positions,
    aes(x = Gene, y = y_pos, label = p_label),
    vjust = 0,
    family = "Helvetica",
    size = 3.2,
    fontface = "italic"
  ) +
  labs(
    y = "% of Positive CD4+ T Cells (>0 values)",
    x = NULL
  ) +
  scale_fill_manual(values = c("Young" = "#296d9e", "Old" = "#cc0b01")) +
  # Leaves 15% headroom dynamically at the top of the plot grid for labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    
    legend.position = "right",
    legend.title = element_blank(),
    legend.justification = c(0, 0.9),
    legend.key.height = unit(0.3, "cm"),
    legend.box.spacing = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(5.5, 30, 5.5, 5.5)
  )

# Output plot display
print(p_positivity_bars)

# Save as vector PDF
ggsave(
  filename = "results/260617/immune_markers_positivity_CD4+ T Cells_bars.pdf", 
  plot = p_positivity_bars, 
  width = 6.2, 
  height = 3.8
)


library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

plot_marker_positivity <- function(obj, genes_to_plot, y_label, save_path = NULL) {
  
  obj_sub <- obj
  
  genes_present <- genes_to_plot[genes_to_plot %in% rownames(obj_sub[["Xenium"]])]
  
  expression_matrix <- GetAssayData(
    obj_sub,
    assay = "Xenium",
    layer = "data"
  )[genes_present, , drop = FALSE]
  
  sample_positivity <- as.data.frame(t(as.matrix(expression_matrix))) %>%
    mutate(
      age = obj_sub$age,
      sample = obj_sub$sample
    ) %>%
    pivot_longer(
      cols = all_of(genes_present),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    mutate(is_positive = Expression > 0) %>%
    group_by(age, sample, Gene) %>%
    summarise(
      total_cells = n(),
      positive_cells = sum(is_positive),
      pct_positive = positive_cells / total_cells * 100,
      .groups = "drop"
    ) %>%
    mutate(
      age = factor(age, levels = c("Young", "Old")),
      Gene = factor(Gene, levels = genes_present)
    )
  
  global_positivity <- as.data.frame(t(as.matrix(expression_matrix))) %>%
    mutate(age = obj_sub$age) %>%
    pivot_longer(
      cols = all_of(genes_present),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    mutate(is_positive = Expression > 0) %>%
    group_by(age, Gene) %>%
    summarise(
      total_cells = n(),
      positive_cells = sum(is_positive),
      pct_positive = positive_cells / total_cells * 100,
      .groups = "drop"
    ) %>%
    mutate(
      age = factor(age, levels = c("Young", "Old")),
      Gene = factor(Gene, levels = genes_present)
    )
  
  p_values <- sample_positivity %>%
    group_by(Gene) %>%
    summarise(
      p_val = if (length(unique(age)) == 2) {
        t.test(pct_positive ~ age)$p.value
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      p_label = case_when(
        is.na(p_val) ~ "",
        p_val < 0.001 ~ "p < 0.001",
        p_val < 0.05  ~ sprintf("p = %.3f *", p_val),
        TRUE          ~ sprintf("p = %.3f", p_val)
      )
    )
  
  label_positions <- sample_positivity %>%
    group_by(Gene) %>%
    summarise(max_dot = max(pct_positive), .groups = "drop") %>%
    left_join(
      global_positivity %>%
        group_by(Gene) %>%
        summarise(max_bar = max(pct_positive), .groups = "drop"),
      by = "Gene"
    ) %>%
    rowwise() %>%
    mutate(upper_bound = max(max_dot, max_bar)) %>%
    ungroup() %>%
    left_join(p_values, by = "Gene") %>%
    mutate(y_pos = upper_bound + 4)
  
  p <- ggplot() +
    geom_col(
      data = global_positivity,
      aes(x = Gene, y = pct_positive, fill = age),
      position = position_dodge(width = 0.75),
      width = 0.65,
      color = "black",
      linewidth = 0.3
    ) +
    geom_point(
      data = sample_positivity,
      aes(x = Gene, y = pct_positive, group = age),
      position = position_dodge(width = 0.75),
      color = "black",
      fill = "white",
      shape = 21,
      size = 1.4,
      stroke = 0.5,
      show.legend = FALSE
    ) +
    geom_text(
      data = label_positions,
      aes(x = Gene, y = y_pos, label = p_label),
      vjust = 0,
      family = "Helvetica",
      size = 3.2,
      fontface = "italic"
    ) +
    labs(y = y_label, x = NULL) +
    scale_fill_manual(values = c("Young" = "#296d9e", "Old" = "#cc0b01")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
    theme_classic(base_size = 12) +
    theme(
      text = element_text(family = "Helvetica"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
      axis.text.y = element_text(color = "black"),
      legend.position = "right",
      legend.title = element_blank()
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 6.2, height = 5)
  }
  
  p
}

genes_to_plot <- c("GZMB", "PRF1", "LAG3", "PDCD1", "HAVCR2", "GZMK")

p_cd4 <- plot_marker_positivity(
  obj = cd4_subset,
  genes_to_plot = genes_to_plot,
  y_label = "% of Positive CD4+ Immune Cells (>0 values)",
  save_path = "results/260617/immune_markers_positivity_CD4pos_immune_cells_bars.pdf"
)

p_cd8 <- plot_marker_positivity(
  obj = cd8_subset,
  genes_to_plot = genes_to_plot,
  y_label = "% of Positive CD8A+ Immune Cells (>0 values)",
  save_path = "results/260617/immune_markers_positivity_CD8Apos_immune_cells_bars.pdf"
)


###################################
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

plot_positivity_simple <- function(obj, celltype_value, genes_to_plot, y_label, save_path = NULL) {
  
  obj_sub <- subset(obj, subset = cell_type_minor == celltype_value)
  
  genes_present <- genes_to_plot[genes_to_plot %in% rownames(obj_sub[["Xenium"]])]
  
  mat <- GetAssayData(obj_sub, assay = "Xenium", layer = "data")[genes_present, , drop = FALSE]
  
  df <- as.data.frame(t(as.matrix(mat))) %>%
    mutate(age = obj_sub$age, sample = obj_sub$sample) %>%
    pivot_longer(
      cols = all_of(genes_present),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    mutate(
      is_positive = Expression > 0,
      age = factor(age, levels = c("Young", "Old")),
      Gene = factor(Gene, levels = genes_present)
    )
  
  sample_positivity <- df %>%
    group_by(age, sample, Gene) %>%
    summarise(
      pct_positive = mean(is_positive) * 100,
      .groups = "drop"
    )
  
  global_positivity <- df %>%
    group_by(age, Gene) %>%
    summarise(
      pct_positive = mean(is_positive) * 100,
      .groups = "drop"
    )
  
  p_values <- sample_positivity %>%
    group_by(Gene) %>%
    summarise(
      p_val = if (length(unique(age)) == 2) {
        t.test(pct_positive ~ age)$p.value
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      p_label = case_when(
        is.na(p_val) ~ "",
        p_val < 0.001 ~ "p < 0.001",
        p_val < 0.05  ~ sprintf("p = %.3f *", p_val),
        TRUE          ~ sprintf("p = %.3f", p_val)
      )
    )
  
  label_positions <- sample_positivity %>%
    group_by(Gene) %>%
    summarise(max_dot = max(pct_positive), .groups = "drop") %>%
    left_join(
      global_positivity %>%
        group_by(Gene) %>%
        summarise(max_bar = max(pct_positive), .groups = "drop"),
      by = "Gene"
    ) %>%
    rowwise() %>%
    mutate(upper_bound = max(max_dot, max_bar)) %>%
    ungroup() %>%
    left_join(p_values, by = "Gene") %>%
    mutate(y_pos = upper_bound + 4)
  
  p <- ggplot() +
    geom_col(
      data = global_positivity,
      aes(x = Gene, y = pct_positive, fill = age),
      position = position_dodge(width = 0.75),
      width = 0.65,
      color = "black",
      linewidth = 0.3
    ) +
    geom_point(
      data = sample_positivity,
      aes(x = Gene, y = pct_positive, group = age),
      position = position_dodge(width = 0.75),
      color = "black",
      fill = "white",
      shape = 21,
      size = 1.5,
      stroke = 0.5,
      show.legend = FALSE
    )+
    geom_text(
      data = label_positions,
      aes(x = Gene, y = y_pos, label = p_label),
      vjust = 0,
      family = "Helvetica",
      size = 3.2,
      fontface = "italic"
    ) +
    labs(y = y_label, x = NULL) +
    scale_fill_manual(values = c("Young" = "#296d9e", "Old" = "#cc0b01")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
      axis.text.y = element_text(color = "black"),
      legend.position = "right",
      legend.title = element_blank()
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 6.2, height = 5)
  }
  
  p
}

p_cd4 <- plot_positivity_simple(
  obj = hsskin,
  celltype_value = "CD4+ T Cells",
  genes_to_plot = genes_to_plot,
  y_label = "% of Positive CD4+ T Cells",
  save_path = "results/260617/immune_markers_CD4_T_cells_simple.pdf"
)

p_cd8 <- plot_positivity_simple(
  obj = hsskin,
  celltype_value = "CD8+ T Cells/NK Cells",
  genes_to_plot = genes_to_plot,
  y_label = "% of Positive CD8+ T Cells/NK Cells",
  save_path = "results/260617/immune_markers_CD8_T_cells_NK_simple.pdf"
)

print(p_cd4)
print(p_cd8)






genes_to_plot <- c("CD68", "ITGAE")
expression_matrix <- GetAssayData(immune_subset, assay = "Xenium", layer = "data")[genes_to_plot, ]

sample_positivity <- as.data.frame(t(as.matrix(expression_matrix))) %>%
  mutate(
    age = immune_subset$age,
    sample = immune_subset$sample
  ) %>%
  pivot_longer(cols = all_of(genes_to_plot), names_to = "Gene", values_to = "Expression") %>%
  mutate(is_positive = Expression > 0) %>%
  group_by(age, sample, Gene) %>%
  summarize(
    total_cells = n(),
    positive_cells = sum(is_positive),
    pct_positive = (positive_cells / total_cells) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    Gene = factor(Gene, levels = genes_to_plot)
  )

# 5. Calculate global pooled positivity percentages (for the bar heights)
global_positivity <- as.data.frame(t(as.matrix(expression_matrix))) %>%
  mutate(age = immune_subset$age) %>%
  pivot_longer(cols = all_of(genes_to_plot), names_to = "Gene", values_to = "Expression") %>%
  mutate(is_positive = Expression > 0) %>%
  group_by(age, Gene) %>%
  summarize(
    total_cells = n(),
    positive_cells = sum(is_positive),
    pct_positive = (positive_cells / total_cells) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    Gene = factor(Gene, levels = genes_to_plot)
  )

# 6. Compute p-values across samples for each gene
p_values <- sample_positivity %>%
  group_by(Gene) %>%
  summarize(
    p_val = if(length(unique(age)) == 2) {
      t.test(pct_positive ~ age)$p.value
    } else {
      NA
    },
    .groups = "drop"
  ) %>%
  mutate(p_label = case_when(
    is.na(p_val) ~ "",
    p_val < 0.001 ~ "p < 0.001",
    p_val < 0.05  ~ sprintf("p = %.3f *", p_val), # Asterisk marks significance
    TRUE          ~ sprintf("p = %.3f", p_val)
  ))

# 7. Coordinate text positions smoothly above the highest elements
label_positions <- sample_positivity %>%
  group_by(Gene) %>%
  summarize(max_dot = max(pct_positive), .groups = "drop") %>%
  left_join(global_positivity %>% group_by(Gene) %>% summarize(max_bar = max(pct_positive), .groups = "drop"), by = "Gene") %>%
  rowwise() %>%
  mutate(upper_bound = max(max_dot, max_bar)) %>%
  ungroup() %>%
  left_join(p_values, by = "Gene") %>%
  mutate(y_pos = upper_bound + 4) # 4% margin headroom above tallest asset

# 8. Plot with global bars, sample dots, and p-value labels
p_positivity_bars <- ggplot() +
  # Global proportion bars
  geom_col(
    data = global_positivity,
    aes(x = Gene, y = pct_positive, fill = age),
    position = position_dodge(width = 0.75), 
    width = 0.65, 
    color = "black", 
    linewidth = 0.3
  ) +
  # Per-sample replicate dots (properly aligned/dodged over respective bars)
  geom_point(
    data = sample_positivity,
    aes(x = Gene, y = pct_positive, group = age),
    position = position_dodge(width = 0.75),
    color = "black",
    fill = "white", # High-contrast white fill with crisp black outline
    shape = 21,
    size = 1.4,
    stroke = 0.5,
    show.legend = FALSE
  ) +
  # P-value text blocks centered over each comparison bracket
  geom_text(
    data = label_positions,
    aes(x = Gene, y = y_pos, label = p_label),
    vjust = 0,
    family = "Helvetica",
    size = 3.2,
    fontface = "italic"
  ) +
  labs(
    y = "% of Positive Immune Cells (>0 values)",
    x = NULL
  ) +
  scale_fill_manual(values = c("Young" = "#296d9e", "Old" = "#cc0b01")) +
  # Leaves 15% headroom dynamically at the top of the plot grid for labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    
    legend.position = "right",
    legend.title = element_blank(),
    legend.justification = c(0, 0.9),
    legend.key.height = unit(0.3, "cm"),
    legend.box.spacing = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(5.5, 30, 5.5, 5.5)
  )

# Output plot display
print(p_positivity_bars)

# Save as vector PDF
ggsave(
  filename = "results/260617/immune_markers_positivity_stat_bars.pdf", 
  plot = p_positivity_bars, 
  width = 6.2, 
  height = 3.8
)



# 1. Master list of human chemokine receptors organized by family
genes_to_plot <- c(
  "CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7",
  "CXCR1", "CXCR2", "CXCR3", "CXCR4", "CXCR5", "CXCR6",
  "CX3CR1"
)

genes_to_plot <- c(
  "IL10RA", "IL12RB2", "IL13RA2", "IL15RA", "IL1R1", 
  "IL1RAP", "IL22RA1", "IL22RA2", "IL4R", "IL2RA", 
  "IL2RB", "IL2RG"
)

expression_matrix <- GetAssayData(immune_subset, assay = "Xenium", layer = "data")[genes_to_plot, ]

sample_positivity <- as.data.frame(t(as.matrix(expression_matrix))) %>%
  mutate(
    age = immune_subset$age,
    sample = immune_subset$sample
  ) %>%
  pivot_longer(cols = all_of(genes_to_plot), names_to = "Gene", values_to = "Expression") %>%
  mutate(is_positive = Expression > 0) %>%
  group_by(age, sample, Gene) %>%
  summarize(
    total_cells = n(),
    positive_cells = sum(is_positive),
    pct_positive = (positive_cells / total_cells) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    Gene = factor(Gene, levels = genes_to_plot)
  )

# 5. Calculate global pooled positivity percentages (for the bar heights)
global_positivity <- as.data.frame(t(as.matrix(expression_matrix))) %>%
  mutate(age = immune_subset$age) %>%
  pivot_longer(cols = all_of(genes_to_plot), names_to = "Gene", values_to = "Expression") %>%
  mutate(is_positive = Expression > 0) %>%
  group_by(age, Gene) %>%
  summarize(
    total_cells = n(),
    positive_cells = sum(is_positive),
    pct_positive = (positive_cells / total_cells) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    Gene = factor(Gene, levels = genes_to_plot)
  )

# 6. Compute p-values across samples for each gene
p_values <- sample_positivity %>%
  group_by(Gene) %>%
  summarize(
    p_val = if(length(unique(age)) == 2) {
      t.test(pct_positive ~ age)$p.value
    } else {
      NA
    },
    .groups = "drop"
  ) %>%
  mutate(p_label = case_when(
    is.na(p_val) ~ "",
    p_val < 0.001 ~ "p < 0.001",
    p_val < 0.05  ~ sprintf("p = %.3f *", p_val), # Asterisk marks significance
    TRUE          ~ sprintf("p = %.3f", p_val)
  ))

# 7. Coordinate text positions smoothly above the highest elements
label_positions <- sample_positivity %>%
  group_by(Gene) %>%
  summarize(max_dot = max(pct_positive), .groups = "drop") %>%
  left_join(global_positivity %>% group_by(Gene) %>% summarize(max_bar = max(pct_positive), .groups = "drop"), by = "Gene") %>%
  rowwise() %>%
  mutate(upper_bound = max(max_dot, max_bar)) %>%
  ungroup() %>%
  left_join(p_values, by = "Gene") %>%
  mutate(y_pos = upper_bound + 4) # 4% margin headroom above tallest asset

# 8. Plot with global bars, sample dots, and p-value labels
p_positivity_bars <- ggplot() +
  # Global proportion bars
  geom_col(
    data = global_positivity,
    aes(x = Gene, y = pct_positive, fill = age),
    position = position_dodge(width = 0.75), 
    width = 0.65, 
    color = "black", 
    linewidth = 0.3
  ) +
  # Per-sample replicate dots (properly aligned/dodged over respective bars)
  geom_point(
    data = sample_positivity,
    aes(x = Gene, y = pct_positive, group = age),
    position = position_dodge(width = 0.75),
    color = "black",
    fill = "white", # High-contrast white fill with crisp black outline
    shape = 21,
    size = 1.4,
    stroke = 0.5,
    show.legend = FALSE
  ) +
  # P-value text blocks centered over each comparison bracket
  geom_text(
    data = label_positions,
    aes(x = Gene, y = y_pos, label = p_label),
    vjust = 0,
    family = "Helvetica",
    size = 3.2,
    fontface = "italic"
  ) +
  labs(
    y = "% of Positive Immune Cells (>0 values)",
    x = NULL
  ) +
  scale_fill_manual(values = c("Young" = "#296d9e", "Old" = "#cc0b01")) +
  # Leaves 15% headroom dynamically at the top of the plot grid for labels
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) + 
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    
    legend.position = "right",
    legend.title = element_blank(),
    legend.justification = c(0, 0.9),
    legend.key.height = unit(0.3, "cm"),
    legend.box.spacing = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(5.5, 30, 5.5, 5.5)
  )

# Output plot display
print(p_positivity_bars)

# Save as vector PDF
ggsave(
  filename = "results/260617/immune_markers_positivity_stat_bars.pdf", 
  plot = p_positivity_bars, 
  width = 6.2, 
  height = 3.8
)


library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# 2. Define the individual markers required for your combinations
target_genes <- c("LAG3", "PDCD1", "HAVCR2", "GZMK")

# 3. Pull the normalized data layer matrix
expression_matrix <- GetAssayData(immune_subset, assay = "Xenium", layer = "data")[target_genes, ]

# 4. Construct a data frame of logical binary true/false vectors per cell (> 0)
cell_combinations <- as.data.frame(t(as.matrix(expression_matrix))) %>%
  mutate(
    age = immune_subset$age,
    sample = immune_subset$sample,
    
    # Evaluate individual positivity
    lag3_pos  = LAG3 > 0,
    pdcd1_pos = PDCD1 > 0,
    havcr2_pos = HAVCR2 > 0,
    gzmk_pos  = GZMK > 0,
    
    # Establish your co-expression state buckets
    `LAG3+PDCD1+HAVCR2+` = lag3_pos & pdcd1_pos & havcr2_pos,
    `LAG3+PDCD1+HAVCR2+GZMK+` = lag3_pos & pdcd1_pos & havcr2_pos & gzmk_pos
  )

# 5. Define combination vectors for plotting tracking
combos_to_plot <- c("LAG3+PDCD1+HAVCR2+", "LAG3+PDCD1+HAVCR2+GZMK+")

# 6. Calculate per-sample positivity percentages (for dots and statistics)
sample_positivity <- cell_combinations %>%
  select(age, sample, all_of(combos_to_plot)) %>%
  pivot_longer(cols = all_of(combos_to_plot), names_to = "Combination", values_to = "is_positive") %>%
  group_by(age, sample, Combination) %>%
  summarize(
    total_cells = n(),
    positive_cells = sum(is_positive),
    pct_positive = (positive_cells / total_cells) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    Combination = factor(Combination, levels = combos_to_plot)
  )

# 7. Calculate global pooled positivity percentages (for the bar heights)
global_positivity <- cell_combinations %>%
  select(age, all_of(combos_to_plot)) %>%
  pivot_longer(cols = all_of(combos_to_plot), names_to = "Combination", values_to = "is_positive") %>%
  group_by(age, Combination) %>%
  summarize(
    total_cells = n(),
    positive_cells = sum(is_positive),
    pct_positive = (positive_cells / total_cells) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    Combination = factor(Combination, levels = combos_to_plot)
  )

# 8. Compute p-values across sample groups for each combination
p_values <- sample_positivity %>%
  group_by(Combination) %>%
  summarize(
    p_val = if(length(unique(age)) == 2) t.test(pct_positive ~ age)$p.value else NA,
    .groups = "drop"
  ) %>%
  mutate(p_label = case_when(
    is.na(p_val) ~ "",
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE          ~ "ns"
  ))

# 9. Calculate label positioning overhead bounds to prevent text clipping
label_positions <- sample_positivity %>%
  group_by(Combination) %>%
  summarize(max_dot = max(pct_positive), .groups = "drop") %>%
  left_join(global_positivity %>% group_by(Combination) %>% summarize(max_bar = max(pct_positive), .groups = "drop"), by = "Combination") %>%
  rowwise() %>%
  mutate(upper_bound = max(max_dot, max_bar)) %>%
  ungroup() %>%
  left_join(p_values, by = "Combination") %>%
  mutate(y_pos = upper_bound + 1) # Margins set slightly lower since combos yield smaller raw % pools

# 10. Build custom styled combinatorial comparison plot
p_combo_bars <- ggplot() +
  # Global proportion bars
  geom_col(
    data = global_positivity,
    aes(x = Combination, y = pct_positive, fill = age),
    position = position_dodge(width = 0.75), 
    width = 0.55, 
    color = "black", 
    linewidth = 0.3
  ) +
  # Per-sample dots aligned cleanly over bars
  geom_point(
    data = sample_positivity,
    aes(x = Combination, y = pct_positive, group = age),
    position = position_dodge(width = 0.75), 
    color = "black", 
    fill = "white", 
    shape = 21, 
    size = 1.5, 
    stroke = 0.5, 
    show.legend = FALSE
  ) +
  # Statistical text significance marks
  geom_text(
    data = label_positions,
    aes(x = Combination, y = y_pos, label = p_label),
    vjust = -0.2, 
    family = "Helvetica", 
    size = 3.8, 
    fontface = "bold"
  ) +
  labs(
    y = "% of Co-expressing Immune Cells", 
    x = NULL
  ) +
  scale_fill_manual(values = c("Young" = "#296d9e", "Old" = "#cc0b01")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) + # Leaves 20% head cushion at top
  theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    
    # Matches your exact legend / plot margin parameters
    legend.position = "right",
    legend.title = element_blank(),
    legend.justification = c(0, 0.9),
    legend.key.height = unit(0.3, "cm"),
    legend.box.spacing = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(5.5, 30, 5.5, 5.5)
  )

# Render layout
print(p_combo_bars)

# Save vector file output
ggsave("results/260617/immune_exhaustion_combinations_bars.pdf", p_combo_bars, width = 6.0, height = 4.2)



#########################################
#proliferation

plot_nonprolif_fraction <- function(obj, celltype_value, save_path = NULL) {
  
  obj_sub <- subset(obj, subset = cell_type_minor == celltype_value)
  
  mat <- GetAssayData(obj_sub, assay = "Xenium", layer = "counts")[c("MKI67", "PCNA"), , drop = FALSE]
  
  df <- obj_sub@meta.data %>%
    mutate(
      MKI67_pos = mat["MKI67", ] > 0,
      PCNA_pos  = mat["PCNA", ] > 0,
      nonprolif = !(MKI67_pos | PCNA_pos)
    ) %>%
    group_by(age, sample) %>%
    summarise(nonprolif_pct = mean(nonprolif) * 100, .groups = "drop") %>%
    mutate(age = factor(age, levels = c("Young", "Old")))
  
  p <- ggplot(df, aes(x = age, y = nonprolif_pct, fill = age)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 2, shape = 21, color = "black") +
    scale_fill_manual(values = c("Young" = "#296d9e", "Old" = "#cc0b01")) +
    theme_classic(base_size = 12) +
    labs(x = NULL, y = "% MKI67-/PCNA- cells") +
    theme(legend.position = "none")
  
  if (!is.null(save_path)) ggsave(save_path, p, width = 3.5, height = 4.5)
  p
}

p_cd4_nonprolif <- plot_nonprolif_fraction(
  hsskin,
  "CD4+ T Cells",
  "results/260617/CD4_T_cells_nonprolif_fraction.pdf"
)

p_cd8_nonprolif <- plot_nonprolif_fraction(
  hsskin,
  "CD8+ T Cells/NK Cells",
  "results/260617/CD8_T_cells_NK_nonprolif_fraction.pdf"
)

plot_nonprolif_fraction <- function(obj, label, save_path = NULL) {
  
  mat <- GetAssayData(
    obj,
    assay = "Xenium",
    layer = "counts"
  )[c("MKI67", "PCNA"), , drop = FALSE]
  
  df <- obj@meta.data %>%
    mutate(
      MKI67_pos = mat["MKI67", ] > 0,
      PCNA_pos  = mat["PCNA", ] > 0,
      nonprolif = !(MKI67_pos | PCNA_pos)
    ) %>%
    group_by(age, sample) %>%
    summarise(
      pct_nonprolif = mean(nonprolif) * 100,
      .groups = "drop"
    ) %>%
    mutate(age = factor(age, levels = c("Young", "Old")))
  
  p <- ggplot(df, aes(x = age, y = pct_nonprolif, fill = age)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.12, size = 2, shape = 21,
                color = "black", fill = "white") +
    scale_fill_manual(values = c("Young" = "#296d9e",
                                 "Old" = "#cc0b01")) +
    theme_classic(base_size = 12) +
    labs(
      title = label,
      x = NULL,
      y = "% MKI67-/PCNA- cells"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 3.8, height = 4.5)
  }
  
  p
}

p_cd4_nonprolif <- plot_nonprolif_fraction(
  cd4_subset,
  label = "CD4+ immune cells",
  save_path = "results/260617/CD4pos_nonprolif.pdf"
)

p_cd8_nonprolif <- plot_nonprolif_fraction(
  cd8_subset,
  label = "CD8A+ immune cells",
  save_path = "results/260617/CD8pos_nonprolif.pdf"
)

print(p_cd4_nonprolif)
print(p_cd8_nonprolif)