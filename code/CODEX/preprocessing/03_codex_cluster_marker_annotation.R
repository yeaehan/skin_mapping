##########################################################
# CODEX cluster marker annotation and sample/age diagnostics
#
# Uses the fixed clusters from 02_codex_integrated_analysis.R. It does not
# recluster cells or apply batch correction.
##########################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(file_arg)) stop("Run this script with Rscript")
script_path <- normalizePath(sub("^--file=", "", file_arg[1]),
                             winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."),
                              winslash = "/", mustWork = TRUE)
setwd(project_root)

input_file <- "data/codex/codex_merged_log1p_uncorrected.rds"
output_dir <- "results/codex/preprocessing/cluster_marker_annotation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Missing ", input_file, "; run preprocessing step 02 first")

obj <- readRDS(input_file)
required_metadata <- c("uncorrected_cluster", "original_sample", "analysis_age")
missing_metadata <- setdiff(required_metadata, colnames(obj@meta.data))
if (length(missing_metadata)) stop("Missing metadata: ", paste(missing_metadata, collapse = ", "))

identity_markers_requested <- c(
  "SMA", "Keratin 14", "Vimentin", "Keratin 5", "CD34", "CD8", "CD14",
  "Filaggrin", "FOXP3", "MelanA", "CD4", "CD68", "HLA-DR", "CD45"
)
identity_markers <- intersect(identity_markers_requested, rownames(obj))
if (!length(identity_markers)) stop("No identity markers found")

raw <- LayerData(obj, assay = "RNA", layer = "counts")[identity_markers, , drop = FALSE]
logx <- LayerData(obj, assay = "CODEXlog1p", layer = "data")[identity_markers, , drop = FALSE]
stopifnot(identical(colnames(raw), colnames(logx)))

cell_metadata <- obj@meta.data[colnames(raw), , drop = FALSE] %>%
  transmute(
    cell_id = rownames(.),
    cluster = factor(as.character(uncorrected_cluster),
                     levels = sort(unique(as.integer(as.character(uncorrected_cluster))))),
    sample = as.character(original_sample),
    age = factor(as.character(analysis_age), levels = c("Young", "Old"))
  )

# Long table: every cell x identity marker. Positive fraction and positive
# magnitude are retained as separate quantities.
raw_cells_markers <- as.matrix(t(raw))
log_cells_markers <- as.matrix(t(logx))
n_cells <- nrow(raw_cells_markers)
n_markers <- ncol(raw_cells_markers)
marker_long <- data.frame(
  cell_id = rep(rownames(raw_cells_markers), times = n_markers),
  marker = rep(colnames(raw_cells_markers), each = n_cells),
  raw_intensity = as.vector(raw_cells_markers),
  log1p_intensity = as.vector(log_cells_markers),
  stringsAsFactors = FALSE
) %>%
  left_join(cell_metadata, by = "cell_id") %>%
  mutate(positive = raw_intensity > 0)

# Pooled cluster marker profiles are used only to name/interpret the already
# fixed clusters, not for age inference.
cluster_marker <- marker_long %>%
  group_by(cluster, marker) %>%
  summarise(
    n_cells = n(),
    n_positive = sum(positive),
    fraction_positive = mean(positive),
    mean_log1p = mean(log1p_intensity),
    median_log1p = median(log1p_intensity),
    median_positive_raw = if (any(positive)) median(raw_intensity[positive]) else NA_real_,
    .groups = "drop"
  ) %>%
  group_by(marker) %>%
  mutate(
    fraction_positive_z = as.numeric(scale(fraction_positive)),
    mean_log1p_z = as.numeric(scale(mean_log1p)),
    combined_marker_score = rowMeans(
      cbind(fraction_positive_z, mean_log1p_z), na.rm = TRUE
    )
  ) %>%
  ungroup()

marker_lineage <- c(
  "SMA" = "Smooth muscle",
  "Keratin 14" = "Basal keratinocyte",
  "Keratin 5" = "Basal keratinocyte",
  "Filaggrin" = "Differentiated keratinocyte",
  "Vimentin" = "Mesenchymal/fibroblast",
  "CD34" = "Vascular/CD34+ stromal",
  "MelanA" = "Melanocyte",
  "CD45" = "Immune",
  "CD4" = "CD4 T cell",
  "CD8" = "CD8 T cell",
  "FOXP3" = "Regulatory T-cell candidate",
  "CD14" = "Monocyte/macrophage",
  "CD68" = "Macrophage",
  "HLA-DR" = "Antigen-presenting immune"
)

enriched_marker <- cluster_marker %>%
  group_by(cluster) %>%
  arrange(desc(combined_marker_score), desc(fraction_positive), .by_group = TRUE) %>%
  summarise(most_enriched_marker = first(marker), .groups = "drop")

cluster_annotation <- cluster_marker %>%
  group_by(cluster) %>%
  arrange(desc(fraction_positive), desc(mean_log1p), .by_group = TRUE) %>%
  summarise(
    top_marker = first(marker),
    top_marker_fraction_positive = first(fraction_positive),
    top_marker_mean_log1p = first(mean_log1p),
    top_3_markers = paste(head(marker, 3), collapse = ";"),
    .groups = "drop"
  ) %>%
  left_join(enriched_marker, by = "cluster") %>%
  mutate(
    candidate_lineage = ifelse(
      top_marker_fraction_positive < 0.10,
      "Low shared-identity signal",
      unname(marker_lineage[top_marker])
    ),
    cluster_marker_label = ifelse(
      top_marker_fraction_positive < 0.10,
      paste0("C", cluster, ": low identity"),
      paste0("C", cluster, ": ", top_marker)
    )
  )

write.csv(cluster_marker, file.path(output_dir, "cluster_identity_marker_profiles.csv"),
          row.names = FALSE)
write.csv(cluster_annotation, file.path(output_dir, "cluster_top_marker_annotations.csv"),
          row.names = FALSE)

# Calculate sample-level marker profiles first. These are the appropriate
# units for evaluating reproducibility and age trends.
cluster_sample_marker <- marker_long %>%
  group_by(cluster, sample, age, marker) %>%
  summarise(
    n_cells = n(),
    n_positive = sum(positive),
    fraction_positive = mean(positive),
    mean_log1p = mean(log1p_intensity),
    median_log1p = median(log1p_intensity),
    median_positive_raw = if (any(positive)) median(raw_intensity[positive]) else NA_real_,
    .groups = "drop"
  )
write.csv(cluster_sample_marker,
          file.path(output_dir, "cluster_sample_identity_marker_profiles.csv"),
          row.names = FALSE)

# Age summaries use sample-level values rather than pooling cells. With only
# two Young samples, these are descriptive and should not be treated as a
# definitive age test.
cluster_age_marker <- cluster_sample_marker %>%
  group_by(cluster, age, marker) %>%
  summarise(
    n_samples = n_distinct(sample),
    mean_fraction_positive_across_samples = mean(fraction_positive),
    sd_fraction_positive_across_samples = sd(fraction_positive),
    sem_fraction_positive_across_samples = sd(fraction_positive) / sqrt(n_samples),
    mean_log1p_across_samples = mean(mean_log1p),
    sd_mean_log1p_across_samples = sd(mean_log1p),
    sem_mean_log1p_across_samples = sd(mean_log1p) / sqrt(n_samples),
    .groups = "drop"
  )
write.csv(cluster_age_marker,
          file.path(output_dir, "cluster_age_identity_marker_profiles.csv"),
          row.names = FALSE)

# Within each cluster, correlate each sample's marker-positive-fraction profile
# with the pooled cluster profile. Low correlations flag inconsistent cluster
# identity and possible technical/sample effects.
pooled_reference <- cluster_marker %>%
  select(cluster, marker, pooled_fraction_positive = fraction_positive,
         pooled_mean_log1p = mean_log1p)
safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3 || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  cor(x[keep], y[keep])
}
sample_profile_correlations <- cluster_sample_marker %>%
  left_join(pooled_reference, by = c("cluster", "marker")) %>%
  group_by(cluster, sample, age) %>%
  summarise(
    n_cells = first(n_cells),
    correlation_fraction_positive = safe_cor(
      fraction_positive, pooled_fraction_positive
    ),
    correlation_mean_log1p = safe_cor(mean_log1p, pooled_mean_log1p),
    .groups = "drop"
  )
write.csv(sample_profile_correlations,
          file.path(output_dir, "cluster_sample_marker_profile_correlations.csv"),
          row.names = FALSE)

# Cluster-level marker heatmaps.
cluster_marker <- cluster_marker %>%
  left_join(cluster_annotation %>% select(cluster, cluster_marker_label), by = "cluster")
marker_order <- identity_markers_requested[identity_markers_requested %in% identity_markers]
cluster_label_order <- cluster_annotation$cluster_marker_label

p_fraction <- ggplot(
  cluster_marker,
  aes(x = factor(marker, levels = marker_order),
      y = factor(cluster_marker_label, levels = rev(cluster_label_order)),
      fill = fraction_positive)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(labels = scales::label_percent(accuracy = 1)) +
  labs(x = NULL, y = NULL, fill = "% positive",
       title = "Identity-marker positivity by fixed cluster") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(output_dir, "cluster_marker_fraction_positive_heatmap.png"),
       p_fraction, width = 11, height = 7, dpi = 300)

p_intensity <- ggplot(
  cluster_marker,
  aes(x = factor(marker, levels = marker_order),
      y = factor(cluster_marker_label, levels = rev(cluster_label_order)),
      fill = mean_log1p_z)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Marker z-score",
       title = "Relative log1p identity-marker intensity by fixed cluster") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(output_dir, "cluster_marker_log1p_z_heatmap.png"),
       p_intensity, width = 11, height = 7, dpi = 300)

# Show whether each cluster's assigned top marker is reproducible by sample and
# whether its distribution differs systematically by age.
top_marker_by_cluster <- cluster_annotation %>% select(cluster, top_marker)
top_marker_sample <- cluster_sample_marker %>%
  inner_join(top_marker_by_cluster, by = "cluster") %>%
  filter(marker == top_marker) %>%
  left_join(cluster_annotation %>% select(cluster, cluster_marker_label), by = "cluster")

p_top_sample <- ggplot(
  top_marker_sample,
  aes(x = cluster_marker_label, y = fraction_positive,
      color = sample, shape = age)
) +
  geom_point(size = 2.5, alpha = 0.85,
             position = position_jitter(width = 0.12, height = 0)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = NULL, y = "Top-marker-positive cells within cluster",
       title = "Assigned cluster marker across individual samples",
       color = "Sample", shape = "Age") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(output_dir, "cluster_top_marker_fraction_by_sample_age.png"),
       p_top_sample, width = 14, height = 7, dpi = 300)

# Put the top-marker annotation onto the existing visualization only; clusters
# and coordinates remain unchanged.
label_map <- setNames(cluster_annotation$cluster_marker_label,
                      as.character(cluster_annotation$cluster))
obj$cluster_top_marker <- factor(
  unname(label_map[as.character(obj$uncorrected_cluster)]),
  levels = cluster_annotation$cluster_marker_label
)
p_umap <- DimPlot(
  obj, reduction = "umap_rounded", group.by = "cluster_top_marker",
  label = TRUE, repel = TRUE
) + ggtitle("Fixed clusters labeled by most enriched identity marker")
ggsave(file.path(output_dir, "umap_cluster_top_marker_annotation.png"),
       p_umap, width = 11, height = 8, dpi = 300)

saveRDS(obj, "data/codex/codex_merged_log1p_uncorrected_marker_annotated.rds")
cat("Saved marker-based cluster annotation diagnostics to:", output_dir, "\n")
