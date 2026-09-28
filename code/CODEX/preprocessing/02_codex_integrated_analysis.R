##########################################################
# 7) CODEX integration + better-looking UMAP
##########################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(patchwork)

# Resolve the project root from this script's location so relative data/results
# paths work even when Rscript is launched from another directory.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(file_arg)) {
  stop("Could not determine this script's path; run it with Rscript")
}
script_path <- normalizePath(
  sub("^--file=", "", file_arg[1]),
  winslash = "/",
  mustWork = TRUE
)
project_root <- normalizePath(
  file.path(dirname(script_path), "..", "..", ".."),
  winslash = "/",
  mustWork = TRUE
)
setwd(project_root)
cat("Project root:", project_root, "\n")

load("data/codex/seurat_list_codex.RData")  # loads seurat_list_codex

##########################################################
# 0) Per-sample QC before merging
#
# The inputs are zero-thresholded, so a cell is removed only when every marker
# available in that sample is zero. Filtering is done before merge so missing
# optional channels in one specimen cannot be interpreted as biological zeros.
##########################################################

qc_output_dir <- "results/codex/preprocessing"
dir.create(qc_output_dir, recursive = TRUE, showWarnings = FALSE)

premerge_zero_qc <- do.call(rbind, lapply(names(seurat_list_codex), function(sample_id) {
  obj <- seurat_list_codex[[sample_id]]
  DefaultAssay(obj) <- "RNA"
  obj$original_sample <- sample_id
  obj$source_age <- as.character(obj$age)
  # Project convention established for CODEX comparisons: only 207UA and 7503
  # are Young; retain source_age separately so the provenance is never lost.
  obj$analysis_age <- ifelse(sample_id %in% c("207UA", "7503"), "Young", "Old")
  obj$analysis_age <- factor(obj$analysis_age, levels = c("Young", "Old"))
  counts <- LayerData(obj, assay = "RNA", layer = "counts")
  n_positive_markers <- Matrix::colSums(counts > 0)
  all_panel_zero <- n_positive_markers == 0
  n_original <- ncol(obj)
  n_removed <- sum(all_panel_zero)

  # Retain per-cell QC in the filtered object for downstream auditing.
  obj$n_panel_markers_positive <- as.integer(n_positive_markers)
  obj$all_panel_zero <- all_panel_zero
  keep_cells <- names(n_positive_markers)[!all_panel_zero]
  if (!length(keep_cells)) stop("All cells are panel-zero in sample ", sample_id)
  seurat_list_codex[[sample_id]] <<- subset(obj, cells = keep_cells)

  data.frame(
    sample = sample_id,
    n_markers_available = nrow(counts),
    n_cells_original = n_original,
    n_all_panel_zero_removed = n_removed,
    fraction_all_panel_zero_removed = n_removed / n_original,
    percent_all_panel_zero_removed = 100 * n_removed / n_original,
    n_cells_retained = n_original - n_removed,
    fraction_cells_retained = (n_original - n_removed) / n_original,
    percent_cells_retained = 100 * (n_original - n_removed) / n_original,
    stringsAsFactors = FALSE
  )
}))

write.csv(
  premerge_zero_qc,
  file.path(qc_output_dir, "codex_premerge_all_panel_zero_qc_by_sample.csv"),
  row.names = FALSE
)
save(
  seurat_list_codex,
  file = "data/codex/seurat_list_codex_nonzero_filtered.RData"
)
print(premerge_zero_qc)
cat("Saved the per-sample all-panel-zero QC audit and filtered list.\n")

# 
# # ---- (A) Minimal per-sample preprocessing for integration ----
# seurat_list_codex <- lapply(seurat_list_codex, function(obj) {
#   DefaultAssay(obj) <- "RNA"
#   
#   # optional: drop super-low-signal cells (prevents distance artifacts)
#   # tune threshold if needed; CODEX intensities can be large, so start ~100
#   tot <- Matrix::colSums(GetAssayData(obj, layer = "counts"))
#   obj <- subset(obj, cells = names(tot)[tot > 100])
#   
#   obj <- NormalizeData(obj, verbose = FALSE)
#   obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = min(200, nrow(obj)), verbose = FALSE)
#   obj
# })
# 
# features <- Reduce(intersect, lapply(seurat_list_codex, rownames))
# 
# # scale + PCA per object (use only selected features)
# seurat_list_codex <- lapply(seurat_list_codex, function(obj) {
#   obj <- ScaleData(obj, features = features, verbose = FALSE)
#   npcs_use <- min(20, length(features) - 1)
#   obj <- RunPCA(obj, features = features, npcs = npcs_use, verbose = FALSE)
#   obj
# })
# 
# # ---- (B) Find anchors + integrate ----
# anchors <- FindIntegrationAnchors(
#   object.list = seurat_list_codex,
#   anchor.features = features,
#   reduction = "rpca",     # usually more stable for small panels
#   dims = 1:15
# )
# 
# codex_int <- IntegrateData(anchorset = anchors, dims = 1:15)
# DefaultAssay(codex_int) <- "integrated"
# 
# codex_int <- ScaleData(codex_int, verbose = FALSE)
# 
# codex_int <- RunPCA(
#   codex_int,
#   npcs = 10,      # very important: keep PCs small
#   verbose = FALSE
# )
# 
# dims_use <- 1:3   # use only first few PCs
# 
# codex_int <- FindNeighbors(
#   codex_int,
#   dims = dims_use,
#   k.param = 100,   # larger neighborhood smooths sparse data
#   verbose = FALSE
# )
# 
# codex_int <- RunUMAP(
#   codex_int,
#   dims = dims_use,
#   umap.method = "uwot",
#   metric = "cosine",
#   n.neighbors = 100,
#   min.dist = 0.6,   # prevents spaghetti structure
#   spread = 1.3,
#   verbose = FALSE
# )
# 
# 
# # ---- (D) Plots ----
# p_major <- DimPlot(codex_int, group.by = "majorskincell", label = TRUE) + NoLegend()
# p_minor <- DimPlot(codex_int, group.by = "minorskincell1", label = TRUE) + NoLegend()
# p_snc   <- DimPlot(codex_int, group.by = "snc") + NoLegend()
# p_age   <- DimPlot(codex_int, group.by = "age") + NoLegend()
# 
# p_major + p_minor + p_snc + p_age
# 
# # ---- (E) Save integrated object ----
# saveRDS(codex_int, file = "data/codex/codex_integrated_umap.rds")


##########################################################
# 1) Merge filtered samples and construct an uncorrected log1p baseline
##########################################################

# Only markers present in every specimen can be compared jointly without
# interpreting an unmeasured channel as zero.
common_markers <- Reduce(intersect, lapply(seurat_list_codex, rownames))
identity_markers_requested <- c(
  "SMA", "Keratin 14", "Vimentin", "Keratin 5", "CD34", "CD8", "CD14",
  "Filaggrin", "FOXP3", "MelanA", "CD4", "CD68", "HLA-DR", "CD45"
)
identity_markers <- intersect(identity_markers_requested, common_markers)
if (length(identity_markers) < 3) {
  stop("Too few shared identity markers for joint clustering")
}

sample_ids <- names(seurat_list_codex)
codex_merged <- merge(
  seurat_list_codex[[1]],
  y = seurat_list_codex[-1],
  add.cell.ids = sample_ids,
  merge.data = FALSE
)
DefaultAssay(codex_merged) <- "RNA"
codex_merged <- JoinLayers(codex_merged, assay = "RNA")

raw_merged <- LayerData(codex_merged, assay = "RNA", layer = "counts")
raw_common <- raw_merged[common_markers, , drop = FALSE]
log1p_common <- log1p(raw_common)

# Store the transformed protein intensities in a separate assay. No CLR,
# library-size normalization or batch correction is applied at this stage.
codex_merged[["CODEXlog1p"]] <- CreateAssay5Object(data = log1p_common)
DefaultAssay(codex_merged) <- "CODEXlog1p"
VariableFeatures(codex_merged) <- identity_markers

# Scaling is across markers for PCA/clustering; the unscaled log1p values remain
# in the data layer for interpretation.
codex_merged <- ScaleData(
  codex_merged,
  features = identity_markers,
  verbose = FALSE
)
npcs_use <- min(10L, length(identity_markers) - 1L)
codex_merged <- RunPCA(
  codex_merged,
  assay = "CODEXlog1p",
  features = identity_markers,
  npcs = npcs_use,
  reduction.name = "pca_uncorrected",
  verbose = FALSE
)
dims_use <- seq_len(min(8L, npcs_use))
codex_merged <- FindNeighbors(
  codex_merged,
  reduction = "pca_uncorrected",
  dims = dims_use,
  k.param = 50,
  graph.name = c("uncorrected_nn", "uncorrected_snn"),
  verbose = FALSE
)
codex_merged <- FindClusters(
  codex_merged,
  graph.name = "uncorrected_snn",
  resolution = 0.5,
  cluster.name = "uncorrected_cluster",
  random.seed = 1234,
  verbose = FALSE
)
codex_merged <- RunUMAP(
  codex_merged,
  reduction = "pca_uncorrected",
  dims = dims_use,
  reduction.name = "umap_uncorrected",
  reduction.key = "UMAPUNCORRECTED_",
  metric = "cosine",
  n.neighbors = 50,
  min.dist = 0.4,
  seed.use = 1234,
  verbose = FALSE
)

# Alternative visualization coordinates only. These calls do not rebuild the
# neighbor graph or modify uncorrected_cluster.
codex_merged <- RunUMAP(
  codex_merged,
  reduction = "pca_uncorrected",
  dims = dims_use,
  reduction.name = "umap_rounded",
  reduction.key = "UMAPROUNDED_",
  metric = "euclidean",
  n.neighbors = 30,
  min.dist = 0.2,
  spread = 1,
  seed.use = 1234,
  verbose = FALSE
)
codex_merged <- RunTSNE(
  codex_merged,
  reduction = "pca_uncorrected",
  dims = dims_use,
  reduction.name = "tsne_visual",
  reduction.key = "TSNEVISUAL_",
  perplexity = 30,
  seed.use = 1234,
  check_duplicates = FALSE,
  verbose = FALSE
)

##########################################################
# 2) Diagnose sample/age structure before any batch correction
##########################################################

diagnostic_dir <- "results/codex/preprocessing/uncorrected_log1p"
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

p_cluster <- DimPlot(
  codex_merged, reduction = "umap_uncorrected",
  group.by = "uncorrected_cluster", label = TRUE, repel = TRUE
) + ggtitle("Uncorrected log1p: clusters")
p_sample <- DimPlot(
  codex_merged, reduction = "umap_uncorrected", group.by = "original_sample"
) + ggtitle("Uncorrected log1p: sample")
p_age <- DimPlot(
  codex_merged, reduction = "umap_uncorrected", group.by = "analysis_age"
) + ggtitle("Uncorrected log1p: analysis age")

p_cluster_rounded <- DimPlot(
  codex_merged, reduction = "umap_rounded",
  group.by = "uncorrected_cluster", label = TRUE, repel = TRUE
) + ggtitle("Compact UMAP: same clusters")
p_sample_rounded <- DimPlot(
  codex_merged, reduction = "umap_rounded", group.by = "original_sample"
) + ggtitle("Compact UMAP: sample")
p_age_rounded <- DimPlot(
  codex_merged, reduction = "umap_rounded", group.by = "analysis_age"
) + ggtitle("Compact UMAP: analysis age")

p_cluster_tsne <- DimPlot(
  codex_merged, reduction = "tsne_visual",
  group.by = "uncorrected_cluster", label = TRUE, repel = TRUE
) + ggtitle("t-SNE: same clusters")
p_sample_tsne <- DimPlot(
  codex_merged, reduction = "tsne_visual", group.by = "original_sample"
) + ggtitle("t-SNE: sample")
p_age_tsne <- DimPlot(
  codex_merged, reduction = "tsne_visual", group.by = "analysis_age"
) + ggtitle("t-SNE: analysis age")

ggsave(file.path(diagnostic_dir, "umap_by_cluster.png"), p_cluster,
       width = 8, height = 6, dpi = 300)
ggsave(file.path(diagnostic_dir, "umap_by_sample.png"), p_sample,
       width = 9, height = 7, dpi = 300)
ggsave(file.path(diagnostic_dir, "umap_by_age.png"), p_age,
       width = 8, height = 6, dpi = 300)

cluster_visualization_panel <- p_cluster + p_cluster_rounded + p_cluster_tsne +
  plot_layout(ncol = 3, guides = "collect")
sample_visualization_panel <- p_sample + p_sample_rounded + p_sample_tsne +
  plot_layout(ncol = 3, guides = "collect")
age_visualization_panel <- p_age + p_age_rounded + p_age_tsne +
  plot_layout(ncol = 3, guides = "collect")

ggsave(
  file.path(diagnostic_dir, "cluster_visualization_umap_compact_tsne.png"),
  cluster_visualization_panel, width = 20, height = 6.5, dpi = 300
)
ggsave(
  file.path(diagnostic_dir, "sample_visualization_umap_compact_tsne.png"),
  sample_visualization_panel, width = 22, height = 7, dpi = 300
)
ggsave(
  file.path(diagnostic_dir, "age_visualization_umap_compact_tsne.png"),
  age_visualization_panel, width = 20, height = 6.5, dpi = 300
)

cluster_sample_counts <- as.data.frame(table(
  cluster = codex_merged$uncorrected_cluster,
  sample = codex_merged$original_sample
))
cluster_sample_counts <- cluster_sample_counts %>%
  group_by(cluster) %>%
  mutate(
    total_cells_in_cluster = sum(Freq),
    fraction_within_cluster = Freq / total_cells_in_cluster,
    percent_within_cluster = 100 * fraction_within_cluster
  ) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(
    total_cells_in_sample = sum(Freq),
    fraction_of_sample_in_cluster = Freq / total_cells_in_sample,
    percent_of_sample_in_cluster = 100 * fraction_of_sample_in_cluster
  ) %>%
  ungroup()
write.csv(cluster_sample_counts,
          file.path(diagnostic_dir, "cluster_composition_by_sample.csv"),
          row.names = FALSE)

cluster_age_counts <- as.data.frame(table(
  cluster = codex_merged$uncorrected_cluster,
  age = codex_merged$analysis_age
))
cluster_age_counts <- cluster_age_counts %>%
  group_by(cluster) %>%
  mutate(
    total_cells_in_cluster = sum(Freq),
    fraction_within_cluster = Freq / total_cells_in_cluster,
    percent_within_cluster = 100 * fraction_within_cluster
  ) %>%
  ungroup() %>%
  group_by(age) %>%
  mutate(
    total_cells_in_age_group = sum(Freq),
    fraction_of_age_group_in_cluster = Freq / total_cells_in_age_group,
    percent_of_age_group_in_cluster = 100 * fraction_of_age_group_in_cluster
  ) %>%
  ungroup()
write.csv(cluster_age_counts,
          file.path(diagnostic_dir, "cluster_composition_by_age.csv"),
          row.names = FALSE)

write.csv(
  data.frame(marker = identity_markers),
  file.path(diagnostic_dir, "identity_markers_used.csv"),
  row.names = FALSE
)
saveRDS(
  codex_merged,
  file = "data/codex/codex_merged_log1p_uncorrected.rds"
)
cat("Saved the uncorrected log1p clustering baseline and diagnostics.\n")

##########################################################
# CODEX merged object: add SnC + major/minor skin cell types
# (works with Seurat v5 layers; uses RNA counts > 0 like 4i/CODEX)
##########################################################

library(dplyr)
library(Seurat)

DefaultAssay(codex_merged) <- "RNA"

# --- v5-safe counts matrix ---
mat_counts <- LayerData(codex_merged, assay = "RNA", layer = "counts")

# helper: positivity for a feature; returns all FALSE if missing
get_pos_merged <- function(feature) {
  if (!feature %in% rownames(mat_counts)) {
    return(setNames(rep(FALSE, ncol(codex_merged)), colnames(codex_merged)))
  }
  v <- mat_counts[feature, ]
  setNames(as.numeric(v) > 0, names(v))
}

# =========================
# 1) SnC annotation
# =========================
p21   <- get_pos_merged("p21")
p16   <- get_pos_merged("p16")
PCNA  <- get_pos_merged("PCNA")
IP10  <- get_pos_merged("IP10")
CXCR3 <- get_pos_merged("CXCR3")

codex_merged$snc <- ifelse(((p21 | p16) & (!PCNA) & (IP10 | CXCR3)), "SnC+", "SnC-")
codex_merged$snc <- factor(codex_merged$snc, levels = c("SnC-", "SnC+"))

# (optional) keep sencelltype name too if you want
# codex_merged$sencelltype <- codex_merged$snc

# =========================
# 2) Major/minor skin cell types (CODEX convention)
# =========================
SMA        <- get_pos_merged("SMA")
MelanA     <- get_pos_merged("MelanA")
Vimentin   <- get_pos_merged("Vimentin")
Keratin5   <- get_pos_merged("Keratin 5")
Keratin14  <- get_pos_merged("Keratin 14")
Filaggrin  <- get_pos_merged("Filaggrin")
PanCyto    <- get_pos_merged("Pan-Cytokeratin")

codex_merged$majorskincell <- dplyr::case_when(
  MelanA ~ "Melanocyte",
  SMA ~ "Smooth Muscle",
  Vimentin ~ "Fibroblast",
  (Keratin5 | Keratin14 | Filaggrin) ~ "Keratinocyte",
  TRUE ~ "none"
)

codex_merged$minorskincell1 <- dplyr::case_when(
  MelanA ~ "Melanocyte",
  SMA ~ "Smooth Muscle",
  Vimentin ~ "Fibroblast",
  (Keratin5 | Keratin14) ~ "Basal Keratinocyte",
  Filaggrin ~ "Corneocyte",
  !Filaggrin & !Keratin5 & !Keratin14 & PanCyto ~ "Suprabasal Keratinocyte",
  TRUE ~ "none"
)

# inverse anatomical factor order (deep -> superficial), "none" first
codex_merged$majorskincell <- factor(
  codex_merged$majorskincell,
  levels = c("none", "Smooth Muscle", "Fibroblast", "Melanocyte", "Keratinocyte")
)
codex_merged$minorskincell1 <- factor(
  codex_merged$minorskincell1,
  levels = c("none", "Smooth Muscle", "Fibroblast", "Melanocyte",
             "Basal Keratinocyte", "Suprabasal Keratinocyte", "Corneocyte")
)

# =========================
# 3) Quick sanity checks
# =========================
table(codex_merged$snc)
table(codex_merged$majorskincell)
table(codex_merged$minorskincell1)

# =========================
# 4) Plot on your UMAP
# =========================
p_major <- DimPlot(codex_merged, group.by = "majorskincell", label = TRUE) + NoLegend()
p_minor <- DimPlot(codex_merged, group.by = "minorskincell1", label = TRUE) + NoLegend()
p_snc   <- DimPlot(codex_merged, group.by = "snc") + NoLegend()
p_major + p_minor + p_snc



p_snc   <- DimPlot(codex_merged, group.by = "snc") + NoLegend()


# SnC was assigned on the merged object above, so summarize that metadata
# rather than the pre-annotation per-sample list.
df_codex <- codex_merged@meta.data %>%
  as.data.frame() %>%
  group_by(original_sample, analysis_age) %>%
  summarise(
    n_cells = n(),
    n_snc = sum(snc == "SnC+", na.rm = TRUE),
    snc_prop = n_snc / n_cells,
    .groups = "drop"
  ) %>%
  transmute(
    sample = original_sample,
    age = factor(as.character(analysis_age), levels = c("Young", "Old")),
    n_cells,
    n_snc,
    snc_prop
  )

blue_col <- "#296d9e"

ymax <- max(df_codex$snc_prop, na.rm = TRUE) * 1.15
if (!is.finite(ymax) || ymax <= 0) ymax <- 0.01

p_age <- ggplot(df_codex, aes(x = age, y = snc_prop)) +
  stat_summary(fun = mean, geom = "bar", width = 0.6,
               fill = blue_col, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15,
               linewidth = 0.8, color = "black") +
  geom_point(position = position_jitter(width = 0.07),
             size = 3, color = "black") +
  scale_y_continuous(
    limits = c(0, ymax),
    breaks = pretty(c(0, ymax), n = 4),
    labels = scales::label_percent(accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 18)

print(p_age)

dir.create("results/260225", recursive = TRUE, showWarnings = FALSE)
ggsave("results/260225/codex_snc_age.pdf", plot = p_age, width = 3.6, height = 4.2, dpi = 500)
