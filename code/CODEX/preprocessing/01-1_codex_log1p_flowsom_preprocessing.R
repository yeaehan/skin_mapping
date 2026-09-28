# FlowSOM overclustering of the cleaned, log1p CODEX data.
# Run after the CODEX SCE has been prepared:
#   Rscript code/CODEX/preprocessing/01-1_codex_log1p_flowsom_preprocessing.R

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(FlowSOM)
  library(dplyr)
  library(ggplot2)
})

project_dir <- normalizePath(
  "//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2",
  winslash = "/", mustWork = TRUE
)
setwd(project_dir)

input_file <- "data/codex/log1p_flowsom_comparison/codex_log1p_preprocessed_sce.rds"
data_dir <- "data/codex/log1p_flowsom_comparison"
results_dir <- "results/codex/log1p_flowsom_comparison"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

sce <- readRDS(input_file)
type_markers <- metadata(sce)$preprocessing$type_markers
eligible <- as.logical(colData(sce)$flowsom_eligible)
stopifnot(length(type_markers) > 2, any(eligible), !"HLA-DR" %in% type_markers)

flow_input <- t(assay(sce, "log1p")[type_markers, eligible, drop = FALSE])
set.seed(1234)
fsom <- FlowSOM(
  flow_input,
  colsToUse = seq_along(type_markers),
  scale = TRUE,
  xdim = 10,
  ydim = 10,
  nClus = 40,
  seed = 1234,
  silent = FALSE
)

som_cluster <- as.integer(fsom$map$mapping[, 1])
metacluster <- as.integer(fsom$metaclustering[som_cluster])
stopifnot(length(metacluster) == sum(eligible), !anyNA(metacluster))

all_som <- all_meta <- rep(NA_integer_, ncol(sce))
all_som[eligible] <- som_cluster
all_meta[eligible] <- metacluster
colData(sce)$flowsom_cluster <- factor(all_som)
colData(sce)$flowsom_metacluster <- factor(all_meta)

assignments <- data.frame(
  cell_id = colnames(sce),
  sample = as.character(colData(sce)$codex_sample),
  age = as.character(colData(sce)$age),
  qc_status = as.character(colData(sce)$qc_status),
  flowsom_cluster = all_som,
  flowsom_metacluster = all_meta
)
stopifnot(all(assignments$age %in% c("Young", "Old")))

assigned <- !is.na(all_meta)
raw_type <- t(assay(sce, "raw_intensity")[type_markers, assigned, drop = FALSE])
log_type <- t(assay(sce, "log1p")[type_markers, assigned, drop = FALSE])
assigned_meta <- all_meta[assigned]

profiles <- bind_rows(lapply(sort(unique(assigned_meta)), function(k) {
  cells <- assigned_meta == k
  data.frame(
    flowsom_metacluster = k,
    marker = type_markers,
    n_cells = sum(cells),
    fraction_positive = colMeans(raw_type[cells, , drop = FALSE] > 0),
    mean_log1p = colMeans(log_type[cells, , drop = FALSE])
  )
}))

assigned_cells <- assignments[assigned, ]
composition_sample <- count(
  assigned_cells, sample, flowsom_metacluster, name = "n_cells"
) %>%
  group_by(sample) %>%
  mutate(fraction_of_sample = n_cells / sum(n_cells)) %>%
  ungroup()
composition_age <- count(
  assigned_cells, age, flowsom_metacluster, name = "n_cells"
) %>%
  group_by(age) %>%
  mutate(fraction_of_age = n_cells / sum(n_cells)) %>%
  ungroup()

write.csv(assignments, file.path(results_dir, "flowsom_assignments.csv"), row.names = FALSE)
write.csv(profiles, file.path(results_dir, "flowsom_metacluster_marker_profiles.csv"), row.names = FALSE)
write.csv(composition_sample, file.path(results_dir, "flowsom_metacluster_composition_by_sample.csv"), row.names = FALSE)
write.csv(composition_age, file.path(results_dir, "flowsom_metacluster_composition_by_age.csv"), row.names = FALSE)
saveRDS(fsom, file.path(data_dir, "codex_log1p_flowsom_model.rds"))
saveRDS(sce, file.path(data_dir, "codex_log1p_flowsom_sce.rds"))

p <- ggplot(
  profiles,
  aes(
    factor(marker, levels = type_markers),
    factor(flowsom_metacluster),
    fill = fraction_positive
  )
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(labels = scales::label_percent(accuracy = 1)) +
  labs(x = NULL, y = "FlowSOM metacluster", fill = "% positive") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

ggsave(
  file.path(results_dir, "flowsom_metacluster_marker_heatmap.png"),
  p, width = 11, height = 9, dpi = 300
)

cat(
  "FlowSOM complete:", sum(eligible), "of", ncol(sce), "cells assigned across",
  length(unique(metacluster)), "occupied metaclusters.\n"
)
