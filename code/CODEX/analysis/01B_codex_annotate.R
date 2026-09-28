##########################################################
# Focused, reproducible CODEX annotation stage
#
# Applies the annotation rules from 01_codex_cell_type_analysis.R to the
# freshly preprocessed list, then saves seurat_list_codex_annotated.RData.
##########################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

setwd("//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2")
load("data/codex/seurat_list_codex.RData")

marker_panel_codex <- c(
  "BCL2", "SMA", "p21", "Keratin 14", "Vimentin", "Keratin 5", "p16", "CD34",
  "IP10", "CD8", "hIL-1alpha", "CD14", "LaminB1", "Ki67", "Filaggrin", "CXCR3",
  "HMGB1", "Sirt1", "IDO1", "FOXP3", "IL6", "MelanA", "IL8", "CD4",
  "CD68", "HLA-DR", "CD45", "PCNA", "ICOS", "Pan-Cytokeratin", "CD3e", "CD11c"
)

get_positive <- function(obj, marker) {
  cells <- colnames(obj)
  if (!marker %in% rownames(obj)) {
    warning("Marker ", marker, " not found in ", unique(obj$sample))
    return(setNames(rep(FALSE, length(cells)), cells))
  }
  values <- LayerData(obj, assay = "RNA", layer = "counts")[marker, cells, drop = TRUE]
  setNames(as.numeric(values) > 0, cells)
}

for (sample_name in names(seurat_list_codex)) {
  cat("Annotating:", sample_name, "\n")
  obj <- seurat_list_codex[[sample_name]]
  cells <- colnames(obj)

  expected_p16_source <- if (tolower(sample_name) == "cams18a") "p16 layer 9" else "p16-BC42"
  actual_p16_source <- unique(as.character(obj$p16_source))
  if (length(actual_p16_source) != 1 || actual_p16_source != expected_p16_source) {
    stop(
      "Unexpected p16 source for ", sample_name, ": ",
      paste(actual_p16_source, collapse = ";"), " (expected ", expected_p16_source, ")"
    )
  }

  p21 <- get_positive(obj, "p21")
  p16 <- get_positive(obj, "p16")
  PCNA <- get_positive(obj, "PCNA")
  IP10 <- get_positive(obj, "IP10")
  IL6 <- get_positive(obj, "IL6")
  IL8 <- get_positive(obj, "IL8")
  CXCR3 <- get_positive(obj, "CXCR3")
  MelanA <- get_positive(obj, "MelanA")
  Vimentin <- get_positive(obj, "Vimentin")
  SMA <- get_positive(obj, "SMA")
  Keratin5 <- get_positive(obj, "Keratin 5")
  Keratin14 <- get_positive(obj, "Keratin 14")
  Filaggrin <- get_positive(obj, "Filaggrin")
  PanCyto <- get_positive(obj, "Pan-Cytokeratin")
  CD4 <- get_positive(obj, "CD4")
  CD68 <- get_positive(obj, "CD68")
  CD8 <- get_positive(obj, "CD8")
  FOXP3 <- get_positive(obj, "FOXP3")
  CD14 <- get_positive(obj, "CD14")
  CD45 <- get_positive(obj, "CD45")

  obj$cdkn_state <- case_when(
    p16 & !p21 ~ "p16+",
    p21 & !p16 ~ "p21+",
    p16 & p21 ~ "p16+p21+",
    TRUE ~ "p16-p21-"
  )
  obj$cdkn_state <- factor(
    obj$cdkn_state,
    levels = c("p16-p21-", "p21+", "p16+", "p16+p21+")
  )

  obj$snc <- ifelse(
    (p21 | p16) & !PCNA & (IP10 | IL6 | IL8),
    "SnC+", "SnC-"
  )
  obj$snc <- factor(obj$snc, levels = c("SnC-", "SnC+"))

  obj$majorskincell <- case_when(
    MelanA ~ "Melanocyte",
    SMA ~ "Smooth Muscle",
    Vimentin ~ "Fibroblast",
    Keratin5 | Keratin14 | Filaggrin ~ "Keratinocyte",
    TRUE ~ "none"
  )
  obj$minorskincell1 <- case_when(
    MelanA ~ "Melanocyte",
    SMA ~ "Smooth Muscle",
    Vimentin ~ "Fibroblast",
    Keratin5 | Keratin14 ~ "Basal Keratinocyte",
    Filaggrin ~ "Corneocyte",
    !Filaggrin & !Keratin5 & !Keratin14 & PanCyto ~ "Suprabasal Keratinocyte",
    TRUE ~ "none"
  )
  obj$majorskincell <- factor(
    obj$majorskincell,
    levels = c("none", "Smooth Muscle", "Fibroblast", "Melanocyte", "Keratinocyte")
  )
  obj$minorskincell1 <- factor(
    obj$minorskincell1,
    levels = c("none", "Smooth Muscle", "Fibroblast", "Melanocyte",
               "Basal Keratinocyte", "Suprabasal Keratinocyte", "Corneocyte")
  )

  obj$immcelltype <- case_when(
    CD4 & FOXP3 ~ "CD4+ FOXP3+",
    CD4 & CXCR3 ~ "CD4+ CXCR3+",
    CD14 ~ "CD14+",
    CD68 ~ "CD68+",
    CD8 ~ "CD8+",
    CD4 ~ "CD4+",
    CD45 ~ "CD45+",
    TRUE ~ "none"
  )

  if ("spatial" %in% names(obj@reductions)) {
    coords <- Embeddings(obj, reduction = "spatial")
    obj$Spatial_X <- coords[cells, 1]
    obj$Spatial_Y <- coords[cells, 2]
  }
  obj$age <- factor(as.character(obj$age), levels = c("Young", "Old"))
  obj$age_major <- paste(obj$age, obj$majorskincell, sep = " ")
  obj$age_minor <- paste(obj$age, obj$minorskincell1, sep = " ")
  seurat_list_codex[[sample_name]] <- obj
}

audit <- bind_rows(lapply(names(seurat_list_codex), function(sample_name) {
  obj <- seurat_list_codex[[sample_name]]
  data.frame(
    sample = sample_name,
    age = as.character(unique(obj$age)),
    p16_source = unique(as.character(obj$p16_source)),
    n_cells = ncol(obj),
    n_p16_positive = sum(obj$cdkn_state %in% c("p16+", "p16+p21+")),
    n_p21_positive = sum(obj$cdkn_state %in% c("p21+", "p16+p21+")),
    n_snc_positive = sum(obj$snc == "SnC+")
  )
}))

write.csv(audit, "data/codex/codex_annotation_audit.csv", row.names = FALSE)
save(seurat_list_codex, file = "data/codex/seurat_list_codex_annotated.RData")
print(audit)
cat("Saved corrected CODEX annotations.\n")
