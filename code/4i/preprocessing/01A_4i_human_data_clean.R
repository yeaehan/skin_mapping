
# =========================
# Excel cleanup functions (4i)
# =========================


##########################################################
#I removed SIRT1
##########################################################
install.packages(c("readxl", "dplyr", "stringr"))

# 2. Install SeuratDisk from GitHub
# First, you need the 'remotes' package to download things from GitHub
if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Next, pull SeuratDisk directly from the creator's repository
remotes::install_github("mojaveazure/seurat-disk")

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(SeuratDisk)
})

setwd("//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2")

load("data/4i/seurat_list_4i.RData")


# Canonical markers you want downstream
marker_panel <- c(
  "p21", "p16", "PCNA", "HMGB1", "LaminB1", "BCL2", "IP10", "CXCR3",
  "Keratin 5", "Keratin 10", "MelanA", "Filaggrin", "Vimentin"
)

# Root-name map: messy variants -> canonical
marker_root_map <- c(
  "p21" = "p21",
  "p16" = "p16",
  "pcna" = "PCNA",
  "hmgb1" = "HMGB1",
  "laminb1" = "LaminB1",
  "laminb" = "LaminB1",
  "bcl2" = "BCL2",
  "ip10" = "IP10",
  "cxcr3" = "CXCR3",
  "k5" = "Keratin 5",
  "keratin5" = "Keratin 5",
  "k10" = "Keratin 10",
  "keratin10" = "Keratin 10",
  "melana" = "MelanA",
  "filaggrin" = "Filaggrin",
  "vimentin" = "Vimentin"
)

get_age_from_path <- function(path) {
  if (grepl("[/\\\\]old[/\\\\]", path, ignore.case = TRUE)) return("old")
  if (grepl("[/\\\\]young[/\\\\]", path, ignore.case = TRUE)) return("young")
  NA_character_
}

normalize_token <- function(x) {
  x %>%
    str_replace(regex("^MeanIntensity[_ ]*", ignore_case = TRUE), "") %>%
    tolower() %>%
    str_replace_all("[^a-z0-9]", "")
}

pick_sheet <- function(path) {
  s <- excel_sheets(path)
  s_use <- s[!grepl("keratinocyte", s, ignore.case = TRUE)]
  i <- grep("^nuclei\\s*props$|nuclei.*props|nucleiprops", s_use, ignore.case = TRUE)
  if (length(i) > 0) return(s_use[i[1]])
  if ("Sheet1" %in% s_use) return("Sheet1")
  stop(sprintf("No nucleiProps-like sheet or Sheet1 in: %s", basename(path)))
}

pick_coord_cols <- function(df_raw) {
  nms <- names(df_raw)
  roots <- normalize_token(nms)
  if (all(c("umap1", "umap2") %in% roots)) return(c(nms[match("umap1", roots)], nms[match("umap2", roots)]))
  if (all(c("x", "y") %in% roots)) return(c(nms[match("x", roots)], nms[match("y", roots)]))
  NULL
}

rename_columns <- function(df) {
  old <- names(df)
  root <- normalize_token(old)
  new <- old

  # metadata names
  new[root == "label"] <- "Label"
  new[root == "x"] <- "X"
  new[root == "y"] <- "Y"

  # marker names
  hit <- root %in% names(marker_root_map)
  new[hit] <- unname(marker_root_map[root[hit]])

  names(df) <- new
  df
}

coalesce_duplicate_columns <- function(df) {
  nms <- names(df)
  out <- list()
  for (nm in unique(nms)) {
    idx <- which(nms == nm)
    if (length(idx) == 1) out[[nm]] <- df[[idx]]
    if (length(idx) > 1) out[[nm]] <- Reduce(dplyr::coalesce, df[idx])
  }
  as.data.frame(out, check.names = FALSE)
}

remove_footer_rows <- function(df, marker_cols) {
  is_footer_label <- grepl(
    "threshold|thresh|cutoff|summary|mean|median|sd|note|background",
    df$Label, ignore.case = TRUE
  )

  marker_num_count <- rowSums(!is.na(as.data.frame(lapply(df[, marker_cols, drop = FALSE], function(z) suppressWarnings(as.numeric(z))))))
  too_few_numeric <- marker_num_count < 2

  df[!(is_footer_label | too_few_numeric), , drop = FALSE]
}

clean_one_excel <- function(path) {
  age <- get_age_from_path(path)
  
  sh <- pick_sheet(path)
  raw <- read_excel(path, sheet = sh)
  raw <- as.data.frame(raw, check.names = FALSE)
  
  coord_cols <- pick_coord_cols(raw)
  df <- rename_columns(raw)
  df <- coalesce_duplicate_columns(df)
  
  if ("keratinocytes" %in% names(df)) df$keratinocytes <- NULL
  if (!"Label" %in% names(df)) stop(sprintf("Missing Label in %s", basename(path)))
  
  present_markers <- intersect(marker_panel, names(df))
  if (length(present_markers) == 0) stop(sprintf("No markers found in %s", basename(path)))
  
  keep <- c("Label", present_markers)
  df <- df[, keep, drop = FALSE]
  
  df$Label <- as.character(df$Label)
  df <- df[!is.na(df$Label) & df$Label != "", , drop = FALSE]
  df <- remove_footer_rows(df, present_markers)
  df$Label <- make.unique(df$Label)
  
  for (m in present_markers) df[[m]] <- suppressWarnings(as.numeric(df[[m]]))
  all_na <- rowSums(is.na(df[, present_markers, drop = FALSE])) == length(present_markers)
  df <- df[!all_na, , drop = FALSE]
  
  # replace remaining NA/NaN/Inf in marker intensities with 0
  for (m in present_markers) {
    v <- suppressWarnings(as.numeric(df[[m]]))
    v[is.na(v) | is.nan(v) | is.infinite(v)] <- 0
    df[[m]] <- v
  }
  
  
  spatial_coords <- NULL
  if (!is.null(coord_cols) && all(c("Label", coord_cols) %in% names(raw))) {
    coord_df <- raw[, c("Label", coord_cols), drop = FALSE]
    names(coord_df) <- c("Label", "spatial_1", "spatial_2")
    coord_df$Label <- as.character(coord_df$Label)
    coord_df <- coord_df[!is.na(coord_df$Label) & coord_df$Label != "", , drop = FALSE]
    coord_df$Label <- make.unique(coord_df$Label)
    coord_df <- coord_df[coord_df$Label %in% df$Label, , drop = FALSE]
    rownames(coord_df) <- coord_df$Label
    spatial_coords <- coord_df[df$Label, c("spatial_1", "spatial_2"), drop = FALSE]
    spatial_coords$spatial_1 <- suppressWarnings(as.numeric(spatial_coords$spatial_1))
    spatial_coords$spatial_2 <- suppressWarnings(as.numeric(spatial_coords$spatial_2))
  }
  
  rownames(df) <- df$Label
  df$Label <- NULL
  
  list(
    data = df,
    spatial_coords = spatial_coords,
    sheet = sh,
    markers = present_markers,
    age = age,
    source_path = path
  )
}


clean_many_excels <- function(paths) {
  out <- vector("list", length(paths))
  names(out) <- tools::file_path_sans_ext(basename(paths))

  log <- data.frame(
    file = basename(paths),
    status = NA_character_,
    sheet = NA_character_,
    n_cells = NA_integer_,
    n_markers = NA_integer_,
    markers = NA_character_,
    message = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(paths)) {
    res <- tryCatch(clean_one_excel(paths[i]), error = function(e) e)
    if (inherits(res, "error")) {
      log$status[i] <- "error"
      log$message[i] <- conditionMessage(res)
      next
    }

    out[[i]] <- res
    log$status[i] <- "ok"
    log$sheet[i] <- res$sheet
    log$n_cells[i] <- nrow(res$data)
    log$n_markers[i] <- length(res$markers)
    log$markers[i] <- paste(res$markers, collapse = ";")
    log$message[i] <- ""
  }

  list(results = out, log = log)
}

# Example:
files <- list.files("data/4i", pattern = "\\.xlsx$", recursive = TRUE, full.names = TRUE)
cleaned <- clean_many_excels(files)

##########################################################
# Turn cleaned data into Seurat objects 
##########################################################


seurat_list <- list()
for (nm in names(cleaned$results)) {
  res <- cleaned$results[[nm]]
  if (is.null(res) || !is.data.frame(res$data) || nrow(res$data) == 0) next
  
  counts <- t(as.matrix(res$data))  # features x cells
  seu <- CreateSeuratObject(counts = counts, project = nm)
  
  # metadata
  seu$source_file <- nm
  seu$age <- if (!is.null(res$age)) res$age else NA_character_
  
  # add spatial reduction if present
  if (!is.null(res$spatial_coords)) {
    common_cells <- intersect(colnames(seu), rownames(res$spatial_coords))
    if (length(common_cells) > 0) {
      seu <- subset(seu, cells = common_cells)
      coords <- as.matrix(res$spatial_coords[common_cells, , drop = FALSE])
      colnames(coords) <- c("spatial_1", "spatial_2")
      seu[["spatial"]] <- CreateDimReducObject(
        embeddings = coords[colnames(seu), , drop = FALSE],
        key = "SPATIAL_",
        assay = DefaultAssay(seu)
      )
    }
  }
  
  seu <- NormalizeData(seu, verbose = FALSE)
  
  # optional PCA only if enough data
  if (nrow(seu) >= 10 && ncol(seu) >= 30) {
    seu <- FindVariableFeatures(seu, verbose = FALSE)
    seu <- ScaleData(seu, verbose = FALSE)
    npcs_use <- min(10, nrow(seu) - 1, ncol(seu) - 1)
    if (npcs_use >= 2) {
      seu <- RunPCA(seu, npcs = npcs_use, approx = FALSE, verbose = FALSE)
    }
  }
  
  seurat_list[[nm]] <- seu
}



conversion_log <- data.frame(
  sample = names(seurat_list),
  status = NA_character_,
  message = NA_character_,
  stringsAsFactors = FALSE
)
library(Matrix)
library(SeuratObject)

for (i in seq_along(seurat_list)) {
  nm <- names(seurat_list)[i]
  seu <- seurat_list[[i]]
  
  out_base <- gsub("[^A-Za-z0-9_\\-]", "_", nm)
  h5s <- file.path("Data", paste0(out_base, ".h5Seurat"))
  h5a <- file.path("Data", paste0(out_base, ".h5ad"))
  
  DefaultAssay(seu) <- "RNA"
  seu_export <- DietSeurat(
    seu,
    assays = "RNA",
    dimreducs = intersect(c("pca", "spatial"), names(seu@reductions)),
    graphs = NULL
  )
  
  # Seurat v5 layer-safe coercion
  rna_layers <- Layers(seu_export[["RNA"]])
  for (lyr in intersect(c("counts", "data"), rna_layers)) {
    m <- LayerData(seu_export, assay = "RNA", layer = lyr)
    if (!inherits(m, "dgCMatrix")) {
      LayerData(seu_export, assay = "RNA", layer = lyr) <- as(m, "dgCMatrix")
    }
  }
  
  if (file.exists(h5s)) file.remove(h5s)
  if (file.exists(h5a)) file.remove(h5a)
  
  ok <- tryCatch({
    SeuratDisk::SaveH5Seurat(seu_export, filename = h5s, overwrite = TRUE)
    SeuratDisk::Convert(h5s, dest = "h5ad", overwrite = TRUE)
    TRUE
  }, error = function(e) {
    conversion_log$status[i] <- "error"
    conversion_log$message[i] <- conditionMessage(e)
    FALSE
  })
  
  if (ok) {
    conversion_log$status[i] <- "ok"
    conversion_log$message[i] <- ""
  }
}

# Extract clean sample IDs
clean_names <- stringr::str_extract(
  names(seurat_list),
  "(CAMS\\d+a|Subject\\d+|\\d{3,4}UA)"
)

# Assign back
names(seurat_list) <- clean_names

save(seurat_list, file = "data/seurat_list_4i.RData")


save(seurat_list, file = "Data/seurat_list_4i.RData")
