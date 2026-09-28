##########################################################
# 4i_mouse cleanup -> Seurat objects
##########################################################

install.packages(c("readxl", "dplyr", "stringr"))
if (!require("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(SeuratDisk)
  library(Matrix)
  library(SeuratObject)
})

setwd("//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2")

marker_panel <- c(
  "p21", "p16", "PCNA", "HMGB1", "LaminB1", "BCL2",
  "IL1a", "IL1b", "RANTES",
  "Keratin 14", "Keratin 10", "MelanA", "DCT",
  "Vimentin", "Involucrin",
  "CD3", "F4_80"
)

marker_root_map <- c(
  "p21" = "p21",
  "p16" = "p16",
  "pcna" = "PCNA",
  "hmgb1" = "HMGB1",
  "laminb1" = "LaminB1",
  "laminb" = "LaminB1",
  "bcl2" = "BCL2",
  "il1a" = "IL1a",
  "il1b" = "IL1b",
  "rantes" = "RANTES",
  "k14" = "Keratin 14",
  "keratin14" = "Keratin 14",
  "k10" = "Keratin 10",
  "keratin10" = "Keratin 10",
  "involucrin" = "Involucrin",
  "involuc" = "Involucrin",
  "melana" = "MelanA",
  "dct" = "DCT",
  "vimentin" = "Vimentin",
  "cd3" = "CD3",
  "f480" = "F4_80"
)

extract_sample_id <- function(path) {
  base <- basename(path)
  id <- stringr::str_extract(base, "DJ\\d+")
  if (is.na(id)) id <- tools::file_path_sans_ext(base)
  id
}

get_age_from_path <- function(path) {
  base <- basename(path)
  sample_id <- stringr::str_extract(base, "DJ\\d+")
  if (is.na(sample_id)) return(NA_character_)
  dj_num <- as.numeric(stringr::str_remove(sample_id, "DJ"))
  if (!is.na(dj_num) && dj_num >= 1020 && dj_num <= 1038) return("Old")
  if (!is.na(dj_num) && dj_num < 1020) return("Young")
  NA_character_
}

normalize_token <- function(x) {
  x %>%
    tolower() %>%
    str_replace(regex("_mean_intensity$", ignore_case = TRUE), "") %>%
    str_replace_all("[^a-z0-9]", "")
}

pick_sheet <- function(path) {
  s <- excel_sheets(path)
  s_use <- s[!grepl("keratinocyte", s, ignore.case = TRUE)]
  i <- grep("^nuclei\\s*props$|nuclei.*props|nucleiprops", s_use, ignore.case = TRUE)
  if (length(i) > 0) return(s_use[i[1]])
  if ("Sheet1" %in% s_use) return("Sheet1")
  if ("Sheet" %in% s_use) return("Sheet")
  if (length(s_use) > 0) return(s_use[1])
  stop(sprintf("No usable sheet found in: %s", basename(path)))
}

pick_coord_cols <- function(df_raw) {
  nms <- names(df_raw)
  roots <- normalize_token(nms)
  
  coord_patterns <- list(
    c("centroidxpixels", "centroidypixels"),
    c("cyto1centroidxpixels", "cyto1centroidypixels"),
    c("nuclear1centroidxpixels", "nuclear1centroidypixels"),
    c("xpixels", "ypixels"),
    c("x", "y"),
    c("umap1", "umap2")
  )
  
  for (pat in coord_patterns) {
    if (all(pat %in% roots)) {
      return(c(
        nms[match(pat[1], roots)],
        nms[match(pat[2], roots)]
      ))
    }
  }
  
  NULL
}

rename_columns <- function(df) {
  old <- names(df)
  root <- normalize_token(old)
  new <- old
  
  new[root %in% c("label", "labelid", "mergelabelid")] <- "Label"
  
  new[root %in% c("x", "xpixels", "centroidxpixels", "nuclear1centroidxpixels", "cyto1centroidxpixels")] <- "X"
  new[root %in% c("y", "ypixels", "centroidypixels", "nuclear1centroidypixels", "cyto1centroidypixels")] <- "Y"
  
  new[root == "includedinepidermis"] <- "included_in_epidermis"
  
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
    if (length(idx) == 1) {
      out[[nm]] <- df[[idx]]
    } else {
      out[[nm]] <- Reduce(dplyr::coalesce, df[idx])
    }
  }
  as.data.frame(out, check.names = FALSE)
}

remove_footer_rows <- function(df, marker_cols) {
  is_footer_label <- grepl(
    "threshold|thresh|cutoff|summary|mean|median|sd|note|background",
    df$Label, ignore.case = TRUE
  )
  
  marker_num_count <- rowSums(!is.na(as.data.frame(
    lapply(df[, marker_cols, drop = FALSE], function(z) suppressWarnings(as.numeric(z)))
  )))
  too_few_numeric <- marker_num_count < 2
  
  df[!(is_footer_label | too_few_numeric), , drop = FALSE]
}

clean_one_excel <- function(path) {
  age <- get_age_from_path(path)
  sample_id <- extract_sample_id(path)
  
  sh <- pick_sheet(path)
  raw <- read_excel(path, sheet = sh)
  raw <- as.data.frame(raw, check.names = FALSE)
  
  coord_cols_raw <- pick_coord_cols(raw)
  
  # keep a copy of original labels + coords BEFORE renaming/subsetting
  coord_df <- NULL
  raw_roots <- normalize_token(names(raw))
  label_col_raw <- names(raw)[match(c("label", "labelid", "mergelabelid")[c("label", "labelid", "mergelabelid") %in% raw_roots][1], raw_roots)]
  
  if (!is.null(coord_cols_raw) && !is.na(label_col_raw)) {
    coord_df <- raw[, c(label_col_raw, coord_cols_raw), drop = FALSE]
    names(coord_df) <- c("Label", "spatial_1", "spatial_2")
  }
  
  df <- rename_columns(raw)
  df <- coalesce_duplicate_columns(df)
  
  if (!"Label" %in% names(df)) stop(sprintf("Missing Label in %s", basename(path)))
  
  present_markers <- intersect(marker_panel, names(df))
  if (length(present_markers) == 0) stop(sprintf("No markers found in %s", basename(path)))
  
  metadata_cols <- intersect(c("included_in_epidermis", "X", "Y"), names(df))
  keep <- c("Label", metadata_cols, present_markers)
  df <- df[, keep, drop = FALSE]
  
  df$Label <- as.character(df$Label)
  df <- df[!is.na(df$Label) & df$Label != "", , drop = FALSE]
  df <- remove_footer_rows(df, present_markers)
  df$Label <- make.unique(df$Label)
  
  if ("included_in_epidermis" %in% names(df)) {
    df$included_in_epidermis <- suppressWarnings(as.numeric(df$included_in_epidermis))
  }
  
  for (m in present_markers) df[[m]] <- suppressWarnings(as.numeric(df[[m]]))
  all_na <- rowSums(is.na(df[, present_markers, drop = FALSE])) == length(present_markers)
  df <- df[!all_na, , drop = FALSE]
  
  for (m in present_markers) {
    v <- suppressWarnings(as.numeric(df[[m]]))
    v[is.na(v) | is.nan(v) | is.infinite(v)] <- 0
    df[[m]] <- v
  }
  
  spatial_coords <- NULL
  
  # first preference: raw coord table captured before renaming
  if (!is.null(coord_df)) {
    coord_df$Label <- as.character(coord_df$Label)
    coord_df <- coord_df[!is.na(coord_df$Label) & coord_df$Label != "", , drop = FALSE]
    coord_df$Label <- make.unique(coord_df$Label)
    coord_df <- coord_df[coord_df$Label %in% df$Label, , drop = FALSE]
    rownames(coord_df) <- coord_df$Label
    spatial_coords <- coord_df[df$Label, c("spatial_1", "spatial_2"), drop = FALSE]
    spatial_coords$spatial_1 <- suppressWarnings(as.numeric(spatial_coords$spatial_1))
    spatial_coords$spatial_2 <- suppressWarnings(as.numeric(spatial_coords$spatial_2))
  } else if (all(c("X", "Y") %in% names(df))) {
    # fallback: use renamed columns directly from df
    spatial_coords <- df[, c("X", "Y"), drop = FALSE]
    names(spatial_coords) <- c("spatial_1", "spatial_2")
    rownames(spatial_coords) <- df$Label
    spatial_coords$spatial_1 <- suppressWarnings(as.numeric(spatial_coords$spatial_1))
    spatial_coords$spatial_2 <- suppressWarnings(as.numeric(spatial_coords$spatial_2))
  }
  
  rownames(df) <- df$Label
  meta_df <- df[, metadata_cols, drop = FALSE]
  expr_df <- df[, present_markers, drop = FALSE]
  
  list(
    data = expr_df,
    metadata = meta_df,
    spatial_coords = spatial_coords,
    sheet = sh,
    markers = present_markers,
    age = age,
    sample_id = sample_id,
    source_path = path
  )
}

clean_many_excels <- function(paths) {
  if (length(paths) == 0) {
    stop("No .xlsx files found in data/4i_mouse. Check the folder path and working directory.")
  }
  
  out <- vector("list", length(paths))
  names(out) <- vapply(paths, extract_sample_id, character(1))
  
  log <- data.frame(
    file = basename(paths),
    sample_id = vapply(paths, extract_sample_id, character(1)),
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

files_mouse <- list.files(
  "data/4i_mouse",
  pattern = "\\.xlsx$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files_mouse) == 0) stop("No .xlsx files found in data/4i_mouse")

cleaned_mouse <- clean_many_excels(files_mouse)
print(cleaned_mouse$log)

seurat_list_mouse <- list()

for (nm in names(cleaned_mouse$results)) {
  res <- cleaned_mouse$results[[nm]]
  if (is.null(res) || !is.data.frame(res$data) || nrow(res$data) == 0) next
  
  counts <- t(as.matrix(res$data))
  seu <- CreateSeuratObject(counts = counts, project = res$sample_id)
  
  seu$source_file <- basename(res$source_path)
  seu$sample_id <- res$sample_id
  seu$age <- if (!is.null(res$age)) res$age else NA_character_
  seu$age <- factor(seu$age, levels = c("Young", "Old"))
  
  if (!is.null(res$metadata) && nrow(res$metadata) > 0) {
    common_cells <- intersect(colnames(seu), rownames(res$metadata))
    if (length(common_cells) > 0) {
      seu <- AddMetaData(seu, metadata = res$metadata[common_cells, , drop = FALSE])
    }
  }
  
  if (!is.null(res$spatial_coords)) {
    common_cells <- intersect(colnames(seu), rownames(res$spatial_coords))
    if (length(common_cells) > 0) {
      coords <- as.matrix(res$spatial_coords[common_cells, , drop = FALSE])
      colnames(coords) <- c("spatial_1", "spatial_2")
      seu[["spatial"]] <- CreateDimReducObject(
        embeddings = coords,
        key = "SPATIAL_",
        assay = DefaultAssay(seu)
      )
      seu$centroid_x_pixels <- coords[colnames(seu), 1]
      seu$centroid_y_pixels <- coords[colnames(seu), 2]
    }
  }
  
  seu <- NormalizeData(seu, verbose = FALSE)
  
  if (nrow(seu) >= 10 && ncol(seu) >= 30) {
    seu <- FindVariableFeatures(seu, verbose = FALSE)
    seu <- ScaleData(seu, verbose = FALSE)
    npcs_use <- min(10, nrow(seu) - 1, ncol(seu) - 1)
    if (npcs_use >= 2) {
      seu <- RunPCA(seu, npcs = npcs_use, approx = FALSE, verbose = FALSE)
    }
  }
  
  seurat_list_mouse[[res$sample_id]] <- seu
}

if (anyDuplicated(names(seurat_list_mouse)) > 0) {
  names(seurat_list_mouse) <- make.unique(names(seurat_list_mouse))
}

conversion_log_mouse <- data.frame(
  sample = names(seurat_list_mouse),
  status = NA_character_,
  message = NA_character_,
  stringsAsFactors = FALSE
)

dir.create("Data/4i_mouse", recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(seurat_list_mouse)) {
  nm <- names(seurat_list_mouse)[i]
  seu <- seurat_list_mouse[[i]]
  
  out_base <- gsub("[^A-Za-z0-9_\\-]", "_", nm)
  h5s <- file.path("Data/4i_mouse", paste0(out_base, ".h5Seurat"))
  h5a <- file.path("Data/4i_mouse", paste0(out_base, ".h5ad"))
  
  DefaultAssay(seu) <- "RNA"
  seu_export <- DietSeurat(
    seu,
    assays = "RNA",
    dimreducs = intersect(c("pca", "spatial"), names(seu@reductions)),
    graphs = NULL
  )
  
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
    conversion_log_mouse$status[i] <- "error"
    conversion_log_mouse$message[i] <- conditionMessage(e)
    FALSE
  })
  
  if (ok) {
    conversion_log_mouse$status[i] <- "ok"
    conversion_log_mouse$message[i] <- ""
  }
}

save(cleaned_mouse, file = "data/4i_mouse/cleaned_4i_mouse.RData")
save(seurat_list_mouse, file = "data/4i_mouse/seurat_list_4i_mouse.RData")
write.csv(cleaned_mouse$log, "data/4i_mouse/cleaned_4i_mouse_log.csv", row.names = FALSE)
write.csv(conversion_log_mouse, "data/4i_mouse/conversion_log_4i_mouse.csv", row.names = FALSE)