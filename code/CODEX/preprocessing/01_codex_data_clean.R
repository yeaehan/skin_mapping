# =========================
# Excel cleanup functions (CODEX) -- 4i-style convention
# =========================

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
})

setwd("//mfad/researchmn/PHYS-BIOMED-ENG/PASSOS/SHARED/Yeaeun/ana/skin_mapping_2")



# Canonical markers available/expected in CODEX
marker_panel <- c(
  "BCL2","SMA","p21","Keratin 14","Vimentin","Keratin 5","p16","CD34","IP10","CD8",
  "hIL-1alpha","CD14","LaminB1","Ki67","Filaggrin","CXCR3","HMGB1","Sirt1","IDO1",
  "FOXP3","IL6","MelanA","IL8","CD4","CD68","HLA-DR","CD45","PCNA","ICOS",
  "Pan-Cytokeratin","CD3e","CD11c"   # optional; only present in some files
)

marker_root_map <- c(
  # senescence / SASP / cycling
  "p21" = "p21",
  # Never allow the generic/alternate p16 channel to become canonical p16.
  "p16" = "p16_non_bc42",
  "p16bc42" = "p16",
  "p16layer9" = "p16_layer9",
  "pcna" = "PCNA",
  "ki67" = "Ki67",
  "hmgb1" = "HMGB1",
  "laminb" = "LaminB1",
  "laminb1" = "LaminB1",
  "ip10" = "IP10",
  "cxcr3" = "CXCR3",
  "hil1a" = "hIL-1alpha",
  "hil1alpha" = "hIL-1alpha",
  "il6" = "IL6",
  "il8" = "IL8",
  "sirt1" = "Sirt1",
  "sirtt1" = "Sirt1",
  "ido1" = "IDO1",
  
  # immune
  "cd3e" = "CD3e",
  "cd4" = "CD4",
  "cd8" = "CD8",
  "cd14" = "CD14",
  "cd68" = "CD68",
  "cd11c" = "CD11c",
  "cd34" = "CD34",
  "cd45" = "CD45",
  "foxp3" = "FOXP3",
  "hladr" = "HLA-DR",
  "icos" = "ICOS",
  
  # skin structure
  "sma" = "SMA",
  "vimentin" = "Vimentin",
  "k5" = "Keratin 5",
  "keratin5" = "Keratin 5",
  "k14" = "Keratin 14",
  "keratin14" = "Keratin 14",
  "mela" = "MelanA",
  "melana" = "MelanA",
  "filagrin" = "Filaggrin",
  "fillagrin" = "Filaggrin",
  "filaggrin" = "Filaggrin",
  
  # pancytokeratin variants -> canonical
  "pancyto" = "Pan-Cytokeratin",
  "pancytokeratin" = "Pan-Cytokeratin",

  
  # BCL2 variants
  "bcl2" = "BCL2",
  "bcl2bc42" = "BCL2",
  "bcl2bc" = "BCL2",
  "bcl02" = "BCL2"
)


# NOTE: The marker_root_map above contains many redundant "bcl2" lines because you pasted a big block.
# You can safely delete duplicates; R will keep the LAST occurrence for the same key.
# If you want, I can give you a cleaned minimal map.

get_age_from_path <- function(path) {
  if (grepl("[/\\\\]old[/\\\\]", path, ignore.case = TRUE)) return("Old")
  if (grepl("[/\\\\]young[/\\\\]", path, ignore.case = TRUE)) return("Young")
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
  # prefer nucleiProps-ish, else Sheet1, else first sheet
  i <- grep("^nuclei\\s*props$|nuclei.*props|nucleiprops", s, ignore.case = TRUE)
  if (length(i) > 0) return(s[i[1]])
  if ("Sheet1" %in% s) return("Sheet1")
  s[1]
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
  
  marker_num_count <- rowSums(
    !is.na(as.data.frame(lapply(df[, marker_cols, drop = FALSE],
                                function(z) suppressWarnings(as.numeric(z)))))
  )
  too_few_numeric <- marker_num_count < 2
  
  df[!(is_footer_label | too_few_numeric), , drop = FALSE]
}

# Sample-specific marker exclusions (your notes)
sample_exclude_markers <- function(sample_id) {
  if (sample_id %in% c("202UA", "207UA")) return(c("CD3e"))
  if (sample_id == "7503") return(c("CD11c"))
  if (sample_id == "4614") return(c("Pan-Cytokeratin"))
  character(0)
}

extract_sample_id <- function(path) {
  fname <- basename(path)
  
  # 1) CAMS / Subject / UA (case-insensitive)
  id <- stringr::str_extract(
    fname,
    regex("(CAMS\\d+a|Subject\\d+|\\d{3,4}UA)", ignore_case = TRUE)
  )
  
  if (!is.na(id)) {
    # standardize casing
    if (grepl("^cams", id, ignore.case = TRUE)) {
      id <- sub("^(?i)cams", "CAMS", id, perl = TRUE)
    }
    if (grepl("^subject", id, ignore.case = TRUE)) {
      id <- sub("^(?i)subject", "Subject", id, perl = TRUE)
    }
    # UA usually already numeric+UA; force UA uppercase
    id <- sub("(?i)ua$", "UA", id, perl = TRUE)
    return(id)
  }
  
  # 2) 4-digit sample id even if followed by -digits (e.g., 5014-21 -> 5014)
  m <- stringr::str_match(fname, "(^|[^0-9])([0-9]{4})(-[0-9]+)?")
  if (!all(is.na(m))) return(m[1, 3])
  
  # fallback
  tools::file_path_sans_ext(fname)
}

detect_label_col <- function(df) {
  nms <- names(df)
  roots <- normalize_token(nms)
  
  # strongest matches
  preferred <- c("label", "cellid", "objectid", "id", "barcode")
  hit <- which(roots %in% preferred)
  if (length(hit) > 0) return(nms[hit[1]])
  
  # weaker: contains "label" / "id" / "cell"
  hit2 <- which(grepl("label|cell|object|barcode|id", roots))
  if (length(hit2) > 0) return(nms[hit2[1]])
  
  # fallback: first column
  nms[1]
}

clean_one_excel <- function(path) {
  age <- get_age_from_path(path)
  sample_id <- extract_sample_id(path)
  
  sh <- pick_sheet(path)
  raw <- read_excel(path, sheet = sh)
  raw <- as.data.frame(raw, check.names = FALSE)
  
  # >>> ADDED
  coord_cols <- pick_coord_cols(raw)
  
  df <- rename_columns(raw)
  df <- coalesce_duplicate_columns(df)

  # p16 provenance rule:
  #   CAMS18a -> p16 layer 9
  #   every other sample -> p16-BC42
  if (identical(tolower(sample_id), "cams18a")) {
    if (!"p16_layer9" %in% names(df)) {
      stop(sprintf("CAMS18a is missing required p16 layer 9 in %s", basename(path)))
    }
    df$p16 <- df$p16_layer9
    p16_source <- "p16 layer 9"
  } else {
    if (!"p16" %in% names(df)) {
      stop(sprintf("Sample %s is missing required p16-BC42 in %s", sample_id, basename(path)))
    }
    p16_source <- "p16-BC42"
  }
  
  if (!"Label" %in% names(df)) stop(sprintf("Missing Label in %s", basename(path)))
  
  present_markers <- intersect(marker_panel, names(df))
  if (length(present_markers) == 0) stop(sprintf("No markers found in %s", basename(path)))
  
  excl <- sample_exclude_markers(sample_id)
  present_markers <- setdiff(present_markers, excl)
  
  keep <- c("Label", present_markers)
  df <- df[, keep, drop = FALSE]
  
  df$Label <- as.character(df$Label)
  df <- df[!is.na(df$Label) & df$Label != "", , drop = FALSE]
  df <- remove_footer_rows(df, present_markers)
  df$Label <- make.unique(df$Label)
  
  for (m in present_markers) df[[m]] <- suppressWarnings(as.numeric(df[[m]]))
  
  all_na <- rowSums(is.na(df[, present_markers, drop = FALSE])) == length(present_markers)
  df <- df[!all_na, , drop = FALSE]
  
  for (m in present_markers) {
    v <- suppressWarnings(as.numeric(df[[m]]))
    v[is.na(v) | is.nan(v) | is.infinite(v)] <- 0
    df[[m]] <- v
  }
  
  # >>> ADDED
  # carry spatial coords through just like 4i
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
    spatial_coords = spatial_coords,   # >>> ADDED
    sheet = sh,
    markers = present_markers,
    age = age,
    sample_id = sample_id,
    p16_source = p16_source,
    source_path = path
  )
}

clean_many_excels <- function(paths) {
  out <- vector("list", length(paths))
  names(out) <- basename(paths)
  
  log <- data.frame(
    file = basename(paths),
    sample_id = NA_character_,
    status = NA_character_,
    sheet = NA_character_,
    n_cells = NA_integer_,
    n_markers = NA_integer_,
    markers = NA_character_,
    p16_source = NA_character_,
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
    log$sample_id[i] <- res$sample_id
    log$status[i] <- "ok"
    log$sheet[i] <- res$sheet
    log$n_cells[i] <- nrow(res$data)
    log$n_markers[i] <- length(res$markers)
    log$markers[i] <- paste(res$markers, collapse = ";")
    log$p16_source[i] <- res$p16_source
    log$message[i] <- ""
  }
  
  list(results = out, log = log)
}

# ============
# Run cleanup
# ============
files_codex <- list.files("data/codex", pattern = "\\.xlsx$", recursive = TRUE, full.names = TRUE)
cleaned_codex <- clean_many_excels(files_codex)

# =========================
# Turn cleaned data into Seurat objects
# =========================
seurat_list_codex <- list()

for (i in seq_along(cleaned_codex$results)) {
  res <- cleaned_codex$results[[i]]
  if (is.null(res) || !is.data.frame(res$data) || nrow(res$data) == 0) next
  
  counts <- t(as.matrix(res$data))  # features x cells
  seu <- CreateSeuratObject(counts = counts, project = res$sample_id)
  
  seu$source_file <- basename(res$source_path)
  seu$sample <- res$sample_id
  seu$age <- res$age
  seu$p16_source <- res$p16_source
  
  # >>> ADDED
  # make CODEX spatial metadata/reduction match 4i structure
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
      
      # >>> ADDED
      seu$Spatial_X <- coords[colnames(seu), "spatial_1"]
      seu$Spatial_Y <- coords[colnames(seu), "spatial_2"]
    }
  }
  
  seu <- NormalizeData(seu, verbose = FALSE)
  
  # >>> ADDED
  # mirror 4i optional PCA block
  if (nrow(seu) >= 10 && ncol(seu) >= 30) {
    seu <- FindVariableFeatures(seu, verbose = FALSE)
    seu <- ScaleData(seu, verbose = FALSE)
    npcs_use <- min(10, nrow(seu) - 1, ncol(seu) - 1)
    if (npcs_use >= 2) {
      seu <- RunPCA(seu, npcs = npcs_use, approx = FALSE, verbose = FALSE)
    }
  }
  
  seurat_list_codex[[res$sample_id]] <- seu
}

# Quick check
names(seurat_list_codex)
cleaned_codex$log
write.csv(cleaned_codex$log, "data/codex/codex_preprocessing_log.csv", row.names = FALSE)

# quick check for spatial columns
lapply(seurat_list_codex, function(x) colnames(x@meta.data)[colnames(x@meta.data) %in% c("Spatial_X", "Spatial_Y")])

save(seurat_list_codex, file = "data/codex/seurat_list_codex.RData")



