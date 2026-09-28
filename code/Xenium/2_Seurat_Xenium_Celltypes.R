library(dplyr)
library(readr)

# ---------------------------------------
# Xenium CSV export: exactly 3 CSV files
# using cell_type_minor restricted to skin cells
# ---------------------------------------

# choose which SnC definition to use
snc_col <- "snc2"   # change to "snc" if needed

# define skin-relevant cell types from cell_type_minor
skin_cell_types <- c(
  "Fibroblast",
  "Keratinocyte/Sebocyte",
  "Differentiated Keratinocyte",
  "Basal/Suprabasal Keratinocyte",
  "Basal Stem Cell",
  "Melanocyte",
  "Smooth Muscle Cell"
)

# metadata table
md <- hsskin.ST@meta.data %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    sample = as.character(sample),
    cell_type_minor = as.character(cell_type_minor)
  )

# ======================================================
# CSV 1) skin cell type distribution young vs old
# ======================================================
xen_skin_celltype_distribution <- md %>%
  filter(!is.na(age), !is.na(cell_type_minor)) %>%
  filter(cell_type_minor %in% skin_cell_types) %>%
  count(age, cell_type_minor, name = "n_cells") %>%
  group_by(age) %>%
  mutate(
    total_cells_in_age = sum(n_cells),
    prop_within_age = n_cells / total_cells_in_age,
    pct_within_age = 100 * prop_within_age
  ) %>%
  ungroup() %>%
  arrange(age, cell_type_minor)

write_csv(
  xen_skin_celltype_distribution,
  "skin_mapping/results/xenium_skin_celltype_distribution_young_vs_old.csv"
)

# ======================================================
# CSV 2) number of SnC in young vs old
#         + individual number of SnCs by sample
# ======================================================
xen_snc_numbers <- bind_rows(
  md %>%
    filter(!is.na(age), !is.na(.data[[snc_col]])) %>%
    group_by(age) %>%
    summarise(
      level = "age",
      sample = NA_character_,
      n_total_cells = n(),
      n_snc_cells = sum(.data[[snc_col]] == "SnC+", na.rm = TRUE),
      snc_prop = n_snc_cells / n_total_cells,
      snc_pct = 100 * snc_prop,
      .groups = "drop"
    ),
  md %>%
    filter(!is.na(age), !is.na(sample), !is.na(.data[[snc_col]])) %>%
    group_by(age, sample) %>%
    summarise(
      level = "sample",
      n_total_cells = n(),
      n_snc_cells = sum(.data[[snc_col]] == "SnC+", na.rm = TRUE),
      snc_prop = n_snc_cells / n_total_cells,
      snc_pct = 100 * snc_prop,
      .groups = "drop"
    )
) %>%
  select(level, age, sample, n_total_cells, n_snc_cells, snc_prop, snc_pct) %>%
  arrange(level, age, sample)

write_csv(
  xen_snc_numbers,
  paste0("skin_mapping/results/xenium_", snc_col, "_numbers_young_vs_old_and_by_sample.csv")
)

# ======================================================
# CSV 3) SnC by cell type young vs old
# ======================================================
xen_snc_by_celltype <- md %>%
  filter(!is.na(age), !is.na(cell_type_minor), !is.na(.data[[snc_col]])) %>%
  filter(cell_type_minor %in% skin_cell_types) %>%
  group_by(age, cell_type_minor) %>%
  summarise(
    n_total_cells = n(),
    n_snc_cells = sum(.data[[snc_col]] == "SnC+", na.rm = TRUE),
    snc_prop_within_celltype = n_snc_cells / n_total_cells,
    snc_pct_within_celltype = 100 * snc_prop_within_celltype,
    .groups = "drop"
  ) %>%
  arrange(age, cell_type_minor)

write_csv(
  xen_snc_by_celltype,
  paste0("skin_mapping/results/xenium_", snc_col, "_by_celltype_young_vs_old.csv")
)





# ======================================================
# CD45+ SnC Proportions
# ======================================================

# ---------------------------------------
# 0) Fetch PTPRC (CD45) expression safely
# ---------------------------------------
# Access raw counts from the Xenium assay
ptprc_data <- GetAssayData(hsskin.ST, assay = "Xenium", layer = "counts")
ptprc_vec <- ptprc_data["PTPRC", ]

# Add CD45 positivity to Seurat object metadata
hsskin.ST$CD45_pos <- ptprc_vec > 0

# ---------------------------------------
# 1) Build metadata
# ---------------------------------------
snc_col <- "snc2" 

md <- hsskin.ST@meta.data %>%
  mutate(
    cell_id = rownames(.),
    age = factor(age, levels = c("Young", "Old")),
    sample = as.character(sample),
    snc_status = as.character(.data[[snc_col]])
  )

# ---------------------------------------
# 2) CD45+ SnC counts (Target Structure)
# ---------------------------------------
cd45_snc_counts_combined <- bind_rows(
  
  # ---- AGE LEVEL ----
  md %>%
    filter(!is.na(age), !is.na(snc_status), !is.na(CD45_pos)) %>%
    group_by(age) %>%
    summarise(
      level = "age",
      sample = NA_character_,
      n_cd45_cells = sum(CD45_pos, na.rm = TRUE),
      n_cd45_snc_cells = sum(CD45_pos & snc_status == "SnC+", na.rm = TRUE),
      cd45_snc_prop = n_cd45_snc_cells / n_cd45_cells,
      cd45_snc_pct = 100 * cd45_snc_prop,
      .groups = "drop"
    ),
  
  # ---- SAMPLE LEVEL ----
  md %>%
    filter(!is.na(age), !is.na(sample), !is.na(snc_status), !is.na(CD45_pos)) %>%
    group_by(age, sample) %>%
    summarise(
      level = "sample",
      n_cd45_cells = sum(CD45_pos, na.rm = TRUE),
      n_cd45_snc_cells = sum(CD45_pos & snc_status == "SnC+", na.rm = TRUE),
      cd45_snc_prop = n_cd45_snc_cells / n_cd45_cells,
      cd45_snc_pct = 100 * cd45_snc_prop,
      .groups = "drop"
    )
) %>%
  dplyr::select(level, age, sample, n_cd45_cells, n_cd45_snc_cells, cd45_snc_prop, cd45_snc_pct) %>%
  arrange(factor(level, levels = c("age", "sample")), age, sample)

# ---------------------------------------
# 3) Save CSV
# ---------------------------------------
# Ensure the directory exists before saving
if (!dir.exists("skin_mapping/results")) { dir.create("skin_mapping/results", recursive = TRUE) }

write_csv(
  cd45_snc_counts_combined,
  paste0("skin_mapping/results/xenium_cd45_", snc_col, "_counts_combined.csv")
)