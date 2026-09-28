##########################################################
# Export compact CSV summaries from CODEX object
##########################################################

library(dplyr)
library(readr)

# choose which SnC label to use
snc_col <- "snc2"   # change to "snc" if you want the BCL2 version

md <- obj_codex@meta.data %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    sample = if ("sample" %in% colnames(.)) as.character(sample) else as.character(source_file),
    source_file = if ("source_file" %in% colnames(.)) as.character(source_file) else as.character(sample),
    snc_status = as.character(.data[[snc_col]])
  )

skin_levels <- c(
  "Fibroblast", "Melanocyte", "Basal Keratinocyte",
  "Suprabasal Keratinocyte", "Corneocyte", "Smooth Muscle"
)

##########################################################
# 1) Skin cell type distribution (out of ALL cells)
##########################################################

celltype_distribution <- md %>%
  filter(!is.na(age), !is.na(minorskincell1)) %>%
  group_by(age) %>%
  mutate(total_cells_in_age = n()) %>%   # <-- ALL cells per age
  ungroup() %>%
  filter(minorskincell1 %in% skin_levels) %>%   # numerator = only skin types
  group_by(age, minorskincell1, total_cells_in_age) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    prop_of_total = n_cells / total_cells_in_age,
    pct_of_total = 100 * prop_of_total
  ) %>%
  arrange(age, minorskincell1)

write_csv(
  celltype_distribution,
  "results/260303/codex_skin_celltype_distribution_minorskincell1.csv"
)

##########################################################
# 2) Number of SnC in young vs old
#    one file with pooled age rows + sample rows together
##########################################################

snc_counts_combined <- bind_rows(
  
  md %>%
    filter(!is.na(age), !is.na(snc_status)) %>%
    group_by(age) %>%
    summarise(
      level = "age",
      sample = NA_character_,
      n_total_cells = n(),
      n_snc_cells = sum(snc_status == "SnC+"),
      snc_prop = n_snc_cells / n_total_cells,
      snc_pct = 100 * snc_prop,
      .groups = "drop"
    ),
  
  md %>%
    filter(!is.na(age), !is.na(sample), !is.na(snc_status)) %>%
    group_by(age, sample) %>%
    summarise(
      level = "sample",
      n_total_cells = n(),
      n_snc_cells = sum(snc_status == "SnC+"),
      snc_prop = n_snc_cells / n_total_cells,
      snc_pct = 100 * snc_prop,
      .groups = "drop"
    )
) %>%
  dplyr::select(level, age, sample, n_total_cells, n_snc_cells, snc_prop, snc_pct) %>%
  arrange(factor(level, levels = c("age", "sample")), age, sample)

write_csv(
  snc_counts_combined,
  paste0("results/260303/codex_", snc_col, "_counts_combined.csv")
)

##########################################################
# 3) SnC by cell type young vs old
#    one file with pooled age rows + sample rows together
#    includes both:
#      - within-cell-type SnC burden
#      - composition of SnC pool
##########################################################

# pooled age-level
snc_by_celltype_age <- md %>%
  filter(!is.na(age), !is.na(minorskincell1), !is.na(snc_status)) %>%
  filter(minorskincell1 %in% skin_levels) %>%
  group_by(age, minorskincell1) %>%
  summarise(
    n_total_cells = n(),
    n_snc_cells = sum(snc_status == "SnC+"),
    snc_prop_within_celltype = n_snc_cells / n_total_cells,
    snc_pct_within_celltype = 100 * snc_prop_within_celltype,
    .groups = "drop"
  ) %>%
  group_by(age) %>%
  mutate(
    total_snc_in_age = sum(n_snc_cells),
    snc_composition_prop = ifelse(total_snc_in_age > 0, n_snc_cells / total_snc_in_age, NA_real_),
    snc_composition_pct = 100 * snc_composition_prop
  ) %>%
  ungroup() %>%
  mutate(
    level = "age",
    sample = NA_character_
  )

# sample-level
snc_by_celltype_sample <- md %>%
  filter(!is.na(age), !is.na(sample), !is.na(minorskincell1), !is.na(snc_status)) %>%
  filter(minorskincell1 %in% skin_levels) %>%
  group_by(age, sample, minorskincell1) %>%
  summarise(
    n_total_cells = n(),
    n_snc_cells = sum(snc_status == "SnC+"),
    snc_prop_within_celltype = n_snc_cells / n_total_cells,
    snc_pct_within_celltype = 100 * snc_prop_within_celltype,
    .groups = "drop"
  ) %>%
  group_by(age, sample) %>%
  mutate(
    total_snc_in_sample = sum(n_snc_cells),
    snc_composition_prop = ifelse(total_snc_in_sample > 0, n_snc_cells / total_snc_in_sample, NA_real_),
    snc_composition_pct = 100 * snc_composition_prop
  ) %>%
  ungroup() %>%
  mutate(
    level = "sample",
    total_snc_in_age = NA_real_
  )

snc_by_celltype_combined <- bind_rows(
  snc_by_celltype_age,
  snc_by_celltype_sample
) %>%
  dplyr::select(
    level, age, sample, minorskincell1,
    n_total_cells, n_snc_cells,
    snc_prop_within_celltype, snc_pct_within_celltype,
    starts_with("total_snc"),
    snc_composition_prop, snc_composition_pct
  ) %>%
  arrange(factor(level, levels = c("age", "sample")), age, sample, minorskincell1)

write_csv(
  snc_by_celltype_combined,
  paste0("results/260303/codex_", snc_col, "_by_celltype_combined_minorskincell1.csv")

)







##########################################################
# CD45+ SnC Proportions
##########################################################

library(dplyr)
library(readr)

snc_col <- "snc2"   # keep consistent with your pipeline

# ---------------------------------------
# 0) Fetch expression safely (Seurat v5)
# ---------------------------------------
expr_df <- FetchData(
  obj_codex,
  vars = c("CD45"),
  assay = "RNA",
  layer = "data"
)

# Add CD45 positivity to metadata
obj_codex$CD45_pos <- expr_df$CD45 > 0

# ---------------------------------------
# 1) Build metadata (same style as yours)
# ---------------------------------------
md <- obj_codex@meta.data %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    sample = if ("sample" %in% colnames(.)) as.character(sample) else as.character(source_file),
    source_file = if ("source_file" %in% colnames(.)) as.character(source_file) else as.character(sample),
    snc_status = as.character(.data[[snc_col]])
  )

# ---------------------------------------
# 2) CD45+ SnC counts (same structure)
# ---------------------------------------
cd45_snc_counts_combined <- bind_rows(
  
  # ---- AGE LEVEL ----
  md %>%
    filter(!is.na(age), !is.na(snc_status), !is.na(CD45_pos)) %>%
    group_by(age) %>%
    summarise(
      level = "age",
      sample = NA_character_,
      n_cd45_cells = sum(CD45_pos),
      n_cd45_snc_cells = sum(CD45_pos & snc_status == "SnC+"),
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
      n_cd45_cells = sum(CD45_pos),
      n_cd45_snc_cells = sum(CD45_pos & snc_status == "SnC+"),
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
write_csv(
  cd45_snc_counts_combined,
  paste0("results/260303/codex_cd45_", snc_col, "_counts_combined.csv")
)

print(cd45_snc_counts_combined)
