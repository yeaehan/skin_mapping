##########################################################
# Export compact CSV summaries from CODEX object
##########################################################

library(dplyr)
library(readr)

load("data/4i/260225.RData")

snc_col <- "snc2"   # or "snc" if you want BCL2 version

md <- obj_4i@meta.data %>%
  mutate(
    age = factor(age, levels = c("Young", "Old")),
    sample = as.character(source_file),
    snc_status = as.character(.data[[snc_col]])
  )
# 
# skin_levels <- c(
#   "Fibroblast", "Melanocyte", "Basal Keratinocyte",
#   "Suprabasal Keratinocyte", "Corneocyte", "Smooth Muscle"
# )

##########################################################
# 1) Skin cell type distribution (out of ALL cells)
##########################################################
celltype_distribution <- md %>%
  filter(!is.na(age), !is.na(minorskincell)) %>%
  group_by(age) %>%
  mutate(total_cells_in_age = n()) %>%   # ALL CELLS
  ungroup() %>%
  group_by(age, minorskincell, total_cells_in_age) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    prop_of_total = n_cells / total_cells_in_age,
    pct_of_total = 100 * prop_of_total
  ) %>%
  arrange(age, minorskincell)

write_csv(
  celltype_distribution,
  "results/260303/4i_skin_celltype_distribution.csv"
)


##########################################################
# 2) Number of SnC in young vs old
#    one file with pooled age rows + sample rows together
##########################################################
snc_counts_combined <- bind_rows(
  
  # age-level
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
  
  # sample-level
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
  arrange(level, age, sample)

write_csv(
  snc_counts_combined,
  "results/260303/4i_snc_counts_combined.csv"
)

##########################################################
# 3) SnC by cell type young vs old
#    one file with pooled age rows + sample rows together
#    includes both:
#      - within-cell-type SnC burden
#      - composition of SnC pool
##########################################################
snc_by_celltype <- bind_rows(
  
  # ---- AGE LEVEL ----
  md %>%
    filter(!is.na(age), !is.na(minorskincell), !is.na(snc_status)) %>%
    group_by(age, minorskincell) %>%
    summarise(
      n_total_cells = n(),
      n_snc_cells = sum(snc_status == "SnC+"),
      snc_prop_within_celltype = n_snc_cells / n_total_cells,
      snc_pct_within_celltype = 100 * snc_prop_within_celltype,
      .groups = "drop"
    ) %>%
    group_by(age) %>%
    mutate(
      total_snc = sum(n_snc_cells),
      snc_composition_prop = ifelse(total_snc > 0, n_snc_cells / total_snc, NA),
      snc_composition_pct = 100 * snc_composition_prop
    ) %>%
    ungroup() %>%
    mutate(
      level = "age",
      sample = NA_character_
    ),
  
  # ---- SAMPLE LEVEL ----
  md %>%
    filter(!is.na(age), !is.na(sample), !is.na(minorskincell), !is.na(snc_status)) %>%
    group_by(age, sample, minorskincell) %>%
    summarise(
      n_total_cells = n(),
      n_snc_cells = sum(snc_status == "SnC+"),
      snc_prop_within_celltype = n_snc_cells / n_total_cells,
      snc_pct_within_celltype = 100 * snc_prop_within_celltype,
      .groups = "drop"
    ) %>%
    group_by(age, sample) %>%
    mutate(
      total_snc = sum(n_snc_cells),
      snc_composition_prop = ifelse(total_snc > 0, n_snc_cells / total_snc, NA),
      snc_composition_pct = 100 * snc_composition_prop
    ) %>%
    ungroup() %>%
    mutate(
      level = "sample"
    )
) %>%
  arrange(level, age, sample, minorskincell)

write_csv(
  snc_by_celltype,
  "results/260303/4i_snc_by_celltype_combined.csv"
)