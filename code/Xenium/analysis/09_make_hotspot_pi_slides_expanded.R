############################################################
# Build expanded PI-facing slides for Xenium SnC hotspot analysis
############################################################

suppressPackageStartupMessages({
  library(officer)
  library(readr)
  library(dplyr)
})

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_file <- normalizePath(sub("^--file=", "", script_file_arg[1]), winslash = "/", mustWork = TRUE)
  project_root <- normalizePath(file.path(dirname(script_file), "../../.."), winslash = "/", mustWork = TRUE)
  setwd(project_root)
}

out_dir <- "results/260727/xenium_snc_hotspot_young_old_comparison"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pptx_out <- file.path(out_dir, "xenium_snc_hotspot_PI_slides_expanded.pptx")
md_out <- file.path(out_dir, "xenium_snc_hotspot_PI_slide_draft_expanded.md")

hotspot_dir <- "results/260727/xenium_snc_hotspots"
neigh_dir <- "results/260727/xenium_snc_hotspot_neighborhoods"
expr_dir <- "results/260727/xenium_snc_hotspot_expression"
scale_dir <- "results/260727/xenium_spatial_scale_diagnostics"

hotspot_summary <- read_csv(file.path(hotspot_dir, "xenium_snc_hotspot_summary.csv"), show_col_types = FALSE)
fraction_intensity <- read_csv(file.path(out_dir, "young_old_hotspot_fraction_intensity_by_sample.csv"), show_col_types = FALSE)
sasp_zone <- read_csv(file.path(out_dir, "young_old_sasp_by_age_zone.csv"), show_col_types = FALSE)
sig_neighbor <- read_csv(file.path(out_dir, "young_old_significant_neighbor_effects_by_age.csv"), show_col_types = FALSE)
neighbor_age <- read_csv(file.path(out_dir, "young_old_neighbor_enrichment_by_age.csv"), show_col_types = FALSE)
comp_age_test <- read_csv(file.path(out_dir, "young_old_celltype_composition_age_tests.csv"), show_col_types = FALSE)
scale_hex <- read_csv(file.path(scale_dir, "xenium_35um_hex_cell_count_summary.csv"), show_col_types = FALSE)
scale_cell <- read_csv(file.path(scale_dir, "xenium_cell_area_spacing_summary.csv"), show_col_types = FALSE)

pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), 100 * x)
num <- function(x, digits = 2) sprintf(paste0("%.", digits, "f"), x)

age_extent <- fraction_intensity %>%
  group_by(age) %>%
  summarise(
    n_samples = n(),
    mean_hot_hex_prop = mean(hot_hex_prop),
    sem_hot_hex_prop = sd(hot_hex_prop) / sqrt(n()),
    mean_snc_prop_in_hot = mean(snc_prop_in_hot_hexes),
    sem_snc_prop_in_hot = sd(snc_prop_in_hot_hexes) / sqrt(n()),
    total_snc_hot = sum(n_snc_in_hot_hexes),
    .groups = "drop"
  )

young_extent <- age_extent %>% filter(age == "Young")
old_extent <- age_extent %>% filter(age == "Old")

sasp_core <- sasp_zone %>%
  filter(zone == "hotspot_core") %>%
  select(age, mean_sasp_num_sample_mean, sem_sasp_num_sample_mean, pooled_snc_prop)
young_sasp <- sasp_core %>% filter(age == "Young")
old_sasp <- sasp_core %>% filter(age == "Old")

hotspot_counts <- hotspot_summary %>%
  filter(metric == "snc_prop") %>%
  group_by(hotspot_class) %>%
  summarise(n_hex = sum(as.numeric(n_hex)), n_cells = sum(as.numeric(n_cells)), n_snc = sum(as.numeric(n_snc)), .groups = "drop")

sample_table <- fraction_intensity %>%
  mutate(
    `Hotspot hexes` = paste0(n_hot_hex, "/", n_hex),
    `Hotspot fraction` = pct(hot_hex_prop),
    `SnC fraction in hotspots` = pct(snc_prop_in_hot_hexes),
    `SnC in hotspots` = n_snc_in_hot_hexes
  ) %>%
  select(Age = age, Sample = sample, `Hotspot hexes`, `Hotspot fraction`, `SnC fraction in hotspots`, `SnC in hotspots`)

sasp_table <- sasp_zone %>%
  filter(zone %in% c("hotspot_core", "neighbor_1", "background")) %>%
  mutate(
    Zone = as.character(zone),
    `Mean SASP genes/cell` = paste0(num(mean_sasp_num_sample_mean), " +/- ", num(sem_sasp_num_sample_mean)),
    `SnC fraction` = pct(pooled_snc_prop)
  ) %>%
  select(Age = age, Zone, `Mean SASP genes/cell`, `SnC fraction`)

immune_age_top <- neighbor_age %>%
  filter(annotation_col == "immune_subclusters1", is.finite(mean_prop_diff)) %>%
  arrange(age, desc(abs(mean_prop_diff))) %>%
  group_by(age) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    `Near %` = pct(mean_near_prop),
    `Reference %` = pct(mean_reference_prop),
    Difference = pct(mean_prop_diff),
    FDR = signif(fdr, 2)
  ) %>%
  select(Age = age, `Immune subtype` = cell_type, `Near %`, `Reference %`, Difference, FDR)

minor_comp_top <- comp_age_test %>%
  filter(annotation_col == "cell_type_minor", zone %in% c("hotspot_core", "neighbor_1")) %>%
  arrange(zone, desc(abs(old_minus_young))) %>%
  group_by(zone) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    `Young %` = pct(mean_prop_young),
    `Old %` = pct(mean_prop_old),
    `Old - Young` = pct(old_minus_young),
    FDR = signif(fdr, 2)
  ) %>%
  select(Zone = zone, `Minor cell type` = cell_type, `Young %`, `Old %`, `Old - Young`, FDR)

sig_text <- if (nrow(sig_neighbor) == 0) {
  "No sample-level age-specific neighbor effects passed FDR < 0.05."
} else {
  paste0(
    sig_neighbor$age[1], " ", sig_neighbor$cell_type[1],
    " was lower near hotspots than reference tissue (",
    pct(sig_neighbor$mean_near_prop[1]), " near vs ",
    pct(sig_neighbor$mean_reference_prop[1]), " reference; FDR ",
    signif(sig_neighbor$fdr[1], 3), ")."
  )
}

fig <- list(
  hotspot_snc = file.path(hotspot_dir, "xenium_snc_hotspots_snc_prop.png"),
  hotspot_sasp = file.path(hotspot_dir, "xenium_snc_hotspots_mean_sasp_num.png"),
  fraction_intensity = file.path(out_dir, "young_old_hotspot_fraction_vs_intensity.png"),
  burden = file.path(out_dir, "young_old_hotspot_burden_by_sample.png"),
  sasp_distance = file.path(out_dir, "young_old_sasp_by_hotspot_distance.png"),
  expr = file.path(expr_dir, "xenium_snc_hotspot_expression_summary.png"),
  neigh_minor = file.path(neigh_dir, "xenium_snc_hotspot_celltype_cell_type_minor.png"),
  neigh_final = file.path(neigh_dir, "xenium_snc_hotspot_celltype_cell_type_final.png"),
  neigh_immune = file.path(neigh_dir, "xenium_snc_hotspot_celltype_immune_subclusters1.png"),
  yo_minor = file.path(out_dir, "young_old_cell_type_minor_composition_hotspot_neighbor.png"),
  yo_final = file.path(out_dir, "young_old_cell_type_final_composition_hotspot_neighbor.png"),
  yo_immune = file.path(out_dir, "young_old_immune_subclusters1_composition_hotspot_neighbor.png")
)

add_blank <- function(ppt, title = NULL) {
  ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
  if (!is.null(title)) {
    ppt <- ph_with(ppt, title, location = ph_location(left = 0.45, top = 0.22, width = 12.4, height = 0.5))
  }
  ppt
}

add_text <- function(ppt, text, left, top, width, height) {
  ph_with(ppt, text, location = ph_location(left = left, top = top, width = width, height = height))
}

add_img <- function(ppt, path, left, top, width, height) {
  if (file.exists(path)) {
    ph_with(ppt, external_img(path, width = width, height = height), location = ph_location(left = left, top = top, width = width, height = height))
  } else {
    add_text(ppt, paste("Missing figure:", path), left, top, width, height)
  }
}

add_tbl <- function(ppt, df, left, top, width, height, font_size = 10) {
  ph_with(ppt, df, location = ph_location(left = left, top = top, width = width, height = height))
}

ppt <- read_pptx()

ppt <- add_slide(ppt, layout = "Title Slide", master = "Office Theme")
ppt <- ph_with(ppt, "Expanded SnC Hotspot Analysis in Xenium Skin", location = ph_location_type(type = "ctrTitle"))
ppt <- ph_with(ppt, "Getis-Ord Gi* hotspots, SASP signal, cell-type neighborhoods, and young-old comparison", location = ph_location_type(type = "subTitle"))

ppt <- add_blank(ppt, "Analysis Roadmap")
ppt <- add_text(ppt, paste(
  "1. Define SnC+ cells from Xenium marker expression.",
  "2. Aggregate cells into 35 um hexagonal spatial bins.",
  "3. Use Getis-Ord Gi* to identify local clusters of high SnC fraction or high SASP signal.",
  "4. Compare hotspot extent and intensity across young vs old samples.",
  "5. Ask which cell types, including immune subtypes, are inside or near hotspots.",
  "6. Keep age comparisons sample-level to avoid pseudo-replication.",
  sep = "\n"
), 0.8, 1.0, 11.8, 5.8)

ppt <- add_blank(ppt, "Definitions and Parameters")
ppt <- add_text(ppt, paste(
  "SnC+ definition:",
  "  (CDKN1A > 0 OR CDKN2A > 0) AND NOT(PCNA > 0 OR MKI67 > 0) AND sasp_num >= 5",
  "",
  "SASP summary:",
  "  sasp_num = number of detected SASP genes per cell; 51 SASP genes present in Xenium panel.",
  "",
  "Hotspot calling:",
  "  35 um hex bins; minimum 5 cells/bin.",
  "  First-order hex neighbors, including the focal hex.",
  "  Hot spot = Gi* z > 0 and BH FDR < 0.05 within each FOV.",
  "",
  paste0("Scale check: retained 35 um hex bins contain median ", scale_hex$median[1], " cells; median equivalent cell diameter is ~", num(scale_cell$median[scale_cell$metric == "equivalent_diameter_um"][1]), " um."),
  sep = "\n"
), 0.8, 0.9, 12.0, 6.0)

ppt <- add_blank(ppt, "What Does Gi* Add?")
ppt <- add_text(ppt, paste(
  "Getis-Ord Gi* asks whether a bin and its immediate neighbors have unusually high values compared with the sample-wide background.",
  "",
  "For this analysis, the value is either:",
  "  SnC fraction per hex bin, or",
  "  mean SASP genes detected per cell per hex bin.",
  "",
  "Important interpretation:",
  "  A hotspot is not simply a bin with one SnC+ cell.",
  "  It is a spatially local cluster whose surrounding region is high relative to that FOV.",
  "",
  "This makes hotspot analysis stricter than simply mapping individual SnC+ cells.",
  sep = "\n"
), 0.8, 0.9, 12.0, 6.0)

ppt <- add_blank(ppt, "SnC Fraction Hotspot Maps")
ppt <- add_img(ppt, fig$hotspot_snc, 0.55, 0.95, 12.2, 6.1)

ppt <- add_blank(ppt, "SASP Hotspot Maps")
ppt <- add_img(ppt, fig$hotspot_sasp, 0.55, 0.95, 12.2, 6.1)

ppt <- add_blank(ppt, "Hotspot Calls by Metric")
ppt <- add_text(ppt, paste(
  "SnC-fraction metric:",
  "  55 hotspot hexes out of 1,673 retained hexes.",
  "  No cold spots passed FDR for SnC fraction.",
  "",
  "Mean-SASP metric:",
  "  47 hotspot hexes and 36 cold-spot hexes.",
  "",
  "Interpretation:",
  "  SnC fraction hotspots highlight clustered SnC-rich regions.",
  "  SASP hotspots highlight regions with elevated SASP burden, which can overlap with but is not identical to SnC fraction.",
  sep = "\n"
), 0.8, 1.0, 5.7, 5.8)
ppt <- add_tbl(ppt, hotspot_counts, 6.8, 1.15, 5.7, 3.0)

ppt <- add_blank(ppt, "Young vs Old: Hotspot Extent and Intensity")
ppt <- add_img(ppt, fig$fraction_intensity, 0.55, 0.95, 7.0, 5.9)
ppt <- add_text(ppt, paste(
  "Extent:",
  paste0("  Young hotspot hex fraction = ", pct(young_extent$mean_hot_hex_prop), "."),
  paste0("  Old hotspot hex fraction = ", pct(old_extent$mean_hot_hex_prop), "."),
  "",
  "Intensity:",
  paste0("  Young SnC fraction inside hotspots = ", pct(young_extent$mean_snc_prop_in_hot), "."),
  paste0("  Old SnC fraction inside hotspots = ", pct(old_extent$mean_snc_prop_in_hot), "."),
  "",
  "Main message:",
  "  Young hotspots are slightly more spatially extensive.",
  "  Old hotspots are more SnC-dense and more intense.",
  sep = "\n"
), 7.85, 1.05, 4.9, 5.7)

ppt <- add_blank(ppt, "Per-Sample Hotspot Burden")
ppt <- add_img(ppt, fig$burden, 0.65, 0.95, 5.2, 5.8)
ppt <- add_tbl(ppt, sample_table, 6.1, 1.0, 6.8, 5.6)

ppt <- add_blank(ppt, "SASP Signal by Distance from Hotspots")
ppt <- add_img(ppt, fig$sasp_distance, 0.65, 0.95, 6.6, 5.8)
ppt <- add_text(ppt, paste(
  "Hotspot core SASP signal:",
  paste0("  Young: ", num(young_sasp$mean_sasp_num_sample_mean), " +/- ", num(young_sasp$sem_sasp_num_sample_mean), " genes/cell."),
  paste0("  Old: ", num(old_sasp$mean_sasp_num_sample_mean), " +/- ", num(old_sasp$sem_sasp_num_sample_mean), " genes/cell."),
  "",
  "SnC fraction in hotspot core:",
  paste0("  Young: ", pct(young_sasp$pooled_snc_prop), "."),
  paste0("  Old: ", pct(old_sasp$pooled_snc_prop), "."),
  "",
  "Interpretation: old hotspot cores carry stronger SnC/SASP burden.",
  sep = "\n"
), 7.55, 1.0, 5.0, 5.8)

ppt <- add_blank(ppt, "SASP and SnC Summary Table")
ppt <- add_tbl(ppt, sasp_table, 0.7, 1.0, 12.0, 5.8)

ppt <- add_blank(ppt, "Raw SnC/SASP Gene Expression by Hotspot Class")
ppt <- add_img(ppt, fig$expr, 0.75, 0.9, 11.8, 6.2)

ppt <- add_blank(ppt, "All Samples: Minor Cell Types Near Hotspots")
ppt <- add_img(ppt, fig$neigh_minor, 0.6, 0.9, 12.1, 6.2)

ppt <- add_blank(ppt, "All Samples: Final Cell Types Near Hotspots")
ppt <- add_img(ppt, fig$neigh_final, 0.6, 0.9, 12.1, 6.2)

ppt <- add_blank(ppt, "All Samples: Immune Subtypes Near Hotspots")
ppt <- add_img(ppt, fig$neigh_immune, 0.6, 0.9, 12.1, 6.2)

ppt <- add_blank(ppt, "Young vs Old: Minor Cell-Type Composition")
ppt <- add_img(ppt, fig$yo_minor, 0.6, 0.9, 12.1, 6.2)

ppt <- add_blank(ppt, "Young vs Old: Final Cell-Type Composition")
ppt <- add_img(ppt, fig$yo_final, 0.6, 0.9, 12.1, 6.2)

ppt <- add_blank(ppt, "Young vs Old: Immune Subtype Composition")
ppt <- add_img(ppt, fig$yo_immune, 0.6, 0.9, 12.1, 6.2)

ppt <- add_blank(ppt, "Immune-Subtype Neighbor Effects")
ppt <- add_text(ppt, paste(
  "These are sample-level effects for immune_subclusters1.",
  "Near = hotspot core + neighbor_1.",
  "Reference = all other cells.",
  "",
  "Because n = 3 per age, these should be interpreted as exploratory.",
  "The table shows the largest near-vs-reference differences, not only FDR-significant hits.",
  sep = "\n"
), 0.7, 1.0, 4.1, 5.9)
ppt <- add_tbl(ppt, immune_age_top, 4.95, 0.95, 7.8, 5.9)

ppt <- add_blank(ppt, "Largest Young-Old Composition Differences")
ppt <- add_text(ppt, paste(
  "Sample-level composition differences for cell_type_minor.",
  "Values are age-group sample means.",
  "FDR values are conservative with n = 3 samples per age.",
  sep = "\n"
), 0.7, 1.0, 4.1, 5.8)
ppt <- add_tbl(ppt, minor_comp_top, 4.95, 0.95, 7.8, 5.9)

ppt <- add_blank(ppt, "Significant Sample-Level Neighbor Effects")
ppt <- add_text(ppt, paste(
  sig_text,
  "",
  "Important distinction:",
  "  Pooled-cell Fisher tests gave many small p-values.",
  "  Sample-level tests are more conservative and better for age comparisons.",
  "  With n = 3 per age group, many large-looking differences remain exploratory.",
  sep = "\n"
), 0.8, 1.0, 11.8, 5.8)

ppt <- add_blank(ppt, "Hotspot Analysis vs Neighborhood Analysis")
ppt <- add_text(ppt, paste(
  "Neighborhood analysis:",
  "  Starts with individual SnC+ cells or a cell cluster of interest.",
  "  Asks which cells are nearby within a radius or nearest-neighbor zone.",
  "  Best for: local context around every SnC+ cell.",
  "",
  "Hotspot analysis:",
  "  First identifies regions where SnC signal is spatially clustered above expectation.",
  "  Then asks what cell types/signals are in the hotspot core and surrounding rings.",
  "  Best for: region-level microenvironments around clustered SnC-rich areas.",
  "",
  "Short version:",
  "  Neighborhood = around individual SnCs.",
  "  Hotspot = around statistically clustered SnC regions.",
  sep = "\n"
), 0.8, 0.9, 11.8, 6.0)

ppt <- add_blank(ppt, "Interpretation")
ppt <- add_text(ppt, paste(
  "1. SnC+ cells form statistically significant spatial hotspots in a subset of tissue regions.",
  "",
  "2. Young samples show slightly higher hotspot spatial extent, but those hotspots contain fewer SnC+ cells.",
  "",
  "3. Old samples show stronger hotspot intensity: higher SnC fraction and higher SASP burden in hotspot cores.",
  "",
  "4. Immune-subtype and minor-cell-type neighbor patterns are visible, but sample-level significance is limited by n = 3 per age.",
  "",
  "5. A useful next step is to inspect old high-intensity hotspots, especially sample 5138-O, against tissue structures and immune states.",
  sep = "\n"
), 0.8, 0.9, 11.8, 6.0)

ppt <- add_blank(ppt, "Suggested PI Discussion Points")
ppt <- add_text(ppt, paste(
  "Do we want to define hotspots by SnC fraction, SASP signal, or both?",
  "",
  "Should old high-intensity hotspots be treated as a distinct phenotype from young low-intensity hotspots?",
  "",
  "Should hotspot-neighbor immune composition be followed up with region-level validation or histology?",
  "",
  "Do we want to repeat with alternative spatial scales, e.g. 25 um and 50 um, as sensitivity analysis?",
  "",
  "Should future statistics model sample as the replicate and hotspot region as nested within sample?",
  sep = "\n"
), 0.8, 0.9, 11.8, 6.0)

print(ppt, target = pptx_out)

md_lines <- c(
  "# Expanded SnC Hotspot Analysis in Xenium Skin",
  "",
  "## Core Story",
  "- SnC+ cells form significant spatial hotspots in a subset of regions.",
  "- Young samples show slightly greater hotspot extent.",
  "- Old samples show greater hotspot intensity and SASP burden.",
  "- Immune-subtype neighbor patterns are included, but sample-level inference is conservative with n = 3 per age.",
  "",
  "## Key Figures Included",
  paste0("- SnC hotspot maps: `", fig$hotspot_snc, "`"),
  paste0("- SASP hotspot maps: `", fig$hotspot_sasp, "`"),
  paste0("- Extent vs intensity: `", fig$fraction_intensity, "`"),
  paste0("- SASP distance plot: `", fig$sasp_distance, "`"),
  paste0("- Expression summary: `", fig$expr, "`"),
  paste0("- All-sample immune subtype neighbors: `", fig$neigh_immune, "`"),
  paste0("- Young-old immune subtype composition: `", fig$yo_immune, "`"),
  "",
  "## Slide-Level Message",
  "Young samples have slightly more hotspot area, while old samples have fewer but stronger SnC/SASP-rich hotspots.",
  "",
  "## Hotspot vs Neighborhood",
  "- Neighborhood analysis asks what is near individual SnC+ cells.",
  "- Hotspot analysis first identifies statistically clustered SnC-rich regions, then asks what cell types/signals occur in and around those regions."
)
writeLines(md_lines, md_out)

message("Wrote ", pptx_out)
message("Wrote ", md_out)
