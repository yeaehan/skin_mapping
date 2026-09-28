############################################################
# Build PI-facing slides for Xenium SnC hotspot analysis
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

pptx_out <- file.path(out_dir, "xenium_snc_hotspot_PI_slides.pptx")
md_out <- file.path(out_dir, "xenium_snc_hotspot_PI_slide_draft.md")

hotspot_summary_csv <- "results/260727/xenium_snc_hotspots/xenium_snc_hotspot_summary.csv"
fraction_intensity_csv <- file.path(out_dir, "young_old_hotspot_fraction_intensity_by_sample.csv")
sasp_zone_csv <- file.path(out_dir, "young_old_sasp_by_age_zone.csv")
sig_neighbor_csv <- file.path(out_dir, "young_old_significant_neighbor_effects_by_age.csv")

spatial_plot <- "results/260727/xenium_snc_hotspots/xenium_snc_hotspots_snc_prop.png"
fraction_intensity_plot <- file.path(out_dir, "young_old_hotspot_fraction_vs_intensity.png")
sasp_plot <- file.path(out_dir, "young_old_sasp_by_hotspot_distance.png")
celltype_plot <- file.path(out_dir, "young_old_cell_type_minor_composition_hotspot_neighbor.png")

pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), 100 * x)
num <- function(x, digits = 2) sprintf(paste0("%.", digits, "f"), x)

hotspot_summary <- read_csv(hotspot_summary_csv, show_col_types = FALSE)
fraction_intensity <- read_csv(fraction_intensity_csv, show_col_types = FALSE)
sasp_zone <- read_csv(sasp_zone_csv, show_col_types = FALSE)
sig_neighbor <- read_csv(sig_neighbor_csv, show_col_types = FALSE)

snc_hotspots <- hotspot_summary %>%
  filter(metric == "snc_prop", hotspot_class == "Hot spot")

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

sasp_core <- sasp_zone %>%
  filter(zone == "hotspot_core") %>%
  select(age, mean_sasp_num_sample_mean, sem_sasp_num_sample_mean, pooled_snc_prop)

young_extent <- age_extent %>% filter(age == "Young")
old_extent <- age_extent %>% filter(age == "Old")
young_sasp <- sasp_core %>% filter(age == "Young")
old_sasp <- sasp_core %>% filter(age == "Old")

slide_size <- c(width = 13.333, height = 7.5)

add_blank <- function(ppt, title = NULL) {
  ppt <- add_slide(ppt, layout = "Blank", master = "Office Theme")
  if (!is.null(title)) {
    ppt <- ph_with(
      ppt,
      value = title,
      location = ph_location(left = 0.55, top = 0.25, width = 12.2, height = 0.55)
    )
  }
  ppt
}

add_text <- function(ppt, text, left, top, width, height, font_size = 20) {
  ph_with(
    ppt,
    value = text,
    location = ph_location(left = left, top = top, width = width, height = height)
  )
}

add_img <- function(ppt, path, left, top, width, height) {
  if (file.exists(path)) {
    ph_with(
      ppt,
      value = external_img(path, width = width, height = height),
      location = ph_location(left = left, top = top, width = width, height = height)
    )
  } else {
    add_text(ppt, paste("Missing figure:", path), left, top, width, height)
  }
}

ppt <- read_pptx()

ppt <- add_slide(ppt, layout = "Title Slide", master = "Office Theme")
ppt <- ph_with(ppt, "SnC Hotspot Analysis in Xenium Skin", location = ph_location_type(type = "ctrTitle"))
ppt <- ph_with(
  ppt,
  paste(
    "Getis-Ord Gi* spatial clustering of senescent-like cells",
    "Young vs old skin samples",
    sep = "\n"
  ),
  location = ph_location_type(type = "subTitle")
)

ppt <- add_blank(ppt, "Study Question")
ppt <- add_text(
  ppt,
  paste(
    "Question: are SnC+ cells randomly scattered, or do they form spatially clustered hotspots?",
    "",
    "Dataset: Xenium skin data, 6 samples: 3 young and 3 old.",
    "Cells analyzed: 32,549 total cells.",
    "SnC+ cells: 505 total cells using the SnC gate.",
    "SASP genes detected in panel: 51.",
    "",
    "Core readout: hotspot extent vs hotspot intensity.",
    sep = "\n"
  ),
  0.8, 1.1, 11.7, 5.5
)

ppt <- add_blank(ppt, "SnC and Hotspot Definitions")
ppt <- add_text(
  ppt,
  paste(
    "SnC+ definition:",
    "  (CDKN1A > 0 OR CDKN2A > 0) AND NOT(PCNA > 0 OR MKI67 > 0) AND sasp_num >= 5",
    "",
    "Spatial binning:",
    "  Cells were assigned to 35 um hexagonal bins within each FOV.",
    "  Each bin was summarized by SnC fraction and mean SASP gene count.",
    "",
    "Hotspot definition:",
    "  Getis-Ord Gi* was calculated within each FOV.",
    "  Hot spot = positive Gi* z-score with BH FDR < 0.05.",
    "  For SnC fraction, 55 of 1,673 retained hex bins were called hotspots.",
    sep = "\n"
  ),
  0.8, 1.0, 11.9, 5.8
)

ppt <- add_blank(ppt, "Spatial Hotspot Maps")
ppt <- add_img(ppt, spatial_plot, 0.7, 1.0, 12.0, 6.0)

ppt <- add_blank(ppt, "Hotspot Extent vs Intensity")
ppt <- add_img(ppt, fraction_intensity_plot, 0.7, 1.0, 7.0, 5.8)
ppt <- add_text(
  ppt,
  paste(
    "Extent:",
    paste0("  Young mean hotspot hex fraction = ", pct(young_extent$mean_hot_hex_prop), "."),
    paste0("  Old mean hotspot hex fraction = ", pct(old_extent$mean_hot_hex_prop), "."),
    "",
    "Intensity:",
    paste0("  Young mean SnC fraction inside hotspots = ", pct(young_extent$mean_snc_prop_in_hot), "."),
    paste0("  Old mean SnC fraction inside hotspots = ", pct(old_extent$mean_snc_prop_in_hot), "."),
    "",
    "Interpretation:",
    "  Young hotspots are slightly more spatially extensive.",
    "  Old hotspots are more SnC-dense, especially sample 5138-O.",
    sep = "\n"
  ),
  8.0, 1.2, 4.8, 5.3
)

ppt <- add_blank(ppt, "SASP Signal Near Hotspots")
ppt <- add_img(ppt, sasp_plot, 0.7, 1.0, 7.4, 5.8)
ppt <- add_text(
  ppt,
  paste(
    "Mean detected SASP genes per cell in hotspot core:",
    paste0("  Young: ", num(young_sasp$mean_sasp_num_sample_mean), " +/- ", num(young_sasp$sem_sasp_num_sample_mean)),
    paste0("  Old: ", num(old_sasp$mean_sasp_num_sample_mean), " +/- ", num(old_sasp$sem_sasp_num_sample_mean)),
    "",
    "SnC fraction in hotspot core:",
    paste0("  Young: ", pct(young_sasp$pooled_snc_prop)),
    paste0("  Old: ", pct(old_sasp$pooled_snc_prop)),
    "",
    "Interpretation: old hotspot cores show stronger SnC/SASP burden.",
    sep = "\n"
  ),
  8.3, 1.2, 4.6, 5.3
)

ppt <- add_blank(ppt, "Cell Types Near Hotspots")
ppt <- add_img(ppt, celltype_plot, 0.7, 1.0, 7.4, 5.8)
neighbor_text <- if (nrow(sig_neighbor) == 0) {
  "No sample-level neighbor effects passed FDR < 0.05."
} else {
  paste0(
    "Significant sample-level effect: ",
    sig_neighbor$age[1], " ",
    sig_neighbor$cell_type[1], " in ",
    sig_neighbor$annotation_col[1], " was lower near hotspots than reference tissue.",
    "\nNear = ", pct(sig_neighbor$mean_near_prop[1]),
    "; reference = ", pct(sig_neighbor$mean_reference_prop[1]),
    "; FDR = ", signif(sig_neighbor$fdr[1], 3), "."
  )
}
ppt <- add_text(
  ppt,
  paste(
    "Cell-type composition was summarized at the sample level first.",
    "This prevents samples with more cells from dominating the result.",
    "",
    neighbor_text,
    "",
    "Treat these as exploratory with n = 3 samples per age group.",
    sep = "\n"
  ),
  8.3, 1.2, 4.6, 5.3
)

ppt <- add_blank(ppt, "Hotspot Analysis vs Neighborhood Analysis")
ppt <- add_text(
  ppt,
  paste(
    "Neighborhood analysis:",
    "  Anchor on individual cells or clusters of interest.",
    "  Ask: what cells are within a fixed radius or nearest-neighbor zone?",
    "  Good for local cell-cell context around all SnC+ cells.",
    "",
    "Hotspot analysis:",
    "  First identifies spatial regions where SnC signal is locally clustered above expectation.",
    "  Ask: what cell types and signals occur inside or around significant SnC-rich regions?",
    "  Higher bar: requires local spatial aggregation, not just proximity to one SnC+ cell.",
    "",
    "Short version:",
    "  Neighborhood = around individual SnCs.",
    "  Hotspot = around statistically clustered SnC regions.",
    sep = "\n"
  ),
  0.75, 1.0, 12.0, 5.9
)

ppt <- add_blank(ppt, "Main Takeaways")
ppt <- add_text(
  ppt,
  paste(
    "1. SnC+ cells are not only individual events; a subset forms spatially significant hotspots.",
    "",
    "2. Young samples show slightly higher hotspot spatial extent, but lower SnC intensity inside hotspots.",
    "",
    "3. Old samples show fewer hotspot regions on average, but those regions are much more SnC-dense and SASP-rich.",
    "",
    "4. Sample-level analysis is more conservative than pooled-cell analysis and better suited for PI-facing age comparisons.",
    "",
    "5. Next step: validate whether old high-intensity hotspots co-localize with specific tissue structures or immune states.",
    sep = "\n"
  ),
  0.8, 1.0, 11.8, 5.8
)

print(ppt, target = pptx_out)

md_lines <- c(
  "# SnC Hotspot Analysis in Xenium Skin",
  "",
  "## Slide 1 - Title",
  "SnC Hotspot Analysis in Xenium Skin",
  "Subtitle: Getis-Ord Gi* spatial clustering of senescent-like cells in young vs old skin.",
  "",
  "Speaker note: This analysis asks whether SnC-like cells are spatially clustered into regions, not just whether individual SnC cells have neighbors.",
  "",
  "## Slide 2 - Study Question",
  "- Are SnC+ cells randomly scattered, or do they form spatially clustered hotspots?",
  "- Dataset: 6 Xenium skin samples, 3 young and 3 old.",
  "- Cells analyzed: 32,549.",
  "- SnC+ cells: 505.",
  "- SASP genes detected in panel: 51.",
  "- Main readout: hotspot extent vs hotspot intensity.",
  "",
  "## Slide 3 - Definitions",
  "- SnC+ = (CDKN1A > 0 OR CDKN2A > 0) AND NOT(PCNA > 0 OR MKI67 > 0) AND sasp_num >= 5.",
  "- Cells were assigned to 35 um hexagonal bins within each FOV.",
  "- Each bin was summarized by SnC fraction and SASP signal.",
  "- Hot spot = positive Getis-Ord Gi* z-score with BH FDR < 0.05 within each FOV.",
  "- For SnC fraction, 55 of 1,673 retained hex bins were called hotspots.",
  "",
  "## Slide 4 - Spatial Hotspot Maps",
  paste0("Figure: `", spatial_plot, "`"),
  "- Red bins indicate significant local clustering of high SnC fraction.",
  "- Gray/non-red bins are not significant.",
  "- This converts single-cell SnC calls into spatial regions for downstream comparison.",
  "",
  "## Slide 5 - Hotspot Extent vs Intensity",
  paste0("Figure: `", fraction_intensity_plot, "`"),
  paste0("- Young mean hotspot hex fraction: ", pct(young_extent$mean_hot_hex_prop), "."),
  paste0("- Old mean hotspot hex fraction: ", pct(old_extent$mean_hot_hex_prop), "."),
  paste0("- Young mean SnC fraction inside hotspot hexes: ", pct(young_extent$mean_snc_prop_in_hot), "."),
  paste0("- Old mean SnC fraction inside hotspot hexes: ", pct(old_extent$mean_snc_prop_in_hot), "."),
  "- Interpretation: young has slightly greater hotspot extent, but old has much greater hotspot intensity.",
  "",
  "## Slide 6 - SASP Signal Near Hotspots",
  paste0("Figure: `", sasp_plot, "`"),
  paste0("- Mean SASP genes detected per cell in hotspot core, young: ", num(young_sasp$mean_sasp_num_sample_mean), " +/- ", num(young_sasp$sem_sasp_num_sample_mean), "."),
  paste0("- Mean SASP genes detected per cell in hotspot core, old: ", num(old_sasp$mean_sasp_num_sample_mean), " +/- ", num(old_sasp$sem_sasp_num_sample_mean), "."),
  paste0("- SnC fraction in hotspot core, young: ", pct(young_sasp$pooled_snc_prop), "."),
  paste0("- SnC fraction in hotspot core, old: ", pct(old_sasp$pooled_snc_prop), "."),
  "- Interpretation: old hotspot cores are more SnC/SASP-rich.",
  "",
  "## Slide 7 - Cell Types Near Hotspots",
  paste0("Figure: `", celltype_plot, "`"),
  "- Cell-type composition was summarized per sample first, then averaged across samples.",
  "- This avoids giving more weight to samples with more cells or larger hotspot regions.",
  paste0("- Significant neighbor effect: ", neighbor_text),
  "- Because n = 3 per age group, treat cell-type findings as exploratory.",
  "",
  "## Slide 8 - Hotspot Analysis vs Neighborhood Analysis",
  "- Neighborhood analysis anchors on individual cells or clusters of interest and asks what is nearby.",
  "- Hotspot analysis first identifies statistically significant spatial regions with high local SnC signal.",
  "- Neighborhood analysis is good for local cell-cell context around all SnCs.",
  "- Hotspot analysis is stricter: it asks what happens in or near clustered SnC-rich regions.",
  "- Short version: neighborhood = around individual SnCs; hotspot = around statistically clustered SnC regions.",
  "",
  "## Slide 9 - Main Takeaways",
  "- SnC+ cells form spatially significant hotspots in a subset of regions.",
  "- Young samples show slightly greater hotspot spatial extent.",
  "- Old samples show greater hotspot intensity: more SnC-rich and SASP-rich hotspot cores.",
  "- Sample-level analysis is more conservative and better for age-group comparison.",
  "- Next step: inspect whether old high-intensity hotspots align with tissue structures or immune states."
)

writeLines(md_lines, md_out)

message("Wrote ", pptx_out)
message("Wrote ", md_out)
