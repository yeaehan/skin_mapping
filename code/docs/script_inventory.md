# Script Inventory

This inventory describes the intended role of each analysis folder and flags scripts that are good candidates for future consolidation.

## Shared Helpers

- `shared/plotfunctions.R` defines shared plotting helpers, including `geom_split_violin()`, `get_marker()`, and `sample_mean_violin()`.
- Prefer adding reusable plotting or marker-access helpers here instead of copying the same function into multiple analysis scripts.

## 4i

The `01A` and `01B` prefixes indicate parallel human and mouse branches. These should stay separate because marker panels and gene names differ by species.

- `4i/preprocessing/01A_4i_human_data_clean.R` cleans human 4i Excel files and creates Seurat objects.
- `4i/preprocessing/01B_4i_mouse_data_clean.R` cleans mouse 4i Excel files and creates Seurat objects.
- `4i/analysis/01A_4i_human_cell_type_analysis.R` annotates and plots human 4i cell types and senescence states.
- `4i/analysis/01B_4i_mouse_cell_type_analysis.R` annotates and plots mouse 4i cell types and senescence states.
- `4i/analysis/02_4i_banksy_spatial_niche.R` runs Banksy spatial niche analysis for 4i data.
- `4i/analysis/03_4i_raw_values_to_doc.R` exports raw values used for documentation or reporting.

Future cleanup: extract shared Excel-cleaning helpers from the two preprocessing scripts into `shared/` or `4i/shared/`, while keeping species-specific marker panels in separate files.

## CODEX

- `CODEX/preprocessing/01_codex_data_clean.R` cleans CODEX source spreadsheets and creates analysis-ready objects.
- `CODEX/analysis/01_codex_cell_type_analysis.R` annotates CODEX cell types and senescence states.
- `CODEX/analysis/02_codex_integrated_analysis.R` creates integrated/merged CODEX objects and UMAP views.
- `CODEX/analysis/03_codex_banksy_spatial_niche.R` runs Banksy spatial niche analysis for CODEX data.
- `CODEX/analysis/04_codex_raw_values_to_doc.R` exports raw values used for documentation or reporting.
- `CODEX/analysis/05_seurat_codex_cell_type_annotation.Rmd` contains Seurat-based CODEX annotation and plotting notes.
- `CODEX/analysis/06_seurat_codex_integration.Rmd` contains Seurat-based CODEX integration notes.
- `CODEX/analysis/08_seurat_codex_cell_type_annotation_from_preprocessing.Rmd` appears related to `05`; review whether both are still needed after the preprocessing workflow is finalized.
- `CODEX/analysis/check_immune_marker_overlaps.R` audits CODEX immune marker overlaps without modifying metadata.
- `CODEX/plotting/01_codex_violinplots.R` builds CODEX violin plots with sample-level dots.
- `CODEX/plotting/02_seurat_codex_plotting.Rmd` contains Seurat CODEX plotting code.

Removed duplicates: `07_seurat_codex_integration_from_preprocessing.Rmd` was identical to `06_seurat_codex_integration.Rmd`, and `03_seurat_codex_plotting.Rmd` duplicated `02_seurat_codex_plotting.Rmd` except for one diagnostic `DimPlot()` chunk.

Future cleanup: merge repeated plotting helpers from `01_codex_cell_type_analysis.R`, `01_codex_violinplots.R`, and 4i plotting scripts into `shared/plotfunctions.R`.

## Xenium

The Xenium folder is still mid-reorganization. Root-level files and nested `preprocessing/` / `analysis/` files overlap in purpose.

- `Xenium/preprocessing/01_seurat_xenium_preprocessing.Rmd` is the preferred location for Xenium preprocessing.
- `Xenium/analysis/01_seurat_xenium_cell_type_annotation.Rmd` is the preferred location for Xenium cell-type annotation.
- `Xenium/analysis/02_seurat_xenium_trajectory.Rmd` is the preferred location for Xenium trajectory analysis.
- Root-level `1_...` and `2_...` Xenium files should be reviewed and either moved into the nested folders or archived once the canonical workflow is confirmed.
- `Xenium/2_Xenium_PrettyPlotting.R`, `2_Xenium_SnCAnalysis.Rmd`, and `2_Xenium_SnCAnalysis_Scoring.Rmd` are substantial analysis/plotting scripts and should not be deleted without validating their outputs.

Future cleanup: decide whether nested folders or root-level numbered files are canonical, then retire the duplicate track.

## Repository-Wide Cleanup Targets

- Move package installation commands out of analysis scripts and into `environment.yml` plus setup notes.
- Replace hard-coded absolute paths with project-relative paths or a small config file.
- Avoid `save.image()` in scripts intended for GitHub; save explicit objects with `saveRDS()` instead.
- Keep generated results and large intermediate data outside Git.
