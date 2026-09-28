# Project organization

This project is organized around a simple analysis lifecycle:

1. Data preprocessing
2. Cell-type annotation
3. Spatial and trajectory analysis
4. Visualization and reporting

## Current folders

- `CODEX/` contains CODEX-specific preprocessing, analysis, plotting, and reports.
- `Xenium/` contains Xenium-specific preprocessing, analysis, and plotting workflows.
- `4i/` contains 4i-specific preprocessing and analysis workflows.
- `shared/` contains reusable helper functions.
- `public_validation/` contains analyses of publicly available datasets used for validation.
- `docs/` contains project organization notes.

## Recommended practice

- Keep source workflows in the assay-specific folders.
- Keep generated outputs and large intermediate files outside Git.
- Keep local input data outside this code directory, for example in `../data`.
- Store reusable documentation in `docs/`.
- Avoid mixing analysis code with generated figures or large data artifacts.
