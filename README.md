# Skin Mapping

Analysis workflows for spatial, single-cell, and bulk studies of human skin aging and senescence.

## Workflows

- `code/Xenium/` — Xenium preprocessing, spatial analyses, hotspot detection, and cell-cell communication.
- `code/CODEX/` — CODEX preprocessing, cell-type annotation, FlowSOM, and spatial analyses.
- `code/4i/` — human and mouse 4i preprocessing and analysis.
- `code/public_validation/` — GTEx, Sinai, Tabula Sapiens, and other external validation workflows.
- `code/shared/` — shared gene sets and analysis helpers.

Raw data, intermediate objects, results, external resources, and presentations are intentionally excluded from version control. See the workflow-specific README files and script headers for expected inputs.

## Python environment

Install the Python dependencies with:

```bash
pip install -r requirements.txt
```

Several workflows also require R/Bioconductor packages documented in their script headers.
