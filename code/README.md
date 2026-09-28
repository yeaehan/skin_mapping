# Skin Mapping 2

This repository contains spatial and single-cell analysis workflows for skin mapping studies across multiple platforms, including CODEX, Xenium, and 4i.

## Repository structure

The project is organized by assay/platform:

- `CODEX/` - CODEX-specific preprocessing, analysis, plotting, and reports
- `Xenium/` - Xenium-specific preprocessing, analysis, and plotting workflows
- `4i/` - 4i-specific preprocessing and analysis workflows
- `shared/` - reusable utilities and cross-platform helper code
- `public_validation/` - analyses of publicly available datasets used for validation
- `docs/` - project organization and workflow notes

See `docs/script_inventory.md` for a script-by-script workflow map and cleanup notes.

## Recommended workflow

1. Start with the assay-specific folder for the dataset you are working on.
2. Run preprocessing steps in the corresponding `preprocessing/` directory when present.
3. Perform downstream analysis in `analysis/` or the assay-level workflow files.
4. Write generated figures, tables, and intermediate data outside version control.
5. Keep reusable helper code in `shared/`.

## Notes

- Raw or sensitive data should be kept outside version control whenever possible.
- Local data are expected outside this code directory, for example in `../data`.
- Generated outputs should generally be excluded from Git unless they are small, curated, and needed for documentation.
