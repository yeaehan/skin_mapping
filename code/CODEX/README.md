# CODEX Workflows

This folder contains CODEX preprocessing, analysis, and plotting code.

Recommended order:

1. `preprocessing/01_codex_data_clean.R`
2. `analysis/01_codex_cell_type_analysis.R`
3. `analysis/02_codex_integrated_analysis.R`
4. `analysis/03_codex_banksy_spatial_niche.R`
5. plotting/report export scripts as needed

Audit scripts that check assumptions or intermediate annotations are in `audits/`.
Generated audit outputs are written under `../results/audits/` and are not tracked by git.

See `docs/script_inventory.md` for merge and cleanup notes.
