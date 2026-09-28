# Public Validation

Current workflows, in order:

1. `01_snc_quartile_scoring_geo_combined.py` — combined GEO SnC scoring.
2. `02_sinai_age_decade_trends_all_sites.py` — Sinai age trends across sites.
3. `03_sinai_snc_criteria_age_site_celltype.py` — Sinai SnC criteria by site and cell type.
4. `04_sinai_shared_sssp_snc_sender_ligrec_enrichment.py` — SnC+ sender communication.
5. `05_sinai_sssp_core_aucell_age_by_celltype_study.py` — Core21 AUCell age analysis.
6. `06_external_gtex_tabula_sssp_validation.py` — external GTEx/Tabula validation.
7. `07_filter_gtex_single_nucleus_skin.py` — GTEx skin single-nucleus subset.
8. `08_merge_restrepo_tabula_harmony.py` — Restrepo/Tabula reference integration.
9. `09_gtex_bulk_skin_qc.Rmd` — GTEx bulk-skin preprocessing and QC.
10. `10_gtex_bulk_secretory_elastic_net_old_young.py` — old-versus-young classifier.
11. `11_gtex_bulk_sssp_core_age.py` — Core21 GTEx age analysis.
12. `12_gtex_bulk_sssp_core_scoring_comparison.py` — Core21 scoring sensitivity analysis.
13. `13_prepare_gtex_skin_cibersortx_reference.py` — CIBERSORTx reference export.
14. `14_prepare_gtex_skin_cibersortx_bulk.py` — CIBERSORTx bulk-mixture export.
15. `15_run_gtex_skin_music_deconvolution.R` — MuSiC deconvolution.
