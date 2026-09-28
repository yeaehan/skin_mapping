#!/usr/bin/env python
"""Subset the public GTEx eight-tissue snRNA-seq atlas to skin nuclei."""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import anndata as ad
import pandas as pd


def parse_args() -> argparse.Namespace:
    project_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=project_root)
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def counts_table(obs: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    available = [column for column in columns if column in obs.columns]
    if not available:
        return pd.DataFrame()
    return obs.groupby(available, observed=True, dropna=False).size().rename("n_cells").reset_index()


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    input_path = (args.input or root / "data" / "external" /
                  "GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad").resolve()
    output_path = (args.output or root / "data" / "external" /
                   "GTEx_skin_snRNAseq_atlas_071421.h5ad").resolve()
    results_dir = root / "results" / "public_validation" / "gtex_skin_snRNA_filter"
    results_dir.mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    if not input_path.exists():
        raise FileNotFoundError(input_path)
    atlas = ad.read_h5ad(input_path, backed="r")
    if "tissue" not in atlas.obs.columns:
        raise KeyError("Expected lowercase 'tissue' metadata column")

    tissue = atlas.obs["tissue"].astype(str).str.strip().str.lower()
    counts_table(atlas.obs.assign(tissue_normalized=tissue), ["tissue_normalized"]).to_csv(
        results_dir / "all_tissues_cell_counts.csv", index=False
    )
    skin_mask = tissue.eq("skin").to_numpy()
    if not skin_mask.any():
        raise ValueError("No cells matched tissue == 'skin'")

    # Materialize only the skin slice; the source atlas remains untouched.
    skin = atlas[skin_mask].to_memory().copy()
    atlas.file.close()
    skin.obs["gtex_skin_filter"] = "tissue == skin"

    # Whole-atlas graphs are no longer valid after tissue subsetting. Existing
    # PCA/UMAP coordinates are retained as reference coordinates, but any skin-
    # specific graph/UMAP should be recomputed downstream.
    for key in list(skin.obsp.keys()):
        del skin.obsp[key]
    for key in ["neighbors", "paga", "leiden"]:
        skin.uns.pop(key, None)
    skin.uns["skin_subset_provenance"] = {
        "source_file": input_path.name,
        "filter": "obs['tissue'].astype(str).str.lower() == 'skin'",
        "n_source_cells": int(len(skin_mask)),
        "n_skin_cells": int(skin.n_obs),
        "neighbor_graph_removed": True,
        "existing_embeddings_retained_as_reference": True,
    }

    audit_specs = {
        "skin_cells_by_site.csv": ["Tissue Site Detail"],
        "skin_cells_by_sample.csv": ["Sample ID", "Participant ID"],
        "skin_cells_by_age_sex.csv": ["Age_bin", "Sex"],
        "skin_cells_by_broad_cell_type.csv": ["Broad cell type"],
        "skin_cells_by_granular_cell_type.csv": ["Granular cell type"],
        "skin_cells_by_sample_cell_type.csv": ["Sample ID", "Broad cell type"],
    }
    for filename, columns in audit_specs.items():
        table = counts_table(skin.obs, columns)
        if not table.empty:
            table.to_csv(results_dir / filename, index=False)

    skin.write_h5ad(output_path, compression="gzip")
    logging.info("Source atlas: %d cells x %d genes", len(skin_mask), skin.n_vars)
    logging.info("Skin subset: %d cells x %d genes", skin.n_obs, skin.n_vars)
    logging.info("Skin sites: %s", skin.obs["Tissue Site Detail"].astype(str).value_counts().to_dict())
    logging.info("Wrote %s", output_path)
    logging.info("Audit tables: %s", results_dir)


if __name__ == "__main__":
    main()
