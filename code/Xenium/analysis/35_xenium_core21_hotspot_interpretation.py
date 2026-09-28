#!/usr/bin/env python
"""Sample-level interpretation of previously computed Core21 MuSpAn hotspots."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "public_validation" / "xenium_core21_muspan_hotspot_snc"
H5AD = ROOT / "data" / "Xenium" / "xenium_full.h5ad"
GENES = "CXCL16 VEGFA IL18 IL1B CSF1 IL4 CCL19 CXCL1 CCL28 CCL26 CXCL11 CCL11 EGF CXCL5 IL6 CXCL6 IL5 MMP3 CCL7 CCL15 CCL16".split()
MIN_SIZE = 5


def main():
    cells = pd.read_csv(OUT / "core21_muspan_cell_results.csv.gz")
    samples = pd.read_csv(OUT / "core21_muspan_sample_summary.csv")
    regions = pd.read_csv(OUT / "core21_muspan_hotspot_summary.csv")
    rows = []
    for threshold in (3, 5, 10):
        for sample in samples.itertuples(index=False):
            keep = regions[(regions["sample"].eq(sample.sample)) & regions["n_cells"].ge(threshold)]
            n_hot_cells = int(keep["n_cells"].sum())
            rows.append({"sample": sample.sample, "age": sample.age, "minimum_hotspot_cells": threshold,
                         "n_substantial_hotspots": len(keep),
                         "hotspots_per_1000_cells": 1000 * len(keep) / sample.n_cells,
                         "substantial_hotspot_cell_fraction": n_hot_cells / sample.n_cells,
                         "median_hotspot_size": keep["n_cells"].median(),
                         "largest_hotspot_size": keep["n_cells"].max(),
                         "mean_hotspot_score": np.average(keep["mean_core21_score"], weights=keep["n_cells"]) if len(keep) else np.nan,
                         "hotspot_to_background_score_ratio":
                             (np.average(keep["mean_core21_score"], weights=keep["n_cells"]) / sample.mean_background_score)
                             if len(keep) else np.nan})
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "core21_substantial_hotspot_sample_metrics.csv", index=False)

    selected = metrics[metrics["minimum_hotspot_cells"].eq(MIN_SIZE)]
    panels = [("substantial_hotspot_cell_fraction", "% cells in hotspots ≥5 cells", 100),
              ("hotspots_per_1000_cells", "Hotspots ≥5 cells per 1,000 cells", 1),
              ("largest_hotspot_size", "Largest hotspot (cells)", 1),
              ("hotspot_to_background_score_ratio", "Hotspot/background Core21 score", 1)]
    colors = {"Young": "#9ecae1", "Old": "#de2d26"}
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    tests = []
    for ax, (metric, label, scale) in zip(axes.flat, panels):
        for x, age in enumerate(("Young", "Old")):
            values = selected.loc[selected["age"].eq(age), metric].dropna() * scale
            ax.scatter(np.full(len(values), x), values, s=55, color=colors[age])
            ax.plot([x - .15, x + .15], [values.mean()] * 2, color="black")
        young = selected.loc[selected["age"].eq("Young"), metric].dropna()
        old = selected.loc[selected["age"].eq("Old"), metric].dropna()
        _, p = mannwhitneyu(old, young, alternative="two-sided")
        tests.append({"metric": metric, "young_mean": young.mean(), "old_mean": old.mean(), "mannwhitney_p": p})
        ax.set_xticks((0, 1), ("Young", "Old")); ax.set_ylabel(label); ax.set_title(f"p={p:.3g}")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Sample-level Core21 hotspot features")
    fig.tight_layout()
    fig.savefig(OUT / "core21_substantial_hotspot_sample_comparison.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(tests).to_csv(OUT / "core21_substantial_hotspot_old_young_tests.csv", index=False)

    adata = sc.read_h5ad(H5AD)
    matrix = adata.X.tocsr() if hasattr(adata.X, "tocsr") else np.asarray(adata.X)
    drivers = []
    for sample in samples.itertuples(index=False):
        eligible = regions[(regions["sample"].eq(sample.sample)) & regions["n_cells"].ge(MIN_SIZE)]["hotspot_id"]
        hot_ids = set(cells.loc[cells["sample"].eq(sample.sample) & cells["hotspot_id"].isin(eligible), "cell_id"])
        sample_mask = adata.obs["sample_concat"].astype(str).eq(sample.sample).to_numpy()
        hot_mask = np.array([cell in hot_ids for cell in adata.obs_names])
        background_mask = sample_mask & ~hot_mask
        for gene in GENES:
            gene_idx = adata.var_names.get_loc(gene)
            hot_gene = float(matrix[hot_mask, gene_idx].sum())
            bg_gene = float(matrix[background_mask, gene_idx].sum())
            hot_total = float(matrix[hot_mask].sum())
            bg_total = float(matrix[background_mask].sum())
            hot_cpm = 1e6 * hot_gene / hot_total if hot_total else 0
            bg_cpm = 1e6 * bg_gene / bg_total if bg_total else 0
            drivers.append({"sample": sample.sample, "age": sample.age, "gene": gene,
                            "hotspot_cpm": hot_cpm, "background_cpm": bg_cpm,
                            "log2_hotspot_background": np.log2((hot_cpm + 1) / (bg_cpm + 1))})
    drivers = pd.DataFrame(drivers)
    drivers.to_csv(OUT / "core21_hotspot_gene_drivers.csv", index=False)
    order = samples.sort_values(["age", "sample"], ascending=[False, True])["sample"]
    table = drivers.pivot(index="gene", columns="sample", values="log2_hotspot_background").reindex(index=GENES, columns=order)
    fig, ax = plt.subplots(figsize=(8, 9))
    limit = max(1, np.nanmax(np.abs(table.to_numpy())))
    image = ax.imshow(table, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    ax.set_xticks(range(len(table.columns)), table.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(table.index)), table.index)
    ax.set_title("Genes driving substantial Core21 hotspots")
    fig.colorbar(image, ax=ax, label="log2 hotspot/background pseudobulk CPM")
    fig.tight_layout()
    fig.savefig(OUT / "core21_hotspot_gene_driver_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    metadata = adata.obs[["sample_concat", "age", "cell_type_major", "cell_type_minor"]].copy()
    metadata.index = adata.obs_names
    metadata.index.name = "cell_id"
    metadata = metadata.reset_index().merge(
        cells[["cell_id", "hotspot_id"]], on="cell_id", how="left", validate="one_to_one"
    )
    eligible_pairs = set(
        map(tuple, regions.loc[regions["n_cells"].ge(MIN_SIZE), ["sample", "hotspot_id"]].to_numpy())
    )
    metadata["in_substantial_hotspot"] = [
        (sample, hotspot) in eligible_pairs
        for sample, hotspot in zip(metadata["sample_concat"], metadata["hotspot_id"])
    ]
    hotspot_cells = metadata[metadata["in_substantial_hotspot"]].copy()
    per_hotspot = (hotspot_cells.groupby(["sample_concat", "age", "hotspot_id", "cell_type_major", "cell_type_minor"], observed=True)
                   .size().rename("n_cells").reset_index())
    per_hotspot["fraction_within_hotspot"] = per_hotspot["n_cells"] / per_hotspot.groupby(
        ["sample_concat", "hotspot_id"], observed=True)["n_cells"].transform("sum")
    per_hotspot.to_csv(OUT / "core21_substantial_hotspot_cell_type_composition.csv", index=False)

    categories = sorted(metadata["cell_type_minor"].astype(str).unique())
    composition = []
    for (sample, age), group in metadata.groupby(["sample_concat", "age"], observed=True):
        hot = group[group["in_substantial_hotspot"]]
        background = group[~group["in_substantial_hotspot"]]
        for cell_type in categories:
            hot_n = hot["cell_type_minor"].astype(str).eq(cell_type).sum()
            bg_n = background["cell_type_minor"].astype(str).eq(cell_type).sum()
            hot_fraction = hot_n / len(hot) if len(hot) else 0
            bg_fraction = bg_n / len(background) if len(background) else 0
            smoothed_hot = (hot_n + .5) / (len(hot) + .5 * len(categories))
            smoothed_bg = (bg_n + .5) / (len(background) + .5 * len(categories))
            composition.append({"sample": sample, "age": age, "cell_type_minor": cell_type,
                                "hotspot_n": hot_n, "hotspot_fraction": hot_fraction,
                                "background_n": bg_n, "background_fraction": bg_fraction,
                                "log2_hotspot_background_enrichment": np.log2(smoothed_hot / smoothed_bg)})
    composition = pd.DataFrame(composition)
    composition.to_csv(OUT / "core21_hotspot_cell_type_composition_by_sample.csv", index=False)

    top_types = (composition.groupby("cell_type_minor")["hotspot_n"].sum().nlargest(12).index)
    plot_data = composition.copy()
    plot_data["plot_type"] = np.where(plot_data["cell_type_minor"].isin(top_types),
                                      plot_data["cell_type_minor"], "Other")
    stacked = plot_data.groupby(["sample", "age", "plot_type"], observed=True)["hotspot_fraction"].sum().reset_index()
    sample_order = samples.sort_values(["age", "sample"], ascending=[False, True])["sample"]
    table = stacked.pivot(index="sample", columns="plot_type", values="hotspot_fraction").fillna(0).reindex(sample_order)
    fig, ax = plt.subplots(figsize=(11, 5.5))
    bottom = np.zeros(len(table))
    for cell_type in table.columns:
        ax.bar(table.index, table[cell_type], bottom=bottom, label=cell_type)
        bottom += table[cell_type].to_numpy()
    ax.set_ylabel("Fraction of substantial-hotspot cells")
    ax.set_title("Cell-type composition of Core21 hotspots")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "core21_hotspot_cell_type_composition_stacked.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    enrichment = composition[composition["cell_type_minor"].isin(top_types)].pivot(
        index="cell_type_minor", columns="sample", values="log2_hotspot_background_enrichment"
    ).reindex(index=top_types, columns=sample_order)
    fig, ax = plt.subplots(figsize=(8, 6.5))
    limit = max(1, np.nanmax(np.abs(enrichment.to_numpy())))
    image = ax.imshow(enrichment, cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    ax.set_xticks(range(len(enrichment.columns)), enrichment.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(enrichment.index)), enrichment.index)
    ax.set_title("Cell types enriched in Core21 hotspots")
    fig.colorbar(image, ax=ax, label="log2 hotspot/background proportion")
    fig.tight_layout()
    fig.savefig(OUT / "core21_hotspot_cell_type_enrichment_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    assert len(metrics) == len(samples) * 3 and len(drivers) == len(samples) * len(GENES)


if __name__ == "__main__":
    main()
