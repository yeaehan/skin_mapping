#!/usr/bin/env python
"""Per-sample MuSpAn hotspots for the 21-gene core AUCell score and SnC proximity."""
from pathlib import Path

import matplotlib.pyplot as plt
import muspan as ms
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import Delaunay, cKDTree

ROOT = Path(__file__).resolve().parents[3]
H5AD = ROOT / "data" / "Xenium" / "xenium_full.h5ad"
OUT = ROOT / "results" / "public_validation" / "xenium_core21_muspan_hotspot_snc"
GENES = ["CXCL16", "VEGFA", "IL18", "IL1B", "CSF1", "IL4", "CCL19", "CXCL1", "CCL28", "CCL26", "CXCL11", "CCL11", "EGF", "CXCL5", "IL6", "CXCL6", "IL5", "MMP3", "CCL7", "CCL15", "CCL16"]


def aucell(matrix: np.ndarray, signature: np.ndarray, fraction: float = 0.05) -> np.ndarray:
    cutoff = max(1, int(matrix.shape[1] * fraction))
    top = np.argpartition(-matrix, cutoff - 1, axis=1)[:, :cutoff]
    top = np.take_along_axis(top, np.argsort(-np.take_along_axis(matrix, top, axis=1), axis=1), axis=1)
    hits = np.isin(top, signature)
    possible = min(len(signature), cutoff)
    maximum = possible * cutoff - possible * (possible - 1) / 2
    return (hits * np.arange(cutoff, 0, -1)).sum(axis=1) / maximum


def getis_ord(coords: np.ndarray, scores: np.ndarray) -> pd.DataFrame:
    domain = ms.domain("Xenium core21")
    domain.add_points(coords, collection_name="cells")
    domain.add_labels("core21 AUCell", scores, add_labels_to=("collection", "cells"), label_type="continuous")
    z, p, indices = ms.spatial_statistics.getis_ord(
        domain, population=("collection", "cells"), label_name="core21 AUCell", alpha=0.05,
        network_kwargs={"network_type": "Delaunay", "min_edge_distance": 0, "max_edge_distance": np.inf},
    )
    result = pd.DataFrame({"getis_z": np.nan, "getis_p_adj": np.nan, "hotspot_class": "ns"}, index=np.arange(len(scores)))
    indices, z, p = np.asarray(indices, int), np.asarray(z, float), np.asarray(p, float)
    result.loc[indices, ["getis_z", "getis_p_adj"]] = np.column_stack([z, p])
    result.loc[indices[(p <= 0.05) & (z > 0)], "hotspot_class"] = "hot"
    result.loc[indices[(p <= 0.05) & (z < 0)], "hotspot_class"] = "cold"
    return result


def hotspot_ids(coords: np.ndarray, hot: np.ndarray) -> np.ndarray:
    """Label contiguous significant cells using the same Delaunay neighborhood concept."""
    labels = np.full(len(coords), -1, dtype=int)
    hot_indices = np.flatnonzero(hot)
    if not len(hot_indices):
        return labels
    edges = Delaunay(coords).simplices[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2)
    hot_edges = edges[hot[edges].all(axis=1)]
    lookup = np.full(len(coords), -1, dtype=int)
    lookup[hot_indices] = np.arange(len(hot_indices))
    graph = csr_matrix((np.ones(len(hot_edges) * 2),
                        (np.r_[lookup[hot_edges[:, 0]], lookup[hot_edges[:, 1]]],
                         np.r_[lookup[hot_edges[:, 1]], lookup[hot_edges[:, 0]]])),
                       shape=(len(hot_indices), len(hot_indices)))
    _, component = connected_components(graph, directed=False)
    labels[hot_indices] = component + 1
    return labels


def old_young_tests(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = ["hotspot_fraction", "mean_hotspot_score", "hotspot_minus_background",
               "snc_hotspot_enrichment", "median_snc_to_hotspot_um", "median_hotspot_to_snc_um"]
    for metric in metrics:
        young = summary.loc[summary["age"].eq("Young"), metric].dropna()
        old = summary.loc[summary["age"].eq("Old"), metric].dropna()
        u, p = stats.mannwhitneyu(old, young, alternative="two-sided") if len(old) and len(young) else (np.nan, np.nan)
        rows.append({"metric": metric, "n_young": len(young), "n_old": len(old),
                     "young_mean": young.mean(), "old_mean": old.mean(), "old_minus_young": old.mean() - young.mean(),
                     "old_over_young": old.mean() / young.mean() if len(young) and young.mean() != 0 else np.nan,
                     "mannwhitney_u": u, "mannwhitney_p": p})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(H5AD)
    missing = sorted(set(GENES) - set(adata.var_names))
    if missing:
        raise ValueError(f"Missing core genes: {missing}")
    matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    signature = adata.var_names.get_indexer(GENES)
    scores = aucell(matrix, signature)
    samples = adata.obs["sample_concat"].astype(str).to_numpy()
    ages = adata.obs["age"].astype(str).to_numpy()
    snc = adata.obs["snc"].astype(str).eq("SnC+").to_numpy()
    coords = np.asarray(adata.obsm["spatial"])
    cell_rows, summaries, hotspot_rows = [], [], []

    for sample in sorted(np.unique(samples)):
        idx = np.flatnonzero(samples == sample)
        result = getis_ord(coords[idx], scores[idx])
        hot = result["hotspot_class"].eq("hot").to_numpy()
        result["hotspot_id"] = hotspot_ids(coords[idx], hot)
        local_snc = snc[idx]
        local_coords = coords[idx]
        distance_snc_to_hot = cKDTree(local_coords[hot]).query(local_coords[local_snc])[0] if hot.any() and local_snc.any() else np.array([])
        distance_hot_to_snc = cKDTree(local_coords[local_snc]).query(local_coords[hot])[0] if hot.any() and local_snc.any() else np.array([])
        background_mean = scores[idx][~hot].mean() if (~hot).any() else np.nan
        hotspot_mean = scores[idx][hot].mean() if hot.any() else np.nan
        expected = hot.mean()
        summaries.append({"sample": sample, "age": ages[idx][0], "n_cells": len(idx), "n_snc": int(local_snc.sum()),
                          "n_hotspots": int(result.loc[hot, "hotspot_id"].nunique()),
                          "n_hotspot_cells": int(hot.sum()), "hotspot_fraction": hot.mean(),
                          "mean_all_score": scores[idx].mean(), "mean_hotspot_score": hotspot_mean,
                          "mean_background_score": background_mean, "hotspot_minus_background": hotspot_mean - background_mean,
                          "snc_in_hotspot_fraction": hot[local_snc].mean() if local_snc.any() else np.nan,
                          "snc_hotspot_enrichment": hot[local_snc].mean() / expected if local_snc.any() and expected > 0 else np.nan,
                          "median_snc_to_hotspot_um": np.median(distance_snc_to_hot) if len(distance_snc_to_hot) else np.nan,
                          "median_hotspot_to_snc_um": np.median(distance_hot_to_snc) if len(distance_hot_to_snc) else np.nan})
        frame = pd.DataFrame({"cell_id": adata.obs_names[idx], "sample": sample, "age": ages[idx], "core21_aucell": scores[idx],
                              "snc": local_snc, "x": local_coords[:, 0], "y": local_coords[:, 1]})
        frame = pd.concat([frame.reset_index(drop=True), result.reset_index(drop=True)], axis=1)
        cell_rows.append(frame)
        hotspot_rows.append(frame.loc[hot].groupby("hotspot_id", as_index=False).agg(
            n_cells=("cell_id", "size"), mean_core21_score=("core21_aucell", "mean"),
            n_snc=("snc", "sum")).assign(sample=sample, age=ages[idx][0]))

    cells = pd.concat(cell_rows, ignore_index=True)
    summary = pd.DataFrame(summaries)
    hotspots = pd.concat(hotspot_rows, ignore_index=True)
    tests = old_young_tests(summary)
    cells.to_csv(OUT / "core21_muspan_cell_results.csv.gz", index=False, compression="gzip")
    summary.to_csv(OUT / "core21_muspan_sample_summary.csv", index=False)
    hotspots.to_csv(OUT / "core21_muspan_hotspot_summary.csv", index=False)
    tests.to_csv(OUT / "core21_muspan_old_vs_young_tests.csv", index=False)

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    sample_order = summary.sort_values(["age", "sample"], ascending=[False, True])["sample"]
    for ax, sample in zip(axes.flat, sample_order):
        frame = cells.loc[cells["sample"].eq(sample)]
        ax.scatter(frame["x"], frame["y"], c=frame["core21_aucell"], s=1, cmap="viridis")
        hot = frame["hotspot_class"].eq("hot")
        ax.scatter(frame.loc[hot, "x"], frame.loc[hot, "y"], facecolors="none", edgecolors="red", s=5, linewidths=.3)
        positive = frame["snc"]
        ax.scatter(frame.loc[positive, "x"], frame.loc[positive, "y"], c="black", marker="x", s=8, linewidths=.5)
        age = frame["age"].iloc[0]
        ax.set_title(f"{sample} ({age}) | hot={hot.mean():.1%}")
        ax.set_xticks([]); ax.set_yticks([])
    fig.suptitle("Core-21 AUCell: MuSpAn hotspots (red) and SnC+ cells (black)")
    fig.tight_layout()
    fig.savefig(OUT / "core21_muspan_hotspots_snc_spatial.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    age_colors = {"Young": "#9ecae1", "Old": "#de2d26"}
    ordered = summary.set_index("sample").loc[sample_order].reset_index()
    axes[0].bar(range(len(ordered)), ordered["n_hotspots"], color=ordered["age"].map(age_colors))
    axes[0].set_xticks(range(len(ordered)), ordered["sample"], rotation=45, ha="right")
    axes[0].set_ylabel("Number of contiguous hotspots")
    for age, frame in hotspots.groupby("age"):
        axes[1].scatter(frame["n_cells"], frame["mean_core21_score"], s=18,
                        alpha=.7, color=age_colors[age], label=age)
    axes[1].set_xlabel("Cells in hotspot")
    axes[1].set_ylabel("Mean Core21 AUCell score")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "core21_hotspot_count_and_mean_score.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    metrics = [("mean_hotspot_score", "Mean AUCell in hotspots"),
               ("hotspot_fraction", "Hotspot-cell fraction"),
               ("median_snc_to_hotspot_um", "Median SnC-to-hotspot distance (µm)")]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for ax, (metric, label) in zip(axes, metrics):
        for x, age in enumerate(["Young", "Old"]):
            values = summary.loc[summary["age"].eq(age), metric].dropna()
            ax.scatter(np.full(len(values), x), values, s=45)
            if len(values): ax.plot([x - .15, x + .15], [values.mean()] * 2, color="black")
        test = tests.loc[tests["metric"].eq(metric)].iloc[0]
        ax.set_xticks([0, 1], ["Young", "Old"]); ax.set_ylabel(label); ax.set_title(f"Old/Young={test.old_over_young:.2f}\np={test.mannwhitney_p:.3g}")
    fig.tight_layout()
    fig.savefig(OUT / "core21_hotspot_old_young_and_snc_proximity.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    assert len(summary) == 6 and len(cells) == adata.n_obs


if __name__ == "__main__":
    main()
