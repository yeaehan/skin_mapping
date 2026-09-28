#!/usr/bin/env python
"""MuSpAn sensitivity, weighted-PCF, and adjacency analyses for Core21 and SnC cells."""
from pathlib import Path

import matplotlib.pyplot as plt
import muspan as ms
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[3]
H5AD = ROOT / "data" / "Xenium" / "xenium_full.h5ad"
OUT = ROOT / "results" / "public_validation" / "xenium_core21_muspan_multiscale"
GENES = ["CXCL16", "VEGFA", "IL18", "IL1B", "CSF1", "IL4", "CCL19", "CXCL1", "CCL28", "CCL26", "CXCL11", "CCL11", "EGF", "CXCL5", "IL6", "CXCL6", "IL5", "MMP3", "CCL7", "CCL15", "CCL16"]
RADII = (25, 50, 100)


def aucell(matrix: np.ndarray, signature: np.ndarray, fraction: float = .05) -> np.ndarray:
    cutoff = max(1, int(matrix.shape[1] * fraction))
    top = np.argpartition(-matrix, cutoff - 1, axis=1)[:, :cutoff]
    top = np.take_along_axis(top, np.argsort(-np.take_along_axis(matrix, top, axis=1), axis=1), axis=1)
    hits = np.isin(top, signature)
    possible = min(len(signature), cutoff)
    maximum = possible * cutoff - possible * (possible - 1) / 2
    return (hits * np.arange(cutoff, 0, -1)).sum(axis=1) / maximum


def make_domain(coords, scores, snc):
    domain = ms.domain("Xenium Core21")
    domain.add_points(coords, collection_name="cells")
    domain.add_labels("Core21 AUCell", scores, add_labels_to=("collection", "cells"), label_type="continuous")
    domain.add_labels("SnC status", np.where(snc, "SnC+", "SnC-"),
                      add_labels_to=("collection", "cells"), label_type="categorical")
    return domain


def hotspots(domain, radius, n):
    z, p, indices = ms.spatial_statistics.getis_ord(
        domain, label_name="Core21 AUCell", population=("collection", "cells"), alpha=.05,
        add_local_value_as_label=False,
        network_kwargs={"network_type": "Proximity", "min_edge_distance": 0, "max_edge_distance": radius},
    )
    hot = np.zeros(n, dtype=bool)
    hot[np.asarray(indices, int)[(np.asarray(p) <= .05) & (np.asarray(z) > 0)]] = True
    return hot


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(H5AD)
    missing = sorted(set(GENES) - set(adata.var_names))
    if missing:
        raise ValueError(f"Missing Core21 genes: {missing}")
    matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    scores = aucell(matrix, adata.var_names.get_indexer(GENES))
    samples = adata.obs["sample_concat"].astype(str).to_numpy()
    ages = adata.obs["age"].astype(str).to_numpy()
    snc = adata.obs["snc"].astype(str).eq("SnC+").to_numpy()
    coords = np.asarray(adata.obsm["spatial"])
    robustness_rows, wpcf_rows, apt_rows, cell_rows = [], [], [], []

    for sample in sorted(np.unique(samples)):
        keep = samples == sample
        xy, score, positive = coords[keep], scores[keep], snc[keep]
        age = ages[keep][0]
        domain = make_domain(xy, score, positive)
        hot_by_radius = {radius: hotspots(domain, radius, len(xy)) for radius in RADII}
        stable = np.column_stack(list(hot_by_radius.values())).all(axis=1)
        for radius, hot in hot_by_radius.items():
            robustness_rows.append({"sample": sample, "age": age, "radius_um": radius,
                                    "n_hotspot_cells": hot.sum(), "hotspot_fraction": hot.mean(),
                                    "mean_hotspot_score": score[hot].mean() if hot.any() else np.nan})
        cell_rows.append(pd.DataFrame({"cell_id": adata.obs_names[keep], "sample": sample, "age": age,
                                      **{f"hot_{r}um": h for r, h in hot_by_radius.items()},
                                      "stable_hotspot_all_radii": stable}))

        if positive.any():
            result = ms.spatial_statistics.weighted_pair_correlation_function(
                domain, population_A=np.flatnonzero(positive), population_B=np.arange(len(xy)),
                mark_pop_B="Core21 AUCell", max_R=100, annulus_step=10, annulus_width=10,
                weighting_function_pop_B=lambda _target, marks: marks,
                mark_min_pop_B=1e-9, mark_max_pop_B=2e-9, mark_step_pop_B=1,
                exclude_zero=True,
            )
            radii, values = np.asarray(result[0]), np.asarray(result[1]).ravel()
            wpcf_rows.extend({"sample": sample, "age": age, "radius_um": r, "weighted_pcf": value}
                              for r, value in zip(radii, values))

        hot50 = hot_by_radius[50]
        state = np.select([positive & hot50, positive, hot50],
                          ["SnC+ hotspot", "SnC+", "Core21 hotspot"], default="Neither")
        domain.add_labels("Spatial state", state, add_labels_to=("collection", "cells"), label_type="categorical")
        ms.networks.generate_network(domain, network_name="50um", objects_as_nodes=("collection", "cells"),
                                     network_type="Proximity", max_edge_distance=50)
        ses, significant_ses, labels = ms.networks.adjacency_permutation_test(
            domain, network_name="50um", label_name="Spatial state", label_shuffle_iterations=1000,
            observation_threshold=3, alpha=.05,
        )
        for i, source in enumerate(labels):
            for j, target in enumerate(labels):
                apt_rows.append({"sample": sample, "age": age, "source": source, "target": target,
                                 "standardized_effect_size": ses[i, j],
                                 "significant_fdr_0.05": np.isfinite(significant_ses[i, j])})

    robustness = pd.DataFrame(robustness_rows)
    wpcf = pd.DataFrame(wpcf_rows)
    apt = pd.DataFrame(apt_rows)
    cells = pd.concat(cell_rows, ignore_index=True)
    comparisons = []
    for analysis, frame, metrics in [
        ("hotspot_sensitivity", robustness, ["hotspot_fraction", "mean_hotspot_score"]),
        ("weighted_pcf", wpcf, ["weighted_pcf"]),
    ]:
        for radius, group in frame.groupby("radius_um"):
            for metric in metrics:
                young = group.loc[group["age"].eq("Young"), metric].dropna()
                old = group.loc[group["age"].eq("Old"), metric].dropna()
                _, p = mannwhitneyu(old, young, alternative="two-sided")
                comparisons.append({"analysis": analysis, "radius_um": radius, "metric": metric,
                                    "young_mean": young.mean(), "old_mean": old.mean(),
                                    "old_minus_young": old.mean() - young.mean(), "mannwhitney_p": p})
    focus = apt[(apt["source"].eq("SnC+")) & (apt["target"].eq("Core21 hotspot"))]
    young = focus.loc[focus["age"].eq("Young"), "standardized_effect_size"]
    old = focus.loc[focus["age"].eq("Old"), "standardized_effect_size"]
    _, p = mannwhitneyu(old, young, alternative="two-sided")
    comparisons.append({"analysis": "adjacency", "radius_um": 50,
                        "metric": "SnC+_to_Core21_hotspot_SES", "young_mean": young.mean(),
                        "old_mean": old.mean(), "old_minus_young": old.mean() - young.mean(),
                        "mannwhitney_p": p})
    robustness.to_csv(OUT / "core21_hotspot_radius_sensitivity.csv", index=False)
    wpcf.to_csv(OUT / "core21_snc_weighted_pcf.csv", index=False)
    apt.to_csv(OUT / "core21_snc_hotspot_adjacency_permutation.csv", index=False)
    pd.DataFrame(comparisons).to_csv(OUT / "core21_old_young_multiscale_tests.csv", index=False)
    cells.to_csv(OUT / "core21_hotspot_stability_by_cell.csv.gz", index=False, compression="gzip")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    colors = {"Young": "#9ecae1", "Old": "#de2d26"}
    for (sample, age), frame in robustness.groupby(["sample", "age"]):
        axes[0].plot(frame["radius_um"], frame["hotspot_fraction"], marker="o", color=colors[age], alpha=.65)
    axes[0].set(xlabel="Proximity radius (µm)", ylabel="Hotspot-cell fraction", title="Network sensitivity")
    for (sample, age), frame in wpcf.groupby(["sample", "age"]):
        axes[1].plot(frame["radius_um"], frame["weighted_pcf"], color=colors[age], alpha=.35)
    for age, frame in wpcf.groupby("age"):
        mean = frame.groupby("radius_um")["weighted_pcf"].mean()
        axes[1].plot(mean.index, mean, color=colors[age], linewidth=3, label=age)
    axes[1].set(xlabel="Distance from SnC+ cells (µm)", ylabel="Core21 weighted PCF", title="Core21 around SnC+")
    axes[1].legend(frameon=False)
    for age, frame in focus.groupby("age"):
        axes[2].scatter(np.repeat(age, len(frame)), frame["standardized_effect_size"], color=colors[age], s=45)
    axes[2].axhline(0, color="black", linewidth=.7)
    axes[2].set(ylabel="Adjacency standardized effect size", title="SnC+ ↔ Core21 hotspot")
    fig.tight_layout()
    fig.savefig(OUT / "core21_muspan_multiscale_summary.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    assert len(cells) == adata.n_obs and set(RADII) == set(robustness["radius_um"])


if __name__ == "__main__":
    main()
