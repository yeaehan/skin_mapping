"""
Squidpy neighborhood enrichment for SnC+ and SSSP-high cells.

This is a package-based complement to the nearest-distance/permutation analysis
in 18_xenium_snc_vs_sssp_score_proximity.py. It asks whether joint SnC/SSSP
categories are overrepresented as spatial graph neighbors within fixed radii.
"""

from __future__ import annotations

from pathlib import Path
import warnings

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import squidpy as sq
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_squidpy_snc_sssp_neighborhood_enrichment"
RADII_UM = [35, 70, 140]
N_PERMUTATIONS = 1000
RANDOM_SEED = 260727

GENE_SET_CONFIGS = [
    {
        "gene_set_id": "main_20",
        "input_cell_table": PROJECT_ROOT
        / "results"
        / "260727"
        / "xenium_snc_vs_sssp_score_proximity"
        / "xenium_snc_sssp_top_quartile_score_proximity_cell_table.csv.gz",
        "score_col": "sssp_20gene_scanpy_score",
        "high_col": "sssp_high_top_quartile",
    },
    {
        "gene_set_id": "diff_filtered_17",
        "input_cell_table": PROJECT_ROOT
        / "results"
        / "260727"
        / "xenium_snc_vs_sssp_17gene_diff_filtered_score_proximity"
        / "xenium_snc_sssp_top_quartile_score_proximity_cell_table.csv.gz",
        "score_col": "sssp_17gene_diff_filtered_scanpy_score",
        "high_col": "sssp_high_top_quartile",
    },
]

JOINT_CATEGORIES = [
    "Neither",
    "SSSP-high only",
    "SnC+ only",
    "SnC+ & SSSP-high",
]
FOCUS_PAIRS = [
    ("SSSP-high only", "SnC+ only"),
    ("SSSP-high only", "SnC+ & SSSP-high"),
    ("SnC+ only", "SnC+ & SSSP-high"),
    ("SnC+ & SSSP-high", "SnC+ & SSSP-high"),
]
PAIR_COLORS = {
    "SSSP-high only vs SnC+ only": "#2CA02C",
    "SSSP-high only vs SnC+ & SSSP-high": "#6A3D9A",
    "SnC+ only vs SnC+ & SSSP-high": "#D62728",
    "SnC+ & SSSP-high vs SnC+ & SSSP-high": "#9467BD",
}


def sample_sort_key(sample: str) -> tuple[int, str]:
    sample = str(sample)
    if sample.upper().endswith("-Y"):
        return 0, sample
    if sample.upper().endswith("-O"):
        return 1, sample
    return 2, sample


def make_joint_state(df: pd.DataFrame, high_col: str) -> pd.Series:
    is_snc = df["snc"].astype(str).eq("SnC+")
    is_high = df[high_col].astype(bool)
    state = np.select(
        [is_snc & is_high, is_snc, is_high],
        ["SnC+ & SSSP-high", "SnC+ only", "SSSP-high only"],
        default="Neither",
    )
    return pd.Categorical(state, categories=JOINT_CATEGORIES, ordered=True)


def load_cell_table(config: dict[str, object]) -> pd.DataFrame:
    path = Path(config["input_cell_table"])
    score_col = str(config["score_col"])
    high_col = str(config["high_col"])
    required_cols = ["cell_id", "sample", "age", "fov", "x", "y", "snc", score_col, high_col]
    df = pd.read_csv(path, usecols=required_cols)
    df = df[df["snc"].isin(["SnC-", "SnC+"])].copy()
    df["joint_state"] = make_joint_state(df, high_col)
    return df


def make_anndata(sub: pd.DataFrame) -> ad.AnnData:
    obs = sub[["cell_id", "sample", "age", "fov", "snc", "joint_state"]].copy()
    obs.index = obs["cell_id"].astype(str)
    obs["joint_state"] = pd.Categorical(obs["joint_state"], categories=JOINT_CATEGORIES, ordered=True)
    adata = ad.AnnData(X=sparse.csr_matrix((len(obs), 1), dtype=np.float32), obs=obs)
    adata.obsm["spatial"] = sub[["x", "y"]].to_numpy(dtype=float)
    return adata


def extract_pair_rows(
    adata: ad.AnnData,
    gene_set_id: str,
    sample: str,
    age: str,
    radius: int,
) -> list[dict[str, object]]:
    result = adata.uns["joint_state_nhood_enrichment"]
    zscores = result["zscore"]
    counts = result["count"]
    categories = list(adata.obs["joint_state"].cat.categories)
    cat_to_idx = {cat: i for i, cat in enumerate(categories)}
    rows = []
    for a, b in FOCUS_PAIRS:
        if a not in cat_to_idx or b not in cat_to_idx:
            continue
        i = cat_to_idx[a]
        j = cat_to_idx[b]
        rows.append(
            {
                "gene_set_id": gene_set_id,
                "sample": sample,
                "age": age,
                "radius_um": radius,
                "category_a": a,
                "category_b": b,
                "pair": f"{a} vs {b}",
                "squidpy_edge_count": float(counts[i, j]),
                "squidpy_zscore": float(zscores[i, j]),
                "n_cells": int(adata.n_obs),
                "n_permutations": N_PERMUTATIONS,
                "method": "squidpy.gr.spatial_neighbors + squidpy.gr.nhood_enrichment",
            }
        )
    return rows


def run_enrichment_for_gene_set(config: dict[str, object]) -> pd.DataFrame:
    gene_set_id = str(config["gene_set_id"])
    df = load_cell_table(config)
    rows = []

    for (sample, fov), sub in df.groupby(["sample", "fov"], observed=True, sort=False):
        sub = sub.reset_index(drop=True)
        age = str(sub["age"].iloc[0])
        if sub["joint_state"].nunique() < 2:
            continue
        for radius in RADII_UM:
            adata = make_anndata(sub)
            sq.gr.spatial_neighbors_radius(
                adata,
                radius=float(radius),
                spatial_key="spatial",
                set_diag=False,
                key_added="spatial",
                n_jobs=1,
            )
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sq.gr.nhood_enrichment(
                    adata,
                    cluster_key="joint_state",
                    connectivity_key="spatial_connectivities",
                    n_perms=N_PERMUTATIONS,
                    seed=RANDOM_SEED,
                    n_jobs=1,
                    show_progress_bar=False,
                )
            rows.extend(extract_pair_rows(adata, gene_set_id, str(sample), age, int(radius)))

    return pd.DataFrame(rows)


def plot_focus_pair(results: pd.DataFrame) -> None:
    focus_pair = "SSSP-high only vs SnC+ only"
    for gene_set_id, sub in results[results["pair"].eq(focus_pair)].groupby("gene_set_id", observed=True):
        sample_order = sorted(sub["sample"].unique(), key=sample_sort_key)
        fig, axes = plt.subplots(1, len(RADII_UM), figsize=(4.4 * len(RADII_UM), 4.0), sharey=True)
        axes = np.atleast_1d(axes)
        for ax, radius in zip(axes, RADII_UM):
            work = sub[sub["radius_um"].eq(radius)].set_index("sample").reindex(sample_order).reset_index()
            colors = ["#4C78A8" if str(sample).upper().endswith("-Y") else "#F58518" for sample in work["sample"]]
            ax.bar(
                np.arange(len(work)),
                work["squidpy_zscore"],
                color=colors,
                edgecolor="#333333",
                linewidth=0.7,
            )
            ax.axhline(0, color="#333333", linewidth=0.8)
            young_n = sum(str(sample).upper().endswith("-Y") for sample in sample_order)
            if 0 < young_n < len(sample_order):
                ax.axvline(young_n - 0.5, color="#777777", linewidth=0.9, linestyle="--", alpha=0.8)
            ax.set_xticks(np.arange(len(work)))
            ax.set_xticklabels(work["sample"], rotation=35, ha="right")
            ax.set_title(f"{radius} um", fontsize=11, fontweight="bold")
            ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
            ax.spines[["top", "right"]].set_visible(False)
        axes[0].set_ylabel("Squidpy neighborhood enrichment z-score")
        fig.suptitle(f"{gene_set_id}: SSSP-high only neighbors SnC+ only", fontsize=14, fontweight="bold")
        fig.tight_layout()
        fig.savefig(OUTDIR / f"{gene_set_id}_squidpy_sssp_high_only_to_snc_only_zscore_by_radius.png", dpi=300, bbox_inches="tight")
        fig.savefig(OUTDIR / f"{gene_set_id}_squidpy_sssp_high_only_to_snc_only_zscore_by_radius.pdf", bbox_inches="tight")
        plt.close(fig)


def write_summary(results: pd.DataFrame) -> pd.DataFrame:
    focus = results[results["pair"].eq("SSSP-high only vs SnC+ only")].copy()
    summary = (
        focus.groupby(["gene_set_id", "radius_um"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            n_samples_positive_z=("squidpy_zscore", lambda x: int((pd.to_numeric(x) > 0).sum())),
            mean_zscore=("squidpy_zscore", "mean"),
            median_zscore=("squidpy_zscore", "median"),
        )
        .reset_index()
    )
    summary.to_csv(OUTDIR / "squidpy_snc_sssp_focus_pair_summary.csv", index=False)
    return summary


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    results = pd.concat([run_enrichment_for_gene_set(config) for config in GENE_SET_CONFIGS], ignore_index=True)
    results.to_csv(OUTDIR / "squidpy_snc_sssp_joint_state_neighborhood_enrichment_by_sample.csv", index=False)
    summary = write_summary(results)
    plot_focus_pair(results)

    run_summary = pd.DataFrame(
        [
            {
                "gene_sets": ", ".join(str(config["gene_set_id"]) for config in GENE_SET_CONFIGS),
                "radii_um": ", ".join(str(radius) for radius in RADII_UM),
                "n_permutations": N_PERMUTATIONS,
                "method": "squidpy.gr.spatial_neighbors(radius) + squidpy.gr.nhood_enrichment",
                "cluster_key": "joint_state: Neither, SSSP-high only, SnC+ only, SnC+ & SSSP-high",
            }
        ]
    )
    run_summary.to_csv(OUTDIR / "squidpy_snc_sssp_neighborhood_enrichment_run_summary.csv", index=False)

    print("Squidpy SnC/SSSP neighborhood enrichment complete.")
    print(summary.to_string(index=False))
    print(f"\nWrote outputs to: {OUTDIR}")


if __name__ == "__main__":
    main()
