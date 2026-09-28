"""
Spatial visualization of Xenium SnC+ and SSSP-high cells.

This script reads the merged cell-level table produced by
18_xenium_snc_vs_sssp_score_proximity.py. It does not read the raw Xenium H5
directly; the SSSP score/high calls in that table come from the Scanpy-scored
H5AD generated in step 16.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
INPUT_CELL_TABLE = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_snc_vs_sssp_score_proximity"
    / "xenium_snc_sssp_top_quartile_score_proximity_cell_table.csv.gz"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_snc_sssp_spatial_plots"

SCORE_COL = "sssp_20gene_scanpy_score"
SSSP_HIGH_COL = "sssp_high_top_quartile"
GENE_SET_CONFIGS = [
    {
        "gene_set_id": "main_20",
        "input_cell_table": INPUT_CELL_TABLE,
        "outdir": OUTDIR,
        "score_col": "sssp_20gene_scanpy_score",
        "high_col": "sssp_high_top_quartile",
        "score_label": "20-gene SSSP score",
        "score_plot_prefix": "xenium_sssp_20gene_score",
    },
    {
        "gene_set_id": "diff_filtered_17",
        "input_cell_table": PROJECT_ROOT
        / "results"
        / "260727"
        / "xenium_snc_vs_sssp_17gene_diff_filtered_score_proximity"
        / "xenium_snc_sssp_top_quartile_score_proximity_cell_table.csv.gz",
        "outdir": PROJECT_ROOT / "results" / "260727" / "xenium_snc_sssp_17gene_diff_filtered_spatial_plots",
        "score_col": "sssp_17gene_diff_filtered_scanpy_score",
        "high_col": "sssp_high_top_quartile",
        "score_label": "17-gene differentiation-filtered SSSP score",
        "score_plot_prefix": "xenium_sssp_17gene_diff_filtered_score",
    },
]
AGE_ORDER = ["Young", "Old"]

CATEGORY_ORDER = [
    "Other",
    "SSSP-high only",
    "SnC+ only",
    "SnC+ & SSSP-high",
]
CATEGORY_COLORS = {
    "Other": "#D0D0D0",
    "SSSP-high only": "#2CA02C",
    "SnC+ only": "#D62728",
    "SnC+ & SSSP-high": "#6A3D9A",
}
CATEGORY_SIZES = {
    "Other": 2.2,
    "SSSP-high only": 4.8,
    "SnC+ only": 8.5,
    "SnC+ & SSSP-high": 10.0,
}
CATEGORY_ALPHAS = {
    "Other": 0.22,
    "SSSP-high only": 0.70,
    "SnC+ only": 0.95,
    "SnC+ & SSSP-high": 1.0,
}


def sample_sort_key(sample: str) -> tuple[int, str]:
    sample = str(sample)
    if sample.upper().endswith("-Y"):
        return 0, sample
    if sample.upper().endswith("-O"):
        return 1, sample
    return 2, sample


def load_data() -> pd.DataFrame:
    required_cols = [
        "cell_id",
        "sample",
        "age",
        "fov",
        "x",
        "y",
        "snc",
        SCORE_COL,
        SSSP_HIGH_COL,
    ]
    df = pd.read_csv(INPUT_CELL_TABLE, usecols=required_cols)
    df["is_snc_plus"] = df["snc"].astype(str).eq("SnC+")
    df["is_sssp_high"] = df[SSSP_HIGH_COL].astype(bool)
    df["spatial_category"] = np.select(
        [
            df["is_snc_plus"] & df["is_sssp_high"],
            df["is_snc_plus"],
            df["is_sssp_high"],
        ],
        [
            "SnC+ & SSSP-high",
            "SnC+ only",
            "SSSP-high only",
        ],
        default="Other",
    )
    df["spatial_category"] = pd.Categorical(
        df["spatial_category"], categories=CATEGORY_ORDER, ordered=True
    )
    return df


def set_spatial_axes(ax: plt.Axes, sub: pd.DataFrame) -> None:
    pad_x = max((sub["x"].max() - sub["x"].min()) * 0.03, 10)
    pad_y = max((sub["y"].max() - sub["y"].min()) * 0.03, 10)
    ax.set_xlim(sub["x"].min() - pad_x, sub["x"].max() + pad_x)
    ax.set_ylim(sub["y"].max() + pad_y, sub["y"].min() - pad_y)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)


def make_grid(n_panels: int) -> tuple[plt.Figure, np.ndarray]:
    ncols = 3
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.0 * ncols, 4.9 * nrows))
    axes = np.asarray(axes).reshape(-1)
    return fig, axes


def plot_overlay_by_sample(df: pd.DataFrame) -> None:
    samples = sorted(df["sample"].unique(), key=sample_sort_key)
    fig, axes = make_grid(len(samples))

    for ax, sample in zip(axes, samples):
        sub = df[df["sample"].eq(sample)]
        age = sub["age"].astype(str).iloc[0]
        for category in CATEGORY_ORDER:
            cat = sub[sub["spatial_category"].astype(str).eq(category)]
            if cat.empty:
                continue
            ax.scatter(
                cat["x"],
                cat["y"],
                s=CATEGORY_SIZES[category],
                c=CATEGORY_COLORS[category],
                alpha=CATEGORY_ALPHAS[category],
                linewidths=0,
                rasterized=True,
            )
        n_snc = int(sub["is_snc_plus"].sum())
        n_high = int(sub["is_sssp_high"].sum())
        n_both = int((sub["is_snc_plus"] & sub["is_sssp_high"]).sum())
        ax.set_title(
            f"{sample} ({age})\nSnC+ {n_snc} | SSSP-high {n_high} | both {n_both}",
            fontsize=10,
            fontweight="bold",
        )
        set_spatial_axes(ax, sub)

    for ax in axes[len(samples) :]:
        ax.axis("off")

    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=CATEGORY_COLORS[category],
            markeredgecolor="none",
            markersize=7,
            label=category,
            alpha=max(CATEGORY_ALPHAS[category], 0.55),
        )
        for category in CATEGORY_ORDER
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False)
    fig.suptitle("Spatial overlap of SnC+ and SSSP-high cells", fontsize=16, fontweight="bold")
    fig.tight_layout(rect=[0, 0.06, 1, 0.95])
    fig.savefig(OUTDIR / "xenium_snc_sssp_high_spatial_overlay_by_sample.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "xenium_snc_sssp_high_spatial_overlay_by_sample.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_score_by_sample(df: pd.DataFrame, score_label: str, score_plot_prefix: str) -> None:
    samples = sorted(df["sample"].unique(), key=sample_sort_key)
    vmin, vmax = df[SCORE_COL].quantile([0.02, 0.98])
    fig, axes = make_grid(len(samples))
    scatter = None

    for ax, sample in zip(axes, samples):
        sub = df[df["sample"].eq(sample)]
        age = sub["age"].astype(str).iloc[0]
        scatter = ax.scatter(
            sub["x"],
            sub["y"],
            s=3.0,
            c=sub[SCORE_COL],
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(f"{sample} ({age})", fontsize=10, fontweight="bold")
        set_spatial_axes(ax, sub)

    for ax in axes[len(samples) :]:
        ax.axis("off")

    if scatter is not None:
        cbar = fig.colorbar(scatter, ax=axes[: len(samples)], shrink=0.72, pad=0.015)
        cbar.set_label(score_label)
    fig.suptitle(f"Spatial distribution of {score_label}", fontsize=16, fontweight="bold")
    fig.savefig(OUTDIR / f"{score_plot_prefix}_spatial_by_sample.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / f"{score_plot_prefix}_spatial_by_sample.pdf", bbox_inches="tight")
    plt.close(fig)


def write_summary(df: pd.DataFrame, gene_set_id: str, score_label: str) -> pd.DataFrame:
    summary = (
        df.groupby(["sample", "age", "spatial_category"], observed=True)
        .agg(n_cells=("cell_id", "nunique"))
        .reset_index()
    )
    totals = df.groupby(["sample", "age"], observed=True)["cell_id"].nunique().reset_index(name="total_cells")
    summary = summary.merge(totals, on=["sample", "age"], how="left")
    summary["fraction_cells"] = summary["n_cells"] / summary["total_cells"]
    summary.to_csv(OUTDIR / "xenium_snc_sssp_high_spatial_category_counts_by_sample.csv", index=False)

    run_summary = pd.DataFrame(
        [
            {
                "gene_set_id": gene_set_id,
                "input_cell_table": str(INPUT_CELL_TABLE),
                "n_cells": int(len(df)),
                "n_samples": int(df["sample"].nunique()),
                "n_snc_plus": int(df["is_snc_plus"].sum()),
                "n_sssp_high": int(df["is_sssp_high"].sum()),
                "n_snc_plus_and_sssp_high": int((df["is_snc_plus"] & df["is_sssp_high"]).sum()),
                "sssp_high_definition": f"top quartile of {score_label}",
            }
        ]
    )
    run_summary.to_csv(OUTDIR / "xenium_snc_sssp_spatial_plot_run_summary.csv", index=False)
    return run_summary


def run_analysis(config: dict[str, object]) -> pd.DataFrame:
    global INPUT_CELL_TABLE, OUTDIR, SCORE_COL, SSSP_HIGH_COL
    INPUT_CELL_TABLE = Path(config["input_cell_table"])
    OUTDIR = Path(config["outdir"])
    SCORE_COL = str(config["score_col"])
    SSSP_HIGH_COL = str(config["high_col"])

    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = load_data()
    plot_overlay_by_sample(df)
    plot_score_by_sample(df, str(config["score_label"]), str(config["score_plot_prefix"]))
    run_summary = write_summary(df, str(config["gene_set_id"]), str(config["score_label"]))
    print("Spatial SnC/SSSP plots complete.")
    print(f"Input: {INPUT_CELL_TABLE}")
    print(f"Output: {OUTDIR}")
    return run_summary


def main() -> None:
    summaries = [run_analysis(config) for config in GENE_SET_CONFIGS]
    pd.concat(summaries, ignore_index=True).to_csv(
        PROJECT_ROOT / "results" / "260727" / "xenium_snc_sssp_spatial_plot_gene_set_run_summary.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
