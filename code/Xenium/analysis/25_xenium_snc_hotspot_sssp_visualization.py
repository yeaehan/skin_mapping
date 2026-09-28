"""
Visualize SSSP signal in and around Xenium SnC hotspots.

This script does not recompute Getis-Ord hotspots. It uses the existing
SnC-hotspot cell neighborhood assignments and joins in the Scanpy SSSP score
from the Xenium h5ad object.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


PROJECT_ROOT = Path(__file__).resolve().parents[3]
H5AD = PROJECT_ROOT / "results" / "260727" / "xenium_sssp_20gene_scanpy_score_genes" / "xenium_20gene_sssp_scanpy_scored.h5ad"
NEIGHBORHOODS = PROJECT_ROOT / "results" / "260727" / "xenium_snc_hotspot_neighborhoods" / "xenium_snc_hotspot_cell_neighborhoods.csv"
HOTSPOT_RESULTS = PROJECT_ROOT / "results" / "260727" / "xenium_snc_hotspots" / "xenium_snc_hotspot_results_long.csv"
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_snc_hotspot_sssp"

SCORE_COL = "sssp_20gene_scanpy_score"
HOTSPOT_METRIC = "snc_prop"
ZONE_ORDER = ["hotspot_core", "neighbor_1", "neighbor_2", "background"]
ZONE_LABELS = {
    "hotspot_core": "Hotspot core",
    "neighbor_1": "Neighbor 1",
    "neighbor_2": "Neighbor 2",
    "background": "Background",
}
AGE_ORDER = ["Young", "Old"]
AGE_COLORS = {"Young": "#2C7FB8", "Old": "#D95F02"}
ZONE_COLORS = {
    "hotspot_core": "#D73027",
    "neighbor_1": "#FC8D59",
    "neighbor_2": "#FEE08B",
    "background": "#C9C9C9",
}


def load_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    adata = ad.read_h5ad(H5AD, backed="r")
    score_cols = [
        "cell_id",
        "sample",
        "age",
        "fov",
        "x",
        "y",
        "snc",
        "cell_type_major",
        "cell_type_minor",
        SCORE_COL,
    ]
    obs = adata.obs[score_cols].copy().reset_index(drop=True)
    adata.file.close()

    neighborhoods = pd.read_csv(NEIGHBORHOODS)
    keep_cols = [
        "cell_id",
        "hex_q",
        "hex_r",
        "hex_id",
        "hex_distance_to_snc_hotspot",
        "zone",
    ]
    df = obs.merge(neighborhoods[keep_cols], on="cell_id", how="left", validate="one_to_one")
    if df["zone"].isna().any():
        missing = int(df["zone"].isna().sum())
        raise ValueError(f"Missing hotspot zone for {missing} cells")

    df["zone"] = pd.Categorical(df["zone"], categories=ZONE_ORDER, ordered=True)
    df["age"] = pd.Categorical(df["age"], categories=AGE_ORDER, ordered=True)

    # Top quartile is defined within sample, matching the proximity-analysis convention.
    df["sssp_high_q4_within_sample"] = False
    for sample, idx in df.groupby("sample", observed=True).groups.items():
        sample_scores = df.loc[idx, SCORE_COL].astype(float)
        cutoff = sample_scores.quantile(0.75)
        df.loc[idx, "sssp_high_q4_within_sample"] = sample_scores.ge(cutoff)

    hotspot_bins = pd.read_csv(HOTSPOT_RESULTS)
    hotspot_bins = hotspot_bins.loc[hotspot_bins["metric"].eq(HOTSPOT_METRIC)].copy()
    hotspot_bins["is_snc_hotspot_hex"] = hotspot_bins["hotspot_class"].eq("Hot spot")

    return df, hotspot_bins


def summarize_by_sample(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    by_sample = (
        df.groupby(["sample", "age", "zone"], observed=True)
        .agg(
            n_cells=("cell_id", "size"),
            mean_sssp=(SCORE_COL, "mean"),
            median_sssp=(SCORE_COL, "median"),
            fraction_sssp_high=("sssp_high_q4_within_sample", "mean"),
            snc_fraction=("snc", lambda x: x.astype(str).eq("SnC+").mean()),
        )
        .reset_index()
    )

    age_summary = (
        by_sample.groupby(["age", "zone"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_sssp_sample_mean=("mean_sssp", "mean"),
            sem_sssp_sample_mean=("mean_sssp", lambda x: x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan),
            mean_fraction_sssp_high=("fraction_sssp_high", "mean"),
            sem_fraction_sssp_high=("fraction_sssp_high", lambda x: x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan),
            mean_snc_fraction=("snc_fraction", "mean"),
        )
        .reset_index()
    )

    test_rows = []
    for zone, zdf in by_sample.groupby("zone", observed=True):
        young = zdf.loc[zdf["age"].astype(str).eq("Young"), "mean_sssp"].dropna()
        old = zdf.loc[zdf["age"].astype(str).eq("Old"), "mean_sssp"].dropna()
        young_high = zdf.loc[zdf["age"].astype(str).eq("Young"), "fraction_sssp_high"].dropna()
        old_high = zdf.loc[zdf["age"].astype(str).eq("Old"), "fraction_sssp_high"].dropna()
        row = {"zone": zone, "n_young_samples": len(young), "n_old_samples": len(old)}
        for metric_name, old_values, young_values in [
            ("mean_sssp", old, young),
            ("fraction_sssp_high", old_high, young_high),
        ]:
            if len(old_values) >= 2 and len(young_values) >= 2:
                t_stat, t_p = stats.ttest_ind(old_values, young_values, equal_var=False)
                u_stat, u_p = stats.mannwhitneyu(old_values, young_values, alternative="two-sided")
            else:
                t_stat = t_p = u_stat = u_p = np.nan
            row.update(
                {
                    f"{metric_name}_young_mean": young_values.mean(),
                    f"{metric_name}_old_mean": old_values.mean(),
                    f"{metric_name}_old_minus_young": old_values.mean() - young_values.mean(),
                    f"{metric_name}_welch_p": t_p,
                    f"{metric_name}_mannwhitney_p": u_p,
                }
            )
        test_rows.append(row)
    tests = pd.DataFrame(test_rows)
    return by_sample, age_summary, tests


def summarize_hexes(df: pd.DataFrame, hotspot_bins: pd.DataFrame) -> pd.DataFrame:
    hex_sssp = (
        df.groupby(["sample", "age", "fov", "hex_id"], observed=True)
        .agg(
            n_cells_with_score=("cell_id", "size"),
            mean_sssp=(SCORE_COL, "mean"),
            median_sssp=(SCORE_COL, "median"),
            fraction_sssp_high=("sssp_high_q4_within_sample", "mean"),
        )
        .reset_index()
    )
    merge_cols = [
        "sample",
        "age",
        "fov",
        "hex_id",
        "hex_q",
        "hex_r",
        "hex_x",
        "hex_y",
        "n_cells",
        "n_snc",
        "snc_prop",
        "mean_sasp_num",
        "gi_star_z",
        "gi_star_fdr",
        "hotspot_class",
        "is_snc_hotspot_hex",
    ]
    return hotspot_bins[merge_cols].merge(hex_sssp, on=["sample", "age", "fov", "hex_id"], how="left")


def plot_zone_sample_points(by_sample: pd.DataFrame, y_col: str, ylabel: str, filename: str) -> None:
    plot_df = by_sample.copy()
    plot_df["zone_label"] = plot_df["zone"].map(ZONE_LABELS)
    plt.figure(figsize=(8.5, 5.0))
    ax = sns.pointplot(
        data=plot_df,
        x="zone_label",
        y=y_col,
        hue="age",
        order=[ZONE_LABELS[z] for z in ZONE_ORDER],
        hue_order=AGE_ORDER,
        palette=AGE_COLORS,
        errorbar="se",
        dodge=0.35,
        join=False,
    )
    sns.stripplot(
        data=plot_df,
        x="zone_label",
        y=y_col,
        hue="age",
        order=[ZONE_LABELS[z] for z in ZONE_ORDER],
        hue_order=AGE_ORDER,
        palette=AGE_COLORS,
        dodge=True,
        size=7,
        alpha=0.8,
        edgecolor="black",
        linewidth=0.4,
        ax=ax,
    )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:2], labels[:2], title="Age", frameon=False)
    ax.set_xlabel("SnC hotspot zone")
    ax.set_ylabel(ylabel)
    ax.set_title("SSSP signal by SnC hotspot zone")
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(OUTDIR / f"{filename}.png", dpi=300)
    plt.savefig(OUTDIR / f"{filename}.pdf")
    plt.close()


def plot_zone_age_bars(age_summary: pd.DataFrame) -> None:
    plot_df = age_summary.copy()
    plot_df["zone_label"] = plot_df["zone"].map(ZONE_LABELS)
    plt.figure(figsize=(8.2, 4.8))
    ax = sns.barplot(
        data=plot_df,
        x="zone_label",
        y="mean_sssp_sample_mean",
        hue="age",
        order=[ZONE_LABELS[z] for z in ZONE_ORDER],
        hue_order=AGE_ORDER,
        palette=AGE_COLORS,
        errorbar=None,
    )
    for container, age in zip(ax.containers, AGE_ORDER):
        sub = plot_df.loc[plot_df["age"].astype(str).eq(age)].set_index("zone")
        for patch, zone in zip(container, ZONE_ORDER):
            sem = sub.loc[zone, "sem_sssp_sample_mean"] if zone in sub.index else np.nan
            if pd.notna(sem):
                x = patch.get_x() + patch.get_width() / 2
                y = patch.get_height()
                ax.errorbar(x, y, yerr=sem, color="black", capsize=3, linewidth=1)
    ax.set_xlabel("SnC hotspot zone")
    ax.set_ylabel("Mean SSSP score across samples")
    ax.set_title("Age-level SSSP signal in SnC hotspot zones")
    ax.legend(title="Age", frameon=False)
    ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    plt.savefig(OUTDIR / "sssp_score_by_hotspot_zone_age_summary.png", dpi=300)
    plt.savefig(OUTDIR / "sssp_score_by_hotspot_zone_age_summary.pdf")
    plt.close()


def plot_cell_spatial_sssp(df: pd.DataFrame) -> None:
    samples = sorted(df["sample"].unique())
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.5), squeeze=False)
    vmin = df[SCORE_COL].quantile(0.02)
    vmax = df[SCORE_COL].quantile(0.98)
    for ax, sample in zip(axes.ravel(), samples):
        sdf = df.loc[df["sample"].eq(sample)].copy()
        sc = ax.scatter(
            sdf["x"],
            sdf["y"],
            c=sdf[SCORE_COL],
            s=2.2,
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            linewidths=0,
            rasterized=True,
        )
        core = sdf.loc[sdf["zone"].astype(str).eq("hotspot_core")]
        if not core.empty:
            ax.scatter(core["x"], core["y"], s=5.5, facecolors="none", edgecolors="#D73027", linewidths=0.35, rasterized=True)
        ax.set_title(f"{sample} ({sdf['age'].iloc[0]})")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes.ravel()[len(samples):]:
        ax.axis("off")
    cbar = fig.colorbar(sc, ax=axes.ravel().tolist(), fraction=0.025, pad=0.015)
    cbar.set_label("SSSP Scanpy score")
    fig.suptitle("Spatial SSSP score with SnC hotspot-core cells outlined", y=0.99)
    plt.tight_layout(rect=[0, 0, 0.98, 0.96])
    fig.savefig(OUTDIR / "spatial_sssp_score_with_snc_hotspot_core_overlay.png", dpi=300)
    fig.savefig(OUTDIR / "spatial_sssp_score_with_snc_hotspot_core_overlay.pdf")
    plt.close(fig)


def plot_hex_spatial_sssp(hex_sssp: pd.DataFrame) -> None:
    samples = sorted(hex_sssp["sample"].unique())
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.5), squeeze=False)
    vmin = hex_sssp["mean_sssp"].quantile(0.02)
    vmax = hex_sssp["mean_sssp"].quantile(0.98)
    for ax, sample in zip(axes.ravel(), samples):
        sdf = hex_sssp.loc[hex_sssp["sample"].eq(sample)].copy()
        sc = ax.scatter(
            sdf["hex_x"],
            sdf["hex_y"],
            c=sdf["mean_sssp"],
            s=34,
            marker="h",
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            linewidths=0,
        )
        hot = sdf.loc[sdf["is_snc_hotspot_hex"]]
        if not hot.empty:
            ax.scatter(hot["hex_x"], hot["hex_y"], s=72, marker="h", facecolors="none", edgecolors="#D73027", linewidths=0.9)
        ax.set_title(f"{sample} ({sdf['age'].iloc[0]})")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes.ravel()[len(samples):]:
        ax.axis("off")
    cbar = fig.colorbar(sc, ax=axes.ravel().tolist(), fraction=0.025, pad=0.015)
    cbar.set_label("Mean SSSP score per hex")
    fig.suptitle("Hex-level SSSP score with SnC hotspot hexes outlined", y=0.99)
    plt.tight_layout(rect=[0, 0, 0.98, 0.96])
    fig.savefig(OUTDIR / "hex_mean_sssp_score_with_snc_hotspot_overlay.png", dpi=300)
    fig.savefig(OUTDIR / "hex_mean_sssp_score_with_snc_hotspot_overlay.pdf")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df, hotspot_bins = load_data()
    by_sample, age_summary, tests = summarize_by_sample(df)
    hex_sssp = summarize_hexes(df, hotspot_bins)

    df.to_csv(OUTDIR / "xenium_snc_hotspot_sssp_cell_table.csv.gz", index=False)
    by_sample.to_csv(OUTDIR / "xenium_snc_hotspot_sssp_by_sample_zone.csv", index=False)
    age_summary.to_csv(OUTDIR / "xenium_snc_hotspot_sssp_by_age_zone.csv", index=False)
    tests.to_csv(OUTDIR / "xenium_snc_hotspot_sssp_old_vs_young_zone_tests.csv", index=False)
    hex_sssp.to_csv(OUTDIR / "xenium_snc_hotspot_hex_sssp_summary.csv", index=False)

    plot_zone_sample_points(
        by_sample,
        "mean_sssp",
        "Mean SSSP Scanpy score",
        "sssp_mean_score_by_hotspot_zone_sample_points",
    )
    plot_zone_sample_points(
        by_sample,
        "fraction_sssp_high",
        "Fraction SSSP-high cells",
        "sssp_high_fraction_by_hotspot_zone_sample_points",
    )
    plot_zone_age_bars(age_summary)
    plot_cell_spatial_sssp(df)
    plot_hex_spatial_sssp(hex_sssp)

    run_summary = pd.DataFrame(
        [
            {
                "h5ad": str(H5AD),
                "neighborhood_input": str(NEIGHBORHOODS),
                "hotspot_input": str(HOTSPOT_RESULTS),
                "score_col": SCORE_COL,
                "sssp_high_definition": "top quartile within sample",
                "hotspot_metric": HOTSPOT_METRIC,
                "n_cells": len(df),
                "n_hexes": len(hex_sssp),
                "output_dir": str(OUTDIR),
            }
        ]
    )
    run_summary.to_csv(OUTDIR / "xenium_snc_hotspot_sssp_run_summary.csv", index=False)
    print(f"Done: {OUTDIR}")


if __name__ == "__main__":
    main()
