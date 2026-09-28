#!/usr/bin/env python
"""
Compare Sinai SnC calls across age bins, anatomical sites, and major cell types.

This uses the cell-level obs export from 02_sinai_age_decade_trends_all_sites.py.
The two SnC definitions compared here are:
  1. Original marker-positive SnC: snc_bool
  2. Senescence-index quartile SnC: snc_SI_Q4_within_lineage
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
INPUT_OBS = PROJECT_ROOT / "results" / "260727" / "sinai_age_decade_trends_all_sites" / "sinai_all_sites_snc_decade_obs.csv.gz"
OUTDIR = PROJECT_ROOT / "results" / "260727" / "sinai_snc_criteria_age_site_celltype"
FILTERED_OUTDIR = OUTDIR / "filtered_exclude_face_inguinal_antfossa"

MIN_CELLS_PER_SAMPLE_CELLTYPE = 20
RANDOM_SEED = 260811
FILTER_EXCLUDE_SITES = {"face", "inguinal fold", "antecubital fossa"}

METHODS = {
    "snc_bool": "Original marker-positive",
    "snc_component_Q4Q4Q1_within_lineage": "Growth Q4 + SASP Q4 + Cycle Q1",
}

OBS_COLS = [
    "sample_for_summary",
    "donor_id_for_summary",
    "age_decade",
    "site_for_summary",
    "broad_cell_type_for_summary",
    "snc_bool",
    "snc_component_Q4Q4Q1_within_lineage",
]


def age_sort_key(label: str) -> int:
    label = str(label)
    if label.endswith("+"):
        return int(label.rstrip("+"))
    return int(label.split("-")[0])


def sem(values: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if len(values) < 2:
        return 0.0
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def coerce_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.astype("string").str.lower().isin({"true", "1", "yes", "snc+", "q4"})


def load_obs() -> pd.DataFrame:
    obs = pd.read_csv(INPUT_OBS, usecols=OBS_COLS)
    for col in METHODS:
        obs[col] = coerce_bool(obs[col])
    for col in ["sample_for_summary", "donor_id_for_summary", "age_decade", "site_for_summary", "broad_cell_type_for_summary"]:
        obs[col] = obs[col].astype("string").fillna("Unknown")
    obs["sample_unit_id"] = (
        obs["sample_for_summary"].astype(str)
        + " | "
        + obs["donor_id_for_summary"].astype(str)
        + " | "
        + obs["site_for_summary"].astype(str)
        + " | "
        + obs["age_decade"].astype(str)
    )
    return obs


def summarize_calls(obs: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    rows = []
    for call_col, method in METHODS.items():
        summary = (
            obs.groupby(group_cols, observed=True)
            .agg(n_cells=(call_col, "size"), n_snc=(call_col, "sum"))
            .reset_index()
        )
        summary["method"] = method
        summary["snc_fraction"] = summary["n_snc"] / summary["n_cells"].replace(0, np.nan)
        summary["snc_percent"] = summary["snc_fraction"] * 100
        rows.append(summary)
    return pd.concat(rows, ignore_index=True)


def summarize_age(sample_summary: pd.DataFrame, extra_cols: list[str] | None = None) -> pd.DataFrame:
    group_cols = ["method", "age_decade"] + (extra_cols or [])
    return (
        sample_summary.groupby(group_cols, observed=True)
        .agg(
            mean_snc_percent=("snc_percent", "mean"),
            median_snc_percent=("snc_percent", "median"),
            sem_snc_percent=("snc_percent", sem),
            n_samples=("sample_unit_id", "nunique"),
            total_cells=("n_cells", "sum"),
            total_snc=("n_snc", "sum"),
        )
        .reset_index()
    )


def site_palette(sites: list[str]) -> dict[str, tuple[float, float, float, float]]:
    cmap = plt.get_cmap("tab20", max(len(sites), 1))
    return {site: cmap(i % cmap.N) for i, site in enumerate(sites)}


def setup_axes_style(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#d8d8d8", linewidth=0.8, alpha=0.65)
    ax.set_axisbelow(True)


def plot_overall_by_age_site(sample_summary: pd.DataFrame, outdir: Path, title_context: str = "") -> None:
    plot_df = sample_summary.copy()
    age_order = sorted(plot_df["age_decade"].dropna().astype(str).unique(), key=age_sort_key)
    sites = sorted(plot_df["site_for_summary"].dropna().astype(str).unique())
    colors = site_palette(sites)
    rng = np.random.default_rng(RANDOM_SEED)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), sharey=True)
    for ax, method in zip(axes, METHODS.values()):
        method_df = plot_df[plot_df["method"] == method].copy()
        x_lookup = {age: i for i, age in enumerate(age_order)}
        method_df = method_df[method_df["age_decade"].astype(str).isin(x_lookup)]
        jitter = rng.normal(0, 0.055, size=len(method_df))

        for site in sites:
            site_df = method_df[method_df["site_for_summary"].astype(str) == site]
            if site_df.empty:
                continue
            x = site_df["age_decade"].astype(str).map(x_lookup).to_numpy(dtype=float)
            positions = method_df.index.get_indexer(site_df.index)
            ax.scatter(
                x + jitter[positions],
                site_df["snc_percent"],
                s=34,
                color=colors[site],
                edgecolor="white",
                linewidth=0.45,
                alpha=0.86,
                label=site,
            )

        age_summary = summarize_age(method_df)
        age_summary = age_summary[age_summary["age_decade"].astype(str).isin(x_lookup)].copy()
        age_summary["x"] = age_summary["age_decade"].astype(str).map(x_lookup)
        ax.errorbar(
            age_summary["x"],
            age_summary["mean_snc_percent"],
            yerr=age_summary["sem_snc_percent"],
            color="black",
            marker="o",
            markersize=5,
            linewidth=1.7,
            capsize=3,
            label="Age-bin mean +/- SEM",
            zorder=5,
        )

        ax.set_title(method)
        ax.set_xticks(range(len(age_order)))
        ax.set_xticklabels(age_order, rotation=45, ha="right")
        ax.set_xlabel("Donor age decade")
        setup_axes_style(ax)

    axes[0].set_ylabel("SnC+ cells per sample (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    fig.legend(
        by_label.values(),
        by_label.keys(),
        title="Anatomical site",
        bbox_to_anchor=(1.01, 0.98),
        loc="upper left",
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    title = "Sinai SnC fraction by age bin and anatomical site"
    if title_context:
        title = f"{title}: {title_context}"
    fig.suptitle(title, y=0.98, fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 0.84, 0.92])
    fig.savefig(outdir / "sinai_snc_original_vs_si_q4_by_age_site.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "sinai_snc_original_vs_si_q4_by_age_site.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_major_celltype_by_age_site(celltype_summary: pd.DataFrame, outdir: Path, title_context: str = "") -> None:
    plot_df = celltype_summary[celltype_summary["n_cells"] >= MIN_CELLS_PER_SAMPLE_CELLTYPE].copy()
    age_order = sorted(plot_df["age_decade"].dropna().astype(str).unique(), key=age_sort_key)
    sites = sorted(plot_df["site_for_summary"].dropna().astype(str).unique())
    celltypes = (
        plot_df.groupby("broad_cell_type_for_summary", observed=True)["n_cells"]
        .sum()
        .sort_values(ascending=False)
        .index.astype(str)
        .tolist()
    )
    colors = site_palette(sites)
    x_lookup = {age: i for i, age in enumerate(age_order)}
    method_offsets = {
        "Original marker-positive": -0.11,
        "Growth Q4 + SASP Q4 + Cycle Q1": 0.11,
    }

    method_markers = {
        "Original marker-positive": "o",
        "Growth Q4 + SASP Q4 + Cycle Q1": "^",
    }

    method_lines = {
        "Original marker-positive": "-",
        "Growth Q4 + SASP Q4 + Cycle Q1": "--",
    }
    rng = np.random.default_rng(RANDOM_SEED + 1)

    n_cols = 4
    n_rows = int(np.ceil(len(celltypes) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 4.1 * n_rows), sharex=True)
    axes = np.asarray(axes).reshape(-1)

    for ax, celltype in zip(axes, celltypes):
        sub = plot_df[plot_df["broad_cell_type_for_summary"].astype(str) == celltype].copy()
        for method, offset in method_offsets.items():
            method_df = sub[sub["method"] == method].copy()
            method_df = method_df[method_df["age_decade"].astype(str).isin(x_lookup)]
            jitter = rng.normal(0, 0.028, size=len(method_df))
            for site in sites:
                site_df = method_df[method_df["site_for_summary"].astype(str) == site]
                if site_df.empty:
                    continue
                positions = method_df.index.get_indexer(site_df.index)
                x = site_df["age_decade"].astype(str).map(x_lookup).to_numpy(dtype=float)
                ax.scatter(
                    x + offset + jitter[positions],
                    site_df["snc_percent"],
                    s=23,
                    color=colors[site],
                    marker=method_markers[method],
                    edgecolor="white",
                    linewidth=0.35,
                    alpha=0.72,
                )

            age_summary = summarize_age(method_df)
            if age_summary.empty:
                continue
            age_summary = age_summary[age_summary["age_decade"].astype(str).isin(x_lookup)].copy()
            age_summary["x"] = age_summary["age_decade"].astype(str).map(x_lookup)
            ax.plot(
                age_summary["x"] + offset,
                age_summary["mean_snc_percent"],
                method_lines[method],
                color="black",
                linewidth=1.55,
                marker=method_markers[method],
                markersize=4,
            )

        ax.set_title(celltype, fontsize=10, fontweight="bold")
        ax.set_xticks(range(len(age_order)))
        ax.set_xticklabels(age_order, rotation=45, ha="right")
        setup_axes_style(ax)

    for ax in axes[len(celltypes) :]:
        ax.axis("off")

    for ax in axes[::n_cols]:
        ax.set_ylabel("SnC+ cells (%)")
    for ax in axes[-n_cols:]:
        ax.set_xlabel("Donor age decade")

    site_handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=colors[site], markeredgecolor="white", markersize=7)
        for site in sites
    ]
    method_handles = [
        plt.Line2D([0], [0], marker=method_markers[m], linestyle=method_lines[m], color="black", markersize=6)
        for m in METHODS.values()
    ]
    fig.legend(
        site_handles,
        sites,
        title="Anatomical site",
        bbox_to_anchor=(1.005, 0.94),
        loc="upper left",
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    fig.legend(
        method_handles,
        list(METHODS.values()),
        title="SnC definition",
        bbox_to_anchor=(1.005, 0.36),
        loc="upper left",
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )
    title = f"Sinai SnC fraction by major cell type (sample-celltype n >= {MIN_CELLS_PER_SAMPLE_CELLTYPE})"
    if title_context:
        title = f"{title}: {title_context}"
    fig.suptitle(title, y=0.985, fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 0.84, 0.94])
    fig.savefig(outdir / "sinai_snc_original_vs_si_q4_by_major_cell_type_age_site.png", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "sinai_snc_original_vs_si_q4_by_major_cell_type_age_site.pdf", bbox_inches="tight")
    plt.close(fig)


def run_analysis(obs: pd.DataFrame, outdir: Path, title_context: str = "") -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    sample_summary = summarize_calls(
        obs,
        ["sample_unit_id", "sample_for_summary", "donor_id_for_summary", "age_decade", "site_for_summary"],
    )
    sample_summary.to_csv(outdir / "sinai_snc_criteria_sample_summary.csv", index=False)
    summarize_age(sample_summary).to_csv(outdir / "sinai_snc_criteria_age_bin_summary.csv", index=False)

    celltype_summary = summarize_calls(
        obs,
        [
            "sample_unit_id",
            "sample_for_summary",
            "donor_id_for_summary",
            "age_decade",
            "site_for_summary",
            "broad_cell_type_for_summary",
        ],
    )
    celltype_summary.to_csv(outdir / "sinai_snc_criteria_major_cell_type_sample_summary.csv", index=False)
    summarize_age(celltype_summary[celltype_summary["n_cells"] >= MIN_CELLS_PER_SAMPLE_CELLTYPE], ["broad_cell_type_for_summary"]).to_csv(
        outdir / "sinai_snc_criteria_major_cell_type_age_bin_summary.csv",
        index=False,
    )

    plot_overall_by_age_site(sample_summary, outdir, title_context)
    plot_major_celltype_by_age_site(celltype_summary, outdir, title_context)

    run_summary = pd.DataFrame(
        [
            {
                "input_obs": str(INPUT_OBS),
                "n_cells": len(obs),
                "n_samples": obs["sample_unit_id"].nunique(),
                "n_sites": obs["site_for_summary"].nunique(),
                "n_age_bins": obs["age_decade"].nunique(),
                "n_major_cell_types": obs["broad_cell_type_for_summary"].nunique(),
                "min_cells_per_sample_celltype": MIN_CELLS_PER_SAMPLE_CELLTYPE,
                "excluded_sites": "; ".join(sorted(FILTER_EXCLUDE_SITES)) if title_context else "",
            }
        ]
    )
    run_summary.to_csv(outdir / "sinai_snc_criteria_run_summary.csv", index=False)


def main() -> None:
    obs = load_obs()
    run_analysis(obs, OUTDIR)

    filtered_obs = obs[
        ~obs["site_for_summary"].astype("string").str.lower().isin(FILTER_EXCLUDE_SITES)
    ].copy()
    run_analysis(
        filtered_obs,
        FILTERED_OUTDIR,
        "excluding face, inguinal fold, and antecubital fossa",
    )
    print(f"Wrote all-sites outputs to {OUTDIR}")
    print(f"Wrote filtered outputs to {FILTERED_OUTDIR}")


if __name__ == "__main__":
    main()
