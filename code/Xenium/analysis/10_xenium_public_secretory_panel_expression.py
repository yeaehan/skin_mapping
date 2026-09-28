"""
Summarize Xenium expression of Sinai-derived age-associated secretory genes.

The public-validation screen produces a broad secretory panel, but only a subset
is measurable in the custom Xenium panel. This script extracts those measurable
genes from the raw Xenium cell_feature_matrix files and summarizes expression by
sample, age, SnC status, and hotspot class.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[3]
XENIUM_DIR = PROJECT_ROOT / "data" / "xenium"
PANEL_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "sinai_age_decade_trends_all_sites"
    / "permissive_age_associated_secretory_panel_annotated.csv"
)
HOTSPOT_CELL_EXPR_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_snc_hotspot_expression"
    / "xenium_snc_hotspot_cell_expression_long.csv.gz"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_public_secretory_panel_expression"


def age_from_sample(sample: str) -> str:
    if sample.upper().endswith("-Y"):
        return "Young"
    if sample.upper().endswith("-O"):
        return "Old"
    return "Unknown"


def decode_array(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def load_panel_genes() -> pd.DataFrame:
    panel = pd.read_csv(PANEL_PATH)
    panel["gene"] = panel["gene"].astype(str).str.upper()
    panel["present_in_xenium_panel"] = panel["present_in_xenium_panel"].astype(bool)
    selected = panel[(panel["n_study_supports"] >= 2) & panel["present_in_xenium_panel"]].copy()
    selected["in_exact_2_study_gene_set"] = selected["n_study_supports"].eq(2)
    selected["in_18_gene_xenium_subset"] = selected["in_exact_2_study_gene_set"]
    selected["in_20_gene_xenium_subset"] = True
    selected = selected.sort_values(
        ["primary_functional_group", "n_study_supports", "gene"],
        ascending=[True, False, True],
    )
    return selected


def load_hotspot_metadata() -> pd.DataFrame:
    columns = [
        "cell_id",
        "sample",
        "age",
        "hex_id",
        "snc",
        "sasp_num",
        "metric",
        "hotspot_class",
        "gi_star_z",
        "gi_star_fdr",
        "snc_prop_bin",
        "mean_sasp_num_bin",
    ]
    meta = pd.read_csv(HOTSPOT_CELL_EXPR_PATH, usecols=columns)
    meta = meta.drop_duplicates("cell_id")
    return meta


def read_cell_totals(sample: str) -> pd.DataFrame:
    cells_path = XENIUM_DIR / sample / "cells.csv.gz"
    with gzip.open(cells_path, "rt") as handle:
        cells = pd.read_csv(handle, usecols=["cell_id", "total_counts", "transcript_counts"])
    cells["cell_id"] = sample + "_" + cells["cell_id"].astype(str)
    cells = cells.rename(columns={"total_counts": "xenium_total_counts"})
    return cells


def load_sample_gene_counts(sample: str, genes: list[str]) -> pd.DataFrame:
    h5_path = XENIUM_DIR / sample / "cell_feature_matrix.h5"
    with h5py.File(h5_path, "r") as handle:
        matrix_group = handle["matrix"]
        feature_names = decode_array(matrix_group["features"]["name"][:])
        barcodes = decode_array(matrix_group["barcodes"][:])
        shape = tuple(int(x) for x in matrix_group["shape"][:])
        matrix = sparse.csc_matrix(
            (
                matrix_group["data"][:],
                matrix_group["indices"][:],
                matrix_group["indptr"][:],
            ),
            shape=shape,
        )

    feature_lookup = {gene.upper(): i for i, gene in enumerate(feature_names)}
    present_genes = [gene for gene in genes if gene in feature_lookup]
    if not present_genes:
        return pd.DataFrame()

    cell_ids = np.array([f"{sample}_{barcode}" for barcode in barcodes], dtype=object)
    rows = []
    for gene in present_genes:
        counts = matrix.getrow(feature_lookup[gene]).toarray().ravel()
        rows.append(
            pd.DataFrame(
                {
                    "cell_id": cell_ids,
                    "sample": sample,
                    "age": age_from_sample(sample),
                    "gene": gene,
                    "count": counts.astype(np.float32),
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def summarize_expression(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    return (
        df.groupby(group_cols, observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            mean_count=("count", "mean"),
            median_count=("count", "median"),
            pct_expressing=("count", lambda x: float((x > 0).mean())),
            mean_log1p_cp10k=("log1p_cp10k", "mean"),
        )
        .reset_index()
    )


def summarize_panel_scores(cell_gene: pd.DataFrame, panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = cell_gene.copy()
    flag_cols = ["in_18_gene_xenium_subset", "in_20_gene_xenium_subset"]
    if not all(col in df.columns for col in flag_cols):
        panel_flags = panel[["gene", *flag_cols]].copy()
        df = df.merge(panel_flags, on="gene", how="left")

    score_rows = []
    for label, flag_col in [
        ("18_gene_exact_2_study_subset", "in_18_gene_xenium_subset"),
        ("20_gene_at_least_2_study_subset", "in_20_gene_xenium_subset"),
    ]:
        sub = df[df[flag_col]].copy()
        scores = (
            sub.groupby(["cell_id", "sample", "age"], observed=True)
            .agg(
                panel_gene_count=("gene", "nunique"),
                panel_total_count=("count", "sum"),
                panel_mean_count=("count", "mean"),
                panel_n_positive_genes=("count", lambda x: int((x > 0).sum())),
                panel_fraction_positive_genes=("count", lambda x: float((x > 0).mean())),
                panel_mean_log1p_cp10k=("log1p_cp10k", "mean"),
            )
            .reset_index()
        )
        scores["panel_set"] = label
        score_rows.append(scores)

    cell_scores = pd.concat(score_rows, ignore_index=True)
    meta_cols = [
        "cell_id",
        "hex_id",
        "snc",
        "sasp_num",
        "metric",
        "hotspot_class",
        "gi_star_z",
        "gi_star_fdr",
        "snc_prop_bin",
        "mean_sasp_num_bin",
    ]
    meta = cell_gene[[c for c in meta_cols if c in cell_gene.columns]].drop_duplicates("cell_id")
    cell_scores = cell_scores.merge(meta, on="cell_id", how="left")

    summary = (
        cell_scores.groupby(["panel_set", "sample", "age"], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            mean_panel_total_count=("panel_total_count", "mean"),
            mean_panel_n_positive_genes=("panel_n_positive_genes", "mean"),
            mean_panel_fraction_positive_genes=("panel_fraction_positive_genes", "mean"),
            mean_panel_log1p_cp10k=("panel_mean_log1p_cp10k", "mean"),
        )
        .reset_index()
    )
    return cell_scores, summary


def write_dotplot(summary: pd.DataFrame, panel: pd.DataFrame, outpath: Path) -> None:
    genes = panel["gene"].tolist()
    samples = summary[["sample", "age"]].drop_duplicates()
    samples["age_order"] = samples["age"].map({"Young": 0, "Old": 1}).fillna(2)
    samples = samples.sort_values(["age_order", "sample"])
    sample_order = samples["sample"].tolist()

    plot_df = summary.copy()
    plot_df["gene"] = pd.Categorical(plot_df["gene"], categories=genes, ordered=True)
    plot_df["sample"] = pd.Categorical(plot_df["sample"], categories=sample_order, ordered=True)
    plot_df = plot_df.sort_values(["sample", "gene"])
    x = plot_df["gene"].cat.codes.to_numpy()
    y = plot_df["sample"].cat.codes.to_numpy()
    sizes = 30 + 450 * plot_df["pct_expressing"].to_numpy()
    colors = np.log1p(plot_df["mean_count"].to_numpy())

    fig_w = max(10, 0.45 * len(genes) + 4)
    fig_h = max(4, 0.45 * len(sample_order) + 2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    scatter = ax.scatter(x, y, s=sizes, c=colors, cmap="viridis", edgecolor="0.35", linewidth=0.25)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=45, ha="right")
    ax.set_yticks(range(len(sample_order)))
    ax.set_yticklabels(sample_order)
    ax.set_xlabel("Xenium-measurable public secretory genes")
    ax.set_ylabel("Xenium sample")
    ax.set_title("Expression of Sinai-derived age-associated secretory genes in Xenium")
    cbar = fig.colorbar(scatter, ax=ax, pad=0.01)
    cbar.set_label("log1p(mean raw count per cell)")

    for pct in [0.01, 0.05, 0.10, 0.25]:
        ax.scatter([], [], s=30 + 450 * pct, c="lightgray", edgecolor="0.35", label=f"{pct:.0%} positive")
    ax.legend(title="Cell fraction", bbox_to_anchor=(1.18, 1.0), loc="upper left", frameon=False)
    ax.grid(axis="both", color="0.9", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    fig.savefig(outpath.with_suffix(".pdf"))
    plt.close(fig)


def write_panel_score_plot(summary: pd.DataFrame, outpath: Path) -> None:
    summary = summary.copy()
    summary["age_order"] = summary["age"].map({"Young": 0, "Old": 1}).fillna(2)
    summary = summary.sort_values(["panel_set", "age_order", "sample"])
    panel_sets = summary["panel_set"].drop_duplicates().tolist()

    fig, axes = plt.subplots(1, len(panel_sets), figsize=(12, 4), sharey=True)
    if len(panel_sets) == 1:
        axes = [axes]
    colors = {"Young": "#4C78A8", "Old": "#F58518", "Unknown": "#9E9E9E"}
    for ax, panel_set in zip(axes, panel_sets):
        sub = summary[summary["panel_set"].eq(panel_set)].copy()
        x = np.arange(len(sub))
        ax.bar(x, sub["mean_panel_n_positive_genes"], color=[colors.get(a, "#9E9E9E") for a in sub["age"]])
        ax.set_xticks(x)
        ax.set_xticklabels(sub["sample"], rotation=45, ha="right")
        ax.set_title(panel_set.replace("_", " "))
        ax.set_ylabel("Mean positive panel genes per cell")
        ax.grid(axis="y", color="0.9", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    fig.savefig(outpath.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    panel = load_panel_genes()
    genes = panel["gene"].tolist()
    panel.to_csv(OUTDIR / "xenium_public_secretory_selected_genes.csv", index=False)

    meta = load_hotspot_metadata()
    sample_dirs = sorted(p.name for p in XENIUM_DIR.iterdir() if (p / "cell_feature_matrix.h5").exists())

    long_frames = []
    for sample in sample_dirs:
        counts = load_sample_gene_counts(sample, genes)
        totals = read_cell_totals(sample)
        counts = counts.merge(totals, on="cell_id", how="left")
        counts["cp10k"] = np.where(
            counts["xenium_total_counts"] > 0,
            counts["count"] / counts["xenium_total_counts"] * 10000,
            0,
        )
        counts["log1p_cp10k"] = np.log1p(counts["cp10k"])
        long_frames.append(counts)

    cell_gene = pd.concat(long_frames, ignore_index=True)
    cell_gene = cell_gene.merge(meta, on=["cell_id", "sample", "age"], how="left")
    cell_gene = cell_gene.merge(
        panel[
            [
                "gene",
                "n_study_supports",
                "supported_studies",
                "primary_functional_group",
                "in_18_gene_xenium_subset",
                "in_20_gene_xenium_subset",
            ]
        ],
        on="gene",
        how="left",
    )
    cell_gene.to_csv(OUTDIR / "xenium_public_secretory_cell_gene_counts.csv.gz", index=False)

    by_sample = summarize_expression(cell_gene, ["sample", "age", "gene"])
    by_sample = by_sample.merge(panel[["gene", "n_study_supports", "primary_functional_group", "in_18_gene_xenium_subset"]], on="gene", how="left")
    by_sample.to_csv(OUTDIR / "xenium_public_secretory_expression_by_sample.csv", index=False)

    by_age = summarize_expression(cell_gene, ["age", "gene"])
    by_age = by_age.merge(panel[["gene", "n_study_supports", "primary_functional_group", "in_18_gene_xenium_subset"]], on="gene", how="left")
    by_age.to_csv(OUTDIR / "xenium_public_secretory_expression_by_age.csv", index=False)

    age_sample_mean = (
        by_sample.groupby(["age", "gene"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_count_sample_mean=("mean_count", "mean"),
            pct_expressing_sample_mean=("pct_expressing", "mean"),
            mean_log1p_cp10k_sample_mean=("mean_log1p_cp10k", "mean"),
        )
        .reset_index()
    )
    age_delta = age_sample_mean.pivot(index="gene", columns="age", values="mean_count_sample_mean").reset_index()
    if "Old" in age_delta.columns and "Young" in age_delta.columns:
        age_delta["old_minus_young_mean_count"] = age_delta["Old"] - age_delta["Young"]
    age_delta = age_delta.merge(
        panel[["gene", "n_study_supports", "primary_functional_group", "in_18_gene_xenium_subset"]],
        on="gene",
        how="left",
    ).sort_values("old_minus_young_mean_count", ascending=False)
    age_sample_mean.to_csv(OUTDIR / "xenium_public_secretory_expression_by_age_sample_mean.csv", index=False)
    age_delta.to_csv(OUTDIR / "xenium_public_secretory_expression_old_young_delta_sample_mean.csv", index=False)

    by_hotspot_class = summarize_expression(
        cell_gene.dropna(subset=["hotspot_class"]),
        ["sample", "age", "hotspot_class", "gene"],
    )
    by_hotspot_class = by_hotspot_class.merge(
        panel[["gene", "n_study_supports", "primary_functional_group", "in_18_gene_xenium_subset"]],
        on="gene",
        how="left",
    )
    by_hotspot_class.to_csv(OUTDIR / "xenium_public_secretory_expression_by_hotspot_class.csv", index=False)

    hotspot_mean = (
        by_hotspot_class.groupby(["hotspot_class", "gene"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_count_sample_mean=("mean_count", "mean"),
            pct_expressing_sample_mean=("pct_expressing", "mean"),
        )
        .reset_index()
    )
    hotspot_delta = hotspot_mean.pivot(index="gene", columns="hotspot_class", values="mean_count_sample_mean").reset_index()
    if "Hot spot" in hotspot_delta.columns and "Not significant" in hotspot_delta.columns:
        hotspot_delta["hotspot_minus_not_significant_mean_count"] = (
            hotspot_delta["Hot spot"] - hotspot_delta["Not significant"]
        )
    hotspot_delta = hotspot_delta.merge(
        panel[["gene", "n_study_supports", "primary_functional_group", "in_18_gene_xenium_subset"]],
        on="gene",
        how="left",
    ).sort_values("hotspot_minus_not_significant_mean_count", ascending=False)
    hotspot_delta.to_csv(OUTDIR / "xenium_public_secretory_expression_hotspot_delta_sample_mean.csv", index=False)

    by_hotspot = summarize_expression(
        cell_gene.dropna(subset=["hotspot_class"]),
        ["sample", "age", "hotspot_class", "snc", "gene"],
    )
    by_hotspot = by_hotspot.merge(
        panel[["gene", "n_study_supports", "primary_functional_group", "in_18_gene_xenium_subset"]],
        on="gene",
        how="left",
    )
    by_hotspot.to_csv(OUTDIR / "xenium_public_secretory_expression_by_hotspot_snc.csv", index=False)

    cell_scores, score_summary = summarize_panel_scores(cell_gene, panel)
    cell_scores.to_csv(OUTDIR / "xenium_public_secretory_cell_panel_scores.csv.gz", index=False)
    score_summary.to_csv(OUTDIR / "xenium_public_secretory_panel_score_by_sample.csv", index=False)

    write_dotplot(by_sample, panel, OUTDIR / "xenium_public_secretory_gene_dotplot_by_sample.png")
    write_panel_score_plot(score_summary, OUTDIR / "xenium_public_secretory_panel_score_by_sample.png")

    print(f"Saved outputs to: {OUTDIR}")
    print(f"20-gene set: {', '.join(panel['gene'])}")
    print(f"18-gene exact-2-study subset: {', '.join(panel.loc[panel['in_18_gene_xenium_subset'], 'gene'])}")
    print(score_summary.to_string(index=False))


if __name__ == "__main__":
    main()
