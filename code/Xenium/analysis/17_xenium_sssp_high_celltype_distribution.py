"""
Explore cell type composition of Xenium cells with Scanpy SSSP score > 0.

Input comes from 16_xenium_sssp_20gene_scanpy_score_genes.py, where cell IDs
and annotations are RDS-derived and raw counts are loaded from Xenium H5 files.

SSSP-high is an exploratory threshold:
    scanpy.tl.score_genes score > 0
meaning the SSSP genes score above matched background/control genes.
"""

from __future__ import annotations

from pathlib import Path
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


PROJECT_ROOT = Path(__file__).resolve().parents[3]
INPUT_CELL_SCORES = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_sssp_20gene_scanpy_score_genes"
    / "xenium_20gene_sssp_scanpy_cell_scores.csv.gz"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_sssp_high_celltype_distribution"
SCORE_COL = "sssp_20gene_scanpy_score"
HIGH_COL = "sssp_high_gt0"
GENE_SET_CONFIGS = [
    {
        "gene_set_id": "main_20",
        "input_cell_scores": INPUT_CELL_SCORES,
        "outdir": OUTDIR,
        "score_col": "sssp_20gene_scanpy_score",
    },
    {
        "gene_set_id": "diff_filtered_17",
        "input_cell_scores": PROJECT_ROOT
        / "results"
        / "260727"
        / "xenium_sssp_17gene_diff_filtered_scanpy_score_genes"
        / "xenium_17gene_diff_filtered_sssp_scanpy_cell_scores.csv.gz",
        "outdir": PROJECT_ROOT / "results" / "260727" / "xenium_sssp_17gene_diff_filtered_high_celltype_distribution",
        "score_col": "sssp_17gene_diff_filtered_scanpy_score",
    },
]
AGE_ORDER = ["Young", "Old"]
AGE_COLORS = {"Young": "#4C78A8", "Old": "#F58518"}
ANNOTATION_COLS = ["cell_type_minor"]
TOP_N = 12


PALETTE = [
    "#4C78A8",
    "#F58518",
    "#54A24B",
    "#E45756",
    "#72B7B2",
    "#B279A2",
    "#FF9DA6",
    "#9D755D",
    "#BAB0AC",
    "#EECA3B",
    "#499894",
    "#D37295",
    "#8CD17D",
]


def sem(values: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if len(values) <= 1:
        return np.nan
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def safe_welch_p(old: pd.Series, young: pd.Series) -> float:
    old = pd.to_numeric(old, errors="coerce").dropna()
    young = pd.to_numeric(young, errors="coerce").dropna()
    if len(old) < 2 or len(young) < 2:
        return np.nan
    if np.isclose(old.std(ddof=1), 0) and np.isclose(young.std(ddof=1), 0):
        return np.nan
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Precision loss occurred.*", category=RuntimeWarning)
        return float(stats.ttest_ind(old, young, equal_var=False).pvalue)


def sample_sort_key(sample: str) -> tuple[int, str]:
    age_rank = 0 if str(sample).upper().endswith("-Y") else 1 if str(sample).upper().endswith("-O") else 2
    return age_rank, str(sample)


def clean_category(series: pd.Series) -> pd.Series:
    cleaned = series.astype("string").fillna("Unknown")
    cleaned = cleaned.replace({"": "Unknown", "nan": "Unknown", "None": "Unknown"})
    return cleaned.astype(str)


def choose_categories(df: pd.DataFrame, col: str) -> list[str]:
    counts = df.loc[df[HIGH_COL], col].value_counts()
    categories = counts.head(TOP_N).index.tolist()
    remaining = sorted(set(df[col].dropna()) - set(categories))
    if remaining:
        categories.append("Other")
    return categories


def collapse_categories(df: pd.DataFrame, col: str, categories: list[str]) -> pd.Series:
    top = set(categories) - {"Other"}
    collapsed = df[col].where(df[col].isin(top), "Other")
    return pd.Categorical(collapsed, categories=categories, ordered=True)


def complete_sample_category_table(
    table: pd.DataFrame,
    sample_info: pd.DataFrame,
    categories: list[str],
    value_cols: list[str],
) -> pd.DataFrame:
    scaffold = pd.MultiIndex.from_product(
        [sample_info["sample"].tolist(), categories],
        names=["sample", "category"],
    ).to_frame(index=False)
    out = scaffold.merge(sample_info, on="sample", how="left").merge(table, on=["sample", "age", "category"], how="left")
    for col in value_cols:
        out[col] = out[col].fillna(0)
    return out


def write_overall_high_fraction(df: pd.DataFrame, sample_order: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_summary = (
        df.groupby(["sample", "age"], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            n_sssp_high=(HIGH_COL, "sum"),
            fraction_sssp_high=(HIGH_COL, "mean"),
            mean_score=(SCORE_COL, "mean"),
            median_score=(SCORE_COL, "median"),
        )
        .reset_index()
    )
    sample_summary["sample"] = pd.Categorical(sample_summary["sample"], categories=sample_order, ordered=True)
    sample_summary = sample_summary.sort_values("sample")

    age_summary = (
        sample_summary.groupby("age", observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_fraction_sssp_high=("fraction_sssp_high", "mean"),
            sem_fraction_sssp_high=("fraction_sssp_high", sem),
            mean_score=("mean_score", "mean"),
            sem_score=("mean_score", sem),
        )
        .reset_index()
    )

    young = sample_summary.loc[sample_summary["age"].eq("Young"), "fraction_sssp_high"]
    old = sample_summary.loc[sample_summary["age"].eq("Old"), "fraction_sssp_high"]
    age_summary["welch_p_old_vs_young_fraction"] = safe_welch_p(old, young)

    sample_summary.to_csv(OUTDIR / "sssp_high_gt0_fraction_by_sample.csv", index=False)
    age_summary.to_csv(OUTDIR / "sssp_high_gt0_fraction_by_age_sample_level.csv", index=False)
    return sample_summary, age_summary


def plot_overall_high_fraction(sample_summary: pd.DataFrame, age_summary: pd.DataFrame, sample_order: list[str]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2), gridspec_kw={"width_ratios": [0.8, 1.6]})

    ax = axes[0]
    for i, age in enumerate(AGE_ORDER):
        row = age_summary[age_summary["age"].eq(age)]
        if row.empty:
            continue
        row = row.iloc[0]
        ax.bar(
            i,
            row["mean_fraction_sssp_high"],
            yerr=row["sem_fraction_sssp_high"],
            color=AGE_COLORS[age],
            edgecolor="#333333",
            linewidth=0.8,
            alpha=0.85,
            width=0.58,
        )
        vals = sample_summary.loc[sample_summary["age"].eq(age), "fraction_sssp_high"].to_numpy()
        jitter = np.linspace(-0.08, 0.08, len(vals)) if len(vals) > 1 else np.array([0.0])
        ax.scatter(i + jitter, vals, color="#222222", s=28, zorder=3)
    ax.set_xticks(range(len(AGE_ORDER)))
    ax.set_xticklabels(AGE_ORDER)
    ax.set_ylabel("Fraction of all cells with SSSP score > 0")
    ax.set_title("Sample-level age summary", fontsize=12, fontweight="bold")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    sample_summary = sample_summary.set_index("sample").loc[sample_order].reset_index()
    colors = [AGE_COLORS.get(age, "#888888") for age in sample_summary["age"]]
    x = np.arange(len(sample_summary))
    ax.bar(x, sample_summary["fraction_sssp_high"], color=colors, edgecolor="#333333", linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(sample_summary["sample"], rotation=35, ha="right")
    ax.set_ylabel("Fraction of all cells with SSSP score > 0")
    ax.set_title("Individual samples", fontsize=12, fontweight="bold")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(OUTDIR / "sssp_high_gt0_fraction_by_age_and_sample.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "sssp_high_gt0_fraction_by_age_and_sample.pdf", bbox_inches="tight")
    plt.close(fig)


def summarize_annotation(df: pd.DataFrame, col: str, sample_order: list[str]) -> None:
    work = df.copy()
    work[col] = clean_category(work[col])
    categories = choose_categories(work, col)
    work["category"] = collapse_categories(work, col, categories)
    sample_info = work[["sample", "age"]].drop_duplicates().sort_values("sample", key=lambda s: s.map({x: i for i, x in enumerate(sample_order)}))

    all_comp = (
        work.groupby(["sample", "age", "category"], observed=True)
        .agg(n_all_type=("cell_id", "nunique"))
        .reset_index()
    )
    total_all = work.groupby(["sample", "age"], observed=True).size().reset_index(name="n_all_cells")
    all_comp = all_comp.merge(total_all, on=["sample", "age"], how="left")
    all_comp["fraction_of_all_cells"] = all_comp["n_all_type"] / all_comp["n_all_cells"]

    high = work[work[HIGH_COL]].copy()
    high_comp = (
        high.groupby(["sample", "age", "category"], observed=True)
        .agg(n_sssp_high_type=("cell_id", "nunique"))
        .reset_index()
    )
    total_high = high.groupby(["sample", "age"], observed=True).size().reset_index(name="n_sssp_high_cells")
    high_comp = high_comp.merge(total_high, on=["sample", "age"], how="left")
    high_comp["fraction_of_sssp_high_cells"] = high_comp["n_sssp_high_type"] / high_comp["n_sssp_high_cells"]
    high_comp = complete_sample_category_table(
        high_comp,
        sample_info,
        categories,
        ["n_sssp_high_type", "n_sssp_high_cells", "fraction_of_sssp_high_cells"],
    )

    within = (
        work.groupby(["sample", "age", "category"], observed=True)
        .agg(
            n_all_type=("cell_id", "nunique"),
            n_sssp_high_type=(HIGH_COL, "sum"),
            mean_score=(SCORE_COL, "mean"),
        )
        .reset_index()
    )
    within["fraction_sssp_high_within_type"] = within["n_sssp_high_type"] / within["n_all_type"]
    within = complete_sample_category_table(
        within,
        sample_info,
        categories,
        ["n_all_type", "n_sssp_high_type", "mean_score", "fraction_sssp_high_within_type"],
    )

    enrichment = high_comp.merge(
        all_comp[["sample", "age", "category", "fraction_of_all_cells"]],
        on=["sample", "age", "category"],
        how="left",
    )
    eps = 1e-4
    enrichment["fraction_of_all_cells"] = enrichment["fraction_of_all_cells"].fillna(0)
    enrichment["log2_high_vs_all_composition"] = np.log2(
        (enrichment["fraction_of_sssp_high_cells"] + eps) / (enrichment["fraction_of_all_cells"] + eps)
    )

    age_comp = (
        high_comp.groupby(["age", "category"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_fraction_of_sssp_high_cells=("fraction_of_sssp_high_cells", "mean"),
            sem_fraction_of_sssp_high_cells=("fraction_of_sssp_high_cells", sem),
        )
        .reset_index()
    )
    within_age = (
        within.groupby(["age", "category"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_fraction_sssp_high_within_type=("fraction_sssp_high_within_type", "mean"),
            sem_fraction_sssp_high_within_type=("fraction_sssp_high_within_type", sem),
            mean_score=("mean_score", "mean"),
            sem_score=("mean_score", sem),
        )
        .reset_index()
    )

    comp_tests = sample_level_tests(high_comp, "fraction_of_sssp_high_cells")
    within_tests = sample_level_tests(within, "fraction_sssp_high_within_type")

    safe_col = col.replace(" ", "_")
    high_comp.to_csv(OUTDIR / f"sssp_high_gt0_composition_by_sample_{safe_col}.csv", index=False)
    all_comp.to_csv(OUTDIR / f"all_cell_composition_by_sample_{safe_col}.csv", index=False)
    within.to_csv(OUTDIR / f"sssp_high_gt0_within_celltype_fraction_by_sample_{safe_col}.csv", index=False)
    enrichment.to_csv(OUTDIR / f"sssp_high_gt0_composition_enrichment_by_sample_{safe_col}.csv", index=False)
    age_comp.to_csv(OUTDIR / f"sssp_high_gt0_composition_by_age_sample_means_{safe_col}.csv", index=False)
    within_age.to_csv(OUTDIR / f"sssp_high_gt0_within_celltype_fraction_by_age_sample_means_{safe_col}.csv", index=False)
    comp_tests.to_csv(OUTDIR / f"sssp_high_gt0_composition_old_vs_young_tests_{safe_col}.csv", index=False)
    within_tests.to_csv(OUTDIR / f"sssp_high_gt0_within_celltype_fraction_old_vs_young_tests_{safe_col}.csv", index=False)
    score_summary = write_mean_score_outputs(work, col, categories, sample_order)

    plot_stacked_by_sample(
        high_comp,
        sample_order,
        categories,
        "fraction_of_sssp_high_cells",
        f"Composition of SSSP-high cells by sample: {col}",
        f"Fraction of SSSP-high cells",
        OUTDIR / f"sssp_high_gt0_composition_by_sample_{safe_col}",
    )
    plot_stacked_by_age(
        age_comp,
        categories,
        "mean_fraction_of_sssp_high_cells",
        f"Mean SSSP-high composition by age: {col}",
        f"Mean fraction of SSSP-high cells",
        OUTDIR / f"sssp_high_gt0_composition_by_age_{safe_col}",
    )
    plot_heatmap(
        within,
        sample_order,
        categories,
        "fraction_sssp_high_within_type",
        f"Within-cell-type SSSP-high fraction: {col}",
        OUTDIR / f"sssp_high_gt0_within_celltype_fraction_heatmap_{safe_col}",
        cmap="YlOrRd",
        vmin=0,
        vmax=min(1.0, max(0.05, within["fraction_sssp_high_within_type"].quantile(0.98))),
    )
    plot_heatmap(
        within,
        sample_order,
        categories,
        "fraction_sssp_high_within_type",
        f"Within-cell-type SSSP-high fraction, row-scaled: {col}",
        OUTDIR / f"sssp_high_gt0_within_celltype_fraction_heatmap_row_scaled_{safe_col}",
        cmap="coolwarm",
        vmin=-2,
        vmax=2,
        row_scale=True,
    )
    plot_heatmap(
        score_summary,
        sample_order,
        categories,
        "mean_sssp_score",
        f"Mean SSSP score by sample: {col}",
        OUTDIR / f"sssp_mean_score_heatmap_by_sample_{safe_col}",
        cmap="coolwarm",
        vmin=-score_summary["mean_sssp_score"].abs().quantile(0.98),
        vmax=score_summary["mean_sssp_score"].abs().quantile(0.98),
    )
    plot_heatmap(
        score_summary,
        sample_order,
        categories,
        "mean_sssp_score",
        f"Mean SSSP score by sample, row-scaled: {col}",
        OUTDIR / f"sssp_mean_score_heatmap_row_scaled_by_sample_{safe_col}",
        cmap="coolwarm",
        vmin=-2,
        vmax=2,
        row_scale=True,
    )
    plot_score_dotplot(
        score_summary,
        sample_order,
        categories,
        f"Mean SSSP score by cell type and sample: {col}",
        OUTDIR / f"sssp_mean_score_dotplot_by_sample_{safe_col}",
    )


def write_mean_score_outputs(
    work: pd.DataFrame,
    annotation_col: str,
    categories: list[str],
    sample_order: list[str],
) -> pd.DataFrame:
    sample_info = work[["sample", "age"]].drop_duplicates()
    summary = (
        work.groupby(["sample", "age", "category"], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            mean_sssp_score=(SCORE_COL, "mean"),
            median_sssp_score=(SCORE_COL, "median"),
            sem_sssp_score=(SCORE_COL, sem),
            fraction_sssp_high=(HIGH_COL, "mean"),
        )
        .reset_index()
    )
    summary = complete_sample_category_table(
        summary,
        sample_info,
        categories,
        ["n_cells"],
    )
    summary["sample"] = pd.Categorical(summary["sample"], categories=sample_order, ordered=True)
    summary["category"] = pd.Categorical(summary["category"], categories=categories, ordered=True)
    summary = summary.sort_values(["sample", "category"])

    age_summary = (
        summary.groupby(["age", "category"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_of_sample_mean_scores=("mean_sssp_score", "mean"),
            sem_of_sample_mean_scores=("mean_sssp_score", sem),
            mean_fraction_sssp_high=("fraction_sssp_high", "mean"),
            sem_fraction_sssp_high=("fraction_sssp_high", sem),
        )
        .reset_index()
    )
    tests = sample_level_tests(summary, "mean_sssp_score")

    safe_col = annotation_col.replace(" ", "_")
    summary.to_csv(OUTDIR / f"sssp_mean_score_by_sample_{safe_col}.csv", index=False)
    age_summary.to_csv(OUTDIR / f"sssp_mean_score_by_age_sample_means_{safe_col}.csv", index=False)
    tests.to_csv(OUTDIR / f"sssp_mean_score_old_vs_young_tests_{safe_col}.csv", index=False)
    return summary
    plot_heatmap(
        enrichment,
        sample_order,
        categories,
        "log2_high_vs_all_composition",
        f"SSSP-high composition enrichment vs all cells: {col}",
        OUTDIR / f"sssp_high_gt0_composition_enrichment_heatmap_{safe_col}",
        cmap="coolwarm",
        vmin=-2,
        vmax=2,
    )
    plot_heatmap(
        enrichment,
        sample_order,
        categories,
        "log2_high_vs_all_composition",
        f"SSSP-high composition enrichment, row-scaled: {col}",
        OUTDIR / f"sssp_high_gt0_composition_enrichment_heatmap_row_scaled_{safe_col}",
        cmap="coolwarm",
        vmin=-2,
        vmax=2,
        row_scale=True,
    )


def sample_level_tests(table: pd.DataFrame, value_col: str) -> pd.DataFrame:
    rows = []
    for category, sub in table.groupby("category", observed=True):
        young = sub.loc[sub["age"].eq("Young"), value_col].dropna()
        old = sub.loc[sub["age"].eq("Old"), value_col].dropna()
        rows.append(
            {
                "category": category,
                "n_young_samples": int(len(young)),
                "n_old_samples": int(len(old)),
                "young_mean": float(young.mean()) if len(young) else np.nan,
                "old_mean": float(old.mean()) if len(old) else np.nan,
                "old_minus_young": float(old.mean() - young.mean()) if len(old) and len(young) else np.nan,
                "welch_p": safe_welch_p(old, young),
            }
        )
    return pd.DataFrame(rows).sort_values("old_minus_young", ascending=False)


def category_colors(categories: list[str]) -> dict[str, str]:
    return {cat: PALETTE[i % len(PALETTE)] for i, cat in enumerate(categories)}


def plot_stacked_by_sample(
    table: pd.DataFrame,
    sample_order: list[str],
    categories: list[str],
    value_col: str,
    title: str,
    ylabel: str,
    outbase: Path,
) -> None:
    colors = category_colors(categories)
    pivot = table.pivot_table(index="sample", columns="category", values=value_col, aggfunc="sum", fill_value=0)
    pivot = pivot.reindex(sample_order).reindex(columns=categories, fill_value=0)

    fig, ax = plt.subplots(figsize=(9.4, 4.8))
    bottom = np.zeros(len(pivot))
    x = np.arange(len(pivot))
    for cat in categories:
        vals = pivot[cat].to_numpy(dtype=float)
        ax.bar(x, vals, bottom=bottom, color=colors[cat], edgecolor="white", linewidth=0.35, label=cat)
        bottom += vals

    young_n = sum(str(sample).upper().endswith("-Y") for sample in sample_order)
    if 0 < young_n < len(sample_order):
        ax.axvline(young_n - 0.5, color="#777777", linewidth=0.9, linestyle="--", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(sample_order, rotation=35, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_ylim(0, max(1.0, float(np.nanmax(bottom)) * 1.05))
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_stacked_by_age(
    table: pd.DataFrame,
    categories: list[str],
    value_col: str,
    title: str,
    ylabel: str,
    outbase: Path,
) -> None:
    colors = category_colors(categories)
    pivot = table.pivot_table(index="age", columns="category", values=value_col, aggfunc="sum", fill_value=0)
    pivot = pivot.reindex(AGE_ORDER).reindex(columns=categories, fill_value=0)

    fig, ax = plt.subplots(figsize=(5.8, 4.5))
    bottom = np.zeros(len(pivot))
    x = np.arange(len(pivot))
    for cat in categories:
        vals = pivot[cat].to_numpy(dtype=float)
        ax.bar(x, vals, bottom=bottom, color=colors[cat], edgecolor="white", linewidth=0.35, label=cat)
        bottom += vals

    ax.set_xticks(x)
    ax.set_xticklabels(AGE_ORDER)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_ylim(0, max(1.0, float(np.nanmax(bottom)) * 1.05))
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(
    table: pd.DataFrame,
    sample_order: list[str],
    categories: list[str],
    value_col: str,
    title: str,
    outbase: Path,
    cmap: str,
    vmin: float,
    vmax: float,
    row_scale: bool = False,
) -> None:
    pivot = table.pivot_table(index="category", columns="sample", values=value_col, aggfunc="mean", fill_value=0, observed=False)
    pivot = pivot.reindex(categories).reindex(columns=sample_order, fill_value=0)
    values = pivot.to_numpy(dtype=float)
    if row_scale:
        row_mean = np.nanmean(values, axis=1, keepdims=True)
        row_sd = np.nanstd(values, axis=1, keepdims=True)
        row_sd[row_sd == 0] = 1
        values = (values - row_mean) / row_sd

    fig_height = max(4.2, 0.34 * len(categories) + 1.5)
    fig, ax = plt.subplots(figsize=(8.8, fig_height))
    im = ax.imshow(values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(np.arange(len(sample_order)))
    ax.set_xticklabels(sample_order, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(categories)))
    ax.set_yticklabels(categories)
    ax.set_title(title, fontsize=12, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    if row_scale:
        cbar.set_label("Row z-score", fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_score_dotplot(
    table: pd.DataFrame,
    sample_order: list[str],
    categories: list[str],
    title: str,
    outbase: Path,
) -> None:
    work = table.copy()
    work["sample"] = pd.Categorical(work["sample"], categories=sample_order, ordered=True)
    work["category"] = pd.Categorical(work["category"], categories=categories, ordered=True)
    work = work.dropna(subset=["mean_sssp_score"])

    if work.empty:
        return

    max_abs = float(work["mean_sssp_score"].abs().quantile(0.98))
    if not np.isfinite(max_abs) or max_abs == 0:
        max_abs = 1.0
    max_n = float(work["n_cells"].max())
    size = 30 + 220 * np.sqrt(work["n_cells"].to_numpy(dtype=float) / max(max_n, 1.0))

    x = work["sample"].cat.codes.to_numpy()
    y = work["category"].cat.codes.to_numpy()

    fig_height = max(4.4, 0.34 * len(categories) + 1.5)
    fig, ax = plt.subplots(figsize=(9.2, fig_height))
    scatter = ax.scatter(
        x,
        y,
        c=work["mean_sssp_score"].to_numpy(dtype=float),
        s=size,
        cmap="coolwarm",
        vmin=-max_abs,
        vmax=max_abs,
        edgecolor="#333333",
        linewidth=0.35,
        alpha=0.9,
    )
    young_n = sum(str(sample).upper().endswith("-Y") for sample in sample_order)
    if 0 < young_n < len(sample_order):
        ax.axvline(young_n - 0.5, color="#777777", linewidth=0.9, linestyle="--", alpha=0.8)
    ax.set_xticks(np.arange(len(sample_order)))
    ax.set_xticklabels(sample_order, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(categories)))
    ax.set_yticklabels(categories)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(color="#E6E6E6", linewidth=0.7)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Mean SSSP score", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    legend_sizes = np.array([0.25, 0.5, 1.0]) * max_n
    handles = [
        ax.scatter([], [], s=30 + 220 * np.sqrt(n / max(max_n, 1.0)), color="#BBBBBB", edgecolor="#333333", linewidth=0.35)
        for n in legend_sizes
    ]
    labels = [f"{int(round(n))} cells" for n in legend_sizes]
    ax.legend(handles, labels, title="Cell count", frameon=False, bbox_to_anchor=(1.17, 0.45), loc="center left", fontsize=8, title_fontsize=9)

    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def run_analysis(config: dict[str, object]) -> pd.DataFrame:
    global INPUT_CELL_SCORES, OUTDIR, SCORE_COL
    INPUT_CELL_SCORES = Path(config["input_cell_scores"])
    OUTDIR = Path(config["outdir"])
    SCORE_COL = str(config["score_col"])

    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(INPUT_CELL_SCORES)
    required = {"cell_id", "sample", "age", SCORE_COL}
    missing = sorted(required - set(df.columns))
    if missing:
        raise KeyError(f"Missing required columns in {INPUT_CELL_SCORES}: {', '.join(missing)}")

    df[HIGH_COL] = df[SCORE_COL] > 0
    sample_order = sorted(df["sample"].dropna().unique().tolist(), key=sample_sort_key)

    sample_summary, age_summary = write_overall_high_fraction(df, sample_order)
    plot_overall_high_fraction(sample_summary, age_summary, sample_order)

    used_cols = []
    for col in ANNOTATION_COLS:
        if col in df.columns:
            used_cols.append(col)
            summarize_annotation(df, col, sample_order)

    run_summary = pd.DataFrame(
        [
            {
                "gene_set_id": config["gene_set_id"],
                "input_cell_scores": str(INPUT_CELL_SCORES),
                "score_column": SCORE_COL,
                "sssp_high_definition": f"{SCORE_COL} > 0",
                "n_cells": int(len(df)),
                "n_sssp_high_cells": int(df[HIGH_COL].sum()),
                "fraction_sssp_high": float(df[HIGH_COL].mean()),
                "annotation_columns_analyzed": ", ".join(used_cols),
            }
        ]
    )
    run_summary.to_csv(OUTDIR / "sssp_high_gt0_celltype_distribution_run_summary.csv", index=False)

    print("SSSP-high cell type distribution complete.")
    print(run_summary.to_string(index=False))
    print("\nOverall SSSP-high fraction by sample:")
    print(sample_summary.to_string(index=False))
    print(f"\nWrote outputs to: {OUTDIR}")
    return run_summary


def main() -> None:
    summaries = [run_analysis(config) for config in GENE_SET_CONFIGS]
    pd.concat(summaries, ignore_index=True).to_csv(
        PROJECT_ROOT / "results" / "260727" / "xenium_sssp_high_celltype_distribution_gene_set_run_summary.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
