"""
Script 29 - Xenium SSSP scoring: 4 methods compared across age cohorts and samples.

Methods:
  A. mean_log1p_cp10k   -- normalize the full 541-gene matrix, then average the SSSP genes
  B. mean_rank          -- mean within-gene percentile rank across cells (non-parametric)
  C. detection_fraction -- fraction of SSSP genes with at least one raw transcript
  D. variance_weighted_mean -- log1p-CP10K weighted by variance across sample means

Outputs (in results/29_sssp_scoring_comparison/):
  sssp_scores_per_cell.csv.gz     -- long-form table with one row per (cell, method)
  per_method_per_age_summary.csv  -- per-method mean +/- sd by age cohort
  per_method_per_sample_summary.csv -- per-method mean +/- sd by sample
  per_method_per_age_sample_heatmap.png -- heatmap (rows = method, cols = sample) split by Old/Young
  per_method_individual_gene_heatmap.png -- per-cell z-scored gene expression heatmap,
                                            faceted by age, top 30 cells per sample
  per_method_old_vs_young_violin.png -- violin plots of scores by age, per method
  per_method_cohort_boxplot.png    -- boxplot of scores by sample, faceted by method
  run_summary.json                 -- run metadata
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.stats import mannwhitneyu
import seaborn as sns

warnings = __import__("warnings")
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s  %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("29_sssp_compare")

ROOT = Path(__file__).resolve().parents[3]
H5AD = ROOT / "data" / "Xenium" / "xenium_full.h5ad"
SSSP_JSON = ROOT / "code" / "shared" / "sssp.json"
OUT_DIR = ROOT / "results" / "260727" / "xenium_sssp_scoring_comparison"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def load_panel() -> list[str]:
    with open(SSSP_JSON) as f:
        panel = json.load(f)
    genes = []
    for e in panel["sssp_xenium"]:
        if isinstance(e, str):
            genes.append(e)
        else:
            genes.append(e.get("gene") or e.get("symbol") or e.get("Gene"))
    genes = [g for g in genes if g]
    log.info("SSSP panel: %d genes", len(genes))
    return genes


def load_adata() -> ad.AnnData:
    log.info("Loading h5ad: %s", H5AD)
    adata = sc.read_h5ad(H5AD)
    log.info("  shape = %s", adata.shape)
    return adata


# ---------------------------------------------------------------------------
# Scoring methods
# ---------------------------------------------------------------------------
def score_mean_log1p_cp10k(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    """Mean SSSP expression; adata was normalized using all Xenium genes."""
    a = adata[:, genes].copy()
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    return X.mean(axis=1)


def score_mean_rank(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    """Mean within-gene percentile rank across all cells."""
    a = adata[:, genes].copy()
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    n_cells = X.shape[0]
    R = np.empty_like(X, dtype=float)
    for j in range(X.shape[1]):
        R[:, j] = (stats.rankdata(X[:, j], method="average") - 1) / max(n_cells - 1, 1)
    return R.mean(axis=1)


def score_detection_fraction(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    """Fraction of SSSP genes with raw count > 0 (AUCell-style detection burden)."""
    X = adata[:, genes].layers["counts"]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return (X > 0).mean(axis=1)


def score_variance_weighted_mean(adata: ad.AnnData, genes: list[str]) -> np.ndarray:
    """Weight genes by variance across the six sample means; age labels are unused."""
    X = adata[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    sample_means = (
        pd.DataFrame(X, columns=genes)
        .assign(sample=adata.obs["sample"].astype(str).to_numpy())
        .groupby("sample", observed=True)[genes]
        .mean()
    )
    variance = sample_means.var(axis=0, ddof=1).to_numpy()
    weights = variance / variance.sum() if variance.sum() > 0 else np.repeat(1 / len(genes), len(genes))
    return X @ weights


METHODS = {
    "mean_log1p_cp10k": score_mean_log1p_cp10k,
    "mean_rank": score_mean_rank,
    "detection_fraction": score_detection_fraction,
    "variance_weighted_mean": score_variance_weighted_mean,
}


# ---------------------------------------------------------------------------
# Per-cohort / per-sample summary
# ---------------------------------------------------------------------------
def summarize(scores_long: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    scores_long has columns: cell_id, method, score, age, sample
    Returns (per_method_per_age_summary, per_method_per_sample_summary).
    """
    g_age = scores_long.groupby(["method", "age"])["score"].agg(["mean", "median", "std", "count"]).reset_index()
    g_age = g_age.rename(columns={"mean": "score_mean", "median": "score_median", "std": "score_sd", "count": "n_cells"})

    g_sample = scores_long.groupby(["method", "age", "sample"])["score"].agg(["mean", "median", "std", "count"]).reset_index()
    g_sample = g_sample.rename(columns={"mean": "score_mean", "median": "score_median", "std": "score_sd", "count": "n_cells"})

    # Old-vs-Young inference uses the six sample means, never cells as replicates.
    pvals = []
    for m in scores_long["method"].unique():
        sub = g_sample[g_sample["method"] == m]
        old = sub.loc[sub["age"] == "Old", "score_mean"].to_numpy()
        young = sub.loc[sub["age"] == "Young", "score_mean"].to_numpy()
        welch_p = stats.ttest_ind(old, young, equal_var=False).pvalue if len(old) >= 2 and len(young) >= 2 else np.nan
        mwu_p = mannwhitneyu(old, young, alternative="two-sided").pvalue if len(old) and len(young) else np.nan
        pvals.append({
            "method": m,
            "inference_unit": "sample",
            "n_old_samples": len(old),
            "n_young_samples": len(young),
            "young_mean_of_sample_means": float(np.mean(young)),
            "old_mean_of_sample_means": float(np.mean(old)),
            "old_minus_young": float(np.mean(old) - np.mean(young)),
            "welch_p_old_vs_young": float(welch_p),
            "mwu_p_old_vs_young": float(mwu_p),
        })
    pvals_df = pd.DataFrame(pvals)
    g_age = g_age.merge(pvals_df, on="method", how="left")
    return g_age, g_sample


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
def plot_per_age_sample_heatmap(per_sample_summary: pd.DataFrame, out_path: Path):
    """Heatmap of mean score per (method, sample), split by age cohort."""
    cohorts = sorted(per_sample_summary["age"].unique().tolist())
    fig, axes = plt.subplots(1, len(cohorts), figsize=(7 * len(cohorts), 1.2 + 0.55 * len(per_sample_summary["method"].unique())),
                             sharey=True)
    if len(cohorts) == 1:
        axes = [axes]

    methods = sorted(per_sample_summary["method"].unique().tolist())

    for ax, cohort in zip(axes, cohorts):
        sub = per_sample_summary[per_sample_summary["age"] == cohort]
        pivot = sub.pivot(index="method", columns="sample", values="score_mean").reindex(methods)
        sns.heatmap(pivot, annot=True, fmt=".3g", cmap="viridis", ax=ax, cbar=True,
                    cbar_kws={"label": "mean score"})
        ax.set_title(f"{cohort}  (n_cells per sample in label)")
        ax.set_xlabel("sample")
        ax.set_ylabel("method" if ax is axes[0] else "")

    fig.suptitle("Per-method SSSP score: mean by age x sample", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved per-age-sample heatmap → %s", out_path)


def plot_old_vs_young_violin(scores_long: pd.DataFrame, out_path: Path):
    methods = sorted(scores_long["method"].unique().tolist())
    fig, axes = plt.subplots(1, len(methods), figsize=(4 * len(methods), 4.5), sharey=False)
    if len(methods) == 1:
        axes = [axes]
    for ax, m in zip(axes, methods):
        sub = scores_long[scores_long["method"] == m]
        sns.violinplot(data=sub, x="age", y="score", ax=ax, inner="quartile", palette={"Old": "#c0392b", "Young": "#2c7fb8"})
        ax.set_title(m, fontsize=11)
        ax.set_xlabel("")
        ax.set_ylabel("score")
    fig.suptitle("SSSP scores by age, per method (Old vs Young)", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved Old-vs-Young violin → %s", out_path)


def plot_per_sample_boxplot(scores_long: pd.DataFrame, out_path: Path):
    methods = sorted(scores_long["method"].unique().tolist())
    n = len(methods)
    ncols = min(n, 2)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 4.5 * nrows))
    axes = np.atleast_2d(axes).reshape(nrows, ncols)
    for ax, m in zip(axes.flatten(), methods):
        sub = scores_long[scores_long["method"] == m]
        sns.boxplot(data=sub, x="sample", y="score", hue="age", ax=ax,
                    palette={"Old": "#c0392b", "Young": "#2c7fb8"})
        ax.set_title(m, fontsize=11)
        ax.set_xlabel("sample")
        ax.set_ylabel("score")
        ax.tick_params(axis="x", rotation=45)
        ax.legend(fontsize=8, title="age")
    # blank trailing axes
    for j in range(n, nrows * ncols):
        axes.flatten()[j].axis("off")
    fig.suptitle("SSSP scores per sample, per method (Old vs Young)", fontsize=13, y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved per-sample boxplot → %s", out_path)


def plot_individual_gene_heatmap(adata: ad.AnnData, genes: list[str], out_path: Path,
                                  cells_per_sample: int = 30):
    """
    Heatmap of per-cell z-scored expression (log1p CP10K) for each SSSP gene,
    faceted by sample and split by age. Top `cells_per_sample` cells (highest
    mean expression across the panel) per sample are shown.
    """
    a = adata[:, genes].copy()
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X).astype(float)
    Xz = (X - X.mean(axis=0, keepdims=True)) / (X.std(axis=0, ddof=1, keepdims=True) + 1e-12)
    df = pd.DataFrame(Xz, columns=genes, index=adata.obs_names)
    df["sample"] = adata.obs["sample"].astype(str).values
    df["age"] = adata.obs["age"].astype(str).values

    samples = sorted(df["sample"].unique().tolist())
    ncols = 3
    nrows = int(np.ceil(len(samples) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 2.0 + 2.5 * nrows))
    axes = np.atleast_2d(axes).reshape(nrows, ncols)

    import matplotlib.colors as mcolors
    cmap = mcolors.LinearSegmentedColormap.from_list("rz", ["#3b4cc0", "#f7f4f9", "#b40426"])

    for i, s in enumerate(samples):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        sub = df[df["sample"] == s]
        # take top cells per sample by mean
        top = sub[genes].mean(axis=1).nlargest(min(cells_per_sample, len(sub))).index
        sub_top = sub.loc[top, genes].T
        # color cells by age along columns
        cell_meta = sub.loc[top, ["age"]]
        col_colors = cell_meta["age"].map({"Old": "#c0392b", "Young": "#2c7fb8"})
        sns.heatmap(sub_top, cmap=cmap, center=0, ax=ax, vmin=-2, vmax=2,
                    cbar=i == 0, cbar_kws={"label": "z-score"} if i == 0 else None,
                    yticklabels=True, xticklabels=False)
        # add age color bar above heatmap
        ax.set_title(f"{s}", fontsize=10)
        for tick, color in zip(ax.get_xticklabels(), col_colors):
            tick.set_color(color)

    for j in range(len(samples), nrows * ncols):
        r, c = divmod(j, ncols)
        axes[r, c].axis("off")

    fig.suptitle(f"Per-cell expression of {len(genes)} SSSP genes (top {cells_per_sample} cells per sample, "
                 f"red=Old, blue=Young)", fontsize=12, y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved individual-gene heatmap → %s", out_path)


def plot_gene_means_by_age(adata: ad.AnnData, genes: list[str], out_path: Path):
    """Bar chart of mean log1p CP10K per gene, split by Old vs Young."""
    a = adata[:, genes].copy()
    X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X).astype(float)
    df = pd.DataFrame(X, columns=genes)
    df["age"] = adata.obs["age"].astype(str).values
    df["sample"] = adata.obs["sample"].astype(str).values

    means = df.groupby("age")[genes].mean().T
    means = means.reindex(sorted(means.columns), axis=1)
    means = means.loc[sorted(means.index)]
    fig, ax = plt.subplots(figsize=(max(10, 0.45 * len(genes) + 2), 5))
    means.plot(kind="bar", ax=ax, color={"Old": "#c0392b", "Young": "#2c7fb8"})
    ax.set_ylabel("mean log1p(CP10K) per cell")
    ax.set_title("SSSP gene expression by age (mean per cell)")
    ax.legend(title="age")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved gene-mean bar chart → %s", out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info("=" * 60)
    log.info("Script 29 - Xenium SSSP scoring: method comparison")
    log.info("=" * 60)

    adata = load_adata()
    genes = load_panel()
    adata.var_names = adata.var_names.astype(str).str.upper()
    if "sample" not in adata.obs:
        for fallback in ["sample_csv", "sample_concat"]:
            if fallback in adata.obs:
                adata.obs["sample"] = adata.obs[fallback].astype(str)
                break
        else:
            raise KeyError("Expected sample, sample_csv, or sample_concat in adata.obs")
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise RuntimeError(f"Missing SSSP genes in adata: {missing}")

    # Preserve raw counts for detection scoring, then normalize the complete
    # targeted panel before selecting SSSP genes.
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    log.info("Computing all 4 scores …")
    records = []
    for method_name, fn in METHODS.items():
        s = fn(adata, genes)
        log.info("  %-22s range = [%.4f, %.4f]", method_name, float(np.min(s)), float(np.max(s)))
        records.append(pd.DataFrame({
            "cell_id": adata.obs_names,
            "method": method_name,
            "score": s,
            "age": adata.obs["age"].astype(str).values,
            "sample": adata.obs["sample"].astype(str).values,
        }))

    scores_long = pd.concat(records, ignore_index=True)
    scores_csv = OUT_DIR / "sssp_scores_per_cell.csv.gz"
    scores_long.to_csv(scores_csv, index=False, compression="gzip")
    log.info("Saved per-cell scores → %s", scores_csv)

    per_age, per_sample = summarize(scores_long)
    per_age.to_csv(OUT_DIR / "per_method_per_age_summary.csv", index=False)
    per_sample.to_csv(OUT_DIR / "per_method_per_sample_summary.csv", index=False)
    test_columns = [
        "method", "inference_unit", "n_old_samples", "n_young_samples",
        "young_mean_of_sample_means", "old_mean_of_sample_means",
        "old_minus_young", "welch_p_old_vs_young", "mwu_p_old_vs_young",
    ]
    per_age[test_columns].drop_duplicates("method").to_csv(
        OUT_DIR / "old_young_sample_level_tests.csv", index=False
    )
    panel_x = adata[:, genes].X
    panel_x = panel_x.toarray() if hasattr(panel_x, "toarray") else np.asarray(panel_x)
    sample_gene_means = (
        pd.DataFrame(panel_x, columns=genes)
        .assign(sample=adata.obs["sample"].astype(str).to_numpy())
        .groupby("sample", observed=True)[genes]
        .mean()
    )
    gene_variance = sample_gene_means.var(axis=0, ddof=1)
    gene_weights = gene_variance / gene_variance.sum()
    pd.DataFrame({
        "gene": genes,
        "sample_mean_variance": gene_variance.reindex(genes).to_numpy(),
        "variance_weight": gene_weights.reindex(genes).to_numpy(),
    }).to_csv(OUT_DIR / "variance_weighted_mean_gene_weights.csv", index=False)
    log.info("Saved per-age / per-sample summaries")

    print("\n--- Per-method per-age summary ---")
    print(per_age.to_string(index=False))
    print("\n--- Per-method per-sample summary ---")
    print(per_sample.to_string(index=False))

    plot_per_age_sample_heatmap(per_sample, OUT_DIR / "per_method_per_age_sample_heatmap.png")
    plot_old_vs_young_violin(scores_long, OUT_DIR / "per_method_old_vs_young_violin.png")
    plot_per_sample_boxplot(scores_long, OUT_DIR / "per_method_per_sample_boxplot.png")
    plot_individual_gene_heatmap(adata, genes, OUT_DIR / "individual_gene_per_sample_heatmap.png", cells_per_sample=30)
    plot_gene_means_by_age(adata, genes, OUT_DIR / "gene_means_old_vs_young.png")

    summary = {
        "n_cells": int(adata.shape[0]),
        "n_panel_genes": len(genes),
        "panel_genes": genes,
        "methods": list(METHODS.keys()),
        "normalization": "CP10K + log1p on the complete Xenium matrix before SSSP subsetting",
        "old_young_inference_unit": "sample",
        "output_dir": str(OUT_DIR),
    }
    with open(OUT_DIR / "run_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    log.info("=" * 60)
    log.info("DONE. Outputs in %s", OUT_DIR)
    log.info("=" * 60)


if __name__ == "__main__":
    main()
