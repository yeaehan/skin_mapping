"""
Score the final GTEx-derived SSSP in Xenium with scanpy.tl.score_genes.

The full SSSP is the HPA-secreted, old-predictive, positive age-trend panel
from the GTEx suprapubic + lower-leg analysis. Xenium only measures a subset
of those genes, so the score here is explicitly the Xenium-measurable SSSP.

Input h5ad is the existing normalized/log1p Xenium object from the earlier
Scanpy scoring workflow; this script adds a new SSSP score and leaves the old
20-gene outputs untouched.
"""

from __future__ import annotations

from pathlib import Path
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats


PROJECT_ROOT = Path(__file__).resolve().parents[3]
INPUT_H5AD = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_sssp_20gene_scanpy_score_genes"
    / "xenium_20gene_sssp_scanpy_scored.h5ad"
)
SSSP_PANEL_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "gtex_bulk_hpa_secreted_elastic_net_old_young_from_16_preprocessed"
    / "suprapubic_plus_lower_leg_elastic_net_stable_old_predictive_positive_spearman_panel.csv"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_gtex_sssp_scanpy_score_analysis"
SCORED_H5AD = OUTDIR / "xenium_gtex_sssp_scanpy_scored.h5ad"

SCORE_COL = "sssp_scanpy_score"
HIGH_COL = "sssp_high_q75_global"
HIGH_QUANTILE = 0.75
AGE_ORDER = ["Young", "Old"]
AGE_COLORS = {"Young": "#4C78A8", "Old": "#F58518"}
CELL_TYPE_COL = "cell_type_minor"
RANDOM_SEED = 260727


def sem(values: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if len(values) <= 1:
        return np.nan
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def sample_sort_key(sample: str) -> tuple[int, str]:
    sample = str(sample)
    age_rank = 0 if sample.upper().endswith("-Y") else 1 if sample.upper().endswith("-O") else 2
    return age_rank, sample


def safe_welch_p(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    if np.isclose(a.std(ddof=1), 0) and np.isclose(b.std(ddof=1), 0):
        return np.nan
    return float(stats.ttest_ind(a, b, equal_var=False).pvalue)


def safe_mannwhitney_p(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 1 or len(b) < 1:
        return np.nan
    return float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)


def load_sssp_panel() -> pd.DataFrame:
    panel = pd.read_csv(SSSP_PANEL_PATH)
    if "gene" not in panel.columns:
        raise KeyError(f"Expected a gene column in {SSSP_PANEL_PATH}")
    panel["gene"] = panel["gene"].astype(str).str.upper()
    return panel.drop_duplicates("gene").copy()


def read_h5ad_compat(path: Path):
    try:
        return sc.read_h5ad(path)
    except Exception:
        OUTDIR.mkdir(parents=True, exist_ok=True)
        tmp = OUTDIR / "_tmp_readable_input_without_null_log1p_base.h5ad"
        shutil.copy2(path, tmp)
        import h5py

        with h5py.File(tmp, "r+") as handle:
            if "uns" in handle and "log1p" in handle["uns"] and "base" in handle["uns/log1p"]:
                del handle["uns/log1p/base"]
        try:
            return sc.read_h5ad(tmp)
        finally:
            tmp.unlink(missing_ok=True)


def add_panel_expression_summaries(adata, genes: list[str]) -> pd.DataFrame:
    panel = adata[:, genes].X
    if sparse.issparse(panel):
        mean_expr = np.asarray(panel.mean(axis=1)).ravel()
        n_positive = np.asarray((panel > 0).sum(axis=1)).ravel()
    else:
        mean_expr = np.asarray(panel).mean(axis=1)
        n_positive = (np.asarray(panel) > 0).sum(axis=1)
    return pd.DataFrame(
        {
            "cell_id": adata.obs["cell_id"].astype(str).to_numpy(),
            "sample": adata.obs["sample"].astype(str).to_numpy(),
            "age": adata.obs["age"].astype(str).to_numpy(),
            SCORE_COL: pd.to_numeric(adata.obs[SCORE_COL], errors="coerce").to_numpy(),
            HIGH_COL: adata.obs[HIGH_COL].astype(bool).to_numpy(),
            "sssp_panel_mean_log1p_cp10k": mean_expr,
            "sssp_n_positive_genes": n_positive.astype(int),
            "sssp_fraction_positive_genes": n_positive / len(genes),
        }
    )


def add_annotations(cell_scores: pd.DataFrame, adata) -> pd.DataFrame:
    wanted = [
        "snc",
        "sasp_num",
        "cell_type_major",
        "cell_type_minor",
        "cell_type_final",
        "immune_subclusters1",
        "immcelltype",
        "fov",
        "x",
        "y",
    ]
    for col in wanted:
        if col in adata.obs.columns:
            cell_scores[col] = adata.obs[col].to_numpy()
    return cell_scores


def summarize_sample_scores(cell_scores: pd.DataFrame) -> pd.DataFrame:
    return (
        cell_scores.groupby(["sample", "age"], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            mean_sssp_score=(SCORE_COL, "mean"),
            median_sssp_score=(SCORE_COL, "median"),
            sem_sssp_score=(SCORE_COL, sem),
            fraction_sssp_high=(HIGH_COL, "mean"),
            mean_sssp_panel_log1p_cp10k=("sssp_panel_mean_log1p_cp10k", "mean"),
            mean_n_positive_sssp_genes=("sssp_n_positive_genes", "mean"),
        )
        .reset_index()
        .sort_values("sample", key=lambda s: s.map(sample_sort_key))
    )


def sample_level_tests(sample_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in [
        "mean_sssp_score",
        "median_sssp_score",
        "fraction_sssp_high",
        "mean_sssp_panel_log1p_cp10k",
        "mean_n_positive_sssp_genes",
    ]:
        young = sample_summary.loc[sample_summary["age"].eq("Young"), metric]
        old = sample_summary.loc[sample_summary["age"].eq("Old"), metric]
        rows.append(
            {
                "comparison": "Old vs Young",
                "unit": "sample",
                "metric": metric,
                "n_young_samples": int(young.notna().sum()),
                "n_old_samples": int(old.notna().sum()),
                "young_mean": float(young.mean()),
                "old_mean": float(old.mean()),
                "old_minus_young": float(old.mean() - young.mean()),
                "welch_ttest_p": safe_welch_p(old, young),
                "mann_whitney_p": safe_mannwhitney_p(old, young),
            }
        )
    return pd.DataFrame(rows)


def summarize_cell_types(cell_scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = cell_scores.dropna(subset=[CELL_TYPE_COL]).copy()
    df[CELL_TYPE_COL] = df[CELL_TYPE_COL].astype(str)

    by_sample_cell_type = (
        df.groupby(["sample", "age", CELL_TYPE_COL], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            mean_sssp_score=(SCORE_COL, "mean"),
            median_sssp_score=(SCORE_COL, "median"),
            fraction_sssp_high=(HIGH_COL, "mean"),
            n_sssp_high=(HIGH_COL, "sum"),
        )
        .reset_index()
    )

    high = df.loc[df[HIGH_COL]].copy()
    comp = (
        high.groupby(["sample", "age", CELL_TYPE_COL], observed=True)
        .agg(n_sssp_high_cells=("cell_id", "nunique"))
        .reset_index()
    )
    totals = high.groupby(["sample", "age"], observed=True)["cell_id"].nunique().rename("total_sssp_high_cells")
    comp = comp.merge(totals.reset_index(), on=["sample", "age"], how="left")
    comp["fraction_of_sssp_high_cells"] = comp["n_sssp_high_cells"] / comp["total_sssp_high_cells"]
    return by_sample_cell_type, comp


def write_score_plots(sample_summary: pd.DataFrame, cell_scores: pd.DataFrame) -> None:
    sample_order = sorted(sample_summary["sample"].astype(str).unique(), key=sample_sort_key)
    colors = [AGE_COLORS.get(str(sample_summary.set_index("sample").loc[s, "age"]), "grey") for s in sample_order]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    ax = axes[0]
    positions = np.arange(len(AGE_ORDER))
    for i, age in enumerate(AGE_ORDER):
        values = sample_summary.loc[sample_summary["age"].eq(age), "mean_sssp_score"].dropna().to_numpy()
        if len(values):
            x = np.full(len(values), positions[i], dtype=float) + np.linspace(-0.06, 0.06, len(values))
            ax.scatter(x, values, s=70, color=AGE_COLORS[age], edgecolor="black", linewidth=0.4, label=age, zorder=3)
            ax.errorbar(
                positions[i],
                values.mean(),
                yerr=values.std(ddof=1) / np.sqrt(len(values)) if len(values) > 1 else np.nan,
                fmt="_",
                color="black",
                markersize=24,
                capsize=4,
                zorder=4,
            )
    ax.set_xticks(positions, AGE_ORDER)
    ax.set_ylabel("Mean SSSP score per sample")
    ax.set_title("Sample-level score")
    ax.grid(axis="y", alpha=0.25)

    ax = axes[1]
    plot_df = sample_summary.set_index("sample").loc[sample_order].reset_index()
    ax.bar(np.arange(len(plot_df)), plot_df["mean_sssp_score"], color=colors, edgecolor="black", linewidth=0.3)
    ax.errorbar(
        np.arange(len(plot_df)),
        plot_df["mean_sssp_score"],
        yerr=plot_df["sem_sssp_score"],
        fmt="none",
        ecolor="black",
        capsize=3,
        linewidth=0.8,
    )
    ax.set_xticks(np.arange(len(plot_df)), plot_df["sample"], rotation=45, ha="right")
    ax.set_ylabel("Mean SSSP score")
    ax.set_title("Individual samples")
    ax.grid(axis="y", alpha=0.25)
    fig.suptitle("Xenium-measurable SSSP by age", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(OUTDIR / "xenium_sssp_score_young_old_and_samples.png", dpi=300)
    fig.savefig(OUTDIR / "xenium_sssp_score_young_old_and_samples.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.5, 4.3))
    grouped = [cell_scores.loc[cell_scores["sample"].eq(sample), SCORE_COL].dropna().to_numpy() for sample in sample_order]
    parts = ax.violinplot(grouped, positions=np.arange(len(sample_order)), widths=0.72, showextrema=False)
    for body, color in zip(parts["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("black")
        body.set_alpha(0.55)
        body.set_linewidth(0.4)
    medians = [np.median(vals) if len(vals) else np.nan for vals in grouped]
    ax.scatter(np.arange(len(sample_order)), medians, color="black", s=18, zorder=3)
    ax.set_xticks(np.arange(len(sample_order)), sample_order, rotation=45, ha="right")
    ax.set_ylabel("Cell-level SSSP score")
    ax.set_title("Cell-level distributions by sample")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUTDIR / "xenium_sssp_score_cell_level_by_sample.png", dpi=300)
    fig.savefig(OUTDIR / "xenium_sssp_score_cell_level_by_sample.pdf")
    plt.close(fig)


def top_cell_types(comp: pd.DataFrame, n: int = 12) -> list[str]:
    totals = comp.groupby(CELL_TYPE_COL)["n_sssp_high_cells"].sum().sort_values(ascending=False)
    return totals.head(n).index.astype(str).tolist()


def write_cell_type_plots(by_sample_cell_type: pd.DataFrame, comp: pd.DataFrame) -> None:
    sample_order = sorted(comp["sample"].astype(str).unique(), key=sample_sort_key)
    keep = top_cell_types(comp, 11)

    comp_plot = comp.copy()
    comp_plot["plot_cell_type"] = np.where(comp_plot[CELL_TYPE_COL].isin(keep), comp_plot[CELL_TYPE_COL], "Other")
    comp_plot = (
        comp_plot.groupby(["sample", "age", "plot_cell_type"], observed=True)
        .agg(fraction_of_sssp_high_cells=("fraction_of_sssp_high_cells", "sum"))
        .reset_index()
    )
    pivot = comp_plot.pivot(index="sample", columns="plot_cell_type", values="fraction_of_sssp_high_cells").fillna(0)
    pivot = pivot.reindex(sample_order)
    col_order = keep + (["Other"] if "Other" in pivot.columns else [])
    pivot = pivot.reindex(columns=col_order, fill_value=0)

    cmap = plt.get_cmap("tab20")
    colors = [cmap(i % cmap.N) for i in range(len(pivot.columns))]
    fig, ax = plt.subplots(figsize=(10.8, 5.0))
    bottom = np.zeros(len(pivot))
    for color, col in zip(colors, pivot.columns):
        vals = pivot[col].to_numpy(dtype=float)
        ax.bar(np.arange(len(pivot)), vals, bottom=bottom, label=col, color=color, edgecolor="white", linewidth=0.25)
        bottom += vals
    ax.axvline(sum(str(s).endswith("-Y") for s in sample_order) - 0.5, color="0.5", linestyle="--", linewidth=1)
    ax.set_xticks(np.arange(len(pivot)), pivot.index, rotation=45, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("Fraction of SSSP-high cells")
    ax.set_title("Composition of globally top-quartile SSSP-high cells")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "xenium_sssp_high_composition_by_sample_cell_type_minor.png", dpi=300)
    fig.savefig(OUTDIR / "xenium_sssp_high_composition_by_sample_cell_type_minor.pdf")
    plt.close(fig)

    frac = by_sample_cell_type.copy()
    totals = frac.groupby(CELL_TYPE_COL)["n_cells"].sum()
    keep_frac = totals.sort_values(ascending=False).head(14).index.tolist()
    frac = frac.loc[frac[CELL_TYPE_COL].isin(keep_frac)].copy()
    age_summary = (
        frac.groupby(["age", CELL_TYPE_COL], observed=True)
        .agg(
            mean_fraction_sssp_high=("fraction_sssp_high", "mean"),
            sem_fraction_sssp_high=("fraction_sssp_high", sem),
            n_samples=("sample", "nunique"),
        )
        .reset_index()
    )
    order = (
        age_summary.groupby(CELL_TYPE_COL)["mean_fraction_sssp_high"]
        .mean()
        .sort_values(ascending=False)
        .index.tolist()
    )

    fig, ax = plt.subplots(figsize=(10.8, 5.2))
    width = 0.36
    x = np.arange(len(order))
    for i, age in enumerate(AGE_ORDER):
        sub = age_summary.loc[age_summary["age"].eq(age)].set_index(CELL_TYPE_COL).reindex(order)
        offset = (i - 0.5) * width
        ax.bar(
            x + offset,
            sub["mean_fraction_sssp_high"],
            width=width,
            yerr=sub["sem_fraction_sssp_high"],
            label=age,
            color=AGE_COLORS[age],
            edgecolor="black",
            linewidth=0.3,
            capsize=2,
        )
    ax.set_xticks(x, order, rotation=50, ha="right")
    ax.set_ylabel("Mean fraction SSSP-high within cell type")
    ax.set_title("SSSP-high enrichment by minor cell type")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUTDIR / "xenium_sssp_high_fraction_by_cell_type_minor_age.png", dpi=300)
    fig.savefig(OUTDIR / "xenium_sssp_high_fraction_by_cell_type_minor_age.pdf")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    panel = load_sssp_panel()
    adata = read_h5ad_compat(INPUT_H5AD)
    adata.var_names = adata.var_names.astype(str).str.upper()

    full_genes = panel["gene"].tolist()
    present_genes = [gene for gene in full_genes if gene in adata.var_names]
    missing_genes = sorted(set(full_genes) - set(present_genes))
    if len(present_genes) < 5:
        raise RuntimeError(f"Only {len(present_genes)} SSSP genes are present in Xenium; aborting score.")

    sc.tl.score_genes(
        adata,
        gene_list=present_genes,
        score_name=SCORE_COL,
        ctrl_size=50,
        n_bins=25,
        random_state=RANDOM_SEED,
        use_raw=False,
    )
    cutoff = float(pd.to_numeric(adata.obs[SCORE_COL], errors="coerce").quantile(HIGH_QUANTILE))
    adata.obs[HIGH_COL] = pd.to_numeric(adata.obs[SCORE_COL], errors="coerce") >= cutoff
    adata.var["in_gtex_sssp_xenium_score"] = adata.var_names.isin(present_genes)
    adata.uns["gtex_sssp_scanpy"] = {
        "score_method": "scanpy.tl.score_genes",
        "score_name": SCORE_COL,
        "high_column": HIGH_COL,
        "high_definition": f"global top {int((1 - HIGH_QUANTILE) * 100)}% of cells by {SCORE_COL}",
        "high_quantile": HIGH_QUANTILE,
        "high_cutoff": cutoff,
        "n_full_sssp_genes": len(full_genes),
        "n_xenium_present_sssp_genes": len(present_genes),
        "gene_list": present_genes,
        "full_panel_source": str(SSSP_PANEL_PATH),
        "input_h5ad": str(INPUT_H5AD),
    }

    used = panel.copy()
    used["in_xenium_score"] = used["gene"].isin(present_genes)
    used.to_csv(OUTDIR / "xenium_gtex_sssp_gene_panel_used.csv", index=False)
    pd.Series(missing_genes, name="gene").to_csv(OUTDIR / "xenium_gtex_sssp_genes_missing_from_xenium.csv", index=False)

    cell_scores = add_panel_expression_summaries(adata, present_genes)
    cell_scores = add_annotations(cell_scores, adata)
    for col in ["sssp_panel_mean_log1p_cp10k", "sssp_n_positive_genes", "sssp_fraction_positive_genes"]:
        adata.obs[col] = cell_scores[col].to_numpy()
    sample_summary = summarize_sample_scores(cell_scores)
    tests = sample_level_tests(sample_summary)
    by_sample_cell_type, comp = summarize_cell_types(cell_scores)

    cell_scores.to_csv(OUTDIR / "xenium_gtex_sssp_scanpy_cell_scores.csv.gz", index=False)
    sample_summary.to_csv(OUTDIR / "xenium_gtex_sssp_scanpy_score_by_sample.csv", index=False)
    tests.to_csv(OUTDIR / "xenium_gtex_sssp_young_old_sample_level_tests.csv", index=False)
    by_sample_cell_type.to_csv(OUTDIR / "xenium_gtex_sssp_score_by_sample_cell_type_minor.csv", index=False)
    comp.to_csv(OUTDIR / "xenium_gtex_sssp_high_composition_by_sample_cell_type_minor.csv", index=False)
    adata.write_h5ad(SCORED_H5AD, compression="gzip")

    write_score_plots(sample_summary, cell_scores)
    write_cell_type_plots(by_sample_cell_type, comp)

    print(f"SSSP full panel genes: {len(full_genes)}")
    print(f"Xenium-measurable SSSP genes scored: {len(present_genes)}")
    print(f"Global SSSP-high cutoff ({HIGH_QUANTILE:.2f} quantile): {cutoff:.4f}")
    print(f"Done. Results written to {OUTDIR}")


if __name__ == "__main__":
    main()
