"""
Correlate Xenium public-secretory/SASP-like scores with CDKN and proliferation markers.

This uses annotated cells from the hotspot-neighborhood analysis and computes
sample-level Spearman correlations so each sample contributes one correlation
coefficient to the final summaries.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats


PROJECT_ROOT = Path(__file__).resolve().parents[3]
XENIUM_DIR = PROJECT_ROOT / "data" / "xenium"
CELL_SCORES_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_public_secretory_panel_scores"
    / "xenium_public_secretory_cell_scores.csv.gz"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_secretory_score_marker_correlations"

MARKER_GENES = ["CDKN1A", "CDKN2A", "PCNA", "MKI67", "STMN1"]
SCORE_COLS = [
    "module_score",
    "panel_mean_log1p_cp10k",
    "auc_like_score",
    "n_positive_panel_genes",
    "sasp_num",
]
MARKER_SCORE_COLS = [
    "CDKN1A",
    "CDKN2A",
    "cdkn_mean",
    "PCNA",
    "MKI67",
    "STMN1",
    "proliferation_mean",
]


def decode_array(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def sem(values: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if len(values) <= 1:
        return np.nan
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def bh_fdr(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if valid.sum() == 0:
        return pd.Series(out, index=pvalues.index)
    valid_idx = np.where(valid)[0]
    order = valid_idx[np.argsort(p[valid])]
    ranked = p[order]
    m = len(ranked)
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out[order] = np.clip(adjusted, 0, 1)
    return pd.Series(out, index=pvalues.index)


def read_marker_log_expression(sample: str, wanted_cells: set[str]) -> pd.DataFrame:
    h5_path = XENIUM_DIR / sample / "cell_feature_matrix.h5"
    with h5py.File(h5_path, "r") as handle:
        matrix_group = handle["matrix"]
        genes = [gene.upper() for gene in decode_array(matrix_group["features"]["name"][:])]
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

    full_cell_ids = np.array([f"{sample}_{barcode}" for barcode in barcodes], dtype=object)
    keep = np.array([cell_id in wanted_cells for cell_id in full_cell_ids])
    gene_lookup = {gene: i for i, gene in enumerate(genes)}
    present = [gene for gene in MARKER_GENES if gene in gene_lookup]
    if not np.any(keep) or not present:
        return pd.DataFrame()

    marker_idx = [gene_lookup[gene] for gene in present]
    marker_counts = matrix[marker_idx, :][:, keep].astype(np.float32).toarray()
    total_counts = matrix[:, keep].sum(axis=0).A1.astype(np.float32)
    norm = np.divide(marker_counts, total_counts, out=np.zeros_like(marker_counts), where=total_counts > 0) * 10000
    log_norm = np.log1p(norm)

    out = pd.DataFrame(log_norm.T, columns=present)
    out.insert(0, "cell_id", full_cell_ids[keep])
    out.insert(1, "sample", sample)
    for gene in MARKER_GENES:
        if gene not in out.columns:
            out[gene] = np.nan
    out["cdkn_mean"] = out[["CDKN1A", "CDKN2A"]].mean(axis=1)
    out["proliferation_mean"] = out[["PCNA", "MKI67", "STMN1"]].mean(axis=1)
    return out[["cell_id", "sample", *MARKER_SCORE_COLS]]


def spearman_or_nan(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    frame = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(frame) < 20 or frame["x"].nunique() < 2 or frame["y"].nunique() < 2:
        return np.nan, np.nan, len(frame)
    rho, p = stats.spearmanr(frame["x"], frame["y"])
    return float(rho), float(p), len(frame)


def sample_correlations(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (panel_set, sample), sub in df.groupby(["panel_set", "sample"], observed=True):
        for score_col in SCORE_COLS:
            for marker_col in MARKER_SCORE_COLS:
                rho, p, n = spearman_or_nan(sub[score_col], sub[marker_col])
                rows.append(
                    {
                        "panel_set": panel_set,
                        "sample": sample,
                        "age": sub["age"].iloc[0],
                        "score": score_col,
                        "marker": marker_col,
                        "n_cells": n,
                        "spearman_rho": rho,
                        "spearman_p_cell_level_descriptive": p,
                    }
                )
    return pd.DataFrame(rows)


def summarize_sample_correlations(corr: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    summary = (
        corr.groupby(group_cols, observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_spearman_rho=("spearman_rho", "mean"),
            sem_spearman_rho=("spearman_rho", sem),
            samples=("sample", lambda x: ", ".join(sorted(x.dropna().unique()))),
        )
        .reset_index()
    )
    p_values = []
    for _, row in summary.iterrows():
        mask = np.ones(len(corr), dtype=bool)
        for col in group_cols:
            mask &= corr[col].eq(row[col]).to_numpy()
        values = corr.loc[mask, "spearman_rho"].dropna()
        p_values.append(stats.ttest_1samp(values, 0).pvalue if len(values) >= 3 else np.nan)
    summary["p_rho_vs_zero_sample_level"] = p_values
    fdr_group_cols = [col for col in group_cols if col not in {"marker"}]
    summary["fdr_rho_vs_zero_sample_level"] = summary.groupby(fdr_group_cols, observed=True)[
        "p_rho_vs_zero_sample_level"
    ].transform(bh_fdr)
    return summary.sort_values(group_cols)


def celltype_correlations(df: pd.DataFrame, min_cells: int = 50) -> pd.DataFrame:
    rows = []
    for (panel_set, sample, cell_type), sub in df.groupby(["panel_set", "sample", "cell_type_minor"], observed=True):
        if len(sub) < min_cells:
            continue
        for score_col in SCORE_COLS:
            for marker_col in MARKER_SCORE_COLS:
                rho, p, n = spearman_or_nan(sub[score_col], sub[marker_col])
                rows.append(
                    {
                        "panel_set": panel_set,
                        "sample": sample,
                        "age": sub["age"].iloc[0],
                        "cell_type_minor": cell_type,
                        "score": score_col,
                        "marker": marker_col,
                        "n_cells": n,
                        "spearman_rho": rho,
                        "spearman_p_cell_level_descriptive": p,
                    }
                )
    return pd.DataFrame(rows)


def plot_heatmap(summary: pd.DataFrame, panel_set: str, score: str, outpath: Path) -> None:
    sub = summary[summary["panel_set"].eq(panel_set) & summary["score"].eq(score)].copy()
    sub["marker"] = pd.Categorical(sub["marker"], categories=MARKER_SCORE_COLS, ordered=True)
    values = sub.pivot(index="score", columns="marker", values="mean_spearman_rho").reindex(columns=MARKER_SCORE_COLS)

    fig, ax = plt.subplots(figsize=(9, 2.5))
    im = ax.imshow(values.to_numpy(), vmin=-0.4, vmax=0.4, cmap="coolwarm", aspect="auto")
    ax.set_xticks(range(len(MARKER_SCORE_COLS)))
    ax.set_xticklabels(MARKER_SCORE_COLS, rotation=35, ha="right")
    ax.set_yticks([0])
    ax.set_yticklabels([score])
    ax.set_title(f"{panel_set.replace('_', ' ')}: sample-level mean Spearman rho")
    for j, marker in enumerate(MARKER_SCORE_COLS):
        value = values.iloc[0, j]
        if pd.notna(value):
            ax.text(j, 0, f"{value:.2f}", ha="center", va="center", fontsize=9)
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("Mean Spearman rho")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    fig.savefig(outpath.with_suffix(".pdf"))
    plt.close(fig)


def plot_celltype_bars(summary: pd.DataFrame, panel_set: str, score: str, marker: str, outpath: Path) -> None:
    sub = summary[
        summary["panel_set"].eq(panel_set) & summary["score"].eq(score) & summary["marker"].eq(marker)
    ].copy()
    sub = sub[sub["n_samples"] >= 3].sort_values("mean_spearman_rho")
    fig, ax = plt.subplots(figsize=(8, max(4, 0.3 * len(sub) + 1.5)))
    ax.barh(sub["cell_type_minor"], sub["mean_spearman_rho"], xerr=sub["sem_spearman_rho"], color="#4C78A8")
    ax.axvline(0, color="0.2", linewidth=1)
    ax.set_xlabel("Mean within-cell-type Spearman rho across samples")
    ax.set_ylabel("")
    ax.set_title(f"{score} vs {marker}, {panel_set.replace('_', ' ')}")
    ax.grid(axis="x", color="0.9", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    fig.savefig(outpath.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    cell_scores = pd.read_csv(CELL_SCORES_PATH)
    wanted_cells = set(cell_scores["cell_id"].unique())

    marker_frames = []
    for sample in sorted(cell_scores["sample"].unique()):
        marker_frames.append(read_marker_log_expression(sample, wanted_cells))
    markers = pd.concat(marker_frames, ignore_index=True)
    markers.to_csv(OUTDIR / "xenium_cdkn_proliferation_marker_log_expression.csv.gz", index=False)

    merged = cell_scores.merge(markers, on=["cell_id", "sample"], how="inner")
    merged.to_csv(OUTDIR / "xenium_secretory_scores_with_cdkn_proliferation_markers.csv.gz", index=False)

    corr = sample_correlations(merged)
    corr.to_csv(OUTDIR / "xenium_secretory_marker_correlations_by_sample.csv", index=False)
    corr_summary = summarize_sample_correlations(corr, ["panel_set", "score", "marker"])
    corr_summary.to_csv(OUTDIR / "xenium_secretory_marker_correlations_sample_summary.csv", index=False)

    ct_corr = celltype_correlations(merged)
    ct_corr.to_csv(OUTDIR / "xenium_secretory_marker_correlations_by_sample_celltype_minor.csv", index=False)
    ct_summary = summarize_sample_correlations(ct_corr, ["panel_set", "cell_type_minor", "score", "marker"])
    ct_summary.to_csv(OUTDIR / "xenium_secretory_marker_correlations_celltype_minor_summary.csv", index=False)

    for panel_set in sorted(merged["panel_set"].unique()):
        plot_heatmap(
            corr_summary,
            panel_set,
            "module_score",
            OUTDIR / f"xenium_{panel_set}_module_score_marker_correlation_heatmap.png",
        )
        plot_heatmap(
            corr_summary,
            panel_set,
            "sasp_num",
            OUTDIR / f"xenium_{panel_set}_sasp_num_marker_correlation_heatmap.png",
        )
        plot_celltype_bars(
            ct_summary,
            panel_set,
            "module_score",
            "cdkn_mean",
            OUTDIR / f"xenium_{panel_set}_module_score_cdkn_mean_correlation_by_celltype.png",
        )
        plot_celltype_bars(
            ct_summary,
            panel_set,
            "module_score",
            "proliferation_mean",
            OUTDIR / f"xenium_{panel_set}_module_score_proliferation_mean_correlation_by_celltype.png",
        )

    print(f"Saved outputs to: {OUTDIR}")
    focus = corr_summary[
        corr_summary["panel_set"].eq("20_gene_at_least_2_study_subset")
        & corr_summary["score"].isin(["module_score", "sasp_num"])
        & corr_summary["marker"].isin(["CDKN1A", "CDKN2A", "cdkn_mean", "PCNA", "MKI67", "STMN1", "proliferation_mean"])
    ].sort_values(["score", "marker"])
    print(focus[["score", "marker", "n_samples", "mean_spearman_rho", "sem_spearman_rho", "p_rho_vs_zero_sample_level", "fdr_rho_vs_zero_sample_level"]].to_string(index=False))


if __name__ == "__main__":
    main()
