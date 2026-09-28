"""
Ask whether public secretory/SASP-like scores are higher near CDKN+ cells.

Distances are computed within each sample and FOV from cell centroid coordinates.
The primary distance is to the nearest *other* CDKN+ cell, so CDKN+ cells do not
trivially have distance zero to themselves.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree


PROJECT_ROOT = Path(__file__).resolve().parents[3]
CELL_METADATA_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_snc_hotspot_neighborhoods"
    / "xenium_snc_hotspot_cell_neighborhoods.csv"
)
MARKER_EXPR_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_secretory_score_marker_correlations"
    / "xenium_cdkn_proliferation_marker_log_expression.csv.gz"
)
CELL_SCORES_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_public_secretory_panel_scores"
    / "xenium_public_secretory_cell_scores.csv.gz"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_secretory_score_cdkn_proximity"

MARKER_GROUPS = {
    "CDKN1A_pos": ["CDKN1A"],
    "CDKN2A_pos": ["CDKN2A"],
    "CDKN1A_or_CDKN2A_pos": ["CDKN1A", "CDKN2A"],
}
RADII_UM = [35, 70, 140]
DISTANCE_BINS = [-0.001, 35, 70, 140, np.inf]
DISTANCE_LABELS = ["<=35 um", "35-70 um", "70-140 um", ">140 um"]
SCORE_COLS = ["module_score", "auc_like_score", "panel_mean_log1p_cp10k", "n_positive_panel_genes", "sasp_num"]


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


def load_cell_table() -> pd.DataFrame:
    meta_cols = [
        "cell_id",
        "fov",
        "x",
        "y",
        "sample",
        "age",
        "snc",
        "sasp_num",
        "cell_type_major",
        "cell_type_minor",
        "cell_type_final",
        "zone",
    ]
    meta = pd.read_csv(CELL_METADATA_PATH, usecols=meta_cols)
    markers = pd.read_csv(MARKER_EXPR_PATH)
    cell = meta.merge(markers, on=["cell_id", "sample"], how="left")
    for group_name, genes in MARKER_GROUPS.items():
        cell[group_name] = cell[genes].gt(0).any(axis=1)
    return cell


def add_nearest_distances(cell: pd.DataFrame) -> pd.DataFrame:
    cell = cell.copy()
    for group_name in MARKER_GROUPS:
        cell[f"nearest_other_{group_name}_um"] = np.nan
        cell[f"nearest_including_self_{group_name}_um"] = np.nan

    for (sample, fov), idx in cell.groupby(["sample", "fov"], observed=True).groups.items():
        sub = cell.loc[idx]
        coords = sub[["x", "y"]].to_numpy(dtype=float)
        cell_ids = sub["cell_id"].to_numpy()
        for group_name in MARKER_GROUPS:
            pos_mask = sub[group_name].to_numpy(dtype=bool)
            if pos_mask.sum() == 0:
                continue

            pos_coords = coords[pos_mask]
            pos_cell_ids = cell_ids[pos_mask]
            tree = cKDTree(pos_coords)

            including_self, nearest_idx = tree.query(coords, k=1)
            cell.loc[idx, f"nearest_including_self_{group_name}_um"] = including_self

            k = min(2, len(pos_coords))
            distances, neighbor_idx = tree.query(coords, k=k)
            if k == 1:
                nearest_other = distances.astype(float)
                is_self = pos_cell_ids[neighbor_idx] == cell_ids
                nearest_other[is_self] = np.nan
            else:
                nearest_other = []
                for cell_id, dists, nidx in zip(cell_ids, distances, neighbor_idx):
                    candidate_dist = np.nan
                    for dist, pos_i in zip(np.atleast_1d(dists), np.atleast_1d(nidx)):
                        if pos_cell_ids[pos_i] != cell_id:
                            candidate_dist = float(dist)
                            break
                    nearest_other.append(candidate_dist)
                nearest_other = np.array(nearest_other, dtype=float)
            cell.loc[idx, f"nearest_other_{group_name}_um"] = nearest_other

    for group_name in MARKER_GROUPS:
        dist_col = f"nearest_other_{group_name}_um"
        cell[f"distance_bin_{group_name}"] = pd.cut(
            cell[dist_col],
            bins=DISTANCE_BINS,
            labels=DISTANCE_LABELS,
            include_lowest=True,
        )
        for radius in RADII_UM:
            cell[f"near_other_{group_name}_{radius}um"] = cell[dist_col].le(radius)
    return cell


def merge_scores(cell: pd.DataFrame) -> pd.DataFrame:
    scores = pd.read_csv(CELL_SCORES_PATH)
    score_cols = [
        "cell_id",
        "sample",
        "panel_set",
        "module_score",
        "auc_like_score",
        "panel_mean_log1p_cp10k",
        "n_positive_panel_genes",
    ]
    merged = scores[score_cols].merge(cell, on=["cell_id", "sample"], how="inner")
    return merged


def marker_counts(cell: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sample, fov), sub in cell.groupby(["sample", "fov"], observed=True):
        row = {"sample": sample, "fov": fov, "n_cells": len(sub)}
        for group_name in MARKER_GROUPS:
            row[f"n_{group_name}"] = int(sub[group_name].sum())
            row[f"pct_{group_name}"] = float(sub[group_name].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def proximity_summaries(scored: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_rows = []
    test_rows = []
    for group_name in MARKER_GROUPS:
        for radius in RADII_UM:
            near_col = f"near_other_{group_name}_{radius}um"
            for (panel_set, sample), sub in scored.groupby(["panel_set", "sample"], observed=True):
                for score_col in SCORE_COLS:
                    near = sub[sub[near_col].eq(True)][score_col].dropna()
                    far = sub[sub[near_col].eq(False)][score_col].dropna()
                    if len(near) == 0 or len(far) == 0:
                        continue
                    sample_rows.append(
                        {
                            "marker_group": group_name,
                            "radius_um": radius,
                            "panel_set": panel_set,
                            "sample": sample,
                            "age": sub["age"].iloc[0],
                            "score": score_col,
                            "n_near": len(near),
                            "n_far": len(far),
                            "mean_near": near.mean(),
                            "mean_far": far.mean(),
                            "diff_near_minus_far": near.mean() - far.mean(),
                        }
                    )
    sample_summary = pd.DataFrame(sample_rows)

    for (marker_group, radius, panel_set, score), sub in sample_summary.groupby(
        ["marker_group", "radius_um", "panel_set", "score"], observed=True
    ):
        diffs = sub["diff_near_minus_far"].dropna()
        test_rows.append(
            {
                "marker_group": marker_group,
                "radius_um": radius,
                "panel_set": panel_set,
                "score": score,
                "n_samples": sub["sample"].nunique(),
                "mean_diff_near_minus_far": diffs.mean(),
                "sem_diff_near_minus_far": sem(diffs),
                "p_diff_vs_zero_sample_level": stats.ttest_1samp(diffs, 0).pvalue if len(diffs) >= 3 else np.nan,
                "samples": ", ".join(sorted(sub["sample"].unique())),
            }
        )
    tests = pd.DataFrame(test_rows)
    if not tests.empty:
        tests["fdr_diff_vs_zero_sample_level"] = tests.groupby(["marker_group", "radius_um", "panel_set"], observed=True)[
            "p_diff_vs_zero_sample_level"
        ].transform(bh_fdr)
    return sample_summary, tests.sort_values(
        ["marker_group", "radius_um", "panel_set", "score"],
        ignore_index=True,
    )


def distance_bin_summary(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group_name in MARKER_GROUPS:
        bin_col = f"distance_bin_{group_name}"
        summary = (
            scored.dropna(subset=[bin_col])
            .groupby(["marker_group_placeholder"], observed=True)
            if False
            else None
        )
        tmp = (
            scored.dropna(subset=[bin_col])
            .groupby(["panel_set", "sample", "age", bin_col], observed=True)
            .agg(
                n_cells=("cell_id", "nunique"),
                mean_module_score=("module_score", "mean"),
                mean_auc_like_score=("auc_like_score", "mean"),
                mean_panel_log1p_cp10k=("panel_mean_log1p_cp10k", "mean"),
                mean_positive_panel_genes=("n_positive_panel_genes", "mean"),
                mean_sasp_num=("sasp_num", "mean"),
            )
            .reset_index()
            .rename(columns={bin_col: "distance_bin"})
        )
        tmp.insert(0, "marker_group", group_name)
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True)


def distance_correlations(scored: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for group_name in MARKER_GROUPS:
        dist_col = f"nearest_other_{group_name}_um"
        for (panel_set, sample), sub in scored.dropna(subset=[dist_col]).groupby(["panel_set", "sample"], observed=True):
            for score_col in SCORE_COLS:
                frame = sub[[dist_col, score_col, "age"]].dropna()
                if len(frame) < 20 or frame[dist_col].nunique() < 2 or frame[score_col].nunique() < 2:
                    rho, p = np.nan, np.nan
                else:
                    rho, p = stats.spearmanr(frame[dist_col], frame[score_col])
                rows.append(
                    {
                        "marker_group": group_name,
                        "panel_set": panel_set,
                        "sample": sample,
                        "age": sub["age"].iloc[0],
                        "score": score_col,
                        "n_cells": len(frame),
                        "spearman_rho_distance_vs_score": rho,
                        "spearman_p_cell_level_descriptive": p,
                    }
                )
    corr = pd.DataFrame(rows)
    summary_rows = []
    for (marker_group, panel_set, score), sub in corr.groupby(["marker_group", "panel_set", "score"], observed=True):
        rhos = sub["spearman_rho_distance_vs_score"].dropna()
        summary_rows.append(
            {
                "marker_group": marker_group,
                "panel_set": panel_set,
                "score": score,
                "n_samples": sub["sample"].nunique(),
                "mean_spearman_rho_distance_vs_score": rhos.mean(),
                "sem_spearman_rho_distance_vs_score": sem(rhos),
                "p_rho_vs_zero_sample_level": stats.ttest_1samp(rhos, 0).pvalue if len(rhos) >= 3 else np.nan,
                "samples": ", ".join(sorted(sub["sample"].unique())),
            }
        )
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary["fdr_rho_vs_zero_sample_level"] = summary.groupby(["marker_group", "panel_set"], observed=True)[
            "p_rho_vs_zero_sample_level"
        ].transform(bh_fdr)
    return corr, summary


def celltype_proximity_tests(scored: pd.DataFrame, radius: int = 35) -> pd.DataFrame:
    rows = []
    for group_name in MARKER_GROUPS:
        near_col = f"near_other_{group_name}_{radius}um"
        for (panel_set, cell_type, sample), sub in scored.groupby(
            ["panel_set", "cell_type_minor", "sample"], observed=True
        ):
            if len(sub) < 50:
                continue
            for score_col in ["module_score", "auc_like_score", "sasp_num"]:
                near = sub[sub[near_col].eq(True)][score_col].dropna()
                far = sub[sub[near_col].eq(False)][score_col].dropna()
                if len(near) < 10 or len(far) < 10:
                    continue
                rows.append(
                    {
                        "marker_group": group_name,
                        "radius_um": radius,
                        "panel_set": panel_set,
                        "cell_type_minor": cell_type,
                        "sample": sample,
                        "age": sub["age"].iloc[0],
                        "score": score_col,
                        "n_near": len(near),
                        "n_far": len(far),
                        "diff_near_minus_far": near.mean() - far.mean(),
                    }
                )
    sample = pd.DataFrame(rows)
    if sample.empty:
        return sample
    tests = (
        sample.groupby(["marker_group", "radius_um", "panel_set", "cell_type_minor", "score"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_diff_near_minus_far=("diff_near_minus_far", "mean"),
            sem_diff_near_minus_far=("diff_near_minus_far", sem),
            samples=("sample", lambda x: ", ".join(sorted(x.unique()))),
        )
        .reset_index()
    )
    p_values = []
    for _, row in tests.iterrows():
        mask = (
            sample["marker_group"].eq(row["marker_group"])
            & sample["panel_set"].eq(row["panel_set"])
            & sample["cell_type_minor"].eq(row["cell_type_minor"])
            & sample["score"].eq(row["score"])
        )
        diffs = sample.loc[mask, "diff_near_minus_far"].dropna()
        p_values.append(stats.ttest_1samp(diffs, 0).pvalue if len(diffs) >= 3 else np.nan)
    tests["p_diff_vs_zero_sample_level"] = p_values
    tests["fdr_diff_vs_zero_sample_level"] = tests.groupby(["marker_group", "radius_um", "panel_set", "score"], observed=True)[
        "p_diff_vs_zero_sample_level"
    ].transform(bh_fdr)
    return tests.sort_values(
        ["marker_group", "panel_set", "score", "mean_diff_near_minus_far"],
        ascending=[True, True, True, False],
    )


def plot_distance_bins(distance_bins: pd.DataFrame, marker_group: str, panel_set: str, outpath: Path) -> None:
    sub = distance_bins[
        distance_bins["marker_group"].eq(marker_group)
        & distance_bins["panel_set"].eq(panel_set)
    ]
    summary = (
        sub.groupby("distance_bin", observed=True)
        .agg(
            mean=("mean_module_score", "mean"),
            error=("mean_module_score", sem),
            n_samples=("sample", "nunique"),
        )
        .reindex(DISTANCE_LABELS)
        .reset_index()
    )
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(summary["distance_bin"].astype(str), summary["mean"], yerr=summary["error"], color="#4C78A8")
    ax.set_xlabel(f"Distance to nearest other {marker_group.replace('_pos', '+').replace('_', ' ')} cell")
    ax.set_ylabel("Mean module score across samples")
    ax.set_title(panel_set.replace("_", " "))
    ax.tick_params(axis="x", rotation=25)
    ax.grid(axis="y", color="0.9", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    fig.savefig(outpath.with_suffix(".pdf"))
    plt.close(fig)


def plot_near_far(sample_summary: pd.DataFrame, marker_group: str, radius: int, panel_set: str, outpath: Path) -> None:
    sub = sample_summary[
        sample_summary["marker_group"].eq(marker_group)
        & sample_summary["radius_um"].eq(radius)
        & sample_summary["panel_set"].eq(panel_set)
        & sample_summary["score"].eq("module_score")
    ].copy()
    plot = pd.DataFrame(
        {
            "group": ["far", "near"],
            "mean": [sub["mean_far"].mean(), sub["mean_near"].mean()],
            "error": [sem(sub["mean_far"]), sem(sub["mean_near"])],
        }
    )
    fig, ax = plt.subplots(figsize=(4.5, 4))
    ax.bar(plot["group"], plot["mean"], yerr=plot["error"], color=["#9E9E9E", "#F58518"])
    ax.set_ylabel("Mean module score across samples")
    ax.set_title(f"{marker_group.replace('_pos', '+').replace('_', ' ')} within {radius} um")
    ax.grid(axis="y", color="0.9", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    fig.savefig(outpath.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    cell = load_cell_table()
    cell = add_nearest_distances(cell)
    cell.to_csv(OUTDIR / "xenium_cdkn_positive_cell_proximity_metadata.csv.gz", index=False)
    marker_counts(cell).to_csv(OUTDIR / "xenium_cdkn_positive_counts_by_sample_fov.csv", index=False)

    scored = merge_scores(cell)
    scored.to_csv(OUTDIR / "xenium_secretory_scores_with_cdkn_proximity.csv.gz", index=False)

    sample_summary, proximity_tests = proximity_summaries(scored)
    sample_summary.to_csv(OUTDIR / "xenium_secretory_score_cdkn_proximity_by_sample.csv", index=False)
    proximity_tests.to_csv(OUTDIR / "xenium_secretory_score_cdkn_proximity_sample_level_tests.csv", index=False)

    bins = distance_bin_summary(scored)
    bins.to_csv(OUTDIR / "xenium_secretory_score_by_cdkn_distance_bin_sample.csv", index=False)

    corr, corr_summary = distance_correlations(scored)
    corr.to_csv(OUTDIR / "xenium_secretory_score_cdkn_distance_correlations_by_sample.csv", index=False)
    corr_summary.to_csv(OUTDIR / "xenium_secretory_score_cdkn_distance_correlation_summary.csv", index=False)

    celltype_tests = celltype_proximity_tests(scored, radius=35)
    celltype_tests.to_csv(OUTDIR / "xenium_secretory_score_cdkn_35um_proximity_by_celltype_tests.csv", index=False)

    panel_set = "20_gene_at_least_2_study_subset"
    for marker_group in MARKER_GROUPS:
        plot_distance_bins(
            bins,
            marker_group,
            panel_set,
            OUTDIR / f"xenium_{panel_set}_{marker_group}_distance_bin_module_score.png",
        )
        plot_near_far(
            sample_summary,
            marker_group,
            35,
            panel_set,
            OUTDIR / f"xenium_{panel_set}_{marker_group}_35um_near_far_module_score.png",
        )

    print(f"Saved outputs to: {OUTDIR}")
    focus = proximity_tests[
        proximity_tests["panel_set"].eq(panel_set)
        & proximity_tests["score"].isin(["module_score", "auc_like_score", "sasp_num"])
        & proximity_tests["radius_um"].isin([35, 70])
    ].sort_values(["marker_group", "radius_um", "score"])
    print(focus[[
        "marker_group",
        "radius_um",
        "score",
        "n_samples",
        "mean_diff_near_minus_far",
        "sem_diff_near_minus_far",
        "p_diff_vs_zero_sample_level",
        "fdr_diff_vs_zero_sample_level",
    ]].to_string(index=False))


if __name__ == "__main__":
    main()
