"""
Score the Xenium-measurable public secretory panel by cell type and hotspot zone.

Scores:
  - Seurat-like module score: mean log-normalized panel expression minus
    expression-matched control genes from the Xenium panel.
  - AUCell-like rank score: normalized recovery of panel genes among the top
    ranked genes in each cell.

The statistical summaries use sample-level means so cells are not treated as
independent biological replicates.
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
PANEL_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_public_secretory_panel_expression"
    / "xenium_public_secretory_selected_genes.csv"
)
CELL_METADATA_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_snc_hotspot_neighborhoods"
    / "xenium_snc_hotspot_cell_neighborhoods.csv"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_public_secretory_panel_scores"

ANNOTATION_COLS = ["cell_type_major", "cell_type_minor", "cell_type_final", "immune_subclusters1"]
PANEL_SETS = {
    "18_gene_exact_2_study_subset": "in_18_gene_xenium_subset",
    "20_gene_at_least_2_study_subset": "in_20_gene_xenium_subset",
}


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


def age_from_sample(sample: str) -> str:
    if sample.upper().endswith("-Y"):
        return "Young"
    if sample.upper().endswith("-O"):
        return "Old"
    return "Unknown"


def load_panel() -> pd.DataFrame:
    panel = pd.read_csv(PANEL_PATH)
    panel["gene"] = panel["gene"].astype(str).str.upper()
    for col in PANEL_SETS.values():
        panel[col] = panel[col].astype(bool)
    return panel


def load_metadata() -> pd.DataFrame:
    cols = [
        "cell_id",
        "sample",
        "age",
        "snc",
        "sasp_num",
        "zone",
        "hex_distance_to_snc_hotspot",
        *ANNOTATION_COLS,
    ]
    meta = pd.read_csv(CELL_METADATA_PATH, usecols=cols)
    meta["age"] = meta["age"].fillna(meta["sample"].map(age_from_sample))
    meta["hotspot_binary"] = np.where(meta["zone"].eq("hotspot_core"), "Hotspot core", "Not hotspot core")
    meta["hotspot_or_neighbor1"] = np.where(
        meta["zone"].isin(["hotspot_core", "neighbor_1"]),
        "Hotspot core + neighbor 1",
        "Background / farther",
    )
    return meta


def read_sample_log_matrix(sample: str, wanted_cells: set[str]) -> tuple[list[str], list[str], np.ndarray]:
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
    if not np.any(keep):
        return genes, [], np.empty((len(genes), 0), dtype=np.float32)

    counts = matrix[:, keep].astype(np.float32).toarray()
    totals = counts.sum(axis=0)
    norm = np.divide(counts, totals, out=np.zeros_like(counts), where=totals > 0) * 10000
    log_norm = np.log1p(norm).astype(np.float32)
    return genes, full_cell_ids[keep].tolist(), log_norm


def matched_control_indices(
    gene_means: np.ndarray,
    target_idx: list[int],
    *,
    n_bins: int = 12,
    ctrl_per_gene: int = 25,
    seed: int = 1,
) -> list[int]:
    rng = np.random.default_rng(seed)
    order = np.argsort(gene_means)
    bins = np.empty(len(gene_means), dtype=int)
    for bin_id, idx in enumerate(np.array_split(order, n_bins)):
        bins[idx] = bin_id

    target_set = set(target_idx)
    controls = []
    for idx in target_idx:
        candidates = np.where(bins == bins[idx])[0]
        candidates = [int(i) for i in candidates if int(i) not in target_set]
        if len(candidates) < 3:
            nearest = np.argsort(np.abs(gene_means - gene_means[idx]))
            candidates = [int(i) for i in nearest if int(i) not in target_set and int(i) != idx]
        if len(candidates) > ctrl_per_gene:
            candidates = rng.choice(candidates, size=ctrl_per_gene, replace=False).tolist()
        controls.extend(candidates)
    return sorted(set(controls))


def auc_like_score(log_norm: np.ndarray, target_idx: list[int], top_fraction: float = 0.10) -> np.ndarray:
    n_genes, n_cells = log_norm.shape
    top_n = max(50, int(np.ceil(n_genes * top_fraction)))
    target_set = set(target_idx)
    scores = np.zeros(n_cells, dtype=np.float32)
    for cell_i in range(n_cells):
        ranked = np.argsort(-log_norm[:, cell_i], kind="mergesort")[:top_n]
        score = 0.0
        for rank_i, gene_i in enumerate(ranked):
            if int(gene_i) in target_set:
                score += (top_n - rank_i) / top_n
        scores[cell_i] = score / max(len(target_idx), 1)
    return scores


def score_sample(sample: str, meta: pd.DataFrame, panel: pd.DataFrame) -> pd.DataFrame:
    sample_meta = meta[meta["sample"].eq(sample)].copy()
    genes, cell_ids, log_norm = read_sample_log_matrix(sample, set(sample_meta["cell_id"]))
    if not cell_ids:
        return pd.DataFrame()

    gene_lookup = {gene: i for i, gene in enumerate(genes)}
    gene_means = log_norm.mean(axis=1)
    rows = []
    for panel_set, flag_col in PANEL_SETS.items():
        panel_genes = [gene for gene in panel.loc[panel[flag_col], "gene"] if gene in gene_lookup]
        target_idx = [gene_lookup[gene] for gene in panel_genes]
        control_idx = matched_control_indices(gene_means, target_idx)

        target_mean = log_norm[target_idx, :].mean(axis=0)
        control_mean = log_norm[control_idx, :].mean(axis=0) if control_idx else np.zeros(len(cell_ids))
        module_score = target_mean - control_mean
        auc_score = auc_like_score(log_norm, target_idx)
        n_positive_genes = (log_norm[target_idx, :] > 0).sum(axis=0)

        rows.append(
            pd.DataFrame(
                {
                    "cell_id": cell_ids,
                    "sample": sample,
                    "age": age_from_sample(sample),
                    "panel_set": panel_set,
                    "n_panel_genes_present": len(panel_genes),
                    "n_control_genes": len(control_idx),
                    "module_score": module_score,
                    "panel_mean_log1p_cp10k": target_mean,
                    "control_mean_log1p_cp10k": control_mean,
                    "auc_like_score": auc_score,
                    "n_positive_panel_genes": n_positive_genes,
                    "fraction_positive_panel_genes": n_positive_genes / max(len(panel_genes), 1),
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def add_cell_type_long(cell_scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    id_cols = [
        "cell_id",
        "sample",
        "age",
        "panel_set",
        "module_score",
        "panel_mean_log1p_cp10k",
        "control_mean_log1p_cp10k",
        "auc_like_score",
        "n_positive_panel_genes",
        "fraction_positive_panel_genes",
        "zone",
        "hotspot_binary",
        "hotspot_or_neighbor1",
        "snc",
        "sasp_num",
    ]
    for annotation_col in ANNOTATION_COLS:
        tmp = cell_scores[id_cols + [annotation_col]].copy()
        tmp = tmp.rename(columns={annotation_col: "cell_type"})
        tmp["annotation_col"] = annotation_col
        tmp = tmp[tmp["cell_type"].notna() & ~tmp["cell_type"].astype(str).eq("none")]
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True)


def summarize_cell_types(cell_type_long: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_summary = (
        cell_type_long.groupby(["annotation_col", "cell_type", "panel_set", "sample", "age"], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            mean_module_score=("module_score", "mean"),
            mean_auc_like_score=("auc_like_score", "mean"),
            mean_panel_log1p_cp10k=("panel_mean_log1p_cp10k", "mean"),
            mean_positive_panel_genes=("n_positive_panel_genes", "mean"),
        )
        .reset_index()
    )
    across_samples = (
        sample_summary.groupby(["annotation_col", "cell_type", "panel_set"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            total_cells=("n_cells", "sum"),
            mean_module_score=("mean_module_score", "mean"),
            sem_module_score=("mean_module_score", sem),
            mean_auc_like_score=("mean_auc_like_score", "mean"),
            sem_auc_like_score=("mean_auc_like_score", sem),
            mean_panel_log1p_cp10k=("mean_panel_log1p_cp10k", "mean"),
            mean_positive_panel_genes=("mean_positive_panel_genes", "mean"),
        )
        .reset_index()
        .sort_values(["annotation_col", "panel_set", "mean_module_score"], ascending=[True, True, False])
    )
    return sample_summary, across_samples


def cell_type_enrichment_tests(cell_type_long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (annotation_col, cell_type, panel_set, sample), sub in cell_type_long.groupby(
        ["annotation_col", "cell_type", "panel_set", "sample"], observed=True
    ):
        sample_all = cell_type_long[
            cell_type_long["annotation_col"].eq(annotation_col)
            & cell_type_long["panel_set"].eq(panel_set)
            & cell_type_long["sample"].eq(sample)
        ]
        in_type = sample_all["cell_type"].eq(cell_type)
        if in_type.sum() < 10 or (~in_type).sum() < 10:
            continue
        rows.append(
            {
                "annotation_col": annotation_col,
                "cell_type": cell_type,
                "panel_set": panel_set,
                "sample": sample,
                "age": sample_all["age"].iloc[0],
                "n_cells_type": int(in_type.sum()),
                "n_cells_other": int((~in_type).sum()),
                "mean_type_module_score": sample_all.loc[in_type, "module_score"].mean(),
                "mean_other_module_score": sample_all.loc[~in_type, "module_score"].mean(),
                "diff_type_minus_other_module_score": sample_all.loc[in_type, "module_score"].mean()
                - sample_all.loc[~in_type, "module_score"].mean(),
                "mean_type_auc_like_score": sample_all.loc[in_type, "auc_like_score"].mean(),
                "mean_other_auc_like_score": sample_all.loc[~in_type, "auc_like_score"].mean(),
                "diff_type_minus_other_auc_like_score": sample_all.loc[in_type, "auc_like_score"].mean()
                - sample_all.loc[~in_type, "auc_like_score"].mean(),
            }
        )
    diffs = pd.DataFrame(rows)
    if diffs.empty:
        return diffs

    tests = []
    for (annotation_col, cell_type, panel_set), sub in diffs.groupby(
        ["annotation_col", "cell_type", "panel_set"], observed=True
    ):
        values = sub["diff_type_minus_other_module_score"].dropna()
        auc_values = sub["diff_type_minus_other_auc_like_score"].dropna()
        p_module = stats.ttest_1samp(values, 0).pvalue if len(values) >= 3 else np.nan
        p_auc = stats.ttest_1samp(auc_values, 0).pvalue if len(auc_values) >= 3 else np.nan
        tests.append(
            {
                "annotation_col": annotation_col,
                "cell_type": cell_type,
                "panel_set": panel_set,
                "n_samples": sub["sample"].nunique(),
                "total_cells_type": sub["n_cells_type"].sum(),
                "mean_diff_module_score": values.mean(),
                "sem_diff_module_score": sem(values),
                "p_module_score_sample_level": p_module,
                "mean_diff_auc_like_score": auc_values.mean(),
                "sem_diff_auc_like_score": sem(auc_values),
                "p_auc_like_score_sample_level": p_auc,
                "samples": ", ".join(sorted(sub["sample"].unique())),
            }
        )
    out = pd.DataFrame(tests)
    out["fdr_module_score"] = out.groupby(["annotation_col", "panel_set"], observed=True)[
        "p_module_score_sample_level"
    ].transform(bh_fdr)
    out["fdr_auc_like_score"] = out.groupby(["annotation_col", "panel_set"], observed=True)[
        "p_auc_like_score_sample_level"
    ].transform(bh_fdr)
    return out.sort_values(["annotation_col", "panel_set", "mean_diff_module_score"], ascending=[True, True, False])


def summarize_hotspots(cell_scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_zone = (
        cell_scores.groupby(["panel_set", "sample", "age", "zone"], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            mean_module_score=("module_score", "mean"),
            mean_auc_like_score=("auc_like_score", "mean"),
            mean_panel_log1p_cp10k=("panel_mean_log1p_cp10k", "mean"),
            mean_positive_panel_genes=("n_positive_panel_genes", "mean"),
        )
        .reset_index()
    )

    comparisons = []
    for panel_set, sub in sample_zone.groupby("panel_set", observed=True):
        for comparator in ["background", "neighbor_1", "neighbor_2"]:
            wide = sub[sub["zone"].isin(["hotspot_core", comparator])].pivot(
                index=["sample", "age"], columns="zone", values="mean_module_score"
            )
            if "hotspot_core" not in wide.columns or comparator not in wide.columns:
                continue
            diffs = (wide["hotspot_core"] - wide[comparator]).dropna()
            p = stats.ttest_1samp(diffs, 0).pvalue if len(diffs) >= 3 else np.nan
            auc_wide = sub[sub["zone"].isin(["hotspot_core", comparator])].pivot(
                index=["sample", "age"], columns="zone", values="mean_auc_like_score"
            )
            auc_diffs = (auc_wide["hotspot_core"] - auc_wide[comparator]).dropna()
            p_auc = stats.ttest_1samp(auc_diffs, 0).pvalue if len(auc_diffs) >= 3 else np.nan
            comparisons.append(
                {
                    "panel_set": panel_set,
                    "comparison": f"hotspot_core_vs_{comparator}",
                    "n_samples": len(diffs),
                    "mean_diff_module_score": diffs.mean(),
                    "sem_diff_module_score": sem(diffs),
                    "p_module_score_sample_level": p,
                    "mean_diff_auc_like_score": auc_diffs.mean(),
                    "sem_diff_auc_like_score": sem(auc_diffs),
                    "p_auc_like_score_sample_level": p_auc,
                    "samples": ", ".join(sorted(diffs.index.get_level_values("sample").unique())),
                }
            )
    comparisons = pd.DataFrame(comparisons)
    if not comparisons.empty:
        comparisons["fdr_module_score"] = comparisons.groupby("panel_set", observed=True)[
            "p_module_score_sample_level"
        ].transform(bh_fdr)
        comparisons["fdr_auc_like_score"] = comparisons.groupby("panel_set", observed=True)[
            "p_auc_like_score_sample_level"
        ].transform(bh_fdr)
    return sample_zone, comparisons


def summarize_hotspots_by_cell_type(cell_type_long: pd.DataFrame) -> pd.DataFrame:
    sample_summary = (
        cell_type_long.groupby(["annotation_col", "cell_type", "panel_set", "sample", "age", "zone"], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            mean_module_score=("module_score", "mean"),
            mean_auc_like_score=("auc_like_score", "mean"),
        )
        .reset_index()
    )

    rows = []
    for (annotation_col, cell_type, panel_set), sub in sample_summary.groupby(
        ["annotation_col", "cell_type", "panel_set"], observed=True
    ):
        wide = sub[sub["zone"].isin(["hotspot_core", "background"])].pivot(
            index=["sample", "age"], columns="zone", values="mean_module_score"
        )
        if "hotspot_core" not in wide.columns or "background" not in wide.columns:
            continue
        diffs = (wide["hotspot_core"] - wide["background"]).dropna()
        if len(diffs) < 3:
            continue
        auc_wide = sub[sub["zone"].isin(["hotspot_core", "background"])].pivot(
            index=["sample", "age"], columns="zone", values="mean_auc_like_score"
        )
        auc_diffs = (auc_wide["hotspot_core"] - auc_wide["background"]).dropna()
        rows.append(
            {
                "annotation_col": annotation_col,
                "cell_type": cell_type,
                "panel_set": panel_set,
                "n_samples": len(diffs),
                "mean_hotspot_minus_background_module_score": diffs.mean(),
                "sem_hotspot_minus_background_module_score": sem(diffs),
                "p_module_score_sample_level": stats.ttest_1samp(diffs, 0).pvalue,
                "mean_hotspot_minus_background_auc_like_score": auc_diffs.mean(),
                "sem_hotspot_minus_background_auc_like_score": sem(auc_diffs),
                "p_auc_like_score_sample_level": stats.ttest_1samp(auc_diffs, 0).pvalue,
                "samples": ", ".join(sorted(diffs.index.get_level_values("sample").unique())),
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_module_score"] = out.groupby(["annotation_col", "panel_set"], observed=True)[
            "p_module_score_sample_level"
        ].transform(bh_fdr)
        out["fdr_auc_like_score"] = out.groupby(["annotation_col", "panel_set"], observed=True)[
            "p_auc_like_score_sample_level"
        ].transform(bh_fdr)
    return out.sort_values(
        ["annotation_col", "panel_set", "mean_hotspot_minus_background_module_score"],
        ascending=[True, True, False],
    )


def plot_top_cell_types(summary: pd.DataFrame, annotation_col: str, panel_set: str, outpath: Path) -> None:
    sub = summary[
        summary["annotation_col"].eq(annotation_col)
        & summary["panel_set"].eq(panel_set)
        & (summary["total_cells"] >= 50)
    ].nlargest(18, "mean_module_score")
    sub = sub.sort_values("mean_module_score")
    fig, ax = plt.subplots(figsize=(8, max(4, 0.28 * len(sub) + 1.5)))
    ax.barh(sub["cell_type"], sub["mean_module_score"], xerr=sub["sem_module_score"], color="#4C78A8")
    ax.set_xlabel("Mean Seurat-like module score across samples")
    ax.set_ylabel("")
    ax.set_title(f"{panel_set.replace('_', ' ')} by {annotation_col}")
    ax.grid(axis="x", color="0.9", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    fig.savefig(outpath.with_suffix(".pdf"))
    plt.close(fig)


def plot_hotspot_zone(sample_zone: pd.DataFrame, panel_set: str, outpath: Path) -> None:
    zone_order = ["background", "neighbor_2", "neighbor_1", "hotspot_core"]
    sub = sample_zone[sample_zone["panel_set"].eq(panel_set)].copy()
    summary = (
        sub.groupby("zone", observed=True)
        .agg(mean=("mean_module_score", "mean"), error=("mean_module_score", sem), n_samples=("sample", "nunique"))
        .reindex(zone_order)
        .reset_index()
    )
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(summary["zone"], summary["mean"], yerr=summary["error"], color="#F58518")
    ax.set_ylabel("Mean Seurat-like module score across samples")
    ax.set_xlabel("")
    ax.set_title(f"{panel_set.replace('_', ' ')} by hotspot zone")
    ax.tick_params(axis="x", rotation=35)
    ax.grid(axis="y", color="0.9", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    fig.savefig(outpath.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    panel = load_panel()
    meta = load_metadata()

    score_frames = []
    for sample in sorted(meta["sample"].unique()):
        score_frames.append(score_sample(sample, meta, panel))
    cell_scores = pd.concat(score_frames, ignore_index=True)
    cell_scores = cell_scores.merge(meta, on=["cell_id", "sample", "age"], how="inner")
    cell_scores.to_csv(OUTDIR / "xenium_public_secretory_cell_scores.csv.gz", index=False)

    cell_type_long = add_cell_type_long(cell_scores)
    cell_type_long.to_csv(OUTDIR / "xenium_public_secretory_cell_scores_celltype_long.csv.gz", index=False)

    celltype_by_sample, celltype_summary = summarize_cell_types(cell_type_long)
    celltype_by_sample.to_csv(OUTDIR / "xenium_public_secretory_scores_by_celltype_sample.csv", index=False)
    celltype_summary.to_csv(OUTDIR / "xenium_public_secretory_scores_by_celltype_summary.csv", index=False)

    celltype_tests = cell_type_enrichment_tests(cell_type_long)
    celltype_tests.to_csv(OUTDIR / "xenium_public_secretory_celltype_enrichment_sample_level_tests.csv", index=False)

    hotspot_by_sample, hotspot_tests = summarize_hotspots(cell_scores)
    hotspot_by_sample.to_csv(OUTDIR / "xenium_public_secretory_scores_by_hotspot_zone_sample.csv", index=False)
    hotspot_tests.to_csv(OUTDIR / "xenium_public_secretory_hotspot_enrichment_sample_level_tests.csv", index=False)

    hotspot_celltype_tests = summarize_hotspots_by_cell_type(cell_type_long)
    hotspot_celltype_tests.to_csv(
        OUTDIR / "xenium_public_secretory_hotspot_enrichment_by_celltype_sample_level_tests.csv",
        index=False,
    )

    for panel_set in PANEL_SETS:
        plot_top_cell_types(
            celltype_summary,
            "cell_type_minor",
            panel_set,
            OUTDIR / f"xenium_public_secretory_{panel_set}_top_cell_type_minor_module_score.png",
        )
        plot_top_cell_types(
            celltype_summary,
            "cell_type_final",
            panel_set,
            OUTDIR / f"xenium_public_secretory_{panel_set}_top_cell_type_final_module_score.png",
        )
        plot_hotspot_zone(
            hotspot_by_sample,
            panel_set,
            OUTDIR / f"xenium_public_secretory_{panel_set}_hotspot_zone_module_score.png",
        )

    print(f"Saved outputs to: {OUTDIR}")
    print("Top cell_type_minor, 20-gene module score:")
    top_minor = celltype_summary[
        celltype_summary["annotation_col"].eq("cell_type_minor")
        & celltype_summary["panel_set"].eq("20_gene_at_least_2_study_subset")
        & (celltype_summary["total_cells"] >= 50)
    ].nlargest(12, "mean_module_score")
    print(top_minor[["cell_type", "n_samples", "total_cells", "mean_module_score", "mean_auc_like_score"]].to_string(index=False))
    print("\nHotspot enrichment tests:")
    print(hotspot_tests.to_string(index=False))


if __name__ == "__main__":
    main()
