"""
Score the Xenium-measurable 20-gene SSSP with scanpy.tl.score_genes.

This is the Python/Scanpy analogue of a Seurat-style module score:
average expression of the target gene set minus average expression of
expression-binned control genes from the same Xenium panel.

Cell IDs and annotations should come from the annotated Seurat RDS metadata
export created by 00_export_xenium_rds_metadata.R. Raw counts come from the
sample-level Xenium cell_feature_matrix.h5 files and are filtered to those
RDS-derived cells.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats


PROJECT_ROOT = Path(__file__).resolve().parents[3]
XENIUM_DIR = PROJECT_ROOT / "data" / "xenium"
RDS_METADATA_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_rds_metadata"
    / "xenium_rds_cell_metadata.csv.gz"
)
SPATIAL_COORDS_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_snc_hotspot_neighborhoods"
    / "xenium_snc_hotspot_cell_neighborhoods.csv"
)
FALLBACK_CELL_SCORE_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_public_secretory_panel_scores"
    / "xenium_public_secretory_cell_scores.csv.gz"
)
GENE_CONFIG_PATH = Path(__file__).resolve().with_name("xenium_sssp_gene_set_config.csv")
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_sssp_20gene_scanpy_score_genes"
SCORED_H5AD_PATH = OUTDIR / "xenium_20gene_sssp_scanpy_scored.h5ad"
PANEL_SET = "20_gene_at_least_2_study_subset"
SCORE_NAME = "sssp_20gene_scanpy_score"
GENE_SET_CONFIGS = [
    {
        "gene_set_id": "main_20",
        "membership_col": "in_xenium_main_20",
        "score_name": "sssp_20gene_scanpy_score",
        "high_gt0_col": "sssp_20gene_scanpy_high_gt0",
        "outdir": OUTDIR,
        "prefix": "xenium_20gene_sssp_scanpy",
        "var_flag": "in_20gene_sssp_panel",
        "uns_key": "sssp_20gene_scanpy",
        "plot_title": "Xenium validation of 20-gene SSSP using scanpy.tl.score_genes",
    },
    {
        "gene_set_id": "diff_filtered_17",
        "membership_col": "in_xenium_diff_filtered",
        "score_name": "sssp_17gene_diff_filtered_scanpy_score",
        "high_gt0_col": "sssp_17gene_diff_filtered_scanpy_high_gt0",
        "outdir": PROJECT_ROOT / "results" / "260727" / "xenium_sssp_17gene_diff_filtered_scanpy_score_genes",
        "prefix": "xenium_17gene_diff_filtered_sssp_scanpy",
        "var_flag": "in_17gene_diff_filtered_sssp_panel",
        "uns_key": "sssp_17gene_diff_filtered_scanpy",
        "plot_title": "Xenium validation of differentiation-filtered SSSP using scanpy.tl.score_genes",
    },
]
AGE_ORDER = ["Young", "Old"]
AGE_COLORS = {"Young": "#4C78A8", "Old": "#F58518"}
ANNOTATION_COLUMNS = [
    "snc",
    "sasp_num",
    "cell_type_major",
    "cell_type_minor",
    "cell_type_final",
    "immune_subclusters1",
    "immcelltype",
]


def decode_array(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def age_from_sample(sample: str) -> str:
    if sample.upper().endswith("-Y"):
        return "Young"
    if sample.upper().endswith("-O"):
        return "Old"
    return "Unknown"


def sem(values: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if len(values) <= 1:
        return np.nan
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def load_gene_list(membership_col: str) -> pd.DataFrame:
    genes = pd.read_csv(GENE_CONFIG_PATH)
    genes["gene"] = genes["gene"].astype(str).str.upper()
    if membership_col not in genes.columns:
        raise KeyError(f"Missing {membership_col} in {GENE_CONFIG_PATH}")
    genes = genes[genes[membership_col].astype(bool)].copy()
    return genes


def load_cell_metadata() -> tuple[pd.DataFrame, str]:
    if RDS_METADATA_PATH.exists():
        meta = pd.read_csv(RDS_METADATA_PATH)
        source = str(RDS_METADATA_PATH)
    else:
        meta = pd.read_csv(FALLBACK_CELL_SCORE_PATH)
        meta = meta[meta["panel_set"].eq(PANEL_SET)].drop(columns="panel_set")
        source = f"{FALLBACK_CELL_SCORE_PATH} (fallback; run 00_export_xenium_rds_metadata.R)"

    if "cell_id" not in meta.columns:
        raise KeyError(f"Expected a cell_id column in metadata source: {source}")

    meta["cell_id"] = meta["cell_id"].astype(str)
    if "sample" not in meta.columns:
        meta["sample"] = meta["cell_id"].str.replace(r"_.*$", "", regex=True)
    meta["sample"] = meta["sample"].astype(str)

    if "age" not in meta.columns:
        meta["age"] = np.nan
    meta["age"] = meta["age"].fillna(meta["sample"].map(age_from_sample))

    h5_samples = {path.parent.name for path in XENIUM_DIR.glob("*/cell_feature_matrix.h5")}
    meta = meta[meta["sample"].isin(h5_samples)].copy()
    meta = meta.drop_duplicates("cell_id")
    return meta, source


def read_sample_counts(sample: str, wanted_cells: set[str], expected_genes: list[str] | None) -> tuple[list[str], pd.DataFrame, sparse.csr_matrix]:
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

    if expected_genes is not None and genes != expected_genes:
        raise ValueError(f"Feature order differs for sample {sample}; cannot safely stack matrices.")

    full_cell_ids = np.array([f"{sample}_{barcode}" for barcode in barcodes], dtype=object)
    keep = np.array([cell_id in wanted_cells for cell_id in full_cell_ids])
    if not np.any(keep):
        return genes, pd.DataFrame(), sparse.csr_matrix((0, len(genes)), dtype=np.float32)

    obs = pd.DataFrame(
        {
            "cell_id": full_cell_ids[keep],
            "sample": sample,
            "age": age_from_sample(sample),
        }
    )
    counts = matrix[:, keep].T.tocsr().astype(np.float32)
    return genes, obs, counts


def build_xenium_anndata(meta: pd.DataFrame) -> tuple[ad.AnnData, pd.DataFrame]:
    pieces = []
    obs_parts = []
    expected_genes = None
    loaded_cell_ids: set[str] = set()
    for sample in sorted(meta["sample"].unique()):
        genes, obs, counts = read_sample_counts(
            sample,
            set(meta.loc[meta["sample"].eq(sample), "cell_id"]),
            expected_genes,
        )
        if expected_genes is None:
            expected_genes = genes
        if len(obs):
            pieces.append(counts)
            obs_parts.append(obs)
            loaded_cell_ids.update(obs["cell_id"].astype(str))

    if expected_genes is None or not pieces:
        raise RuntimeError("No Xenium matrices/cells were loaded.")

    obs = pd.concat(obs_parts, ignore_index=True)
    obs = obs.merge(meta, on=["cell_id", "sample", "age"], how="left")
    obs.index = obs["cell_id"].astype(str)
    var = pd.DataFrame(index=pd.Index(expected_genes, name="gene"))
    matrix = sparse.vstack(pieces, format="csr")
    metadata_cell_ids = set(meta["cell_id"].astype(str))
    audit = pd.DataFrame(
        [
            {
                "n_metadata_cells": int(len(metadata_cell_ids)),
                "n_h5_loaded_cells": int(len(loaded_cell_ids)),
                "n_metadata_cells_not_in_h5": int(len(metadata_cell_ids - loaded_cell_ids)),
                "metadata_only_examples": "; ".join(sorted(metadata_cell_ids - loaded_cell_ids)[:10]),
            }
        ]
    )
    return ad.AnnData(X=matrix, obs=obs, var=var), audit


def add_spatial_coordinates(adata: ad.AnnData, audit: pd.DataFrame) -> tuple[ad.AnnData, pd.DataFrame]:
    if not SPATIAL_COORDS_PATH.exists():
        audit["spatial_coordinate_source"] = str(SPATIAL_COORDS_PATH)
        audit["n_cells_with_spatial_coordinates"] = 0
        audit["n_cells_missing_spatial_coordinates"] = int(adata.n_obs)
        return adata, audit

    coords = pd.read_csv(SPATIAL_COORDS_PATH, usecols=["cell_id", "fov", "x", "y"]).drop_duplicates("cell_id")
    coords["cell_id"] = coords["cell_id"].astype(str)

    obs = adata.obs.reset_index(drop=False).rename(columns={"index": "_obs_index"})
    obs = obs.merge(coords, on="cell_id", how="left")
    has_coords = obs[["x", "y"]].notna().all(axis=1)

    adata = adata[obs.loc[has_coords, "_obs_index"].astype(str)].copy()
    obs = obs.loc[has_coords].set_index("_obs_index")
    adata.obs = obs
    adata.obs_names = adata.obs.index.astype(str)
    adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy(dtype=float)

    audit["spatial_coordinate_source"] = str(SPATIAL_COORDS_PATH)
    audit["n_cells_with_spatial_coordinates"] = int(has_coords.sum())
    audit["n_cells_missing_spatial_coordinates"] = int((~has_coords).sum())
    return adata, audit


def add_panel_expression_summaries(adata: ad.AnnData, panel_genes: list[str], score_name: str) -> pd.DataFrame:
    panel = adata[:, panel_genes].X
    if sparse.issparse(panel):
        panel_mean = np.asarray(panel.mean(axis=1)).ravel()
        n_positive = np.asarray((panel > 0).sum(axis=1)).ravel()
    else:
        panel_mean = panel.mean(axis=1)
        n_positive = (panel > 0).sum(axis=1)
    return pd.DataFrame(
        {
            "cell_id": adata.obs["cell_id"].to_numpy(),
            "sample": adata.obs["sample"].to_numpy(),
            "age": adata.obs["age"].to_numpy(),
            score_name: adata.obs[score_name].to_numpy(),
            "panel_mean_log1p_cp10k": panel_mean,
            "n_positive_panel_genes": n_positive.astype(int),
            "fraction_positive_panel_genes": n_positive / len(panel_genes),
        }
    )


def add_available_annotations(cell_scores: pd.DataFrame, adata: ad.AnnData) -> pd.DataFrame:
    for col in ANNOTATION_COLUMNS:
        if col in adata.obs.columns:
            cell_scores[col] = adata.obs[col].to_numpy()
    return cell_scores


def add_score_summaries_to_adata_obs(
    adata: ad.AnnData,
    cell_scores: pd.DataFrame,
    score_name: str,
    high_gt0_col: str,
) -> None:
    summaries = cell_scores.set_index("cell_id")
    obs_cell_ids = adata.obs["cell_id"].astype(str)
    adata.obs[high_gt0_col] = adata.obs[score_name] > 0
    for col in ["panel_mean_log1p_cp10k", "n_positive_panel_genes", "fraction_positive_panel_genes"]:
        adata.obs[col] = summaries.loc[obs_cell_ids, col].to_numpy()


def write_plots(
    sample_summary: pd.DataFrame,
    cell_scores: pd.DataFrame,
    outbase: Path,
    score_name: str,
    plot_title: str,
) -> None:
    sample_summary = sample_summary.copy()
    sample_summary["age"] = pd.Categorical(sample_summary["age"], categories=AGE_ORDER, ordered=True)
    sample_summary = sample_summary.sort_values(["age", "sample"])
    sample_order = sample_summary["sample"].tolist()

    age_summary = (
        sample_summary.groupby("age", observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_score=("mean_scanpy_score", "mean"),
            sem_score=("mean_scanpy_score", sem),
        )
        .reset_index()
    )

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.3), gridspec_kw={"width_ratios": [0.8, 1.7]})
    ax = axes[0]
    for i, age in enumerate(AGE_ORDER):
        sub = sample_summary[sample_summary["age"].astype(str).eq(age)]
        row = age_summary[age_summary["age"].astype(str).eq(age)].iloc[0]
        ax.bar(
            i,
            row["mean_score"],
            yerr=row["sem_score"],
            color=AGE_COLORS[age],
            alpha=0.82,
            width=0.58,
            edgecolor="#333333",
            linewidth=0.8,
            capsize=4,
        )
        jitter = np.linspace(-0.08, 0.08, len(sub)) if len(sub) > 1 else np.array([0])
        ax.scatter(
            np.full(len(sub), i) + jitter,
            sub["mean_scanpy_score"],
            s=58,
            color="white",
            edgecolor="#222222",
            linewidth=0.9,
            zorder=4,
        )
    ax.set_xticks(range(len(AGE_ORDER)))
    ax.set_xticklabels(AGE_ORDER)
    ax.set_ylabel("Mean Scanpy SSSP score")
    ax.set_title("Composite young vs old", fontsize=12, fontweight="bold")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    x = np.arange(len(sample_summary))
    colors = [AGE_COLORS[str(age)] for age in sample_summary["age"]]
    ax.bar(
        x,
        sample_summary["mean_scanpy_score"],
        color=colors,
        alpha=0.85,
        edgecolor="#333333",
        linewidth=0.7,
        width=0.72,
    )
    ax.errorbar(
        x,
        sample_summary["mean_scanpy_score"],
        yerr=sample_summary["sem_scanpy_score"],
        fmt="none",
        ecolor="#333333",
        elinewidth=0.8,
        capsize=2.5,
        zorder=3,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(sample_order, rotation=35, ha="right")
    ax.set_ylabel("Mean Scanpy SSSP score")
    ax.set_title("Individual samples", fontsize=12, fontweight="bold")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    young_n = int((sample_summary["age"].astype(str) == "Young").sum())
    if 0 < young_n < len(sample_summary):
        ax.axvline(young_n - 0.5, color="#777777", linewidth=0.9, linestyle="--", alpha=0.7)
    handles = [
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=AGE_COLORS[age], markersize=10, label=age)
        for age in AGE_ORDER
    ]
    ax.legend(handles=handles, frameon=False, loc="upper left")

    fig.suptitle(plot_title, fontsize=15, fontweight="bold", y=1.02)
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    grouped = [
        cell_scores.loc[cell_scores["sample"].eq(sample), score_name].dropna().to_numpy()
        for sample in sample_order
    ]
    parts = ax.violinplot(grouped, positions=np.arange(len(sample_order)), widths=0.74, showmeans=False, showextrema=False)
    for body, color in zip(parts["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("#333333")
        body.set_alpha(0.55)
        body.set_linewidth(0.7)
    medians = [np.median(vals) if len(vals) else np.nan for vals in grouped]
    ax.scatter(np.arange(len(sample_order)), medians, color="#222222", s=18, zorder=3, label="cell median")
    ax.set_xticks(np.arange(len(sample_order)))
    ax.set_xticklabels(sample_order, rotation=35, ha="right")
    ax.set_ylabel("Cell-level Scanpy SSSP score")
    ax.set_title("Cell-level distributions by sample", fontsize=13, fontweight="bold")
    if 0 < young_n < len(sample_summary):
        ax.axvline(young_n - 0.5, color="#777777", linewidth=0.9, linestyle="--", alpha=0.7)
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outbase.with_name(outbase.name + "_cell_distributions").with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_name(outbase.name + "_cell_distributions").with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def score_gene_set(
    base_adata: ad.AnnData,
    input_audit: pd.DataFrame,
    metadata_source: str,
    config: dict[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    outdir = Path(config["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    genes = load_gene_list(str(config["membership_col"]))
    panel_genes = genes["gene"].tolist()
    score_name = str(config["score_name"])
    prefix = str(config["prefix"])
    adata = base_adata.copy()

    present_panel_genes = [gene for gene in panel_genes if gene in adata.var_names]
    missing = sorted(set(panel_genes) - set(present_panel_genes))
    if missing:
        raise RuntimeError(f"Missing expected panel genes from Xenium matrix: {', '.join(missing)}")

    sc.tl.score_genes(
        adata,
        gene_list=present_panel_genes,
        score_name=score_name,
        ctrl_size=50,
        n_bins=25,
        random_state=1,
        use_raw=False,
    )

    cell_scores = add_panel_expression_summaries(adata, present_panel_genes, score_name)
    cell_scores = add_available_annotations(cell_scores, adata)
    sample_summary = (
        cell_scores.groupby(["sample", "age"], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            mean_scanpy_score=(score_name, "mean"),
            median_scanpy_score=(score_name, "median"),
            sem_scanpy_score=(score_name, sem),
            mean_panel_log1p_cp10k=("panel_mean_log1p_cp10k", "mean"),
            mean_n_positive_panel_genes=("n_positive_panel_genes", "mean"),
        )
        .reset_index()
    )
    sample_summary["age_order"] = sample_summary["age"].map({"Young": 0, "Old": 1}).fillna(2)
    sample_summary = sample_summary.sort_values(["age_order", "sample"]).drop(columns="age_order")

    young = sample_summary.loc[sample_summary["age"].eq("Young"), "mean_scanpy_score"]
    old = sample_summary.loc[sample_summary["age"].eq("Old"), "mean_scanpy_score"]
    test = pd.DataFrame(
        [
            {
                "score_method": "scanpy.tl.score_genes",
                "gene_set_id": config["gene_set_id"],
                "score_name": score_name,
                "n_genes": len(present_panel_genes),
                "ctrl_size": 50,
                "n_bins": 25,
                "n_young_samples": int(len(young)),
                "n_old_samples": int(len(old)),
                "young_mean_of_sample_means": float(young.mean()),
                "old_mean_of_sample_means": float(old.mean()),
                "old_minus_young": float(old.mean() - young.mean()),
                "welch_ttest_p_sample_level": float(stats.ttest_ind(old, young, equal_var=False).pvalue)
                if len(old) >= 2 and len(young) >= 2
                else np.nan,
                "mannwhitney_p_sample_level": float(stats.mannwhitneyu(old, young, alternative="two-sided").pvalue)
                if len(old) >= 1 and len(young) >= 1
                else np.nan,
            }
        ]
    )

    cell_score_path = outdir / f"{prefix}_cell_scores.csv.gz"
    cell_scores.to_csv(cell_score_path, index=False)
    add_score_summaries_to_adata_obs(adata, cell_scores, score_name, str(config["high_gt0_col"]))
    adata.var[str(config["var_flag"])] = adata.var_names.isin(present_panel_genes)
    adata.uns[str(config["uns_key"])] = {
        "score_method": "scanpy.tl.score_genes",
        "gene_set_id": config["gene_set_id"],
        "score_name": score_name,
        "gene_list": present_panel_genes,
        "ctrl_size": 50,
        "n_bins": 25,
        "random_state": 1,
        "definition_high_gt0": f"{score_name} > 0",
        "metadata_source": metadata_source,
        "gene_config": str(GENE_CONFIG_PATH),
    }
    adata.write_h5ad(outdir / f"{prefix}_scored.h5ad", compression="gzip")
    sample_summary.to_csv(outdir / f"{prefix}_score_by_sample.csv", index=False)
    test.to_csv(outdir / f"{prefix}_young_old_sample_level_test.csv", index=False)
    genes.to_csv(outdir / f"{prefix}_genes.csv", index=False)
    audit = input_audit.copy()
    audit.insert(0, "metadata_source", metadata_source)
    audit.insert(1, "gene_set_id", config["gene_set_id"])
    audit.to_csv(outdir / f"{prefix}_input_audit.csv", index=False)
    write_plots(
        sample_summary,
        cell_scores,
        outdir / f"{prefix}_young_old_and_samples",
        score_name,
        str(config["plot_title"]),
    )

    print(f"\nScanpy SSSP score by sample: {config['gene_set_id']}")
    print(sample_summary.to_string(index=False))
    print("\nSample-level young vs old test:")
    print(test.to_string(index=False))
    print("\nInput audit:")
    print(audit.to_string(index=False))
    print(f"\nWrote outputs to: {outdir}")
    return sample_summary, test


def main() -> None:
    meta, metadata_source = load_cell_metadata()
    adata, input_audit = build_xenium_anndata(meta)
    adata, input_audit = add_spatial_coordinates(adata, input_audit)
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)

    all_tests = []
    for config in GENE_SET_CONFIGS:
        _, test = score_gene_set(adata, input_audit, metadata_source, config)
        all_tests.append(test)

    combined = pd.concat(all_tests, ignore_index=True)
    combined.to_csv(
        PROJECT_ROOT / "results" / "260727" / "xenium_sssp_scanpy_gene_set_comparison_tests.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
