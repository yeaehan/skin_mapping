"""
Plot 0-35 um Squidpy co-occurrence heatmaps for SnC status + minor cell type.

This uses the scored Xenium h5ad from step 16. The h5ad is expected to contain
spatial coordinates in ``adata.obsm["spatial"]``. Squidpy computes the
co-occurrence probability ratio:

    P(neighbor label within 35 um | anchor label) / P(neighbor label overall)

Values > 1 mean the row label is enriched near the column anchor label.
Values < 1 mean the row label is depleted near the column anchor label.

Co-occurrence is computed separately per sample, then averaged by age, so one
sample with more cells does not dominate the Young/Old heatmaps.
"""

from __future__ import annotations

from pathlib import Path
import sys
import types
import warnings

import anndata as ad
import anndata._core.views as anndata_views
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import csc_matrix, csr_matrix

# Squidpy 1.5 imports spatialdata at package import time, but co_occurrence only
# needs the AnnData path. This avoids unrelated spatialdata image-stack imports.
if "spatialdata" not in sys.modules:
    spatialdata_stub = types.ModuleType("spatialdata")
    spatialdata_stub.SpatialData = type("SpatialData", (), {})
    sys.modules["spatialdata"] = spatialdata_stub

# Squidpy 1.5 also imports sparse view classes that newer AnnData no longer
# exposes under these names. Co-occurrence does not touch expression matrices,
# so scipy sparse classes are enough for this compatibility path.
if not hasattr(anndata_views, "SparseCSRView"):
    anndata_views.SparseCSRView = csr_matrix
if not hasattr(anndata_views, "SparseCSCView"):
    anndata_views.SparseCSCView = csc_matrix

import squidpy as sq


PROJECT_ROOT = Path(__file__).resolve().parents[3]
SCORED_H5AD = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_sssp_20gene_scanpy_score_genes"
    / "xenium_20gene_sssp_scanpy_scored.h5ad"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_snc_minor_cooccurrence_0_35um"

CLUSTER_KEY = "snc_minor_label"
SPATIAL_KEY = "spatial"
RADIUS_UM = 35.0
INTERVAL_UM = np.array([0.0, RADIUS_UM])
AGE_ORDER = ["Young", "Old"]
MIN_ANCHOR_CELLS_TOTAL = 10
MIN_ROW_CELLS_TOTAL = 20
MAX_ROW_LABELS = 28


def load_spatial_adata() -> ad.AnnData:
    required = {"sample", "age", "snc", "cell_type_minor", "fov", "cell_id"}
    with h5py.File(SCORED_H5AD, "r") as handle:
        if SPATIAL_KEY not in handle["obsm"]:
            raise KeyError(
                f"{SCORED_H5AD} does not contain adata.obsm[{SPATIAL_KEY!r}]. "
                "Rerun step 16 after the spatial-coordinate update."
            )

        obs_group = handle["obs"]
        missing = sorted(required - set(obs_group.keys()))
        if missing:
            raise KeyError(f"Missing required obs columns in h5ad: {', '.join(missing)}")

        obs = pd.DataFrame({col: read_h5ad_obs_column(obs_group, col) for col in sorted(required)})
        spatial = handle["obsm"][SPATIAL_KEY][:]

    obs.index = obs["cell_id"].astype(str)
    adata = ad.AnnData(
        X=sparse.csr_matrix((len(obs), 0), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=[]),
    )
    adata.obsm[SPATIAL_KEY] = spatial.astype(float)
    adata.obs["snc"] = adata.obs["snc"].astype(str)
    adata.obs["cell_type_minor"] = adata.obs["cell_type_minor"].astype(str)
    adata = adata[adata.obs["snc"].isin(["SnC-", "SnC+"])].copy()
    adata.obs[CLUSTER_KEY] = adata.obs["snc"] + " | " + adata.obs["cell_type_minor"]
    adata.obs["age"] = pd.Categorical(adata.obs["age"], categories=AGE_ORDER, ordered=True)

    label_order = sorted(adata.obs[CLUSTER_KEY].unique())
    adata.obs[CLUSTER_KEY] = pd.Categorical(adata.obs[CLUSTER_KEY], categories=label_order)
    return adata


def decode_h5_values(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def read_h5ad_obs_column(obs_group: h5py.Group, col: str) -> pd.Series:
    obj = obs_group[col]
    if isinstance(obj, h5py.Dataset):
        values = obj[:]
        if values.dtype.kind in {"S", "O", "U"}:
            return pd.Series(decode_h5_values(values), name=col)
        return pd.Series(values, name=col)

    if isinstance(obj, h5py.Group) and obj.attrs.get("encoding-type") == "categorical":
        categories = decode_h5_values(obj["categories"][:])
        codes = obj["codes"][:]
        values = [categories[code] if code >= 0 else np.nan for code in codes]
        return pd.Series(pd.Categorical(values, categories=categories), name=col)

    raise TypeError(f"Unsupported h5ad obs column encoding for {col!r}")


def choose_plot_labels(adata: ad.AnnData) -> tuple[list[str], list[str]]:
    total_counts = adata.obs[CLUSTER_KEY].value_counts()

    anchor_labels = [
        label
        for label, count in total_counts.items()
        if str(label).startswith("SnC+ | ") and count >= MIN_ANCHOR_CELLS_TOTAL
    ]
    anchor_labels = sorted(anchor_labels, key=lambda label: (-total_counts[label], str(label)))

    row_labels = [
        label
        for label, count in total_counts.items()
        if count >= MIN_ROW_CELLS_TOTAL or label in anchor_labels
    ]
    row_labels = sorted(row_labels, key=lambda label: (-total_counts[label], str(label)))[:MAX_ROW_LABELS]
    for label in anchor_labels:
        if label not in row_labels:
            row_labels.append(label)
    return [str(x) for x in row_labels], [str(x) for x in anchor_labels]


def run_sample_cooccurrence(adata: ad.AnnData) -> pd.DataFrame:
    rows = []
    all_labels = list(adata.obs[CLUSTER_KEY].cat.categories.astype(str))

    for sample, sample_obs in adata.obs.groupby("sample", observed=True, sort=False):
        sample_idx = sample_obs.index
        sample_adata = adata[sample_idx].copy()
        sample_adata.obs[CLUSTER_KEY] = sample_adata.obs[CLUSTER_KEY].cat.remove_unused_categories()
        categories = list(sample_adata.obs[CLUSTER_KEY].cat.categories.astype(str))

        if sample_adata.n_obs < 2 or len(categories) < 2:
            continue

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            occ, interval = sq.gr.co_occurrence(
                sample_adata,
                cluster_key=CLUSTER_KEY,
                spatial_key=SPATIAL_KEY,
                interval=INTERVAL_UM,
                copy=True,
                n_splits=1,
                n_jobs=1,
                show_progress_bar=False,
            )

        if occ.shape[2] != 1:
            raise ValueError(f"Expected one 0-35 um interval, got occ shape {occ.shape} and interval {interval}")

        counts = sample_adata.obs[CLUSTER_KEY].value_counts().reindex(all_labels, fill_value=0)
        age = str(sample_adata.obs["age"].iloc[0])
        for anchor_label in all_labels:
            for neighbor_label in all_labels:
                if anchor_label in categories and neighbor_label in categories:
                    anchor_idx = categories.index(anchor_label)
                    neighbor_idx = categories.index(neighbor_label)
                    ratio = float(occ[anchor_idx, neighbor_idx, 0])
                else:
                    ratio = np.nan
                rows.append(
                    {
                        "sample": sample,
                        "age": age,
                        "neighbor_label": neighbor_label,
                        "anchor_label": anchor_label,
                        "radius_um": RADIUS_UM,
                        "cooccurrence_ratio": ratio,
                        "log2_cooccurrence_ratio": np.log2(ratio) if np.isfinite(ratio) and ratio > 0 else np.nan,
                        "neighbor_cells_in_sample": int(counts.get(neighbor_label, 0)),
                        "anchor_cells_in_sample": int(counts.get(anchor_label, 0)),
                    }
                )

    return pd.DataFrame(rows)


def summarize_by_age(sample_pairs: pd.DataFrame) -> pd.DataFrame:
    return (
        sample_pairs.groupby(["age", "neighbor_label", "anchor_label", "radius_um"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_cooccurrence_ratio=("cooccurrence_ratio", "mean"),
            mean_log2_cooccurrence_ratio=("log2_cooccurrence_ratio", "mean"),
            total_anchor_cells=("anchor_cells_in_sample", "sum"),
            total_neighbor_cells=("neighbor_cells_in_sample", "sum"),
        )
        .reset_index()
    )


def matrix_from_age_summary(
    age_summary: pd.DataFrame,
    age: str,
    row_labels: list[str],
    anchor_labels: list[str],
    value_col: str = "mean_log2_cooccurrence_ratio",
) -> pd.DataFrame:
    sub = age_summary[age_summary["age"].astype(str).eq(age)]
    mat = sub.pivot(index="neighbor_label", columns="anchor_label", values=value_col)
    return mat.reindex(index=row_labels, columns=anchor_labels)


def plot_heatmap(matrix: pd.DataFrame, title: str, outbase: Path) -> None:
    data = matrix.to_numpy(dtype=float)
    if data.size == 0:
        return
    vmax = np.nanmax(np.abs(data))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    vmax = min(max(vmax, 1.0), 3.0)

    fig_w = max(7.0, 0.55 * len(matrix.columns) + 4.5)
    fig_h = max(7.0, 0.32 * len(matrix.index) + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(data, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")

    ax.set_xticks(np.arange(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=8)
    ax.set_xlabel("Anchor cell label")
    ax.set_ylabel("Neighbor cell label within 35 um")
    ax.set_title(title, fontsize=12, fontweight="bold")

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(matrix.columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(matrix.index), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.6)
    ax.tick_params(which="minor", bottom=False, left=False)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Mean log2 co-occurrence ratio")
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    adata = load_spatial_adata()
    row_labels, anchor_labels = choose_plot_labels(adata)

    sample_pairs = run_sample_cooccurrence(adata)
    age_summary = summarize_by_age(sample_pairs)
    sample_pairs.to_csv(OUTDIR / "squidpy_snc_minor_cooccurrence_ratio_0_35um_by_sample.csv", index=False)
    age_summary.to_csv(OUTDIR / "squidpy_snc_minor_cooccurrence_ratio_0_35um_by_age.csv", index=False)

    label_counts = (
        adata.obs.groupby(["age", CLUSTER_KEY], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["age", "n_cells"], ascending=[True, False])
    )
    label_counts.to_csv(OUTDIR / "squidpy_snc_minor_label_counts_by_age.csv", index=False)

    young_mat = matrix_from_age_summary(age_summary, "Young", row_labels, anchor_labels)
    old_mat = matrix_from_age_summary(age_summary, "Old", row_labels, anchor_labels)
    plot_heatmap(
        young_mat,
        "Young: Squidpy 0-35 um co-occurrence near SnC+ minor cell anchors",
        OUTDIR / "young_squidpy_sncplus_anchor_snc_minor_cooccurrence_log2ratio_0_35um",
    )
    plot_heatmap(
        old_mat,
        "Old: Squidpy 0-35 um co-occurrence near SnC+ minor cell anchors",
        OUTDIR / "old_squidpy_sncplus_anchor_snc_minor_cooccurrence_log2ratio_0_35um",
    )

    run_summary = pd.DataFrame(
        [
            {
                "input_h5ad": str(SCORED_H5AD),
                "cluster_key": CLUSTER_KEY,
                "squidpy_function": "sq.gr.co_occurrence",
                "interval_um": "0,35",
                "n_cells": int(adata.n_obs),
                "n_samples": int(adata.obs["sample"].nunique()),
                "n_fovs": int(adata.obs[["sample", "fov"]].drop_duplicates().shape[0]),
                "n_combined_labels": int(adata.obs[CLUSTER_KEY].nunique()),
                "n_anchor_labels_plotted": int(len(anchor_labels)),
                "n_neighbor_labels_plotted": int(len(row_labels)),
                "anchor_label_filter": f"SnC+ labels with >= {MIN_ANCHOR_CELLS_TOTAL} cells total",
                "neighbor_label_filter": f"labels with >= {MIN_ROW_CELLS_TOTAL} cells total, plus plotted anchors",
                "summary_level": "Squidpy co-occurrence computed per sample, then averaged by age",
            }
        ]
    )
    run_summary.to_csv(OUTDIR / "squidpy_snc_minor_cooccurrence_0_35um_run_summary.csv", index=False)

    print("Wrote Squidpy 0-35 um SnC/minor co-occurrence heatmaps to:")
    print(OUTDIR)
    print(run_summary.to_string(index=False))


if __name__ == "__main__":
    main()
