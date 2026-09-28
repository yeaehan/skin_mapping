"""
Sinai public validation: SnC trends across 10-year age bins, all sites.

This workflow intentionally does not subset anatomical sites. It starts from the
full integrated Sinai AnnData object, recomputes rule-based SnC calls, bins
donor age by decade, and summarizes trends at the donor level before averaging
within each age bin.

Two SnC definitions are retained for comparison:

1. Original marker-positive SnC rule:
       SnC+ = (CDKN1A+ OR CDKN2A+) AND NOT (PCNA+ OR MKI67+)
              AND sasp_num_user >= SASP_THRESHOLD

2. Experimentalist component-quartile gate, calculated within cell lineage:
       growth arrest = Q4 (top 25%)
       SASP          = Q4 (top 25%)
       cell cycle    = Q1 (bottom 25%; least proliferative)
       SnC+          = growth Q4 AND SASP Q4 AND cell-cycle Q1

The continuous senescence index is then retained as:
       SI = Z_growth + Z_SASP - Z_cycle
where each component is z-scored within lineage. SI is descriptive/ranking
information and is NOT itself quartile-thresholded for the experimentalist call.

Expression positivity means expression > 0 in the raw count layer when a
`counts` layer is present. If no count layer is available, the script falls back
to X, which preserves zero/non-zero positivity for log-normalized data.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy import stats


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_H5AD = PROJECT_ROOT / "data" / "scrna.integrated_annotated_sinai.h5ad"
DEFAULT_OUTDIR = PROJECT_ROOT / "results" / "260727" / "sinai_age_decade_trends_all_sites"

SASP_THRESHOLD = 10
MIN_CELLS_PER_DONOR_GROUP = 20
MIN_CELLS_PER_LINEAGE = 50
TOP_N_CELL_TYPES = 12
UMAP_MAX_POINTS = 80000
TOP_N_SAMPLE_COMPOSITION_CELL_TYPES = 14
TOP_N_SNC_SAMPLE_CELL_TYPES = 12
SITE_FILTER_EXCLUDE_FOR_AGE_FOCUSED = {"face", "inguinal fold", "antecubital fossa"}

SNC_MARKER_GENES = ["CDKN1A", "CDKN2A", "PCNA", "MKI67", "CXCL10", "CXCR3", "BCL2", "HMGB1", "LMNB1"]

SASP_ORIGINAL = [
    "CCL11", "CCL13", "CCL14", "CCL15", "CCL16", "CCL19", "CCL2", "CCL20",
    "CCL21", "CCL26", "CCL28", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8",
    "CXCL1", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL14", "CXCL16", "CXCL2",
    "CXCL3", "CXCL5", "CXCL6", "CXCL9",
    "IL12B", "IL13", "IL15", "IL18", "IL1A", "IL1B", "IL2", "IL23A", "IL4", "IL5", "IL6",
    "MMP12", "MMP9", "MMP3",
    "VEGFA", "TNF", "TGFB1", "CSF1", "EGF", "FGF2", "IFNB1", "CMKLR1",
]

# User-supplied extended SASP panel from the Sinai SnC annotation notebook.
SASP_USER_EXTENDED = [
    "CXCL1", "CXCL10", "CXCL11", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL5",
    "CXCL6", "CXCL8", "CCL2", "CCL20", "CCL28", "CCL7", "IL11", "IL12A", "IL15",
    "IL16", "IL17D", "IL18", "IL18BP", "IL1A", "IL1B", "IL23A", "IL24", "IL32",
    "IL33", "IL6", "IL6ST", "IL7", "MMP1", "MMP10", "MMP11", "MMP12", "MMP14",
    "MMP15", "MMP16", "MMP17", "MMP19", "MMP22", "MMP24", "MMP24OS", "MMP25-AS1",
    "MMP3", "TNFAIP1", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFAIP8L1",
    "TNFAIP8L3", "TNFSF10", "TNFSF12", "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15",
    "TNFSF18", "TNFSF4", "VEGFA", "VEGFB", "VEGFC", "EGF", "PDGFA", "PDGFB",
    "GDF11", "GDF15", "GDF5", "GDF6", "FGF1", "FGF10", "FGF13", "FGF14", "FGF2",
    "FGF5", "FGF7", "FGF9",
]

GROWTH_ARREST_GENES = ["CDKN1A", "CDKN2A"]

S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1",
    "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1",
    "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]

G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
    "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
    "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1",
    "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5",
    "CENPA",
]


def canonical_gene_name(name: str) -> str:
    return str(name).upper()


def decode_h5_values(values) -> np.ndarray:
    values = np.asarray(values)
    if values.dtype.kind in {"S", "O"}:
        return np.array([
            x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else x
            for x in values
        ], dtype=object)
    return values


def read_h5ad_column(obj) -> np.ndarray:
    if isinstance(obj, h5py.Group) and {"categories", "codes"}.issubset(set(obj.keys())):
        categories = decode_h5_values(obj["categories"][:])
        codes = np.asarray(obj["codes"][:], dtype=int)
        out = np.full(codes.shape, pd.NA, dtype=object)
        ok = codes >= 0
        out[ok] = categories[codes[ok]]
        return out
    return decode_h5_values(obj[:])


def read_h5ad_obs(obs_group: h5py.Group) -> pd.DataFrame:
    index = read_h5ad_column(obs_group["_index"]) if "_index" in obs_group else np.arange(len(next(iter(obs_group.values()))))
    n_obs = len(index)
    obs = pd.DataFrame(index=pd.Index(index.astype(str), name="cell_id"))
    for key in obs_group.keys():
        if key == "_index":
            continue
        values = read_h5ad_column(obs_group[key])
        if len(values) == n_obs:
            obs[key] = values
    return obs


def read_var_names(var_group: h5py.Group) -> list[str]:
    if "_index" not in var_group:
        raise KeyError("Expected var/_index in the h5ad file.")
    return [str(x) for x in decode_h5_values(var_group["_index"][:])]


def build_gene_lookup(var_names: Iterable[str]) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for gene in var_names:
        lookup.setdefault(canonical_gene_name(gene), str(gene))
    return lookup


def present_genes(genes: Iterable[str], lookup: dict[str, str]) -> list[str]:
    present = []
    for gene in genes:
        mapped = lookup.get(canonical_gene_name(gene))
        if mapped is not None and mapped not in present:
            present.append(mapped)
    return present


def unique_preserve_order(values: Iterable[str]) -> list[str]:
    seen = set()
    out = []
    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def get_sparse_group(h5: h5py.File, preferred_path: str) -> tuple[h5py.Group, str]:
    if preferred_path in h5:
        return h5[preferred_path], preferred_path
    return h5["X"], "X"


def read_csr_columns_dense(
    csr_group: h5py.Group,
    col_indices: list[int],
    *,
    chunk_rows: int = 1000,
    dtype=np.float32,
) -> np.ndarray:
    shape = tuple(int(x) for x in csr_group.attrs["shape"])
    n_obs, n_vars = shape
    out = np.zeros((n_obs, len(col_indices)), dtype=dtype)
    if not col_indices:
        return out

    col_map = np.full(n_vars, -1, dtype=np.int32)
    col_map[np.asarray(col_indices, dtype=np.int32)] = np.arange(len(col_indices), dtype=np.int32)

    indptr = csr_group["indptr"][:]
    indices_ds = csr_group["indices"]
    data_ds = csr_group["data"]

    for row_start in range(0, n_obs, chunk_rows):
        row_end = min(row_start + chunk_rows, n_obs)
        ptr_start = int(indptr[row_start])
        ptr_end = int(indptr[row_end])
        if ptr_end <= ptr_start:
            continue

        indices = indices_ds[ptr_start:ptr_end]
        mapped_cols = col_map[indices]
        hit_mask = mapped_cols >= 0
        if not np.any(hit_mask):
            continue

        row_lengths = np.diff(indptr[row_start : row_end + 1]).astype(np.int64)
        local_rows = np.repeat(np.arange(row_end - row_start, dtype=np.int32), row_lengths)[hit_mask]
        values = data_ds[ptr_start:ptr_end][hit_mask].astype(dtype, copy=False)
        out[row_start + local_rows, mapped_cols[hit_mask]] = values

    return out


def add_age_decade(obs: pd.DataFrame) -> pd.DataFrame:
    obs = obs.copy()
    obs["donor_age_num"] = pd.to_numeric(obs["donor_age"], errors="coerce")
    bins = list(range(0, 101, 10)) + [np.inf]
    labels = [f"{start}-{start + 9}" for start in range(0, 100, 10)] + ["100+"]
    obs["age_decade"] = pd.cut(obs["donor_age_num"], bins=bins, labels=labels, right=False)
    return obs


def choose_first_existing(columns: list[str], candidates: list[str], label: str) -> str:
    for candidate in candidates:
        if candidate in columns:
            return candidate
    raise KeyError(f"Could not find a {label} column. Tried: {', '.join(candidates)}")


def sem(values: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if len(values) <= 1:
        return np.nan
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def zscore_within_group(values: pd.Series, groups: pd.Series, min_cells: int = MIN_CELLS_PER_LINEAGE) -> pd.Series:
    values = pd.Series(values, index=groups.index, dtype="float64")
    groups = groups.astype("string").fillna("Unknown")
    out = pd.Series(np.nan, index=values.index, dtype="float64")

    global_mean = values.mean(skipna=True)
    global_sd = values.std(skipna=True, ddof=0)
    if not np.isfinite(global_sd) or global_sd == 0:
        global_sd = 1.0

    for _, idx in groups.groupby(groups, observed=True).groups.items():
        idx = list(idx)
        subset = values.loc[idx]
        if len(idx) < min_cells:
            out.loc[idx] = (subset - global_mean) / global_sd
            continue
        mean = subset.mean(skipna=True)
        sd = subset.std(skipna=True, ddof=0)
        if not np.isfinite(sd) or sd == 0:
            out.loc[idx] = (subset - global_mean) / global_sd
        else:
            out.loc[idx] = (subset - mean) / sd
    return out


def percentile_call_within_group(
    values: pd.Series,
    groups: pd.Series,
    percentile: float,
    min_cells: int = MIN_CELLS_PER_LINEAGE,
) -> tuple[pd.Series, pd.DataFrame]:
    values = pd.Series(values, index=groups.index, dtype="float64")
    groups = groups.astype("string").fillna("Unknown")
    out = pd.Series(False, index=values.index)
    rows = []
    global_threshold = values.quantile(percentile)

    for group, idx in groups.groupby(groups, observed=True).groups.items():
        idx = list(idx)
        subset = values.loc[idx]
        threshold = global_threshold if len(idx) < min_cells else subset.quantile(percentile)
        out.loc[idx] = subset >= threshold
        rows.append({"quartile_lineage": group, "percentile": percentile, "threshold": threshold, "n_cells": len(idx)})
    return out, pd.DataFrame(rows)


def quartile_within_group(
    values: pd.Series,
    groups: pd.Series,
    min_cells: int = MIN_CELLS_PER_LINEAGE,
) -> pd.Categorical:

    values = pd.Series(
        values,
        index=groups.index,
        dtype="float64",
    )

    groups = (
        groups
        .astype("string")
        .fillna("Unknown")
    )

    out = pd.Series(
        pd.NA,
        index=values.index,
        dtype="object",
    )

    global_edges = values.quantile(
        [0.25, 0.50, 0.75]
    ).to_numpy()


    def assign(
        vals: pd.Series,
        edges: np.ndarray,
    ) -> pd.Series:

        q25, q50, q75 = edges

        result = pd.Series(
            pd.NA,
            index=vals.index,
            dtype="object",
        )

        valid = vals.notna()

        result.loc[
            valid & (vals <= q25)
        ] = "Q1"

        result.loc[
            valid & (vals > q25) & (vals <= q50)
        ] = "Q2"

        result.loc[
            valid & (vals > q50) & (vals <= q75)
        ] = "Q3"

        result.loc[
            valid & (vals > q75)
        ] = "Q4"

        return result


    for _, idx in groups.groupby(
        groups,
        observed=True,
    ).groups.items():

        idx = list(idx)

        subset = values.loc[idx]

        if len(idx) < min_cells:
            edges = global_edges
        else:
            edges = subset.quantile(
                [0.25, 0.50, 0.75]
            ).to_numpy()

        out.loc[idx] = assign(
            subset,
            edges,
        )


    return pd.Categorical(
        out,
        categories=[
            "Q1",
            "Q2",
            "Q3",
            "Q4",
        ],
        ordered=True,
    )


def summarize_by_decade(donor_summary: pd.DataFrame, value_cols: list[str]) -> pd.DataFrame:
    pieces = []
    for value_col in value_cols:
        tmp = (
            donor_summary.groupby("age_decade", observed=True)[value_col]
            .agg(mean="mean", median="median", sd="std", n_donors="count", sem=sem)
            .reset_index()
        )
        tmp.insert(1, "metric", value_col)
        pieces.append(tmp)
    return pd.concat(pieces, ignore_index=True)


def save_line_with_points(
    donor_summary: pd.DataFrame,
    decade_summary: pd.DataFrame,
    metric: str,
    ylabel: str,
    outpath: Path,
) -> None:
    plot_df = decade_summary.query("metric == @metric").copy()
    order = list(plot_df["age_decade"].astype(str))
    x_lookup = {label: idx for idx, label in enumerate(order)}

    fig, ax = plt.subplots(figsize=(8, 4.8))
    x = np.arange(len(order))
    ax.errorbar(x, plot_df["mean"], yerr=plot_df["sem"], marker="o", color="#2f5d7c", linewidth=2)

    points = donor_summary.dropna(subset=[metric]).copy()
    points = points[points["age_decade"].astype(str).isin(x_lookup)]
    jitter = np.random.default_rng(11).normal(0, 0.045, size=len(points))
    ax.scatter(
        points["age_decade"].astype(str).map(x_lookup).to_numpy() + jitter,
        points[metric],
        s=28,
        color="#f28e2b",
        alpha=0.7,
        edgecolor="white",
        linewidth=0.4,
        zorder=3,
    )

    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=35, ha="right")
    ax.set_xlabel("Donor age decade")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def summarize_call_methods_by_donor(obs: pd.DataFrame, call_cols: dict[str, str]) -> pd.DataFrame:
    rows = []
    for method, col in call_cols.items():
        tmp = (
            obs.groupby(["donor_id_for_summary", "age_decade"], observed=True)
            .agg(
                donor_age=("donor_age_num", "first"),
                n_cells=(col, "size"),
                n_called=(col, "sum"),
                pct_called=(col, lambda x: 100 * x.mean()),
                mean_sasp_num_user=("sasp_num_user", "mean"),
                mean_senescence_index=("senescence_index", "mean"),
            )
            .reset_index()
        )
        tmp.insert(2, "method", method)
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True)


def summarize_call_methods_by_sample(obs: pd.DataFrame, call_cols: dict[str, str]) -> pd.DataFrame:
    rows = []
    group_cols = ["sample_for_summary", "donor_id_for_summary", "age_decade", "site_for_summary"]
    for method, col in call_cols.items():
        tmp = (
            obs.groupby(group_cols, observed=True)
            .agg(
                donor_age=("donor_age_num", "first"),
                n_cells=(col, "size"),
                n_called=(col, "sum"),
                pct_called=(col, lambda x: 100 * x.mean()),
                mean_sasp_num_user=("sasp_num_user", "mean"),
                mean_senescence_index=("senescence_index", "mean"),
            )
            .reset_index()
        )
        tmp.insert(len(group_cols), "method", method)
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True)


def summarize_method_decade_from_donor(donor_method_summary: pd.DataFrame) -> pd.DataFrame:
    return (
        donor_method_summary.groupby(["method", "age_decade"], observed=True)
        .agg(
            mean_pct_called=("pct_called", "mean"),
            median_pct_called=("pct_called", "median"),
            sem_pct_called=("pct_called", sem),
            n_donors=("donor_id_for_summary", "nunique"),
            total_cells=("n_cells", "sum"),
            total_called=("n_called", "sum"),
        )
        .reset_index()
    )


def method_age_correlations(sample_method_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for method, frame in sample_method_summary.groupby("method", observed=True):
        valid = frame.dropna(subset=["donor_age", "pct_called"])
        if len(valid) >= 3:
            pearson_result = stats.pearsonr(valid["donor_age"], valid["pct_called"])
            spearman_result = stats.spearmanr(valid["donor_age"], valid["pct_called"])
            slope_result = stats.linregress(valid["donor_age"], valid["pct_called"])
            pearson = pearson_result.statistic
            pearson_p = pearson_result.pvalue
            spearman = spearman_result.statistic
            spearman_p = spearman_result.pvalue
            slope_pct_per_year = slope_result.slope
            slope_p = slope_result.pvalue
        else:
            pearson = np.nan
            pearson_p = np.nan
            spearman = np.nan
            spearman_p = np.nan
            slope_pct_per_year = np.nan
            slope_p = np.nan
        rows.append(
            {
                "method": method,
                "n_samples": len(valid),
                "pearson_age_vs_pct_called": pearson,
                "pearson_p": pearson_p,
                "spearman_age_vs_pct_called": spearman,
                "spearman_p": spearman_p,
                "linear_slope_pct_points_per_year": slope_pct_per_year,
                "linear_slope_p": slope_p,
            }
        )
    return pd.DataFrame(rows)


def save_method_comparison_plot(
    donor_method_summary: pd.DataFrame,
    decade_method_summary: pd.DataFrame,
    outpath: Path,
    title: str,
) -> None:
    methods = list(decade_method_summary["method"].drop_duplicates())
    decades = list(decade_method_summary["age_decade"].astype(str).drop_duplicates())
    decade_order = sorted(decades, key=lambda x: (int(x.split("-")[0]) if x[0].isdigit() else 999, x))
    x_lookup = {label: i for i, label in enumerate(decade_order)}
    offsets = np.linspace(-0.22, 0.22, max(len(methods), 1))
    colors = plt.get_cmap("tab10", max(len(methods), 1))

    fig, ax = plt.subplots(figsize=(10, 5.6))
    for idx, method in enumerate(methods):
        summary = decade_method_summary[decade_method_summary["method"] == method].copy()
        summary["x"] = summary["age_decade"].astype(str).map(x_lookup)
        ax.errorbar(
            summary["x"] + offsets[idx],
            summary["mean_pct_called"],
            yerr=summary["sem_pct_called"],
            marker="o",
            linewidth=2,
            color=colors(idx),
            label=method,
        )
        points = donor_method_summary[donor_method_summary["method"] == method].copy()
        points = points[points["age_decade"].astype(str).isin(x_lookup)]
        jitter = np.random.default_rng(23 + idx).normal(0, 0.025, size=len(points))
        ax.scatter(
            points["age_decade"].astype(str).map(x_lookup).to_numpy() + offsets[idx] + jitter,
            points["pct_called"],
            s=18,
            alpha=0.45,
            color=colors(idx),
            edgecolor="white",
            linewidth=0.3,
        )

    ax.set_xticks(range(len(decade_order)))
    ax.set_xticklabels(decade_order, rotation=35, ha="right")
    ax.set_xlabel("Donor age decade")
    ax.set_ylabel("Called SnC / SI-high cells (%)")
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def save_bar_plot(df: pd.DataFrame, x_col: str, y_col: str, fill_col: str, ylabel: str, outpath: Path) -> None:
    pivot = df.pivot_table(index=x_col, columns=fill_col, values=y_col, fill_value=0, observed=True)
    pivot.index = pivot.index.astype(str)
    pivot = pivot.loc[sorted(pivot.index, key=lambda x: (int(x.split("-")[0]) if x[0].isdigit() else 999, x))]

    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    pivot.plot(kind="bar", stacked=True, ax=ax, width=0.8, colormap="tab20")
    ax.set_xlabel("Donor age decade")
    ax.set_ylabel(ylabel)
    ax.legend(title=fill_col, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    ax.grid(axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def save_heatmap(df: pd.DataFrame, value_col: str, row_col: str, col_col: str, outpath: Path, title: str) -> None:
    pivot = df.pivot_table(index=row_col, columns=col_col, values=value_col, observed=True)
    if pivot.empty:
        return
    fig_height = max(4.5, 0.35 * len(pivot))
    fig, ax = plt.subplots(figsize=(9.5, fig_height))
    im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns.astype(str), rotation=35, ha="right")
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index.astype(str))
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label(value_col)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def save_umap_plot(
    obs: pd.DataFrame,
    color_col: str,
    title: str,
    outpath: Path,
    *,
    max_points: int = UMAP_MAX_POINTS,
) -> None:
    plot_df = obs.dropna(subset=["umap_1", "umap_2", color_col]).copy()
    if plot_df.empty:
        return
    if len(plot_df) > max_points:
        plot_df = plot_df.sample(n=max_points, random_state=17)

    categories = pd.Series(plot_df[color_col].astype("string")).fillna("Unknown")
    levels = sorted(categories.unique(), key=lambda x: str(x))
    cmap = plt.get_cmap("tab20", max(len(levels), 1))
    color_lookup = {level: cmap(i % cmap.N) for i, level in enumerate(levels)}

    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    for level in levels:
        mask = categories == level
        ax.scatter(
            plot_df.loc[mask, "umap_1"],
            plot_df.loc[mask, "umap_2"],
            s=2,
            alpha=0.55,
            linewidths=0,
            color=color_lookup[level],
            label=str(level),
        )
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        title=color_col,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        markerscale=4,
        fontsize=8,
        title_fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def sanitize_filename(value: str) -> str:
    keep = []
    for char in str(value):
        if char.isalnum():
            keep.append(char)
        elif char in {" ", "_", "-"}:
            keep.append("_")
    out = "".join(keep).strip("_")
    return out or "unknown"


def prepare_sample_feature_matrix(
    sample_summary: pd.DataFrame,
    sample_celltype_composition: pd.DataFrame,
    top_n_celltypes: int = TOP_N_SAMPLE_COMPOSITION_CELL_TYPES,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    top_celltypes = (
        sample_celltype_composition.groupby("cell_type_for_summary", observed=True)["n_cells_celltype"]
        .sum()
        .sort_values(ascending=False)
        .head(top_n_celltypes)
        .index
    )
    comp = sample_celltype_composition.copy()
    comp["cell_type_feature"] = np.where(
        comp["cell_type_for_summary"].isin(top_celltypes),
        comp["cell_type_for_summary"].astype(str),
        "Other",
    )
    comp = (
        comp.groupby(["sample_unit_id", "cell_type_feature"], observed=True)["cell_type_fraction"]
        .sum()
        .reset_index()
    )
    wide = comp.pivot_table(
        index="sample_unit_id",
        columns="cell_type_feature",
        values="cell_type_fraction",
        fill_value=0,
        observed=True,
    )
    metric_cols = ["snc_prop_percent", "mean_sasp_num_user", "mean_sasp_avg_expr_user"]
    metrics = sample_summary.set_index("sample_unit_id")[metric_cols].copy()
    metrics = metrics.apply(pd.to_numeric, errors="coerce")
    metrics = (metrics - metrics.mean()) / metrics.std(ddof=0).replace(0, 1)
    features = pd.concat([wide, metrics.add_prefix("metric_")], axis=1).fillna(0)

    metadata_cols = [
        "sample_unit_id",
        "sample_for_summary",
        "donor_id_for_summary",
        "age_decade",
        "site_for_summary",
        "n_cells",
        "n_snc",
        "snc_prop_percent",
        "mean_sasp_num_user",
    ]
    metadata = sample_summary[metadata_cols].drop_duplicates("sample_unit_id").set_index("sample_unit_id")
    metadata["sample_label"] = (
        metadata["sample_for_summary"].astype(str)
        + " | "
        + metadata["site_for_summary"].astype(str)
        + " | "
        + metadata["age_decade"].astype(str)
    )
    return features, metadata


def compute_sample_similarity_matrix(features: pd.DataFrame) -> pd.DataFrame:
    if features.shape[0] < 2:
        return pd.DataFrame()

    distances = pdist(features.to_numpy(dtype=float), metric="cosine")
    distances = np.nan_to_num(distances, nan=1.0, posinf=1.0, neginf=1.0)
    similarity = 1 - squareform(distances)
    np.fill_diagonal(similarity, 1.0)
    return pd.DataFrame(similarity, index=features.index, columns=features.index)


def save_sample_pca_plot(
    features: pd.DataFrame,
    metadata: pd.DataFrame,
    color_col: str,
    outpath: Path,
    title: str,
) -> pd.DataFrame:
    if features.shape[0] < 3:
        return pd.DataFrame()

    x = features.to_numpy(dtype=float)
    x = x - x.mean(axis=0, keepdims=True)
    sd = x.std(axis=0, keepdims=True)
    sd[sd == 0] = 1
    x = x / sd
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    pcs = u[:, :2] * s[:2]
    var = (s ** 2) / np.sum(s ** 2)
    pc_df = metadata.copy()
    pc_df["PC1"] = pcs[:, 0]
    pc_df["PC2"] = pcs[:, 1]

    categories = pc_df[color_col].astype("string").fillna("Unknown")
    levels = sorted(categories.unique(), key=lambda x: str(x))
    cmap = plt.get_cmap("tab20", max(len(levels), 1))

    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    for i, level in enumerate(levels):
        mask = categories == level
        ax.scatter(
            pc_df.loc[mask, "PC1"],
            pc_df.loc[mask, "PC2"],
            s=58,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.5,
            color=cmap(i % cmap.N),
            label=str(level),
        )
    for sample_id, row in pc_df.iterrows():
        ax.text(row["PC1"], row["PC2"], str(row["sample_for_summary"]), fontsize=5, alpha=0.75)
    ax.set_title(title)
    ax.set_xlabel(f"PC1 ({100 * var[0]:.1f}% feature variance)")
    ax.set_ylabel(f"PC2 ({100 * var[1]:.1f}% feature variance)")
    ax.grid(alpha=0.2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(title=color_col, bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)
    return pc_df.reset_index()


def save_sample_composition_stacked(
    sample_celltype_composition: pd.DataFrame,
    sample_summary: pd.DataFrame,
    outpath: Path,
    top_n_celltypes: int = TOP_N_SAMPLE_COMPOSITION_CELL_TYPES,
) -> None:
    top_celltypes = (
        sample_celltype_composition.groupby("cell_type_for_summary", observed=True)["n_cells_celltype"]
        .sum()
        .sort_values(ascending=False)
        .head(top_n_celltypes)
        .index
    )
    comp = sample_celltype_composition.copy()
    comp["cell_type_feature"] = np.where(
        comp["cell_type_for_summary"].isin(top_celltypes),
        comp["cell_type_for_summary"].astype(str),
        "Other",
    )
    comp = (
        comp.groupby(["sample_unit_id", "cell_type_feature"], observed=True)["cell_type_fraction"]
        .sum()
        .reset_index()
    )
    wide = comp.pivot_table(
        index="sample_unit_id",
        columns="cell_type_feature",
        values="cell_type_fraction",
        fill_value=0,
        observed=True,
    )
    metadata = sample_summary.drop_duplicates("sample_unit_id").set_index("sample_unit_id")
    metadata["sort_key"] = (
        metadata["site_for_summary"].astype(str)
        + " | "
        + metadata["age_decade"].astype(str)
        + " | "
        + metadata["sample_for_summary"].astype(str)
    )
    ordered_ids = metadata.sort_values("sort_key").index
    wide = wide.reindex(ordered_ids).fillna(0)
    labels = (
        metadata.loc[ordered_ids, "sample_for_summary"].astype(str)
        + "\n"
        + metadata.loc[ordered_ids, "site_for_summary"].astype(str)
        + " / "
        + metadata.loc[ordered_ids, "age_decade"].astype(str)
    )

    fig_width = max(12, 0.34 * len(ordered_ids))
    fig, ax = plt.subplots(figsize=(fig_width, 6.2))
    wide.plot(kind="bar", stacked=True, ax=ax, width=0.82, colormap="tab20")
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_xlabel("Sample / anatomical site / age decade")
    ax.set_ylabel("Cell-type fraction")
    ax.legend(title="cell type", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def save_snc_sample_composition_stacked(
    sample_snc_celltype_composition: pd.DataFrame,
    sample_summary: pd.DataFrame,
    outpath: Path,
    *,
    top_n_celltypes: int = TOP_N_SNC_SAMPLE_CELL_TYPES,
) -> None:
    plot_source = sample_snc_celltype_composition.copy()
    if plot_source.empty:
        return

    top_celltypes = (
        plot_source.groupby("cell_type_for_summary", observed=True)["n_snc_celltype"]
        .sum()
        .sort_values(ascending=False)
        .head(top_n_celltypes)
        .index
    )
    plot_source["cell_type_feature"] = np.where(
        plot_source["cell_type_for_summary"].isin(top_celltypes),
        plot_source["cell_type_for_summary"].astype(str),
        "Other",
    )
    plot_source = (
        plot_source.groupby(["sample_unit_id", "cell_type_feature"], observed=True)["snc_celltype_fraction"]
        .sum()
        .reset_index()
    )
    wide = plot_source.pivot_table(
        index="sample_unit_id",
        columns="cell_type_feature",
        values="snc_celltype_fraction",
        fill_value=0,
        observed=True,
    )

    metadata = sample_summary.drop_duplicates("sample_unit_id").set_index("sample_unit_id")
    metadata = metadata.loc[metadata.index.intersection(wide.index)].copy()
    metadata["sort_key"] = (
        metadata["age_decade"].astype(str)
        + " | "
        + metadata["site_for_summary"].astype(str)
        + " | "
        + metadata["sample_for_summary"].astype(str)
    )
    ordered_ids = metadata.sort_values("sort_key").index
    wide = wide.reindex(ordered_ids).fillna(0)
    labels = (
        metadata.loc[ordered_ids, "sample_for_summary"].astype(str)
        + "\n"
        + metadata.loc[ordered_ids, "site_for_summary"].astype(str)
        + " / "
        + metadata.loc[ordered_ids, "age_decade"].astype(str)
        + "\nSnC n="
        + metadata.loc[ordered_ids, "n_snc"].astype(int).astype(str)
    )

    fig_width = max(14, 0.45 * len(ordered_ids))
    fig, ax = plt.subplots(figsize=(fig_width, 6.8))
    wide.plot(kind="bar", stacked=True, ax=ax, width=0.82, colormap="tab20")
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_xlabel("Sample / anatomical site / age decade")
    ax.set_ylabel("Fraction of SnC+ cells within sample")
    ax.set_ylim(0, 1.02)
    ax.legend(title="SnC+ cell type", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False, fontsize=8)
    ax.grid(axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def save_sample_snc_proportion_bar(
    sample_summary: pd.DataFrame,
    outpath: Path,
    title: str,
) -> None:
    plot_df = sample_summary.copy()
    if plot_df.empty:
        return
    plot_df["sort_key"] = (
        plot_df["age_decade"].astype(str)
        + " | "
        + plot_df["site_for_summary"].astype(str)
        + " | "
        + plot_df["sample_for_summary"].astype(str)
    )
    plot_df = plot_df.sort_values("sort_key").reset_index(drop=True)

    sites = sorted(plot_df["site_for_summary"].astype(str).unique())
    cmap = plt.get_cmap("tab20", max(len(sites), 1))
    color_lookup = {site: cmap(i % cmap.N) for i, site in enumerate(sites)}
    colors = plot_df["site_for_summary"].astype(str).map(color_lookup)
    labels = (
        plot_df["sample_for_summary"].astype(str)
        + "\n"
        + plot_df["site_for_summary"].astype(str)
        + " / "
        + plot_df["age_decade"].astype(str)
        + "\nSnC n="
        + plot_df["n_snc"].astype(int).astype(str)
    )

    fig_width = max(13, 0.42 * len(plot_df))
    fig, ax = plt.subplots(figsize=(fig_width, 6.2))
    ax.bar(np.arange(len(plot_df)), plot_df["snc_prop_percent"], color=colors, width=0.82)
    ax.set_xticks(np.arange(len(plot_df)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6)
    ax.set_xlabel("Sample / anatomical site / age decade")
    ax.set_ylabel("SnC+ cells (%)")
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)

    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", color=color_lookup[site], label=site, markersize=8)
        for site in sites
    ]
    ax.legend(handles=handles, title="Anatomical site", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", type=Path, default=DEFAULT_INPUT_H5AD)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--sasp-threshold", type=int, default=SASP_THRESHOLD)
    parser.add_argument("--cell-type-col", default=None)
    args = parser.parse_args()

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Reading: {args.input_h5ad}")
    with h5py.File(args.input_h5ad, "r") as h5:
        var_names = read_var_names(h5["var"])
        obs = read_h5ad_obs(h5["obs"])
        n_obs = obs.shape[0]
        print(f"AnnData-like h5ad with n_obs x n_vars = {n_obs} x {len(var_names)}")
        if "obsm" in h5 and "umap" in h5["obsm"]:
            umap = np.asarray(h5["obsm"]["umap"][:], dtype=np.float32)
            obs["umap_1"] = umap[:, 0]
            obs["umap_2"] = umap[:, 1]
        else:
            print("No obsm/umap coordinates found; UMAP plots will be skipped.")

        count_group, count_layer_name = get_sparse_group(h5, "layers/counts")
        expr_group, expr_layer_name = get_sparse_group(h5, "layers/log1p")
        print(f"Positivity layer: {count_layer_name}")
        print(f"SASP average-expression layer: {expr_layer_name}")

        lookup = build_gene_lookup(var_names)
        gene_to_position = {gene: i for i, gene in enumerate(var_names)}
        sasp_original_present = present_genes(SASP_ORIGINAL, lookup)
        sasp_user_present = present_genes(SASP_USER_EXTENDED, lookup)
        marker_present = present_genes(SNC_MARKER_GENES, lookup)
        growth_present = present_genes(GROWTH_ARREST_GENES, lookup)
        s_phase_present = present_genes(S_GENES, lookup)
        g2m_present = present_genes(G2M_GENES, lookup)

        print(f"Original SASP genes present: {len(sasp_original_present)}/{len(SASP_ORIGINAL)}")
        print(f"User-extended SASP genes present: {len(sasp_user_present)}/{len(SASP_USER_EXTENDED)}")
        print(f"SnC marker genes present: {len(marker_present)}/{len(SNC_MARKER_GENES)}")
        print(f"Quartile growth genes present: {len(growth_present)}/{len(GROWTH_ARREST_GENES)}")
        print(f"Quartile S-phase genes present: {len(s_phase_present)}/{len(S_GENES)}")
        print(f"Quartile G2M genes present: {len(g2m_present)}/{len(G2M_GENES)}")

        count_genes = unique_preserve_order(marker_present + sasp_original_present + sasp_user_present)
        count_gene_to_matrix_col = {gene: i for i, gene in enumerate(count_genes)}
        count_cols = [gene_to_position[gene] for gene in count_genes]
        print(f"Reading {len(count_cols)} selected columns from {count_layer_name}")
        count_matrix = read_csr_columns_dense(count_group, count_cols)

        expr_genes = unique_preserve_order(
            sasp_original_present
            + sasp_user_present
            + growth_present
            + s_phase_present
            + g2m_present
        )
        expr_gene_to_matrix_col = {gene: i for i, gene in enumerate(expr_genes)}
        expr_cols = [gene_to_position[gene] for gene in expr_genes]
        print(f"Reading {len(expr_cols)} selected columns from {expr_layer_name}")
        expr_matrix = read_csr_columns_dense(expr_group, expr_cols)

    obs = add_age_decade(obs)
    donor_col = choose_first_existing(list(obs.columns), ["donor_id", "sample_barcode", "orig.ident"], "donor")
    sample_col = choose_first_existing(list(obs.columns), ["sample_barcode", "orig.ident", "study_id"], "sample")
    site_col = choose_first_existing(list(obs.columns), ["anatomic_site", "anatomic_site.detailed"], "anatomic site")

    if args.cell_type_col is not None:
        cell_type_col = args.cell_type_col
    else:
        cell_type_col = choose_first_existing(
            list(obs.columns),
            ["cell_type.detailed", "cell_type.reclustered", "cell_type", "cell_type.broad", "reported.cell_type"],
            "cell type",
        )
    broad_cell_type_col = "cell_type.broad" if "cell_type.broad" in obs.columns else cell_type_col

    for gene in SNC_MARKER_GENES:
        mapped = lookup.get(canonical_gene_name(gene))
        if mapped is None:
            obs[f"{gene}_pos"] = False
        else:
            obs[f"{gene}_pos"] = count_matrix[:, count_gene_to_matrix_col[mapped]] > 0

    original_count_cols = [count_gene_to_matrix_col[gene] for gene in sasp_original_present]
    user_count_cols = [count_gene_to_matrix_col[gene] for gene in sasp_user_present]
    original_expr_cols = [expr_gene_to_matrix_col[gene] for gene in sasp_original_present]
    user_expr_cols = [expr_gene_to_matrix_col[gene] for gene in sasp_user_present]

    obs["sasp_num_original"] = (count_matrix[:, original_count_cols] > 0).sum(axis=1).astype("int32")
    obs["sasp_num_user"] = (count_matrix[:, user_count_cols] > 0).sum(axis=1).astype("int32")
    obs["sasp_avg_expr_original"] = expr_matrix[:, original_expr_cols].mean(axis=1) if original_expr_cols else np.nan
    obs["sasp_avg_expr_user"] = expr_matrix[:, user_expr_cols].mean(axis=1) if user_expr_cols else np.nan

    base_snc = (
        (obs["CDKN1A_pos"] | obs["CDKN2A_pos"])
        & ~(obs["PCNA_pos"] | obs["MKI67_pos"])
    )
    obs["snc"] = np.where(base_snc & (obs["sasp_num_user"] >= args.sasp_threshold), "SnC+", "SnC-")
    obs["snc_bool"] = obs["snc"].eq("SnC+")

    obs["donor_id_for_summary"] = obs[donor_col].astype("string").fillna("Unknown")
    obs["sample_for_summary"] = obs[sample_col].astype("string").fillna("Unknown")
    obs["site_for_summary"] = obs[site_col].astype("string").fillna("Unknown")
    obs["cell_type_for_summary"] = obs[cell_type_col].astype("string").fillna("Unknown")
    obs["broad_cell_type_for_summary"] = obs[broad_cell_type_col].astype("string").fillna("Unknown")
    obs["age_decade_for_plot"] = obs["age_decade"].astype("string").fillna("Unknown")
    obs["quartile_lineage"] = obs["cell_type_for_summary"].astype("string").fillna("Unknown")

    growth_expr_cols = [expr_gene_to_matrix_col[gene] for gene in growth_present]
    s_phase_expr_cols = [expr_gene_to_matrix_col[gene] for gene in s_phase_present]
    g2m_expr_cols = [expr_gene_to_matrix_col[gene] for gene in g2m_present]
    obs["growth_raw"] = expr_matrix[:, growth_expr_cols].mean(axis=1) if growth_expr_cols else np.nan
    obs["S_score_custom"] = expr_matrix[:, s_phase_expr_cols].mean(axis=1) if s_phase_expr_cols else np.nan
    obs["G2M_score_custom"] = expr_matrix[:, g2m_expr_cols].mean(axis=1) if g2m_expr_cols else np.nan
    obs["cycle_raw"] = obs[["S_score_custom", "G2M_score_custom"]].mean(axis=1)
    obs["sasp_raw"] = obs["sasp_avg_expr_user"]

    # ----------------------------------------------------------------------
    # Experimentalist component-quartile gate
    #
    # Quartiles are calculated separately within each cell lineage.
    # A candidate senescent cell must simultaneously have:
    #   - top 25% growth-arrest score      (growth Q4)
    #   - top 25% SASP score               (SASP Q4)
    #   - bottom 25% cell-cycle score       (cycle Q1)
    #
    # We deliberately do NOT first sum the components and then take Q4 of SI.
    # ----------------------------------------------------------------------
    obs["growth_quartile_within_lineage"] = quartile_within_group(
        obs["growth_raw"], obs["quartile_lineage"]
    )
    obs["sasp_quartile_within_lineage"] = quartile_within_group(
        obs["sasp_raw"], obs["quartile_lineage"]
    )
    obs["cycle_quartile_within_lineage"] = quartile_within_group(
        obs["cycle_raw"], obs["quartile_lineage"]
    )

    obs["growth_Q4_within_lineage"] = (
        obs["growth_quartile_within_lineage"].astype("string").eq("Q4")
    )
    obs["sasp_Q4_within_lineage"] = (
        obs["sasp_quartile_within_lineage"].astype("string").eq("Q4")
    )
    obs["cycle_Q1_within_lineage"] = (
        obs["cycle_quartile_within_lineage"].astype("string").eq("Q1")
    )

    obs["snc_component_Q4Q4Q1_within_lineage"] = (
        obs["growth_Q4_within_lineage"]
        & obs["sasp_Q4_within_lineage"]
        & obs["cycle_Q1_within_lineage"]
    )

    # Continuous SI is calculated after defining the component-wise gate.
    # It remains useful for ranking/description, but it does not determine
    # whether a cell passes the experimentalist Q4/Q4/Q1 SnC call.
    obs["Z_growth"] = zscore_within_group(obs["growth_raw"], obs["quartile_lineage"])
    obs["Z_cycle"] = zscore_within_group(obs["cycle_raw"], obs["quartile_lineage"])
    obs["Z_cycle_suppression"] = -obs["Z_cycle"]
    obs["Z_SASP"] = zscore_within_group(obs["sasp_raw"], obs["quartile_lineage"])
    obs["senescence_index"] = obs["Z_growth"] + obs["Z_SASP"] - obs["Z_cycle"]

    # SI restricted to cells that pass all three component quartile criteria.
    # This makes it easy to inspect/rank the selected cells without changing
    # the gate itself.
    obs["senescence_index_component_gate"] = obs["senescence_index"].where(
        obs["snc_component_Q4Q4Q1_within_lineage"]
    )

    obs = obs.dropna(subset=["donor_age_num", "age_decade"])

    call_method_cols = {
        "Original marker-positive SnC": "snc_bool",
        "Component Q4 growth + Q4 SASP + Q1 cycle": "snc_component_Q4Q4Q1_within_lineage",
    }
    obs_age_focused = obs[
        ~obs["site_for_summary"].astype("string").str.lower().isin(SITE_FILTER_EXCLUDE_FOR_AGE_FOCUSED)
    ].copy()

    method_donor_all = summarize_call_methods_by_donor(obs, call_method_cols)
    method_sample_all = summarize_call_methods_by_sample(obs, call_method_cols)
    method_decade_all = summarize_method_decade_from_donor(method_donor_all)
    method_correlations_all = method_age_correlations(method_sample_all)

    method_donor_age_focused = summarize_call_methods_by_donor(obs_age_focused, call_method_cols)
    method_sample_age_focused = summarize_call_methods_by_sample(obs_age_focused, call_method_cols)
    method_decade_age_focused = summarize_method_decade_from_donor(method_donor_age_focused)
    method_correlations_age_focused = method_age_correlations(method_sample_age_focused)

    donor_summary = (
        obs.groupby(["donor_id_for_summary", "age_decade"], observed=True)
        .agg(
            donor_age=("donor_age_num", "first"),
            n_cells=("snc_bool", "size"),
            n_snc=("snc_bool", "sum"),
            snc_prop=("snc_bool", "mean"),
            mean_sasp_num_user=("sasp_num_user", "mean"),
            mean_sasp_avg_expr_user=("sasp_avg_expr_user", "mean"),
            snc_mean_sasp_num_user=("sasp_num_user", lambda x: x[obs.loc[x.index, "snc_bool"]].mean()),
            cdkn1a_pos_rate=("CDKN1A_pos", "mean"),
            cdkn2a_pos_rate=("CDKN2A_pos", "mean"),
            pcna_pos_rate=("PCNA_pos", "mean"),
            mki67_pos_rate=("MKI67_pos", "mean"),
            n_sites=("site_for_summary", "nunique"),
        )
        .reset_index()
    )
    donor_summary["snc_prop_percent"] = 100 * donor_summary["snc_prop"]

    decade_summary = summarize_by_decade(
        donor_summary,
        [
            "snc_prop_percent",
            "mean_sasp_num_user",
            "mean_sasp_avg_expr_user",
            "snc_mean_sasp_num_user",
            "cdkn1a_pos_rate",
            "cdkn2a_pos_rate",
            "pcna_pos_rate",
            "mki67_pos_rate",
        ],
    )

    donor_site_summary = (
        obs.groupby(["donor_id_for_summary", "age_decade", "site_for_summary"], observed=True)
        .agg(
            donor_age=("donor_age_num", "first"),
            n_cells=("snc_bool", "size"),
            n_snc=("snc_bool", "sum"),
            snc_prop=("snc_bool", "mean"),
            mean_sasp_num_user=("sasp_num_user", "mean"),
        )
        .reset_index()
    )
    donor_site_summary["snc_prop_percent"] = 100 * donor_site_summary["snc_prop"]

    sample_group_cols = ["sample_for_summary", "donor_id_for_summary", "age_decade", "site_for_summary"]
    sample_summary = (
        obs.groupby(sample_group_cols, observed=True)
        .agg(
            donor_age=("donor_age_num", "first"),
            n_cells=("snc_bool", "size"),
            n_snc=("snc_bool", "sum"),
            snc_prop=("snc_bool", "mean"),
            mean_sasp_num_user=("sasp_num_user", "mean"),
            mean_sasp_avg_expr_user=("sasp_avg_expr_user", "mean"),
            cdkn1a_pos_rate=("CDKN1A_pos", "mean"),
            cdkn2a_pos_rate=("CDKN2A_pos", "mean"),
            pcna_pos_rate=("PCNA_pos", "mean"),
            mki67_pos_rate=("MKI67_pos", "mean"),
            n_celltypes=("cell_type_for_summary", "nunique"),
        )
        .reset_index()
    )
    sample_summary["snc_prop_percent"] = 100 * sample_summary["snc_prop"]
    sample_summary["sample_unit_id"] = (
        sample_summary["sample_for_summary"].astype(str)
        + " | "
        + sample_summary["donor_id_for_summary"].astype(str)
        + " | "
        + sample_summary["site_for_summary"].astype(str)
        + " | "
        + sample_summary["age_decade"].astype(str)
    )

    sample_celltype_composition = (
        obs.groupby(sample_group_cols + ["cell_type_for_summary"], observed=True)
        .size()
        .rename("n_cells_celltype")
        .reset_index()
    )
    sample_celltype_composition = sample_celltype_composition.merge(
        sample_summary[sample_group_cols + ["sample_unit_id", "n_cells"]],
        on=sample_group_cols,
        how="left",
    )
    sample_celltype_composition["cell_type_fraction"] = (
        sample_celltype_composition["n_cells_celltype"] / sample_celltype_composition["n_cells"]
    )

    sample_snc_celltype_composition = (
        obs[obs["snc_bool"]]
        .groupby(sample_group_cols + ["cell_type_for_summary"], observed=True)
        .size()
        .rename("n_snc_celltype")
        .reset_index()
    )
    sample_snc_celltype_composition = sample_snc_celltype_composition.merge(
        sample_summary[sample_group_cols + ["sample_unit_id", "n_snc"]],
        on=sample_group_cols,
        how="left",
    )
    sample_snc_celltype_composition = sample_snc_celltype_composition[sample_snc_celltype_composition["n_snc"] > 0].copy()
    sample_snc_celltype_composition["snc_celltype_fraction"] = (
        sample_snc_celltype_composition["n_snc_celltype"] / sample_snc_celltype_composition["n_snc"]
    )

    celltype_donor_summary = (
        obs.groupby(["donor_id_for_summary", "age_decade", "cell_type_for_summary"], observed=True)
        .agg(
            donor_age=("donor_age_num", "first"),
            n_cells=("snc_bool", "size"),
            n_snc=("snc_bool", "sum"),
            snc_prop=("snc_bool", "mean"),
            mean_sasp_num_user=("sasp_num_user", "mean"),
        )
        .reset_index()
    )
    celltype_donor_summary = celltype_donor_summary[celltype_donor_summary["n_cells"] >= MIN_CELLS_PER_DONOR_GROUP]
    celltype_donor_summary["snc_prop_percent"] = 100 * celltype_donor_summary["snc_prop"]

    top_celltypes = (
        celltype_donor_summary.groupby("cell_type_for_summary", observed=True)["n_cells"]
        .sum()
        .sort_values(ascending=False)
        .head(TOP_N_CELL_TYPES)
        .index
    )
    celltype_decade_summary = (
        celltype_donor_summary[celltype_donor_summary["cell_type_for_summary"].isin(top_celltypes)]
        .groupby(["age_decade", "cell_type_for_summary"], observed=True)
        .agg(
            mean_snc_prop_percent=("snc_prop_percent", "mean"),
            sem_snc_prop_percent=("snc_prop_percent", sem),
            mean_sasp_num_user=("mean_sasp_num_user", "mean"),
            n_donors=("donor_id_for_summary", "nunique"),
        )
        .reset_index()
    )

    snc_obs = obs[obs["snc_bool"]].copy()
    snc_composition_donor = (
        snc_obs.groupby(["donor_id_for_summary", "age_decade", "cell_type_for_summary"], observed=True)
        .size()
        .rename("n_snc_celltype")
        .reset_index()
    )
    totals = (
        snc_obs.groupby(["donor_id_for_summary", "age_decade"], observed=True)
        .size()
        .rename("n_snc_total")
        .reset_index()
    )
    snc_composition_donor = snc_composition_donor.merge(totals, on=["donor_id_for_summary", "age_decade"], how="left")
    snc_composition_donor["snc_composition_fraction"] = (
        snc_composition_donor["n_snc_celltype"] / snc_composition_donor["n_snc_total"]
    )
    top_snc_celltypes = (
        snc_composition_donor.groupby("cell_type_for_summary", observed=True)["n_snc_celltype"]
        .sum()
        .sort_values(ascending=False)
        .head(TOP_N_CELL_TYPES)
        .index
    )
    snc_composition_decade = (
        snc_composition_donor[snc_composition_donor["cell_type_for_summary"].isin(top_snc_celltypes)]
        .groupby(["age_decade", "cell_type_for_summary"], observed=True)
        .agg(
            mean_snc_composition_fraction=("snc_composition_fraction", "mean"),
            sem_snc_composition_fraction=("snc_composition_fraction", sem),
            n_donors=("donor_id_for_summary", "nunique"),
        )
        .reset_index()
    )

    site_fraction_donor = (
        obs.groupby(["donor_id_for_summary", "age_decade", "site_for_summary"], observed=True)
        .size()
        .rename("n_cells_site")
        .reset_index()
    )
    donor_totals = (
        obs.groupby(["donor_id_for_summary", "age_decade"], observed=True)
        .size()
        .rename("n_cells_total")
        .reset_index()
    )
    site_fraction_donor = site_fraction_donor.merge(donor_totals, on=["donor_id_for_summary", "age_decade"], how="left")
    site_fraction_donor["site_cell_fraction"] = site_fraction_donor["n_cells_site"] / site_fraction_donor["n_cells_total"]
    site_fraction_decade = (
        site_fraction_donor.groupby(["age_decade", "site_for_summary"], observed=True)
        .agg(mean_site_cell_fraction=("site_cell_fraction", "mean"), n_donors=("donor_id_for_summary", "nunique"))
        .reset_index()
    )

    obs.to_csv(outdir / "sinai_all_sites_snc_decade_obs.csv.gz", index=True, compression="gzip")
    donor_summary.to_csv(outdir / "snc_summary_by_donor_decade.csv", index=False)
    decade_summary.to_csv(outdir / "snc_trends_by_age_decade_donor_level.csv", index=False)
    donor_site_summary.to_csv(outdir / "snc_summary_by_donor_site_decade.csv", index=False)
    sample_summary.to_csv(outdir / "snc_summary_by_sample_site_decade.csv", index=False)
    sample_celltype_composition.to_csv(outdir / "sample_celltype_composition_long.csv", index=False)
    sample_snc_celltype_composition.to_csv(outdir / "snc_celltype_composition_by_sample_long.csv", index=False)
    celltype_donor_summary.to_csv(outdir / "snc_summary_by_donor_celltype_decade.csv", index=False)
    celltype_decade_summary.to_csv(outdir / "snc_celltype_trends_by_age_decade.csv", index=False)
    snc_composition_donor.to_csv(outdir / "snc_celltype_composition_by_donor_decade.csv", index=False)
    snc_composition_decade.to_csv(outdir / "snc_celltype_composition_by_age_decade.csv", index=False)
    site_fraction_decade.to_csv(outdir / "site_composition_by_age_decade.csv", index=False)
    method_donor_all.to_csv(outdir / "snc_method_comparison_by_donor_decade_all_sites.csv", index=False)
    method_sample_all.to_csv(outdir / "snc_method_comparison_by_sample_all_sites.csv", index=False)
    method_decade_all.to_csv(outdir / "snc_method_comparison_by_age_decade_all_sites.csv", index=False)
    method_correlations_all.to_csv(outdir / "snc_method_age_correlations_all_sites.csv", index=False)
    method_donor_age_focused.to_csv(outdir / "snc_method_comparison_by_donor_decade_age_focused_filtered_sites.csv", index=False)
    method_sample_age_focused.to_csv(outdir / "snc_method_comparison_by_sample_age_focused_filtered_sites.csv", index=False)
    method_decade_age_focused.to_csv(outdir / "snc_method_comparison_by_age_decade_age_focused_filtered_sites.csv", index=False)
    method_correlations_age_focused.to_csv(outdir / "snc_method_age_correlations_age_focused_filtered_sites.csv", index=False)

    params = {
        "input_h5ad": str(args.input_h5ad),
        "sasp_threshold": args.sasp_threshold,
        "snc_rule": "(CDKN1A+ OR CDKN2A+) AND NOT (PCNA+ OR MKI67+) AND sasp_num_user >= threshold",
        "age_binning": "10-year bins using donor_age, right-open intervals [20,30), [30,40), etc.",
        "site_filter": "none; all anatomical sites retained",
        "donor_col": donor_col,
        "sample_col": sample_col,
        "site_col": site_col,
        "cell_type_col": cell_type_col,
        "broad_cell_type_col": broad_cell_type_col,
        "count_layer_for_positivity": count_layer_name,
        "expr_layer_for_sasp_average": expr_layer_name,
        "sasp_original_requested": SASP_ORIGINAL,
        "sasp_original_present": sasp_original_present,
        "sasp_user_requested": SASP_USER_EXTENDED,
        "sasp_user_present": sasp_user_present,
        "component_quartile_rule": "within cell_type_for_summary lineage: growth Q4 AND SASP Q4 AND cell-cycle Q1",
        "senescence_index_formula": "SI = Z_growth + Z_SASP - Z_cycle; component z-scores calculated within lineage; SI is not used as the quartile gate",
        "quartile_growth_genes_present": growth_present,
        "quartile_s_phase_genes_present": s_phase_present,
        "quartile_g2m_genes_present": g2m_present,
        "age_focused_excluded_sites": sorted(SITE_FILTER_EXCLUDE_FOR_AGE_FOCUSED),
    }
    (outdir / "sinai_age_decade_snc_parameters.json").write_text(json.dumps(params, indent=2))

    save_line_with_points(
        donor_summary,
        decade_summary,
        "snc_prop_percent",
        "SnC+ cells (%)",
        outdir / "snc_prop_by_age_decade_donor_level.png",
    )
    save_line_with_points(
        donor_summary,
        decade_summary,
        "mean_sasp_num_user",
        "Mean detected SASP genes per cell",
        outdir / "mean_sasp_num_by_age_decade_donor_level.png",
    )
    save_line_with_points(
        donor_summary,
        decade_summary,
        "mean_sasp_avg_expr_user",
        "Mean SASP expression per cell",
        outdir / "mean_sasp_expr_by_age_decade_donor_level.png",
    )
    save_heatmap(
        celltype_decade_summary,
        "mean_snc_prop_percent",
        "cell_type_for_summary",
        "age_decade",
        outdir / "celltype_snc_prop_heatmap_by_age_decade.png",
        "SnC+ percentage by cell type and age decade",
    )
    save_heatmap(
        celltype_decade_summary,
        "mean_sasp_num_user",
        "cell_type_for_summary",
        "age_decade",
        outdir / "celltype_sasp_num_heatmap_by_age_decade.png",
        "Mean SASP gene count by cell type and age decade",
    )
    save_bar_plot(
        snc_composition_decade,
        "age_decade",
        "mean_snc_composition_fraction",
        "cell_type_for_summary",
        "Mean fraction of SnC+ cells",
        outdir / "snc_celltype_composition_by_age_decade.png",
    )
    save_bar_plot(
        site_fraction_decade,
        "age_decade",
        "mean_site_cell_fraction",
        "site_for_summary",
        "Mean fraction of cells",
        outdir / "site_composition_by_age_decade_all_sites.png",
    )
    save_method_comparison_plot(
        method_donor_all,
        method_decade_all,
        outdir / "snc_method_comparison_by_age_decade_all_sites.png",
        "Original marker-positive SnC vs component Q4/Q4/Q1 gate, all sites",
    )
    save_method_comparison_plot(
        method_donor_age_focused,
        method_decade_age_focused,
        outdir / "snc_method_comparison_by_age_decade_age_focused_filtered_sites.png",
        "Original marker-positive SnC vs component Q4/Q4/Q1 gate, site-filtered age-focused set",
    )

    sample_features, sample_metadata = prepare_sample_feature_matrix(sample_summary, sample_celltype_composition)
    sample_features.to_csv(outdir / "sample_similarity_features.csv")
    sample_similarity = compute_sample_similarity_matrix(sample_features)
    if not sample_similarity.empty:
        sample_similarity.to_csv(outdir / "sample_similarity_cosine_matrix.csv")
        nearest = []
        for sample_id in sample_similarity.index:
            ranked = sample_similarity.loc[sample_id].drop(index=sample_id).sort_values(ascending=False)
            source = sample_metadata.loc[sample_id]
            for rank, (neighbor_id, similarity) in enumerate(ranked.head(8).items(), start=1):
                neighbor = sample_metadata.loc[neighbor_id]
                nearest.append(
                    {
                        "sample_unit_id": sample_id,
                        "sample": source["sample_for_summary"],
                        "site": source["site_for_summary"],
                        "age_decade": source["age_decade"],
                        "snc_prop_percent": source["snc_prop_percent"],
                        "mean_sasp_num_user": source["mean_sasp_num_user"],
                        "rank": rank,
                        "neighbor_sample_unit_id": neighbor_id,
                        "neighbor_sample": neighbor["sample_for_summary"],
                        "neighbor_site": neighbor["site_for_summary"],
                        "neighbor_age_decade": neighbor["age_decade"],
                        "neighbor_snc_prop_percent": neighbor["snc_prop_percent"],
                        "neighbor_mean_sasp_num_user": neighbor["mean_sasp_num_user"],
                        "cosine_similarity": similarity,
                        "same_site": source["site_for_summary"] == neighbor["site_for_summary"],
                        "same_age_decade": source["age_decade"] == neighbor["age_decade"],
                    }
                )
        pd.DataFrame(nearest).to_csv(outdir / "sample_nearest_neighbors_by_similarity.csv", index=False)

    pca_site = save_sample_pca_plot(
        sample_features,
        sample_metadata,
        "site_for_summary",
        outdir / "sample_feature_pca_by_anatomical_site.png",
        "Sample-level similarity features colored by anatomical site",
    )
    if not pca_site.empty:
        pca_site.to_csv(outdir / "sample_feature_pca_coordinates.csv", index=False)
    save_sample_pca_plot(
        sample_features,
        sample_metadata,
        "age_decade",
        outdir / "sample_feature_pca_by_age_decade.png",
        "Sample-level similarity features colored by age decade",
    )
    save_sample_composition_stacked(
        sample_celltype_composition,
        sample_summary,
        outdir / "sample_celltype_composition_stacked_by_site_age.png",
    )
    save_snc_sample_composition_stacked(
        sample_snc_celltype_composition,
        sample_summary,
        outdir / "snc_celltype_composition_by_sample_stacked_all_sites.png",
    )
    age_focused_sample_ids = sample_summary.loc[
        ~sample_summary["site_for_summary"].astype("string").str.lower().isin(SITE_FILTER_EXCLUDE_FOR_AGE_FOCUSED),
        "sample_unit_id",
    ]
    save_sample_snc_proportion_bar(
        sample_summary,
        outdir / "snc_proportion_by_sample_all_sites.png",
        "Original marker-positive SnC proportion by sample, all sites",
    )
    save_sample_snc_proportion_bar(
        sample_summary[sample_summary["sample_unit_id"].isin(age_focused_sample_ids)],
        outdir / "snc_proportion_by_sample_age_focused_filtered_sites.png",
        "Original marker-positive SnC proportion by sample, site-filtered age-focused set",
    )
    save_snc_sample_composition_stacked(
        sample_snc_celltype_composition[
            sample_snc_celltype_composition["sample_unit_id"].isin(age_focused_sample_ids)
        ],
        sample_summary[sample_summary["sample_unit_id"].isin(age_focused_sample_ids)],
        outdir / "snc_celltype_composition_by_sample_stacked_age_focused_filtered_sites.png",
    )
    save_umap_plot(
        obs,
        "site_for_summary",
        "Sinai UMAP colored by anatomical site",
        outdir / "umap_by_anatomical_site.png",
    )
    save_umap_plot(
        obs,
        "age_decade_for_plot",
        "Sinai UMAP colored by donor age decade",
        outdir / "umap_by_age_decade.png",
    )

    print(f"Saved outputs to: {outdir}")


if __name__ == "__main__":
    main()
