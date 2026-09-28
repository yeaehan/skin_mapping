from __future__ import annotations

import json
import logging
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scanpy as sc
import muspan as ms
from scipy import stats
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("28_sssp_muspan")

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parents[3]          # skin_mapping_2/
H5AD = ROOT / "data" / "Xenium" / "xenium_full.h5ad"
SSSP_JSON = ROOT / "code" / "shared" / "sssp.json"
OUT_DIR = ROOT / "results" / "28_sssp_muspan"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
GETIS_ALPHA = 0.05        # MuSpAn significance level
MAX_EDGE_DISTANCE = np.inf # Delaunay network maximum edge distance
SNC_STATUS_COLUMN = None    # set explicitly if known; otherwise auto-detect from obs
SNC_STATUS_CANDIDATES = [
    "is_snc", "snc", "SnC", "snc_status", "SnC_status",
    "senescent", "senescence_status", "cell_senescence_status",
]
N_BOOTSTRAP = 2000        # permutations for empirical p-value


def load_data() -> sc.AnnData:
    log.info("Loading h5ad: %s", H5AD)
    adata = sc.read_h5ad(H5AD)
    log.info("  shape = %s", adata.shape)
    log.info("  obs columns: %s", list(adata.obs.columns)[:20])
    return adata


def get_sssp_genes() -> list[str]:
    """Read the Xenium-specific SSSP gene list from shared/sssp.json."""
    log.info("Reading SSSP gene panel: %s", SSSP_JSON)
    with open(SSSP_JSON) as f:
        panel = json.load(f)

    if not isinstance(panel, dict) or "sssp_xenium" not in panel:
        raise KeyError(
            "Expected shared/sssp.json to contain a top-level 'sssp_xenium' key."
        )

    entries = panel["sssp_xenium"]
    genes = []
    for entry in entries:
        if isinstance(entry, str):
            genes.append(entry)
        elif isinstance(entry, dict):
            gene = entry.get("gene") or entry.get("symbol") or entry.get("Gene")
            if gene:
                genes.append(gene)

    genes = [g for g in genes if g]
    log.info("  SSSP Xenium panel: %d genes", len(genes))
    return genes


def score_cells(adata: sc.AnnData, genes: list[str]) -> np.ndarray:
    """
    Rank-based z-score per cell across the SSSP gene set.
    Returns a 1-D array (n_obs,) of SSSP scores.
    """
    present = [g for g in genes if g in adata.var.index]
    missing = [g for g in genes if g not in adata.var.index]
    if missing:
        log.warning("  %d SSSP genes not in var: %s", len(missing), missing[:10])
    log.info("  Scoring with %d / %d genes", len(present), len(genes))

    # Work in a copy to avoid mutating adata
    adata_sscp = adata[:, present].copy()
    sc.pp.normalize_total(adata_sscp, target_sum=1e4)
    sc.pp.log1p(adata_sscp)

    # Rank each gene across cells, average ranks → score
    X = adata_sscp.X.toarray() if hasattr(adata_sscp.X, "toarray") else np.asarray(adata_sscp.X)
    ranks = np.zeros_like(X, dtype=float)
    for j in range(X.shape[1]):
        ranks[:, j] = stats.rankdata(X[:, j])
    scores = ranks.mean(axis=1)

    # z-score
    scores = (scores - scores.mean()) / (scores.std(ddof=1) + 1e-12)
    adata.obsm["X_sscp_score"] = scores.reshape(-1, 1)
    log.info("  SSSP score range: [%.3f, %.3f]", scores.min(), scores.max())
    return scores


def run_muspan_getis_ord(coords: np.ndarray, scores: np.ndarray) -> pd.DataFrame:
    """
    Build a MuSpAn domain from cell centroids and run MuSpAn's native
    local Getis-Ord G* implementation on the continuous SSSP score.

    MuSpAn returns local z-scores, Bonferroni-adjusted p-values, and the
    object indices used in the analysis.
    """
    log.info("Building MuSpAn domain on %d cells", len(coords))

    domain = ms.domain("Xenium SSSP")
    domain.add_points(coords, collection_name="cells")
    domain.add_labels(
        label_name="SSSP score",
        labels=np.asarray(scores, dtype=float),
        add_labels_to=("collection", "cells"),
        label_type="continuous",
    )

    log.info("Running MuSpAn Getis-Ord G* on SSSP score")
    z_G, p_adj, object_indices = ms.spatial_statistics.getis_ord(
        domain,
        population=("collection", "cells"),
        label_name="SSSP score",
        alpha=GETIS_ALPHA,
        network_kwargs={
            "network_type": "Delaunay",
            "min_edge_distance": 0,
            "max_edge_distance": MAX_EDGE_DISTANCE,
        },
        add_local_value_as_label=True,
        local_getis_label_name="Local Getis Ord G*: SSSP",
    )

    # Preserve cell order in the AnnData object.
    result = pd.DataFrame(
        {
            "z": np.nan,
            "p_adj": np.nan,
            "significance": "ns",
        },
        index=np.arange(len(coords)),
    )

    object_indices = np.asarray(object_indices, dtype=int)
    z_G = np.asarray(z_G, dtype=float)
    p_adj = np.asarray(p_adj, dtype=float)

    result.loc[object_indices, "z"] = z_G
    result.loc[object_indices, "p_adj"] = p_adj

    significant = p_adj <= GETIS_ALPHA
    result.loc[object_indices[significant & (z_G > 0)], "significance"] = "hot"
    result.loc[object_indices[significant & (z_G < 0)], "significance"] = "cold"

    result.index.name = "cell"
    return result


def identify_snc_cells(adata: sc.AnnData) -> np.ndarray:
    """
    Read SnC status directly from adata.obs.

    If SNC_STATUS_COLUMN is set, that exact column is used. Otherwise the
    function searches common SnC-status column names. Boolean, numeric, and
    categorical/string encodings are supported.
    """
    if SNC_STATUS_COLUMN is not None:
        if SNC_STATUS_COLUMN not in adata.obs.columns:
            raise KeyError(
                f"SNC_STATUS_COLUMN={SNC_STATUS_COLUMN!r} not found in adata.obs. "
                f"Available columns include: {list(adata.obs.columns)[:50]}"
            )
        col = SNC_STATUS_COLUMN
    else:
        col = next((c for c in SNC_STATUS_CANDIDATES if c in adata.obs.columns), None)
        if col is None:
            raise KeyError(
                "Could not auto-detect the SnC status column in adata.obs. "
                "Set SNC_STATUS_COLUMN near the top of the script to the correct obs column. "
                f"Available columns include: {list(adata.obs.columns)[:50]}"
            )

    values = adata.obs[col]
    log.info("Using existing SnC status from adata.obs[%r]", col)

    if pd.api.types.is_bool_dtype(values):
        positive = values.fillna(False).to_numpy(dtype=bool)

    elif pd.api.types.is_numeric_dtype(values):
        positive = values.fillna(0).to_numpy() > 0

    else:
        s = values.astype(str).str.strip().str.lower()
        positive_labels = {
            "1", "true", "yes", "y", "snc", "snc+", "positive", "pos",
            "senescent", "senescence", "senescence+", "senescent cell",
        }
        negative_labels = {
            "0", "false", "no", "n", "snc-", "negative", "neg",
            "non-snc", "non_snc", "non-senescent", "nonsenescent",
        }

        observed = set(s.dropna().unique())
        unknown = observed - positive_labels - negative_labels - {"nan", "none", ""}
        if unknown:
            log.warning(
                "Unrecognized values in SnC status column %r: %s. "
                "They will be treated as SnC-negative unless they contain 'snc+' or 'senesc'.",
                col, sorted(unknown)[:20]
            )

        positive = (
            s.isin(positive_labels)
            | s.str.contains(r"snc\+", regex=True, na=False)
            | s.str.contains("senesc", regex=False, na=False)
        ).to_numpy()

    log.info(
        "SnC+ cells from obs metadata: %d / %d (%.1f%%)",
        int(positive.sum()), len(positive), 100 * positive.mean()
    )
    return positive


def bootstrap_enrichment(snc_mask: np.ndarray, hotspot_mask: np.ndarray, n_boot: int) -> float:
    """
    Permutation test: randomise SnC label, compute hotspot fraction each time.
    Returns empirical p-value.
    """
    n_total = len(snc_mask)
    n_hot = hotspot_mask.sum()
    n_snc = snc_mask.sum()
    if n_hot == 0 or n_snc == 0:
        return 1.0

    obs_enrich = snc_mask[hotspot_mask].sum() / n_hot

    rng = np.random.default_rng(42)
    count_ge = 0
    for _ in range(n_boot):
        perm = rng.permutation(snc_mask)
        boot_enrich = perm[hotspot_mask].sum() / n_hot
        if boot_enrich >= obs_enrich:
            count_ge += 1
    p = (count_ge + 1) / (n_boot + 1)
    return p


def distance_to_hotspots(coords: np.ndarray, snc_mask: np.ndarray, hotspot_mask: np.ndarray) -> np.ndarray:
    """
    For each SnC+ cell, compute Euclidean distance to the nearest hotspot cell.
    """
    hot_coords = coords[hotspot_mask]
    snc_coords = coords[snc_mask]
    if len(hot_coords) == 0:
        return np.full(snc_mask.sum(), np.nan)
    tree = cKDTree(hot_coords)
    dists, _ = tree.query(snc_coords, k=1)
    return dists


def make_summary(
    snc_mask: np.ndarray,
    hotspot_mask: np.ndarray,
    coldspot_mask: np.ndarray,
    ns_mask: np.ndarray,
    p_boot: float,
) -> pd.DataFrame:
    total = len(snc_mask)
    n_snc = snc_mask.sum()

    def frac(region):
        if region.sum() == 0:
            return np.nan
        return snc_mask[region].sum() / region.sum()

    rows = {
        "region": ["hotspot", "coldspot", "non-significant", "all"],
        "n_cells": [hotspot_mask.sum(), coldspot_mask.sum(), ns_mask.sum(), total],
        "n_snc+": [
            snc_mask[hotspot_mask].sum(),
            snc_mask[coldspot_mask].sum(),
            snc_mask[ns_mask].sum(),
            n_snc,
        ],
        "snc_fraction": [frac(hotspot_mask), frac(coldspot_mask), frac(ns_mask), n_snc / total],
    }
    df = pd.DataFrame(rows)
    df.loc[df.index[0], "bootstrap_p"] = p_boot

    # Enrichment ratio (hotspot vs background)
    bg_frac = n_snc / total
    hot_frac = rows["snc_fraction"][0]
    if hot_frac is not np.nan and bg_frac > 0:
        df.loc[0, "enrichment_ratio"] = hot_frac / bg_frac
    else:
        df.loc[0, "enrichment_ratio"] = np.nan

    return df


# ---------------------------------------------------------------------------
# Age-stratified analysis
# ---------------------------------------------------------------------------
def compute_age_stats(
    adata: sc.AnnData,
    hotspot_mask: np.ndarray,
    snc_mask: np.ndarray,
    coords: np.ndarray,
    n_boot: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Stratify hotspot stats by `age` (Old / Young).
    Returns:
      age_summary       – per-cohort per-region SnC fraction, enrichment, bootstrap p, MW-U
      age_distance_mwu  – MW-U comparing SnC+ vs SnC- distances to hotspots, per cohort
      age_distance_long – long-form table of per-SnC+ cell distance-to-nearest-hotspot
                         (for plotting and downstream stats).
    """
    ages = adata.obs["age"].astype(str).values
    samples = adata.obs["sample"].astype(str).values if "sample" in adata.obs.columns else np.array(["all"] * len(ages))
    n = len(ages)

    # ---- Per-cohort SnC+ fractions and enrichment ----
    cohorts = sorted(set(ages))
    rows = []
    distance_rows = []
    distance_mwu_rows = []

    hot_idx_global = np.where(hotspot_mask)[0]

    for cohort in cohorts:
        cohort_mask = (ages == cohort)
        n_c = int(cohort_mask.sum())
        n_snc_c = int(snc_mask[cohort_mask].sum())
        bg_frac = n_snc_c / n_c if n_c else np.nan

        for region_name, region_mask in [
            ("hotspot",  hotspot_mask & cohort_mask),
            ("coldspot", ~hotspot_mask & cohort_mask),    # we don't recompute coldspot per cohort; keep simple
            ("non_sig",  ~hotspot_mask & cohort_mask),     # alias for "all non-hotspot" cells in this cohort
            ("all",      cohort_mask),
        ]:
            n_r = int(region_mask.sum())
            n_snc_r = int(snc_mask[region_mask].sum())
            snc_frac = n_snc_r / n_r if n_r else np.nan

            if region_name == "hotspot":
                enrich = (snc_frac / bg_frac) if (snc_frac and bg_frac) else np.nan
                p_boot = bootstrap_enrichment(
                    snc_mask[cohort_mask],          # restrict permutation to the cohort
                    hotspot_mask[cohort_mask],
                    n_boot,
                )
            else:
                enrich = np.nan
                p_boot = np.nan

            rows.append({
                "age": cohort,
                "region": region_name,
                "n_cells": n_r,
                "n_snc+": n_snc_r,
                "snc_fraction": snc_frac,
                "background_snc_fraction": bg_frac,
                "enrichment_ratio_vs_cohort_bg": enrich,
                "bootstrap_p": p_boot,
            })

        # ---- Distance-to-hotspot (per cohort, SnC+ vs SnC-) ----
        if hot_idx_global.size == 0:
            d_snc = np.array([])
            d_non = np.array([])
        else:
            tree = cKDTree(coords[hot_idx_global])
            d_all, _ = tree.query(coords, k=1)
        d_snc = d_all[snc_mask & cohort_mask]
        d_non = d_all[(~snc_mask) & cohort_mask]

        for label, d in [("SnC+", d_snc), ("SnC-", d_non)]:
            for v in d:
                distance_rows.append({"age": cohort, "snc_status": label, "distance_px": float(v)})

        if len(d_snc) > 0 and len(d_non) > 0:
            U, p = stats.mannwhitneyu(d_snc, d_non, alternative="two-sided")
        else:
            U, p = np.nan, np.nan
        distance_mwu_rows.append({
            "age": cohort,
            "n_snc+_with_distance": int(len(d_snc)),
            "n_snc-_with_distance": int(len(d_non)),
            "median_dist_snc+_px": float(np.median(d_snc)) if len(d_snc) else np.nan,
            "median_dist_snc-_px": float(np.median(d_non)) if len(d_non) else np.nan,
            "mannwhitney_U": float(U) if not np.isnan(U) else np.nan,
            "mannwhitney_p": float(p) if not np.isnan(p) else np.nan,
        })

    age_summary = pd.DataFrame(rows)
    age_distance_mwu = pd.DataFrame(distance_mwu_rows)
    age_distance_long = pd.DataFrame(distance_rows)
    return age_summary, age_distance_mwu, age_distance_long


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_spatial(
    coords: np.ndarray,
    scores: np.ndarray,
    sig: pd.Series,
    snc_mask: np.ndarray,
    out_path: Path,
    samples: np.ndarray | None = None,
):
    """
    4-panel spatial plot for the whole-cohort view:
      (a) SSSP score, (b) Getis-Ord significance, (c) SnC+ overlay,
      (d) hotspot-density-per-sample bar chart (only if `samples` provided).
    """
    fig, axes = plt.subplots(1, 4, figsize=(30, 6.5))

    # (a) SSSP score
    im0 = axes[0].scatter(coords[:, 0], coords[:, 1], c=scores, s=0.6, cmap="viridis", vmin=np.percentile(scores, 1), vmax=np.percentile(scores, 99))
    axes[0].set_title("SSSP gene score (24 genes)", fontsize=12)
    axes[0].set_xticks([]); axes[0].set_yticks([])
    plt.colorbar(im0, ax=axes[0], fraction=0.04)

    # (b) Getis-Ord significance
    cmap = mcolors.ListedColormap(["#87CEEB", "#FFFFFF", "#FF6347"])
    vals = np.array([2 if s == "hot" else (0 if s == "cold" else 1) for s in sig])
    axes[1].scatter(coords[:, 0], coords[:, 1], c=vals, s=0.6, cmap=cmap, vmin=0, vmax=2)
    axes[1].set_title("Getis-Ord G* (SSSP) hot/cold spots", fontsize=12)
    axes[1].set_xticks([]); axes[1].set_yticks([])
    import matplotlib.patches as mpatches
    handles = [mpatches.Patch(color="#FF6347", label="Hot"),
               mpatches.Patch(color="#FFFFFF", label="ns"),
               mpatches.Patch(color="#87CEEB", label="Cold")]
    axes[1].legend(handles=handles, loc="upper right", fontsize=9)

    # (c) SnC+ overlay
    axes[2].scatter(coords[~snc_mask, 0], coords[~snc_mask, 1], c="lightgrey", s=0.3)
    axes[2].scatter(coords[snc_mask, 0], coords[snc_mask, 1], c="crimson", s=1.2, alpha=0.7)
    hot_idx = np.where((sig == "hot") & snc_mask)[0]
    axes[2].scatter(coords[hot_idx, 0], coords[hot_idx, 1], c="gold", s=4, edgecolors="k", linewidths=0.3, zorder=5)
    axes[2].set_title("SnC+ cells (red) in SSSP hotspots (gold)", fontsize=12)
    axes[2].set_xticks([]); axes[2].set_yticks([])

    # (d) per-sample hotspot density
    if samples is not None:
        per_sample = {}
        for s in np.unique(samples):
            m = (samples == s)
            n_m = int(m.sum())
            n_hot = int(((sig == "hot") & m).sum())
            n_snc_p = int((snc_mask & m).sum())
            per_sample[s] = (n_m, n_hot, n_hot / max(n_m, 1), n_snc_p)
        labels = list(per_sample.keys())
        hot_fracs = [per_sample[s][2] * 100 for s in labels]
        axes[3].bar(range(len(labels)), hot_fracs, color="#FF6347", alpha=0.85)
        axes[3].set_xticks(range(len(labels)))
        axes[3].set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
        axes[3].set_ylabel("% hotspot cells")
        axes[3].set_title("Hotspot density per sample", fontsize=12)
        for i, s in enumerate(labels):
            axes[3].text(i, hot_fracs[i] + max(hot_fracs) * 0.02, f"{per_sample[s][1]}/{per_sample[s][0]}", ha="center", fontsize=8)

    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved spatial plot → %s", out_path)


def plot_spatial_per_sample(
    coords: np.ndarray,
    sig: pd.Series,
    snc_mask: np.ndarray,
    samples: np.ndarray,
    out_path: Path,
):
    """
    One panel per sample (2 cols × N rows).
      grey: non-hotspot cells
      gold edge: SSSP hotspot
      red: SnC+
      gold edge + red: SnC+ in hotspot
    """
    import matplotlib.patches as mpatches
    unique_samples = sorted(np.unique(samples).tolist())
    n = len(unique_samples)
    ncols = 2
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 5.5 * nrows))
    axes = np.atleast_2d(axes).reshape(nrows, ncols)

    sig_arr = np.asarray(sig)
    is_hot = (sig_arr == "hot")
    is_snc = snc_mask.astype(bool)
    is_hot_snc = is_hot & is_snc

    for i, s in enumerate(unique_samples):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        m = (samples == s)
        if m.sum() == 0:
            ax.set_title(f"{s} (empty)", fontsize=11); ax.set_xticks([]); ax.set_yticks([])
            continue
        xy = coords[m]
        snc_in = is_snc[m]
        hot_in = is_hot[m]
        hot_snc_in = is_hot_snc[m]

        ax.scatter(xy[~snc_in & ~hot_in, 0], xy[~snc_in & ~hot_in, 1], c="lightgrey", s=2, alpha=0.6)
        ax.scatter(xy[hot_in & ~snc_in, 0], xy[hot_in & ~snc_in, 1], c="gold", s=6, edgecolors="k", linewidths=0.3, label="Hotspot cell")
        ax.scatter(xy[snc_in & ~hot_in, 0], xy[snc_in & ~hot_in, 1], c="crimson", s=8, alpha=0.85, label="SnC+ cell")
        ax.scatter(xy[hot_snc_in, 0], xy[hot_snc_in, 1], c="gold", s=18, edgecolors="k", linewidths=0.6, zorder=5, label="SnC+ in hotspot")

        n_total = int(m.sum())
        n_snc_s = int(snc_in.sum())
        n_hot_s = int(hot_in.sum())
        n_hot_snc_s = int(hot_snc_in.sum())
        title = f"{s}  (n={n_total:,}, SnC+={n_snc_s}, hot={n_hot_s}, SnC+\u2229hot={n_hot_snc_s})"
        ax.set_title(title, fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
        if i == 0:
            ax.legend(fontsize=8, loc="upper right")

    for j in range(n, nrows * ncols):
        r, c = divmod(j, ncols)
        axes[r, c].axis("off")

    fig.suptitle("SnC+ \u00d7 SSSP hotspots per sample", fontsize=14, y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved per-sample map → %s", out_path)


def plot_age_comparison(
    age_summary: pd.DataFrame,
    age_distance_long: pd.DataFrame,
    out_path: Path,
):
    """Two-panel figure: (a) per-cohort SnC+ fraction by region, (b) distance density."""
    cohorts_sorted = sorted(age_summary["age"].unique().tolist())

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (a) SnC+ fraction by cohort x region
    region_order = ["hotspot", "coldspot", "non_sig", "all"]
    pivot = (age_summary
             .pivot(index="age", columns="region", values="snc_fraction")
             .reindex(cohorts_sorted)[region_order])
    pivot.plot(kind="bar", ax=axes[0],
               color=["#FF6347", "#87CEEB", "#bbbbbb", "#444444"],
               edgecolor="k", linewidth=0.4)
    axes[0].set_ylabel("SnC+ fraction")
    axes[0].set_title("SnC+ fraction by age x region")
    axes[0].set_xticklabels(cohorts_sorted, rotation=0)
    axes[0].legend(title="region", fontsize=8)

    hot_only = age_summary[age_summary["region"] == "hotspot"].set_index("age")
    ymax = axes[0].get_ylim()[1] if axes[0].get_ylim()[1] > 0 else 1
    for i, cohort in enumerate(cohorts_sorted):
        if cohort in hot_only.index:
            e = hot_only.loc[cohort, "enrichment_ratio_vs_cohort_bg"]
            p = hot_only.loc[cohort, "bootstrap_p"]
            label = f"hot x{e:.2f}"
            if not (p is None or np.isnan(p)):
                label += f"\np={p:.3g}"
            axes[0].text(i - 0.27, ymax * 0.95, label, fontsize=8, ha="left", va="top", color="#a00")

    # (b) Distance density per cohort x SnC status
    for cohort in cohorts_sorted:
        for label, color, ls in [("SnC-", "steelblue", "-"), ("SnC+", "crimson", "--")]:
            d = age_distance_long.loc[(age_distance_long["age"] == cohort) & (age_distance_long["snc_status"] == label), "distance_px"].to_numpy()
            if len(d) == 0:
                continue
            axes[1].hist(d, bins=60, alpha=0.35, density=True, linestyle=ls,
                         label=f"{cohort} {label}",
                         color=color)
    axes[1].set_xlabel("Distance to nearest SSSP-hotspot cell (px)")
    axes[1].set_ylabel("Density")
    axes[1].set_title("Distance to SSSP hotspots by age x SnC")
    axes[1].legend(fontsize=8)

    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved age-comparison plot → %s", out_path)


def plot_scatter(
    scores: np.ndarray,
    sig: pd.Series,
    snc_mask: np.ndarray,
    dist_snc: np.ndarray,
    dist_non: np.ndarray,
    out_path: Path,
):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # (a) SSSP score by SnC status & hotspot
    groups = [
        ("SnC+  hot",  snc_mask & (sig == "hot")),
        ("SnC+  ns",   snc_mask & (sig == "ns")),
        ("SnC+  cold", snc_mask & (sig == "cold")),
        ("−SnC  hot",  ~snc_mask & (sig == "hot")),
        ("−SnC  ns",   ~snc_mask & (sig == "ns")),
        ("−SnC  cold", ~snc_mask & (sig == "cold")),
    ]
    data_list, labels = [], []
    for label, mask in groups:
        if mask.sum() > 0:
            data_list.append(np.percentile(scores[mask], [5, 50, 95]))
            labels.append(label)
    bp = axes[0].boxplot(data_list, labels=labels, patch_artist=True,
                         boxprops=dict(facecolor="#ddd"),
                         medianprops=dict(color="red", lw=1.5))
    axes[0].set_ylabel("SSSP score (z)")
    axes[0].set_title("SSSP score by SnC status × hotspot")
    axes[0].tick_params(axis="x", rotation=30)

    # (b) Distance to nearest hotspot cell
    axes[1].hist(dist_snc[~np.isnan(dist_snc)], bins=60, alpha=0.6, label="SnC+", density=True, color="crimson")
    axes[1].hist(dist_non[~np.isnan(dist_non)], bins=60, alpha=0.6, label="−SnC", density=True, color="steelblue")
    axes[1].set_xlabel("Distance to nearest SSSP-hotspot cell (px)")
    axes[1].set_ylabel("Density")
    axes[1].set_title("Proximity to SSSP hotspots")
    axes[1].legend()

    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved scatter/summary plot → %s", out_path)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info("=" * 60)
    log.info("Script 28 – Xenium SSSP MuSpAn hotspot × SnC+ enrichment")
    log.info("=" * 60)

    # 1. Load
    adata = load_data()

    # 2. SSSP score
    sssp_genes = get_sssp_genes()
    scores = score_cells(adata, sssp_genes)

    # 3. Spatial coords
    if "spatial" in adata.obsm:
        coords = adata.obsm["spatial"]
    elif "X_embedding" in adata.obsm:
        coords = adata.obsm["X_embedding"]
    elif "X_pca" in adata.obsm:
        log.warning("No spatial coords; falling back to PCA embedding.")
        coords = adata.obsm["X_pca"]
    else:
        raise ValueError("Cannot find spatial coordinates in obsm.")
    log.info("Using %s for spatial coordinates", "obsm['spatial']" if "spatial" in adata.obsm else "fallback")

    # 4–5. MuSpAn domain + native Getis-Ord G*
    muspan = run_muspan_getis_ord(coords, scores)
    adata.obsm["X_muspan_sscp"] = muspan[["z", "p_adj"]].values
    adata.obs["muspan_sscp"] = muspan["significance"].values

    n_hot = (muspan["significance"] == "hot").sum()
    n_cold = (muspan["significance"] == "cold").sum()
    log.info("  Hotspots: %d (%.1f%%)   Coldspots: %d (%.1f%%)",
             n_hot, 100 * n_hot / len(muspan),
             n_cold, 100 * n_cold / len(muspan))

    # 6. SnC+ status from Xenium obs metadata
    snc_mask = identify_snc_cells(adata)
    adata.obs["is_snc"] = snc_mask

    # 7. Enrichment + bootstrap
    hotspot_mask = (muspan["significance"] == "hot").values
    coldspot_mask = (muspan["significance"] == "cold").values
    ns_mask = (muspan["significance"] == "ns").values

    log.info("Bootstrap permutation test (%d iters) …", N_BOOTSTRAP)
    p_boot = bootstrap_enrichment(snc_mask, hotspot_mask, N_BOOTSTRAP)
    log.info("  SnC+ in hotspots: %.4f   (bootstrap p = %.4f)",
             snc_mask[hotspot_mask].mean() if hotspot_mask.sum() > 0 else float("nan"),
             p_boot)

    # 8. Distance to hotspots
    dist_snc = distance_to_hotspots(coords, snc_mask, hotspot_mask)
    non_snc_mask = ~snc_mask
    dist_non = distance_to_hotspots(coords, non_snc_mask, hotspot_mask)

    # Mann-Whitney U test on distances
    d_snc_valid = dist_snc[~np.isnan(dist_snc)]
    d_non_valid = dist_non[~np.isnan(dist_non)]
    if len(d_snc_valid) > 0 and len(d_non_valid) > 0:
        U_stat, p_mwu = stats.mannwhitneyu(d_snc_valid, d_non_valid, alternative="two-sided")
        log.info("  MW-U distance test: U=%.1f  p=%.4g", U_stat, p_mwu)
    else:
        U_stat, p_mwu = np.nan, np.nan

    # 9. Summary table
    summary = make_summary(snc_mask, hotspot_mask, coldspot_mask, ns_mask, p_boot)
    summary.loc[0, "mannwhitney_U"] = U_stat
    summary.loc[0, "mannwhitney_p"] = p_mwu
    summary_csv = OUT_DIR / "snc_enrichment_summary.csv"
    summary.to_csv(summary_csv, index=False)
    log.info("Saved summary → %s", summary_csv)
    print("\n" + summary.to_string(index=False) + "\n")

    # 10. Save muspan results
    out_muspan = OUT_DIR / "muspan_getisord.parquet"
    muspan.to_parquet(out_muspan)
    log.info("Saved MuSpAn results → %s", out_muspan)

    out_scores = OUT_DIR / "sssp_score_matrix.parquet"
    pd.Series(scores, name="sssp_score").to_frame().to_parquet(out_scores)
    log.info("Saved SSSP scores → %s", out_scores)

    # ---- Per-sample + age stratified analysis ----
    if "sample" in adata.obs.columns:
        samples_arr = adata.obs["sample"].astype(str).to_numpy()
    else:
        # fall back to a single 'all' bucket
        samples_arr = np.array(["all"] * len(adata), dtype=object)

    if "age" in adata.obs.columns:
        log.info("Computing age-stratified stats (age column found)")
        age_summary, age_distance_mwu, age_distance_long = compute_age_stats(
            adata, hotspot_mask, snc_mask, coords, N_BOOTSTRAP,
        )
        age_csv = OUT_DIR / "age_stratified_summary.csv"
        mwu_csv = OUT_DIR / "age_distance_mwu.csv"
        dlong_csv = OUT_DIR / "age_distance_per_cell.csv.gz"
        age_summary.to_csv(age_csv, index=False)
        age_distance_mwu.to_csv(mwu_csv, index=False)
        age_distance_long.to_csv(dlong_csv, index=False, compression="gzip")
        log.info("Saved age_stratified_summary → %s", age_csv)
        log.info("Saved age_distance_mwu → %s", mwu_csv)
        log.info("Saved age_distance_per_cell → %s", dlong_csv)
        print("\n--- Age-stratified hotspot enrichment ---")
        print(age_summary.to_string(index=False))
        print("\n--- Age-stratified MW-U distance test ---")
        print(age_distance_mwu.to_string(index=False))
    else:
        log.warning("No 'age' column in adata.obs; skipping age-stratified analysis.")
        age_summary = age_distance_long = None

    # 11. Plots
    plot_spatial(coords, scores, muspan["significance"], snc_mask,
                 OUT_DIR / "sssp_hotspot_snc_map.png",
                 samples=samples_arr)
    plot_spatial_per_sample(
        coords, muspan["significance"], snc_mask, samples_arr,
        OUT_DIR / "sssp_hotspot_snc_map_per_sample.png",
    )
    if age_summary is not None:
        plot_age_comparison(
            age_summary, age_distance_long,
            OUT_DIR / "sssp_hotspot_snc_age_comparison.png",
        )
    plot_scatter(scores, muspan["significance"], snc_mask, dist_snc, dist_non,
                 OUT_DIR / "sssp_hotspot_snc_scatter.png")

    # 12. Save annotated h5ad (optional)
    adata_out = OUT_DIR / "xenium_sscp_muspan_annotated.h5ad"
    adata.write_h5ad(adata_out)
    log.info("Saved annotated h5ad → %s", adata_out)

    log.info("=" * 60)
    log.info("DONE. All outputs in: %s", OUT_DIR)
    log.info("=" * 60)


if __name__ == "__main__":
    main()