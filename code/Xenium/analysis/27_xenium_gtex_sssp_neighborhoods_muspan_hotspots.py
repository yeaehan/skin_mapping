"""
Neighborhood and MuSpAn hotspot analysis for Xenium-measurable SSSP.

Inputs:
  results/260727/xenium_gtex_sssp_scanpy_score_analysis/
    xenium_gtex_sssp_scanpy_scored.h5ad

Definitions:
  SSSP-high = global top quartile of scanpy.tl.score_genes SSSP score.
  Neighborhood = all cells within 35 um of SSSP-high anchor cells, computed
    within each sample/FOV.
  Hotspots = 35 um hex bins tested with MuSpAn Getis-Ord G* using either mean
    SSSP score or SSSP-high fraction per bin.
"""

from __future__ import annotations

from pathlib import Path
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import cKDTree


PROJECT_ROOT = Path(__file__).resolve().parents[3]
INPUT_H5AD = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_gtex_sssp_scanpy_score_analysis"
    / "xenium_gtex_sssp_scanpy_scored.h5ad"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_gtex_sssp_neighborhoods_muspan_hotspots"

SCORE_COL = "sssp_scanpy_score"
HIGH_COL = "sssp_high_q75_global"
CELL_TYPE_COL = "cell_type_minor"
RADIUS_UM = 35.0
N_PERMUTATIONS = 1000
RANDOM_SEED = 260727

HEX_SIZE_UM = 35.0
HEX_NEIGHBOR_ORDER = 1
MIN_CELLS_PER_HEX = 5
ALPHA = 0.05
HOTSPOT_METRICS = ("mean_sssp_score", "sssp_high_prop")
HEX_ADJACENCY_DISTANCE_UM = np.sqrt(3.0) * HEX_SIZE_UM * HEX_NEIGHBOR_ORDER + 1e-6
AGE_ORDER = ["Young", "Old"]


def require_muspan():
    try:
        import muspan as ms
    except Exception as exc:
        raise RuntimeError(
            "Real MuSpAn is required for Getis-Ord hotspot testing. "
            "Install the authors' package with the credentials supplied by MuSpAn, e.g. "
            "python -m pip install -U https://docs.muspan.co.uk/code/latest.zip"
        ) from exc

    missing = [name for name in ("domain", "spatial_statistics") if not hasattr(ms, name)]
    if missing:
        raise RuntimeError(f"Installed muspan is missing expected API: {missing}")
    if not hasattr(ms.spatial_statistics, "getis_ord"):
        raise RuntimeError("Installed muspan.spatial_statistics has no getis_ord function.")
    return ms


def sample_sort_key(sample: str) -> tuple[int, str]:
    sample = str(sample)
    age_rank = 0 if sample.upper().endswith("-Y") else 1 if sample.upper().endswith("-O") else 2
    return age_rank, sample


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


def load_cell_table() -> pd.DataFrame:
    adata = read_h5ad_compat(INPUT_H5AD)
    obs = adata.obs.copy()
    if "x" not in obs.columns or "y" not in obs.columns:
        if "spatial" not in adata.obsm:
            raise KeyError("Expected x/y in adata.obs or adata.obsm['spatial'].")
        obs["x"] = adata.obsm["spatial"][:, 0]
        obs["y"] = adata.obsm["spatial"][:, 1]

    required = {"cell_id", "sample", "age", "fov", "x", "y", SCORE_COL, HIGH_COL, CELL_TYPE_COL}
    missing = sorted(required - set(obs.columns))
    if missing:
        raise KeyError(f"Missing required h5ad obs columns: {missing}")

    optional = [col for col in ["snc", "sasp_num", "cell_type_major", "immcelltype"] if col in obs.columns]
    df = obs[list(required) + optional].copy()
    df["cell_id"] = df["cell_id"].astype(str)
    df["sample"] = df["sample"].astype(str)
    df["age"] = df["age"].astype(str)
    df["fov"] = df["fov"].astype(str)
    df[CELL_TYPE_COL] = df[CELL_TYPE_COL].astype(str).fillna("Unknown")
    df[SCORE_COL] = pd.to_numeric(df[SCORE_COL], errors="coerce")
    df[HIGH_COL] = df[HIGH_COL].astype(bool)
    df["x"] = pd.to_numeric(df["x"], errors="coerce")
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    if "snc" in df.columns:
        df["is_snc_plus"] = df["snc"].astype(str).eq("SnC+")
    else:
        df["is_snc_plus"] = False
    return df.dropna(subset=["x", "y", SCORE_COL]).copy()


def count_neighbors_from_lists(
    neighbor_lists: list[list[int]],
    anchor_idx: np.ndarray,
    cell_label_codes: np.ndarray,
    n_labels: int,
) -> tuple[np.ndarray, int, int]:
    counts = np.zeros(n_labels, dtype=float)
    total_neighbors = 0
    if len(anchor_idx) == 0:
        return counts, 0, 0
    for anchor_position in anchor_idx:
        neighbors = [idx for idx in neighbor_lists[anchor_position] if idx != anchor_position]
        if not neighbors:
            continue
        counts += np.bincount(cell_label_codes[neighbors], minlength=n_labels)
        total_neighbors += len(neighbors)
    return counts, total_neighbors, int(len(anchor_idx))


def run_neighborhood_analysis(cell_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(RANDOM_SEED)
    labels = sorted(cell_df[CELL_TYPE_COL].dropna().astype(str).unique())
    label_to_idx = {label: idx for idx, label in enumerate(labels)}
    rows = []
    observed_long = []

    for sample, sample_df in cell_df.groupby("sample", sort=False):
        age = str(sample_df["age"].iloc[0])
        observed_counts = np.zeros(len(labels), dtype=float)
        observed_total = 0
        n_anchors_total = 0
        expected_counts = np.zeros((N_PERMUTATIONS, len(labels)), dtype=float)
        expected_totals = np.zeros(N_PERMUTATIONS, dtype=float)

        for fov, group in sample_df.groupby("fov", sort=False):
            group = group.reset_index(drop=True)
            coords = group[["x", "y"]].to_numpy(dtype=float)
            tree = cKDTree(coords)
            neighbor_lists = tree.query_ball_point(coords, r=RADIUS_UM)
            label_codes = group[CELL_TYPE_COL].astype(str).map(label_to_idx).to_numpy(dtype=int)
            high_mask = group[HIGH_COL].to_numpy(dtype=bool)
            n_high = int(high_mask.sum())
            obs_counts, obs_total, n_anchors = count_neighbors_from_lists(
                neighbor_lists,
                np.flatnonzero(high_mask),
                label_codes,
                len(labels),
            )
            observed_counts += obs_counts
            observed_total += obs_total
            n_anchors_total += n_anchors

            for i in range(N_PERMUTATIONS):
                if n_high == 0:
                    random_mask = np.zeros(len(group), dtype=bool)
                else:
                    random_anchor_idx = rng.choice(len(group), size=n_high, replace=False)
                    random_mask = np.zeros(len(group), dtype=bool)
                    random_mask[random_anchor_idx] = True
                perm_counts, perm_total, _ = count_neighbors_from_lists(
                    neighbor_lists,
                    np.flatnonzero(random_mask),
                    label_codes,
                    len(labels),
                )
                expected_counts[i, :] += perm_counts
                expected_totals[i] += perm_total

        observed_fraction = observed_counts / observed_total if observed_total else np.full(len(labels), np.nan)
        expected_fraction = np.divide(
            expected_counts,
            expected_totals[:, None],
            out=np.full_like(expected_counts, np.nan, dtype=float),
            where=expected_totals[:, None] > 0,
        )
        expected_mean_fraction = np.nanmean(expected_fraction, axis=0)
        log2_enrichment = np.log2((observed_fraction + 1e-9) / (expected_mean_fraction + 1e-9))
        empirical_p_high = np.array(
            [
                (np.nansum(expected_fraction[:, j] >= observed_fraction[j]) + 1) / (np.sum(np.isfinite(expected_fraction[:, j])) + 1)
                for j in range(len(labels))
            ]
        )

        for label, obs_count, obs_frac, exp_frac, l2e, pval in zip(
            labels, observed_counts, observed_fraction, expected_mean_fraction, log2_enrichment, empirical_p_high
        ):
            rows.append(
                {
                    "sample": sample,
                    "age": age,
                    "anchor_definition": HIGH_COL,
                    "radius_um": RADIUS_UM,
                    "n_sssp_high_anchors": n_anchors_total,
                    "neighbor_cell_type_minor": label,
                    "observed_neighbor_cells": int(obs_count),
                    "observed_total_neighbor_cells": int(observed_total),
                    "observed_fraction": obs_frac,
                    "expected_mean_fraction_random_anchors": exp_frac,
                    "log2_enrichment_vs_random_anchors": l2e,
                    "empirical_p_enriched": pval,
                    "n_permutations": N_PERMUTATIONS,
                }
            )
            if obs_count > 0:
                observed_long.append(
                    {
                        "sample": sample,
                        "age": age,
                        "neighbor_cell_type_minor": label,
                        "observed_neighbor_cells": int(obs_count),
                    }
                )

    result = pd.DataFrame(rows)
    observed = pd.DataFrame(observed_long)
    result.to_csv(OUTDIR / "xenium_gtex_sssp_high_35um_neighbor_cell_type_enrichment_by_sample.csv", index=False)
    observed.to_csv(OUTDIR / "xenium_gtex_sssp_high_35um_observed_neighbor_counts_by_sample.csv", index=False)
    return result, observed


def plot_neighborhood_heatmap(result: pd.DataFrame) -> None:
    top = (
        result.groupby("neighbor_cell_type_minor")["observed_neighbor_cells"]
        .sum()
        .sort_values(ascending=False)
        .head(18)
        .index.tolist()
    )
    plot = result.loc[result["neighbor_cell_type_minor"].isin(top)].copy()
    sample_order = sorted(plot["sample"].unique(), key=sample_sort_key)
    pivot = plot.pivot(index="neighbor_cell_type_minor", columns="sample", values="log2_enrichment_vs_random_anchors")
    pivot = pivot.reindex(top).reindex(columns=sample_order)

    fig, ax = plt.subplots(figsize=(9.5, 6.8))
    im = ax.imshow(pivot.to_numpy(dtype=float), cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
    ax.set_xticks(np.arange(len(sample_order)), sample_order, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(top)), top)
    ax.set_title("Cell types within 35 um of SSSP-high anchors")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Neighbor minor cell type")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("log2 enrichment vs random anchors")
    fig.tight_layout()
    fig.savefig(OUTDIR / "xenium_gtex_sssp_high_35um_neighbor_enrichment_heatmap.png", dpi=300)
    fig.savefig(OUTDIR / "xenium_gtex_sssp_high_35um_neighbor_enrichment_heatmap.pdf")
    plt.close(fig)


def axial_round(q: np.ndarray, r: np.ndarray) -> pd.DataFrame:
    x = q
    z = r
    y = -x - z
    rx = np.round(x)
    ry = np.round(y)
    rz = np.round(z)
    x_diff = np.abs(rx - x)
    y_diff = np.abs(ry - y)
    z_diff = np.abs(rz - z)
    use_x = (x_diff > y_diff) & (x_diff > z_diff)
    use_y = (~use_x) & (y_diff > z_diff)
    use_z = (~use_x) & (~use_y)
    rx[use_x] = -ry[use_x] - rz[use_x]
    ry[use_y] = -rx[use_y] - rz[use_y]
    rz[use_z] = -rx[use_z] - ry[use_z]
    return pd.DataFrame({"hex_q": rx.astype(int), "hex_r": rz.astype(int)})


def assign_hexes(df: pd.DataFrame, size_um: float) -> pd.DataFrame:
    out = df.copy()
    x0 = out["x"].min()
    y0 = out["y"].min()
    x_shift = out["x"].to_numpy(dtype=float) - x0
    y_shift = out["y"].to_numpy(dtype=float) - y0
    q_float = (np.sqrt(3.0) / 3.0 * x_shift - y_shift / 3.0) / size_um
    r_float = (2.0 / 3.0 * y_shift) / size_um
    axial = axial_round(q_float, r_float)
    out["hex_q"] = axial["hex_q"].to_numpy()
    out["hex_r"] = axial["hex_r"].to_numpy()
    out["hex_x"] = size_um * np.sqrt(3.0) * (out["hex_q"] + out["hex_r"] / 2.0) + x0
    out["hex_y"] = size_um * 1.5 * out["hex_r"] + y0
    out["hex_id"] = out["fov"].astype(str) + "__" + out["hex_q"].astype(str) + "__" + out["hex_r"].astype(str)
    return out


def build_hex_bins(cell_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cell_hex = pd.concat(
        [assign_hexes(group, HEX_SIZE_UM) for _, group in cell_df.groupby("fov", sort=False)],
        ignore_index=True,
    )
    hex_bins = (
        cell_hex.groupby(["fov", "sample", "age", "hex_id", "hex_q", "hex_r", "hex_x", "hex_y"], observed=True)
        .agg(
            n_cells=("cell_id", "size"),
            n_sssp_high=(HIGH_COL, "sum"),
            mean_sssp_score=(SCORE_COL, "mean"),
            median_sssp_score=(SCORE_COL, "median"),
            n_snc=("is_snc_plus", "sum"),
        )
        .reset_index()
    )
    hex_bins["sssp_high_prop"] = hex_bins["n_sssp_high"] / hex_bins["n_cells"]
    hex_bins["snc_prop"] = hex_bins["n_snc"] / hex_bins["n_cells"]
    hex_bins = hex_bins.loc[hex_bins["n_cells"] >= MIN_CELLS_PER_HEX].copy()
    return cell_hex, hex_bins


def run_muspan_getis_ord_for_group(ms, df: pd.DataFrame, metric: str) -> pd.DataFrame:
    out = df.sort_values(["hex_q", "hex_r"]).reset_index(drop=True).copy()
    out["gi_star_z"] = np.nan
    out["gi_star_p_adj_muspan"] = np.nan
    values = pd.to_numeric(out[metric], errors="coerce").to_numpy(dtype=float)
    valid = np.isfinite(values)
    if valid.sum() < 30 or np.nanstd(values) == 0:
        return out

    domain = ms.domain(f"{out.loc[0, 'fov']} {metric}")
    domain.add_points(out.loc[valid, ["hex_x", "hex_y"]].to_numpy(dtype=float), "Hex bins")
    domain.add_labels(metric, values[valid])
    zscores, pvals_adj, object_indices = ms.spatial_statistics.getis_ord(
        domain,
        population=("Collection", "Hex bins"),
        label_name=metric,
        alpha=ALPHA,
        network_kwargs={
            "network_type": "Proximity",
            "min_edge_distance": 0,
            "max_edge_distance": HEX_ADJACENCY_DISTANCE_UM,
        },
        add_local_value_as_label=False,
        local_getis_label_name=f"Local Getis Ord - {metric}",
    )
    valid_positions = np.flatnonzero(valid)
    out_positions = valid_positions[np.asarray(object_indices, dtype=int)]
    out.loc[out_positions, "gi_star_z"] = np.asarray(zscores, dtype=float)
    out.loc[out_positions, "gi_star_p_adj_muspan"] = np.asarray(pvals_adj, dtype=float)
    return out


def classify_hotspots(results: pd.DataFrame) -> pd.DataFrame:
    out = results.copy()
    hot = (out["gi_star_p_adj_muspan"] < ALPHA) & (out["gi_star_z"] > 0)
    cold = (out["gi_star_p_adj_muspan"] < ALPHA) & (out["gi_star_z"] < 0)
    out["hotspot_class"] = "Not significant"
    out.loc[hot, "hotspot_class"] = "Hot spot"
    out.loc[cold, "hotspot_class"] = "Cold spot"
    out["gi_star_fdr"] = out["gi_star_p_adj_muspan"]
    return out


def run_hotspots(cell_hex: pd.DataFrame, hex_bins: pd.DataFrame) -> pd.DataFrame | None:
    try:
        ms = require_muspan()
    except RuntimeError as exc:
        (OUTDIR / "muspan_hotspots_not_run.txt").write_text(str(exc) + "\n", encoding="utf-8")
        print(str(exc))
        return None

    results = []
    for metric in HOTSPOT_METRICS:
        print(f"Running MuSpAn Getis-Ord for {metric}")
        for fov, group in hex_bins.groupby("fov", sort=False):
            g = run_muspan_getis_ord_for_group(ms, group, metric)
            g["metric"] = metric
            results.append(g)
    results_long = classify_hotspots(pd.concat(results, ignore_index=True))
    results_long.to_csv(OUTDIR / "xenium_gtex_sssp_muspan_hotspot_results_long.csv", index=False)

    summary = (
        results_long.groupby(["metric", "sample", "age", "hotspot_class"], observed=True)
        .agg(
            n_hex=("hex_id", "size"),
            n_cells=("n_cells", "sum"),
            n_sssp_high=("n_sssp_high", "sum"),
            median_gi_star_z=("gi_star_z", "median"),
            mean_sssp_score=("mean_sssp_score", "mean"),
            mean_sssp_high_prop=("sssp_high_prop", "mean"),
        )
        .reset_index()
    )
    sample_totals = results_long.groupby(["metric", "sample"], observed=True)["hex_id"].nunique().rename("total_hex")
    summary = summary.merge(sample_totals.reset_index(), on=["metric", "sample"], how="left")
    summary["hex_fraction"] = summary["n_hex"] / summary["total_hex"]
    summary.to_csv(OUTDIR / "xenium_gtex_sssp_muspan_hotspot_summary.csv", index=False)

    comp = cell_hex.merge(
        results_long[["metric", "fov", "hex_id", "hotspot_class"]],
        on=["fov", "hex_id"],
        how="inner",
    )
    comp = (
        comp.loc[comp["hotspot_class"].eq("Hot spot")]
        .groupby(["metric", "sample", "age", CELL_TYPE_COL], observed=True)
        .agg(n_cells=("cell_id", "nunique"), n_sssp_high=(HIGH_COL, "sum"), mean_sssp_score=(SCORE_COL, "mean"))
        .reset_index()
    )
    totals = comp.groupby(["metric", "sample", "age"], observed=True)["n_cells"].sum().rename("total_hotspot_cells")
    comp = comp.merge(totals.reset_index(), on=["metric", "sample", "age"], how="left")
    comp["fraction_hotspot_cells"] = comp["n_cells"] / comp["total_hotspot_cells"]
    comp.to_csv(OUTDIR / "xenium_gtex_sssp_muspan_hotspot_cell_type_minor_composition.csv", index=False)

    for metric in HOTSPOT_METRICS:
        plot_hotspots(results_long, metric)
    return results_long


def plot_hotspots(results: pd.DataFrame, metric: str) -> None:
    df = results.loc[results["metric"].eq(metric)].copy()
    samples = sorted(df["sample"].unique(), key=sample_sort_key)
    ncols = 3
    nrows = int(np.ceil(len(samples) / ncols))
    fig, axes = plt.subplots(nrows * 2, ncols, figsize=(4.1 * ncols, 3.7 * nrows * 2))
    axes = np.asarray(axes).reshape(nrows * 2, ncols)
    vmax = float(np.nanmax(np.abs(df["gi_star_z"]))) if df["gi_star_z"].notna().any() else 1.0
    vmax = max(vmax, 1.0)
    colors = {"Cold spot": "#2C7BB6", "Not significant": "grey84", "Hot spot": "#D7191C"}

    for idx, sample in enumerate(samples):
        row = (idx // ncols) * 2
        col = idx % ncols
        sdf = df.loc[df["sample"].eq(sample)]
        ax = axes[row, col]
        sca = ax.scatter(sdf["hex_x"], sdf["hex_y"], c=sdf["gi_star_z"], s=5, cmap="RdBu_r", vmin=-vmax, vmax=vmax, linewidths=0)
        ax.set_title(sample, fontsize=10, fontweight="bold")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.axis("off")

        ax = axes[row + 1, col]
        for cls in ["Cold spot", "Not significant", "Hot spot"]:
            cdf = sdf.loc[sdf["hotspot_class"].eq(cls)]
            if len(cdf):
                ax.scatter(cdf["hex_x"], cdf["hex_y"], c=colors[cls], s=5, linewidths=0)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.axis("off")

    for ax in axes.ravel()[len(samples) * 2 :]:
        ax.axis("off")

    fig.suptitle(f"MuSpAn Getis-Ord G*: {metric}", fontsize=15, fontweight="bold")
    cbar = fig.colorbar(sca, ax=axes[::2, :].ravel().tolist(), shrink=0.55)
    cbar.set_label("Gi* z")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUTDIR / f"xenium_gtex_sssp_muspan_hotspots_{metric}.png", dpi=300)
    fig.savefig(OUTDIR / f"xenium_gtex_sssp_muspan_hotspots_{metric}.pdf")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    cell_df = load_cell_table()
    cell_df.to_csv(OUTDIR / "xenium_gtex_sssp_spatial_cell_table.csv.gz", index=False)

    print(f"Cells available for spatial analysis: {len(cell_df):,}")
    neighborhood, observed = run_neighborhood_analysis(cell_df)
    plot_neighborhood_heatmap(neighborhood)

    cell_hex, hex_bins = build_hex_bins(cell_df)
    cell_hex.to_csv(OUTDIR / "xenium_gtex_sssp_cell_hex_assignments.csv.gz", index=False)
    hex_bins.to_csv(OUTDIR / "xenium_gtex_sssp_hotspot_hex_bins.csv", index=False)
    print(f"Hex bins retained for hotspots: {len(hex_bins):,}")
    run_hotspots(cell_hex, hex_bins)
    print(f"Done. Results written to {OUTDIR}")


if __name__ == "__main__":
    main()
