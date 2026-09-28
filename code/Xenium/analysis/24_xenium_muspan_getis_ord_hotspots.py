"""
Run Xenium SnC/SASP hotspot analysis with MuSpAn Getis-Ord G*.

This is the h5ad/Python analogue of 03_xenium_snc_getis_ord_hotspots.R.
It keeps the same biological/statistical unit used in the earlier analysis:
35 um hex bins within each FOV, with one continuous value per hex bin.

Inputs:
  results/260727/xenium_sssp_20gene_scanpy_score_genes/
    xenium_20gene_sssp_scanpy_scored.h5ad

Outputs:
  results/260727/xenium_muspan_snc_hotspots/
    xenium_muspan_snc_hotspot_bins.csv
    xenium_muspan_snc_hotspot_results_long.csv
    xenium_muspan_snc_hotspot_summary.csv
    xenium_muspan_snc_hotspots_<metric>.png/.pdf

Note:
  MuSpAn is not the small PyPI bootstrap package. The real package is
  credential-protected by the MuSpAn authors. Install it in senoquant-dev with:

    python -m pip install -U https://docs.muspan.co.uk/code/latest.zip
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


PROJECT_ROOT = Path(__file__).resolve().parents[3]
INPUT_H5AD = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_sssp_20gene_scanpy_score_genes"
    / "xenium_20gene_sssp_scanpy_scored.h5ad"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_muspan_snc_hotspots"

HEX_SIZE_UM = 35.0
HEX_NEIGHBOR_ORDER = 1
MIN_CELLS_PER_HEX = 5
ALPHA = 0.05
METRICS_TO_TEST = ("snc_prop", "mean_sasp_num")
HEX_ADJACENCY_DISTANCE_UM = np.sqrt(3.0) * HEX_SIZE_UM * HEX_NEIGHBOR_ORDER + 1e-6


def require_muspan():
    try:
        import muspan as ms
    except Exception as exc:
        raise RuntimeError(
            "\nReal MuSpAn is required for this script.\n"
            "The official docs install command is:\n"
            "  python -m pip install -U https://docs.muspan.co.uk/code/latest.zip\n"
            "That installer asks for the username/password supplied by MuSpAn.\n"
            "The tiny PyPI package is only a bootstrapper and is not enough.\n"
        ) from exc

    required = ("domain", "spatial_statistics")
    missing = [name for name in required if not hasattr(ms, name)]
    if missing:
        raise RuntimeError(f"Installed muspan is missing expected API: {missing}")
    if not hasattr(ms.spatial_statistics, "getis_ord"):
        raise RuntimeError("Installed muspan.spatial_statistics has no getis_ord function.")
    return ms


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
    out["hex_id"] = (
        out["fov"].astype(str)
        + "__"
        + out["hex_q"].astype(str)
        + "__"
        + out["hex_r"].astype(str)
    )
    return out


def load_cell_table() -> pd.DataFrame:
    adata = sc.read_h5ad(INPUT_H5AD)
    obs = adata.obs.copy()

    if "x" not in obs.columns or "y" not in obs.columns:
        if "spatial" not in adata.obsm:
            raise KeyError("Expected x/y in adata.obs or adata.obsm['spatial'].")
        obs["x"] = adata.obsm["spatial"][:, 0]
        obs["y"] = adata.obsm["spatial"][:, 1]

    required = {"cell_id", "sample", "age", "fov", "x", "y", "snc", "sasp_num"}
    missing = sorted(required - set(obs.columns))
    if missing:
        raise KeyError(f"Missing required h5ad obs columns: {missing}")

    cell_df = obs[list(required)].copy()
    cell_df["cell_id"] = cell_df["cell_id"].astype(str)
    cell_df["sample"] = cell_df["sample"].astype(str)
    cell_df["age"] = cell_df["age"].astype(str)
    cell_df["fov"] = cell_df["fov"].astype(str)
    cell_df["snc"] = cell_df["snc"].astype(str)
    cell_df["x"] = pd.to_numeric(cell_df["x"], errors="coerce")
    cell_df["y"] = pd.to_numeric(cell_df["y"], errors="coerce")
    cell_df["sasp_num"] = pd.to_numeric(cell_df["sasp_num"], errors="coerce")
    cell_df = cell_df.dropna(subset=["x", "y", "sasp_num", "snc"])
    return cell_df


def build_hex_bins(cell_df: pd.DataFrame) -> pd.DataFrame:
    cell_hex = pd.concat(
        [assign_hexes(group, HEX_SIZE_UM) for _, group in cell_df.groupby("fov", sort=False)],
        ignore_index=True,
    )

    hex_bins = (
        cell_hex.groupby(["fov", "sample", "age", "hex_id", "hex_q", "hex_r", "hex_x", "hex_y"])
        .agg(
            n_cells=("cell_id", "size"),
            n_snc=("snc", lambda x: int(np.sum(np.asarray(x) == "SnC+"))),
            mean_sasp_num=("sasp_num", "mean"),
            sum_sasp_num=("sasp_num", "sum"),
        )
        .reset_index()
    )
    hex_bins["snc_prop"] = hex_bins["n_snc"] / hex_bins["n_cells"]
    hex_bins = hex_bins.loc[hex_bins["n_cells"] >= MIN_CELLS_PER_HEX].copy()
    return hex_bins


def run_muspan_getis_ord_for_group(ms, df: pd.DataFrame, metric: str) -> pd.DataFrame:
    out = df.sort_values(["hex_q", "hex_r"]).reset_index(drop=True).copy()
    out["gi_star_z"] = np.nan
    out["gi_star_p_adj_muspan"] = np.nan
    out["n_neighbors_expected"] = np.nan

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
    object_indices = np.asarray(object_indices, dtype=int)
    if len(object_indices) != len(zscores):
        raise RuntimeError("MuSpAn returned mismatched object_indices and z-score lengths.")

    out_positions = valid_positions[object_indices]
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

    # Compatibility alias for older downstream scripts; this is Bonferroni-adjusted
    # by MuSpAn, not BH FDR from the older R implementation.
    out["gi_star_fdr"] = out["gi_star_p_adj_muspan"]
    return out


def plot_hotspots(results: pd.DataFrame, metric: str) -> None:
    df = results.loc[results["metric"].eq(metric)].copy()
    classes = ["Cold spot", "Not significant", "Hot spot"]
    colors = {"Cold spot": "#2C7BB6", "Not significant": "grey84", "Hot spot": "#D7191C"}

    samples = sorted(df["sample"].unique())
    ncols = 3
    nrows = int(np.ceil(len(samples) / ncols))
    fig, axes = plt.subplots(nrows * 2, ncols, figsize=(4.1 * ncols, 3.7 * nrows * 2))
    axes = np.asarray(axes).reshape(nrows * 2, ncols)

    vmax = float(np.nanmax(np.abs(df["gi_star_z"]))) if df["gi_star_z"].notna().any() else 1.0
    vmax = max(vmax, 1.0)

    for idx, sample in enumerate(samples):
        row = (idx // ncols) * 2
        col = idx % ncols
        sdf = df.loc[df["sample"].eq(sample)]

        ax = axes[row, col]
        sca = ax.scatter(
            sdf["hex_x"],
            sdf["hex_y"],
            c=sdf["gi_star_z"],
            s=5,
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            linewidths=0,
        )
        ax.set_title(sample, fontsize=10, fontweight="bold")
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.axis("off")

        ax = axes[row + 1, col]
        for cls in classes:
            cdf = sdf.loc[sdf["hotspot_class"].eq(cls)]
            if len(cdf):
                ax.scatter(cdf["hex_x"], cdf["hex_y"], c=colors[cls], s=5, label=cls, linewidths=0)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.axis("off")

    for ax in axes.ravel()[len(samples) * 2 :]:
        ax.axis("off")

    fig.suptitle(f"MuSpAn Getis-Ord G*: {metric}", fontsize=15, fontweight="bold")
    cbar = fig.colorbar(sca, ax=axes[::2, :].ravel().tolist(), shrink=0.55)
    cbar.set_label("Gi* z")
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", label=cls, markerfacecolor=colors[cls], markersize=6)
        for cls in classes
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))

    safe_metric = "".join(ch if ch.isalnum() or ch == "_" else "_" for ch in metric)
    fig.savefig(OUTDIR / f"xenium_muspan_snc_hotspots_{safe_metric}.png", dpi=300)
    fig.savefig(OUTDIR / f"xenium_muspan_snc_hotspots_{safe_metric}.pdf")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    ms = require_muspan()

    print(f"Reading {INPUT_H5AD}")
    cell_df = load_cell_table()
    print(f"Cells with coordinates and metadata: {len(cell_df):,}")

    hex_bins = build_hex_bins(cell_df)
    hex_bins.to_csv(OUTDIR / "xenium_muspan_snc_hotspot_bins.csv", index=False)
    print(f"Hex bins retained: {len(hex_bins):,}")

    results = []
    for metric in METRICS_TO_TEST:
        print(f"Running MuSpAn Getis-Ord for {metric}")
        for fov, group in hex_bins.groupby("fov", sort=False):
            print(f"  {fov}: {len(group):,} hex bins")
            g = run_muspan_getis_ord_for_group(ms, group, metric)
            g["metric"] = metric
            results.append(g)

    results_long = classify_hotspots(pd.concat(results, ignore_index=True))
    results_long.to_csv(OUTDIR / "xenium_muspan_snc_hotspot_results_long.csv", index=False)

    summary = (
        results_long.groupby(["metric", "fov", "sample", "age", "hotspot_class"], dropna=False)
        .agg(
            n_hex=("hex_id", "size"),
            n_cells=("n_cells", "sum"),
            n_snc=("n_snc", "sum"),
            median_gi_star_z=("gi_star_z", "median"),
        )
        .reset_index()
        .sort_values(["metric", "fov", "hotspot_class"])
    )
    summary.to_csv(OUTDIR / "xenium_muspan_snc_hotspot_summary.csv", index=False)

    for metric in METRICS_TO_TEST:
        plot_hotspots(results_long, metric)

    print(f"Done. Results written to {OUTDIR}")


if __name__ == "__main__":
    main()
