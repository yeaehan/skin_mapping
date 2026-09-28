"""
Compare Xenium SnC calls with the 20-gene Scanpy SSSP score.

Questions:
1. Is SSSP score higher in SnC+ cells than SnC- cells?
2. Are SnC+ cells spatially close to SSSP-high cells?
3. Are SSSP-high cells spatially close to SnC+ cells?

SSSP-high is defined as the top quartile of scanpy.tl.score_genes scores.
Spatial distances are computed within each sample/FOV only.
Permutation tests shuffle SSSP-high labels within each sample/FOV, preserving
the number of SSSP-high cells per local field of view.
"""

from __future__ import annotations

from pathlib import Path
import warnings

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import squidpy as sq
from scipy import sparse, stats
from scipy.spatial import cKDTree


PROJECT_ROOT = Path(__file__).resolve().parents[3]
SSSP_CELL_SCORES = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_sssp_20gene_scanpy_score_genes"
    / "xenium_20gene_sssp_scanpy_cell_scores.csv.gz"
)
COORDS_WITH_SNC = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_snc_hotspot_neighborhoods"
    / "xenium_snc_hotspot_cell_neighborhoods.csv"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_snc_vs_sssp_score_proximity"
SCORE_COL = "sssp_20gene_scanpy_score"
SSSP_HIGH_COL = "sssp_high_top_quartile"
SSSP_HIGH_QUANTILE = 0.75
GENE_SET_CONFIGS = [
    {
        "gene_set_id": "main_20",
        "cell_scores": SSSP_CELL_SCORES,
        "outdir": OUTDIR,
        "score_col": "sssp_20gene_scanpy_score",
        "high_col": "sssp_high_top_quartile",
        "score_label": "20-gene SSSP score",
    },
]
AGE_ORDER = ["Young", "Old"]
SNC_ORDER = ["SnC-", "SnC+"]
AGE_COLORS = {"Young": "#4C78A8", "Old": "#F58518"}
SNC_COLORS = {"SnC-": "#9AA0A6", "SnC+": "#D62728"}
RADII_UM = [35, 70, 140]
DISTANCE_BINS_UM = [0, 35, 70, 140, 280, 560, np.inf]
DISTANCE_LABELS = ["0-35", "35-70", "70-140", "140-280", "280-560", ">560"]
N_PERMUTATIONS = 1000
RANDOM_SEED = 260727


def sem(values: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if len(values) <= 1:
        return np.nan
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def safe_welch_p(a: pd.Series, b: pd.Series) -> float:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) < 2 or len(b) < 2:
        return np.nan
    if np.isclose(a.std(ddof=1), 0) and np.isclose(b.std(ddof=1), 0):
        return np.nan
    return float(stats.ttest_ind(a, b, equal_var=False).pvalue)


def sample_sort_key(sample: str) -> tuple[int, str]:
    sample = str(sample)
    age_rank = 0 if sample.upper().endswith("-Y") else 1 if sample.upper().endswith("-O") else 2
    return age_rank, sample


def load_data() -> pd.DataFrame:
    scores = pd.read_csv(SSSP_CELL_SCORES)
    cutoff = float(scores[SCORE_COL].quantile(SSSP_HIGH_QUANTILE))
    scores[SSSP_HIGH_COL] = scores[SCORE_COL] >= cutoff
    keep_score_cols = [
        "cell_id",
        SCORE_COL,
        SSSP_HIGH_COL,
        "panel_mean_log1p_cp10k",
        "n_positive_panel_genes",
        "fraction_positive_panel_genes",
    ]
    scores = scores[keep_score_cols].drop_duplicates("cell_id")

    coords = pd.read_csv(COORDS_WITH_SNC)
    required = {"cell_id", "sample", "age", "fov", "x", "y", "snc"}
    missing = sorted(required - set(coords.columns))
    if missing:
        raise KeyError(f"Missing required columns in {COORDS_WITH_SNC}: {', '.join(missing)}")

    df = coords.merge(scores, on="cell_id", how="inner")
    df = df[df["snc"].isin(SNC_ORDER)].copy()
    df["snc"] = pd.Categorical(df["snc"], categories=SNC_ORDER, ordered=True)
    df["age"] = pd.Categorical(df["age"], categories=AGE_ORDER, ordered=True)
    df.attrs["sssp_high_cutoff"] = cutoff
    return df


def nearest_other_target_distance(group: pd.DataFrame, target_col: str) -> np.ndarray:
    coords = group[["x", "y"]].to_numpy(dtype=float)
    target_mask = group[target_col].to_numpy(dtype=bool)
    target_positions = np.flatnonzero(target_mask)
    out = np.full(len(group), np.nan, dtype=float)
    if len(target_positions) == 0:
        return out

    tree = cKDTree(coords[target_positions])
    k = min(2, len(target_positions))
    distances, target_rank = tree.query(coords, k=k)
    if k == 1:
        distances = distances.reshape(-1, 1)
        target_rank = target_rank.reshape(-1, 1)

    for i in range(len(group)):
        for rank in range(k):
            candidate_original_index = target_positions[target_rank[i, rank]]
            if candidate_original_index != i:
                out[i] = distances[i, rank]
                break
    return out


def add_proximity(df: pd.DataFrame) -> pd.DataFrame:
    pieces = []
    df = df.copy()
    df["is_snc_plus"] = df["snc"].astype(str).eq("SnC+")
    for _, group in df.groupby(["sample", "fov"], observed=True, sort=False):
        group = group.copy()
        group["nearest_other_sssp_high_um"] = nearest_other_target_distance(group, SSSP_HIGH_COL)
        group["nearest_other_snc_plus_um"] = nearest_other_target_distance(group, "is_snc_plus")
        for radius in RADII_UM:
            group[f"near_other_sssp_high_{radius}um"] = group["nearest_other_sssp_high_um"].le(radius)
            group[f"near_other_snc_plus_{radius}um"] = group["nearest_other_snc_plus_um"].le(radius)
        pieces.append(group)
    return pd.concat(pieces, ignore_index=True)


def summarize_score_by_snc(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    by_sample = (
        df.groupby(["sample", "age", "snc"], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            mean_sssp_score=(SCORE_COL, "mean"),
            median_sssp_score=(SCORE_COL, "median"),
            sem_sssp_score=(SCORE_COL, sem),
            fraction_sssp_high=(SSSP_HIGH_COL, "mean"),
            mean_n_positive_panel_genes=("n_positive_panel_genes", "mean"),
        )
        .reset_index()
    )
    wide = by_sample.pivot(index=["sample", "age"], columns="snc", values=["mean_sssp_score", "fraction_sssp_high"])
    wide.columns = [f"{metric}_{snc}" for metric, snc in wide.columns]
    wide = wide.reset_index()
    wide["delta_mean_sssp_score_snc_plus_minus_snc_minus"] = wide["mean_sssp_score_SnC+"] - wide["mean_sssp_score_SnC-"]
    wide["delta_fraction_sssp_high_snc_plus_minus_snc_minus"] = wide["fraction_sssp_high_SnC+"] - wide["fraction_sssp_high_SnC-"]

    age_summary = (
        wide.groupby("age", observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_delta_score=("delta_mean_sssp_score_snc_plus_minus_snc_minus", "mean"),
            sem_delta_score=("delta_mean_sssp_score_snc_plus_minus_snc_minus", sem),
            mean_delta_fraction_high=("delta_fraction_sssp_high_snc_plus_minus_snc_minus", "mean"),
            sem_delta_fraction_high=("delta_fraction_sssp_high_snc_plus_minus_snc_minus", sem),
        )
        .reset_index()
    )
    by_sample.to_csv(OUTDIR / "snc_status_sssp_score_by_sample.csv", index=False)
    wide.to_csv(OUTDIR / "snc_status_sssp_score_sample_deltas.csv", index=False)
    age_summary.to_csv(OUTDIR / "snc_status_sssp_score_delta_by_age.csv", index=False)
    return by_sample, wide, age_summary


def summarize_proximity(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for (sample, age, snc), sub in df.groupby(["sample", "age", "snc"], observed=True):
        row = {
            "sample": sample,
            "age": age,
            "snc": snc,
            "n_cells": sub["cell_id"].nunique(),
            "median_nearest_other_sssp_high_um": sub["nearest_other_sssp_high_um"].median(),
            "mean_nearest_other_sssp_high_um": sub["nearest_other_sssp_high_um"].mean(),
        }
        for radius in RADII_UM:
            row[f"fraction_near_other_sssp_high_{radius}um"] = sub[f"near_other_sssp_high_{radius}um"].mean()
        rows.append(row)
    snc_to_sssp = pd.DataFrame(rows)

    rows = []
    for (sample, age, sssp_high), sub in df.groupby(["sample", "age", SSSP_HIGH_COL], observed=True):
        row = {
            "sample": sample,
            "age": age,
            SSSP_HIGH_COL: bool(sssp_high),
            "n_cells": sub["cell_id"].nunique(),
            "median_nearest_other_snc_plus_um": sub["nearest_other_snc_plus_um"].median(),
            "mean_nearest_other_snc_plus_um": sub["nearest_other_snc_plus_um"].mean(),
        }
        for radius in RADII_UM:
            row[f"fraction_near_other_snc_plus_{radius}um"] = sub[f"near_other_snc_plus_{radius}um"].mean()
        rows.append(row)
    sssp_to_snc = pd.DataFrame(rows)

    snc_to_sssp.to_csv(OUTDIR / "snc_status_proximity_to_sssp_high_top_quartile_by_sample.csv", index=False)
    sssp_to_snc.to_csv(OUTDIR / "sssp_high_top_quartile_status_proximity_to_snc_plus_by_sample.csv", index=False)
    return snc_to_sssp, sssp_to_snc


def sample_level_tests(by_sample: pd.DataFrame, proximity: pd.DataFrame) -> pd.DataFrame:
    rows = []

    score_wide = by_sample.pivot(index=["sample", "age"], columns="snc", values=["mean_sssp_score", "fraction_sssp_high"])
    score_wide.columns = [f"{metric}_{snc}" for metric, snc in score_wide.columns]
    score_wide = score_wide.reset_index()
    for metric in ["mean_sssp_score", "fraction_sssp_high"]:
        plus = score_wide[f"{metric}_SnC+"]
        minus = score_wide[f"{metric}_SnC-"]
        rows.append(
            {
                "comparison": f"SnC+ vs SnC- {metric}",
                "n_samples": int(len(score_wide)),
                "snc_minus_mean": float(minus.mean()),
                "snc_plus_mean": float(plus.mean()),
                "snc_plus_minus_snc_minus": float((plus - minus).mean()),
                "paired_ttest_p": float(stats.ttest_rel(plus, minus).pvalue) if len(score_wide) >= 2 else np.nan,
                "wilcoxon_p": float(stats.wilcoxon(plus, minus).pvalue) if len(score_wide) >= 2 and not np.allclose(plus, minus) else np.nan,
            }
        )

    prox_wide = proximity.pivot(index=["sample", "age"], columns="snc")
    for metric in ["median_nearest_other_sssp_high_um"] + [f"fraction_near_other_sssp_high_{r}um" for r in RADII_UM]:
        plus = prox_wide[(metric, "SnC+")]
        minus = prox_wide[(metric, "SnC-")]
        rows.append(
            {
                "comparison": f"SnC+ vs SnC- {metric}",
                "n_samples": int(len(plus.dropna())),
                "snc_minus_mean": float(minus.mean()),
                "snc_plus_mean": float(plus.mean()),
                "snc_plus_minus_snc_minus": float((plus - minus).mean()),
                "paired_ttest_p": float(stats.ttest_rel(plus, minus, nan_policy="omit").pvalue) if len(plus.dropna()) >= 2 else np.nan,
                "wilcoxon_p": float(stats.wilcoxon(plus, minus).pvalue) if len(plus.dropna()) >= 2 and not np.allclose(plus, minus, equal_nan=True) else np.nan,
            }
        )
    tests = pd.DataFrame(rows)
    tests.to_csv(OUTDIR / "snc_vs_sssp_top_quartile_sample_level_tests.csv", index=False)
    return tests


def permutation_test_sssp_to_snc(
    df: pd.DataFrame,
    n_permutations: int = N_PERMUTATIONS,
    random_seed: int = RANDOM_SEED,
) -> pd.DataFrame:
    rng = np.random.default_rng(random_seed)
    rows = []
    metrics = ["median_nearest_other_snc_plus_um"] + [
        f"fraction_near_other_snc_plus_{radius}um" for radius in RADII_UM
    ]

    for sample, sample_df in df.groupby("sample", observed=True, sort=False):
        sample_df = sample_df.reset_index(drop=True)
        age = sample_df["age"].astype(str).iloc[0]
        observed_high = sample_df[SSSP_HIGH_COL].to_numpy(dtype=bool)
        n_high = int(observed_high.sum())
        if n_high == 0:
            continue

        distance_to_snc = sample_df["nearest_other_snc_plus_um"].to_numpy(dtype=float)
        radius_masks = {
            radius: sample_df[f"near_other_snc_plus_{radius}um"].to_numpy(dtype=bool)
            for radius in RADII_UM
        }

        observed = {
            "median_nearest_other_snc_plus_um": float(np.nanmedian(distance_to_snc[observed_high])),
        }
        for radius, mask in radius_masks.items():
            observed[f"fraction_near_other_snc_plus_{radius}um"] = float(np.nanmean(mask[observed_high]))

        groups = []
        for _, group in sample_df.groupby("fov", observed=True, sort=False):
            idx = group.index.to_numpy()
            n_group_high = int(group[SSSP_HIGH_COL].sum())
            groups.append((idx, n_group_high))

        null_values = {metric: np.full(n_permutations, np.nan, dtype=float) for metric in metrics}
        for perm_i in range(n_permutations):
            perm_high = np.zeros(len(sample_df), dtype=bool)
            for idx, n_group_high in groups:
                if n_group_high <= 0:
                    continue
                perm_high[rng.choice(idx, size=n_group_high, replace=False)] = True

            null_values["median_nearest_other_snc_plus_um"][perm_i] = np.nanmedian(distance_to_snc[perm_high])
            for radius, mask in radius_masks.items():
                null_values[f"fraction_near_other_snc_plus_{radius}um"][perm_i] = np.nanmean(mask[perm_high])

        for metric in metrics:
            null = pd.Series(null_values[metric]).dropna().to_numpy()
            obs = observed[metric]
            if len(null) == 0 or np.isnan(obs):
                p_lower = np.nan
                p_greater = np.nan
                p_two_sided = np.nan
                null_mean = np.nan
                null_sd = np.nan
                z_score = np.nan
            else:
                p_lower = float((np.sum(null <= obs) + 1) / (len(null) + 1))
                p_greater = float((np.sum(null >= obs) + 1) / (len(null) + 1))
                p_two_sided = float(min(1.0, 2 * min(p_lower, p_greater)))
                null_mean = float(np.mean(null))
                null_sd = float(np.std(null, ddof=1)) if len(null) > 1 else np.nan
                z_score = float((obs - null_mean) / null_sd) if null_sd and not np.isclose(null_sd, 0) else np.nan

            rows.append(
                {
                    "sample": sample,
                    "age": age,
                    "metric": metric,
                    "n_cells": int(len(sample_df)),
                    "n_sssp_high": n_high,
                    "n_permutations": n_permutations,
                    "permutation_scope": "SSSP-high labels shuffled within sample/FOV",
                    "observed": obs,
                    "null_mean": null_mean,
                    "observed_minus_null_mean": obs - null_mean if not np.isnan(null_mean) else np.nan,
                    "null_sd": null_sd,
                    "z_score": z_score,
                    "empirical_p_lower_or_equal": p_lower,
                    "empirical_p_greater_or_equal": p_greater,
                    "empirical_p_two_sided": p_two_sided,
                }
            )

    result = pd.DataFrame(rows)
    result.to_csv(OUTDIR / "sssp_high_top_quartile_to_snc_plus_permutation_tests_by_sample.csv", index=False)
    return result


def make_distance_anndata(df: pd.DataFrame) -> ad.AnnData:
    work = df.copy()
    work["library_id"] = work["sample"].astype(str) + "__" + work["fov"].astype(str)
    work["SnC_plus"] = work["snc"].astype(str).eq("SnC+").astype(float)
    work["SSSP_high"] = work[SSSP_HIGH_COL].astype(bool).astype(float)
    work["SSSP_score"] = pd.to_numeric(work[SCORE_COL], errors="coerce")
    work["snc_anchor"] = np.where(work["SnC_plus"].eq(1), "SnC+", "Other")
    work["sssp_anchor"] = np.where(work["SSSP_high"].eq(1), "SSSP-high", "Other")

    var_names = ["SSSP_score", "SSSP_high", "SnC_plus"]
    obs_cols = ["cell_id", "sample", "age", "fov", "library_id", "snc", "snc_anchor", "sssp_anchor"]
    obs = work[obs_cols].copy()
    obs.index = work["cell_id"].astype(str)
    obs["age"] = pd.Categorical(obs["age"], categories=AGE_ORDER, ordered=True)
    obs["snc_anchor"] = pd.Categorical(obs["snc_anchor"], categories=["Other", "SnC+"], ordered=True)
    obs["sssp_anchor"] = pd.Categorical(obs["sssp_anchor"], categories=["Other", "SSSP-high"], ordered=True)

    adata = ad.AnnData(
        X=sparse.csr_matrix(work[var_names].to_numpy(dtype=np.float32)),
        obs=obs,
        var=pd.DataFrame(index=var_names),
    )
    adata.obsm["spatial"] = work[["x", "y"]].to_numpy(dtype=float)
    return adata


def run_var_by_distance(adata: ad.AnnData) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sq.tl.var_by_distance(
            adata=adata,
            groups="SnC+",
            cluster_key="snc_anchor",
            library_key="library_id",
            covariates=["age", "sample", "fov", "snc_anchor", "sssp_anchor"],
            design_matrix_key="distance_to_snc_plus",
            spatial_key="spatial",
        )
        sq.tl.var_by_distance(
            adata=adata,
            groups="SSSP-high",
            cluster_key="sssp_anchor",
            library_key="library_id",
            covariates=["age", "sample", "fov", "snc_anchor", "sssp_anchor"],
            design_matrix_key="distance_to_sssp_high",
            spatial_key="spatial",
        )


def extract_distance_design(adata: ad.AnnData, key: str, anchor: str) -> pd.DataFrame:
    design = adata.obsm[key].copy()
    values = pd.DataFrame(
        adata.X.toarray() if sparse.issparse(adata.X) else adata.X,
        index=adata.obs_names,
        columns=adata.var_names,
    )
    design = design.join(values)
    raw_col = f"{anchor}_raw"
    norm_col = anchor
    if raw_col not in design.columns:
        raise KeyError(f"Expected {raw_col} in adata.obsm[{key!r}]")
    design = design.rename(columns={raw_col: "distance_um", norm_col: "distance_normalized"})
    design["distance_bin_um"] = pd.cut(
        design["distance_um"],
        bins=DISTANCE_BINS_UM,
        labels=DISTANCE_LABELS,
        right=False,
        include_lowest=True,
    )
    return design.reset_index(names="cell_id")


def summarize_distance_bins(design: pd.DataFrame, metric_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_bin = (
        design.dropna(subset=["distance_bin_um"])
        .groupby(["sample", "age", "distance_bin_um"], observed=True)
        .agg(
            n_cells=("cell_id", "size"),
            mean_value=(metric_col, "mean"),
            median_distance_um=("distance_um", "median"),
        )
        .reset_index()
    )
    age_bin = (
        sample_bin.groupby(["age", "distance_bin_um"], observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            total_cells=("n_cells", "sum"),
            mean_value=("mean_value", "mean"),
            sem_value=("mean_value", sem),
        )
        .reset_index()
    )
    return sample_bin, age_bin


def plot_binned_age_trend(age_bin: pd.DataFrame, ylabel: str, title: str, outbase: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.1))
    x = np.arange(len(DISTANCE_LABELS))
    for age in AGE_ORDER:
        sub = age_bin[age_bin["age"].astype(str).eq(age)].set_index("distance_bin_um").reindex(DISTANCE_LABELS)
        ax.errorbar(
            x,
            sub["mean_value"],
            yerr=sub["sem_value"],
            marker="o",
            linewidth=1.8,
            capsize=3,
            color=AGE_COLORS[age],
            label=age,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(DISTANCE_LABELS, rotation=35, ha="right")
    ax.set_xlabel("Distance bin (um)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def write_squidpy_distance_outputs(df: pd.DataFrame, score_label: str) -> dict[str, int]:
    adata = make_distance_anndata(df)
    run_var_by_distance(adata)

    snc_design = extract_distance_design(adata, "distance_to_snc_plus", "SnC+")
    sssp_design = extract_distance_design(adata, "distance_to_sssp_high", "SSSP-high")
    snc_design.to_csv(OUTDIR / "squidpy_distance_to_snc_plus_design_matrix.csv.gz", index=False)
    sssp_design.to_csv(OUTDIR / "squidpy_distance_to_sssp_high_design_matrix.csv.gz", index=False)

    distance_outputs = [
        (
            snc_design,
            "SSSP_score",
            "sssp_score_by_distance_to_snc_plus",
            score_label,
            f"{score_label} by distance to SnC+ cells",
        ),
        (
            snc_design,
            "SSSP_high",
            "sssp_high_fraction_by_distance_to_snc_plus",
            "Fraction SSSP-high",
            "SSSP-high fraction by distance to SnC+ cells",
        ),
        (
            sssp_design,
            "SnC_plus",
            "snc_fraction_by_distance_to_sssp_high",
            "Fraction SnC+",
            "SnC+ fraction by distance to SSSP-high cells",
        ),
    ]
    for design, metric_col, suffix, ylabel, title in distance_outputs:
        sample_bin, age_bin = summarize_distance_bins(design, metric_col)
        sample_bin.to_csv(OUTDIR / f"{suffix}_sample_bins.csv", index=False)
        age_bin.to_csv(OUTDIR / f"{suffix}_age_bins.csv", index=False)
        plot_binned_age_trend(age_bin, ylabel, title, OUTDIR / f"{suffix}_young_vs_old")

    return {
        "squidpy_distance_to_snc_plus_rows": int(len(snc_design)),
        "squidpy_distance_to_sssp_high_rows": int(len(sssp_design)),
    }


def run_analysis(config: dict[str, object]) -> pd.DataFrame:
    global SSSP_CELL_SCORES, OUTDIR, SCORE_COL, SSSP_HIGH_COL
    SSSP_CELL_SCORES = Path(config["cell_scores"])
    OUTDIR = Path(config["outdir"])
    SCORE_COL = str(config["score_col"])
    SSSP_HIGH_COL = str(config["high_col"])

    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = load_data()
    df = add_proximity(df)
    cutoff = df.attrs.get("sssp_high_cutoff", np.nan)
    df.to_csv(OUTDIR / "xenium_snc_sssp_top_quartile_score_proximity_cell_table.csv.gz", index=False)

    by_sample, deltas, age_summary = summarize_score_by_snc(df)
    snc_to_sssp, sssp_to_snc = summarize_proximity(df)
    tests = sample_level_tests(by_sample, snc_to_sssp)
    permutation = permutation_test_sssp_to_snc(df)
    distance_summary = write_squidpy_distance_outputs(df, str(config["score_label"]))

    run_summary = pd.DataFrame(
        [
            {
                "gene_set_id": config["gene_set_id"],
                "n_cells": int(len(df)),
                "n_snc_plus": int(df["snc"].astype(str).eq("SnC+").sum()),
                "n_sssp_high": int(df[SSSP_HIGH_COL].sum()),
                "sssp_high_definition": f"{SCORE_COL} >= top quartile cutoff",
                "sssp_high_quantile": SSSP_HIGH_QUANTILE,
                "sssp_high_score_cutoff": cutoff,
                "distance_scope": "within sample and FOV",
                "n_permutations": N_PERMUTATIONS,
                "permutation_scope": "SSSP-high labels shuffled within sample/FOV",
                **distance_summary,
            }
        ]
    )
    run_summary.to_csv(OUTDIR / "snc_vs_sssp_top_quartile_run_summary.csv", index=False)

    print("SnC vs SSSP analysis complete.")
    print(run_summary.to_string(index=False))
    print("\nSnC+ minus SnC- score deltas by sample:")
    print(deltas.to_string(index=False))
    print("\nSample-level tests:")
    print(tests.to_string(index=False))
    print("\nSSSP-high to SnC+ permutation tests:")
    print(permutation.to_string(index=False))
    print(f"\nWrote outputs to: {OUTDIR}")
    return run_summary


def main() -> None:
    summaries = [run_analysis(config) for config in GENE_SET_CONFIGS]
    pd.concat(summaries, ignore_index=True).to_csv(
        PROJECT_ROOT / "results" / "260727" / "xenium_snc_vs_sssp_score_proximity_gene_set_run_summary.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
