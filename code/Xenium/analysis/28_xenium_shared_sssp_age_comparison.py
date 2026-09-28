"""Compare the shared SSSP Xenium overlap between individual young and old samples."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats


PROJECT_ROOT = Path(__file__).resolve().parents[3]
INPUT_H5AD = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_sssp_20gene_scanpy_score_genes"
    / "xenium_20gene_sssp_scanpy_scored.h5ad"
)
CONFIG_PATH = PROJECT_ROOT / "code" / "shared" / "sssp.json"
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_shared_sssp_age_comparison"
SCORE_NAME = "sssp_xenium_control_adjusted_score"
RANDOM_SEED = 260727


def read_h5ad_compat(path: Path):
    try:
        return sc.read_h5ad(path)
    except Exception:
        OUTDIR.mkdir(parents=True, exist_ok=True)
        temporary = OUTDIR / "_temporary_input_without_null_log1p_base.h5ad"
        shutil.copy2(path, temporary)
        import h5py

        with h5py.File(temporary, "r+") as handle:
            if "uns" in handle and "log1p" in handle["uns"] and "base" in handle["uns/log1p"]:
                del handle["uns/log1p/base"]
        try:
            return sc.read_h5ad(temporary)
        finally:
            temporary.unlink(missing_ok=True)


def age_from_sample(sample: str) -> str:
    sample = str(sample).upper()
    if sample.endswith("-Y"):
        return "Young"
    if sample.endswith("-O"):
        return "Old"
    return "Unknown"


def sem(values: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce").dropna()
    return float(values.std(ddof=1) / np.sqrt(len(values))) if len(values) > 1 else np.nan


def load_genes() -> list[str]:
    config = json.loads(CONFIG_PATH.read_text())
    genes = [str(g).upper() for g in config["sssp_xenium"]]
    if len(genes) != len(set(genes)):
        raise ValueError("sssp_xenium contains duplicate genes")
    return genes


def safe_test(old: pd.Series, young: pd.Series) -> tuple[float, float]:
    old = pd.to_numeric(old, errors="coerce").dropna()
    young = pd.to_numeric(young, errors="coerce").dropna()
    if len(old) < 2 or len(young) < 2:
        return np.nan, np.nan
    welch = float(stats.ttest_ind(old, young, equal_var=False).pvalue)
    mann_whitney = float(stats.mannwhitneyu(old, young, alternative="two-sided").pvalue)
    return welch, mann_whitney


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    genes = load_genes()
    adata = read_h5ad_compat(INPUT_H5AD)
    adata.var_names = adata.var_names.astype(str).str.upper()
    adata.var_names_make_unique()
    present = [gene for gene in genes if gene in adata.var_names]
    missing = sorted(set(genes) - set(present))
    if not present:
        raise RuntimeError("None of the configured sssp_xenium genes are present in Xenium.")

    expression = adata[:, present].X
    if sparse.issparse(expression):
        expression = expression.tocsr()
        mean_expression = np.asarray(expression.mean(axis=1)).ravel()
        detected = np.asarray((expression > 0).sum(axis=1)).ravel()
    else:
        expression = np.asarray(expression)
        mean_expression = expression.mean(axis=1)
        detected = (expression > 0).sum(axis=1)

    sc.tl.score_genes(
        adata,
        gene_list=present,
        score_name=SCORE_NAME,
        ctrl_size=min(50, max(1, adata.n_vars - len(present))),
        n_bins=25,
        random_state=RANDOM_SEED,
        use_raw=False,
    )

    cell_scores = pd.DataFrame(
        {
            "cell_id": adata.obs["cell_id"].astype(str).to_numpy(),
            "sample": adata.obs["sample"].astype(str).to_numpy(),
            "age": adata.obs["age"].astype(str).to_numpy(),
            "mean_log1p_panel": mean_expression,
            "fraction_detected": detected / len(present),
            "n_detected": detected.astype(int),
            SCORE_NAME: pd.to_numeric(adata.obs[SCORE_NAME], errors="coerce").to_numpy(),
        }
    )
    cell_scores["age"] = cell_scores["age"].where(
        cell_scores["age"].isin(["Young", "Old"]),
        cell_scores["sample"].map(age_from_sample),
    )
    metrics = ["mean_log1p_panel", "fraction_detected", "n_detected", SCORE_NAME]
    sample_summary = (
        cell_scores.groupby(["sample", "age"], observed=True)[metrics]
        .agg(["mean", "median", sem, "count"])
    )
    sample_summary.columns = [f"{metric}_{statistic}" for metric, statistic in sample_summary.columns]
    sample_summary = sample_summary.reset_index().sort_values(["age", "sample"])

    test_rows = []
    for metric in metrics:
        old = sample_summary.loc[sample_summary["age"].eq("Old"), f"{metric}_mean"]
        young = sample_summary.loc[sample_summary["age"].eq("Young"), f"{metric}_mean"]
        welch, mann_whitney = safe_test(old, young)
        test_rows.append(
            {
                "metric": metric,
                "n_old_samples": int(old.notna().sum()),
                "n_young_samples": int(young.notna().sum()),
                "old_mean_across_samples": float(old.mean()),
                "young_mean_across_samples": float(young.mean()),
                "old_minus_young": float(old.mean() - young.mean()),
                "n_old_samples_above_young_mean": int((old > young.mean()).sum()),
                "n_young_samples_below_old_mean": int((young < old.mean()).sum()),
                "welch_p": welch,
                "mann_whitney_p": mann_whitney,
            }
        )
    tests = pd.DataFrame(test_rows)

    cell_scores.to_csv(OUTDIR / "cell_scores.csv.gz", index=False)
    sample_summary.to_csv(OUTDIR / "sample_summary.csv", index=False)
    tests.to_csv(OUTDIR / "old_young_sample_tests.csv", index=False)
    pd.DataFrame({"gene": genes, "present_in_xenium": [gene in present for gene in genes]}).to_csv(
        OUTDIR / "gene_panel_used.csv", index=False
    )
    (OUTDIR / "metadata.json").write_text(
        json.dumps(
            {
                "input_h5ad": str(INPUT_H5AD),
                "config_path": str(CONFIG_PATH),
                "n_configured_genes": len(genes),
                "n_present_genes": len(present),
                "missing_genes": missing,
                "methods": metrics,
                "n_cells": len(cell_scores),
            },
            indent=2,
        )
    )
    print(f"Configured genes: {len(genes)}")
    print(f"Present genes scored: {len(present)}")
    print(f"Samples: {sample_summary['sample'].nunique()} ({(sample_summary['age'] == 'Young').sum()} young, {(sample_summary['age'] == 'Old').sum()} old)")
    print(tests.to_string(index=False))
    print(f"Results written to {OUTDIR}")


if __name__ == "__main__":
    main()