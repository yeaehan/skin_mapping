"""
Assess whether GTEx and Tabula Sapiens skin data can validate the SSSP.

This script intentionally does not integrate these data with Sinai. Instead it
summarizes metadata coverage, SSSP gene overlap, and within-dataset SSSP scores
by donor/age/cell type. That keeps the external datasets as independent
validation/triangulation cohorts rather than adding more batch/site confounding
to the discovery analysis.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats


PROJECT_ROOT = Path(__file__).resolve().parents[2]
GTEX_H5AD = PROJECT_ROOT / "data" / "external" / "GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"
TABULA_H5AD = PROJECT_ROOT / "data" / "external" / "Skin_TSP1_30_version2d_10X_smartseq_scvi_Nov122024_updated.h5ad"
GENE_CONFIG = PROJECT_ROOT / "code" / "Xenium" / "analysis" / "xenium_sssp_gene_set_config.csv"
OUTDIR = PROJECT_ROOT / "results" / "260727" / "external_gtex_tabula_sssp_validation"

GENE_SETS = {
    "SSSP_164": "in_skin_public_secretory_164",
    "SSSP_20": "in_xenium_main_20",
    "SSSP_17_diff_filtered": "in_xenium_diff_filtered",
}
CONVENTIONAL_SASP_GENES = [
    "CCL11", "CCL13", "CCL14", "CCL15", "CCL16", "CCL19", "CCL2", "CCL20", "CCL21",
    "CCL26", "CCL28", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8",
    "CXCL1", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL14", "CXCL16", "CXCL2",
    "CXCL3", "CXCL5", "CXCL6", "CXCL9",
    "IL12B", "IL13", "IL15", "IL18", "IL1A", "IL1B", "IL2", "IL23A", "IL4", "IL5", "IL6",
    "MMP12", "MMP9", "MMP3",
    "VEGFA", "TNF", "TGFB1", "CSF1", "EGF", "FGF2", "IFNB1", "CMKLR1",
]
PRIMARY_SCORE_SET = "SSSP_164"
GTEX_PREP_FILTER = None
SENESCENCE_MARKERS = {
    "CDKN1A_p21": "CDKN1A",
    "CDKN2A_p16": "CDKN2A",
}
PROLIFERATION_MARKERS = {
    "PCNA": "PCNA",
    "MKI67": "MKI67",
}
MIN_CELLS_PER_GROUP = 20
MIN_CELLS_PER_QUARTILE_GROUP = 50


def age_mid_from_gtex_bin(value: str) -> float:
    value = str(value)
    if "-" in value:
        lo, hi = value.split("-", 1)
        return (float(lo) + float(hi)) / 2.0
    return np.nan


def make_scanpy_score_input(adata: ad.AnnData, layer: str | None, symbol_col: str | None) -> ad.AnnData:
    matrix = adata.layers[layer] if layer else adata.X
    matrix = sparse.csr_matrix(matrix, dtype=np.float32) if sparse.issparse(matrix) else np.asarray(matrix, dtype=np.float32)
    work = ad.AnnData(X=matrix.copy(), obs=adata.obs.copy(), var=adata.var.copy())
    symbols = get_var_symbols(work, symbol_col).astype(str).str.upper()
    fallback = pd.Series(work.var_names.astype(str), index=work.var_names).str.upper()
    symbols = symbols.where(~symbols.isin(["", "NAN", "NONE"]), fallback)
    work.var_names = pd.Index(symbols)
    work.var_names_make_unique()
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    return work


def load_gene_sets() -> dict[str, list[str]]:
    cfg = pd.read_csv(GENE_CONFIG)
    cfg["gene"] = cfg["gene"].astype(str).str.upper()
    out = {}
    for name, col in GENE_SETS.items():
        out[name] = sorted(cfg.loc[cfg[col].astype(bool), "gene"].unique())
    return out


def get_var_symbols(adata: ad.AnnData, preferred_col: str | None = None) -> pd.Series:
    if preferred_col and preferred_col in adata.var.columns:
        symbols = adata.var[preferred_col].astype(str)
    elif "gene_symbol" in adata.var.columns:
        symbols = adata.var["gene_symbol"].astype(str)
    elif "Approved symbol" in adata.var.columns:
        symbols = adata.var["Approved symbol"].astype(str)
    elif "gene_name" in adata.var.columns:
        symbols = adata.var["gene_name"].astype(str)
    else:
        symbols = pd.Series(adata.var_names.astype(str), index=adata.var_names)
    symbols = symbols.str.upper()
    symbols.index = adata.var_names
    return symbols


def summarize_metadata(obs: pd.DataFrame, dataset: str) -> pd.DataFrame:
    rows = []
    cols = [
        "dataset",
        "n_cells",
        "n_donors",
        "age_values",
        "site_values",
        "method_values",
        "major_cell_type_values",
    ]
    donor_col = "donor" if "donor" in obs.columns else "Participant ID"
    age_col = "age" if "age" in obs.columns else "Age_bin"
    site_col = "anatomical_position" if "anatomical_position" in obs.columns else "Tissue Site Detail"
    method_col = "method" if "method" in obs.columns else "prep"
    cell_col = "broad_cell_class" if "broad_cell_class" in obs.columns else "Broad cell type"
    rows.append(
        {
            "dataset": dataset,
            "n_cells": int(len(obs)),
            "n_donors": int(obs[donor_col].astype(str).nunique()),
            "age_values": "; ".join(map(str, sorted(obs[age_col].astype(str).unique()))),
            "site_values": "; ".join(map(str, sorted(obs[site_col].astype(str).unique()))),
            "method_values": "; ".join(map(str, sorted(obs[method_col].astype(str).unique()))),
            "major_cell_type_values": "; ".join(map(str, obs[cell_col].astype(str).value_counts().head(12).index)),
        }
    )
    return pd.DataFrame(rows, columns=cols)


def summarize_gene_overlap(dataset: str, var_symbols: pd.Series, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    present = set(var_symbols.dropna().astype(str).str.upper())
    rows = []
    for set_name, genes in gene_sets.items():
        overlap = sorted(set(genes) & present)
        rows.append(
            {
                "dataset": dataset,
                "gene_set": set_name,
                "n_genes_in_set": len(genes),
                "n_present": len(overlap),
                "fraction_present": len(overlap) / len(genes) if genes else np.nan,
                "present_genes": ";".join(overlap),
                "missing_genes": ";".join(sorted(set(genes) - present)),
            }
        )
    return pd.DataFrame(rows)


def add_scores(adata: ad.AnnData, dataset: str, gene_sets: dict[str, list[str]], layer: str | None, symbol_col: str | None) -> pd.DataFrame:
    score_input = make_scanpy_score_input(adata, layer, symbol_col)
    scores = adata.obs.copy()
    scores["dataset"] = dataset
    for marker_name, gene in {**SENESCENCE_MARKERS, **PROLIFERATION_MARKERS}.items():
        expr_col = f"{marker_name}_log1p_norm"
        if gene in score_input.var_names:
            values = score_input[:, gene].X
            scores[expr_col] = np.asarray(values.toarray() if sparse.issparse(values) else values).ravel()
        else:
            scores[expr_col] = np.nan
    sasp_genes_present = [gene for gene in CONVENTIONAL_SASP_GENES if gene in score_input.var_names]
    if sasp_genes_present:
        sc.tl.score_genes(
            score_input,
            gene_list=sasp_genes_present,
            score_name="conventional_sasp_scanpy_score",
            random_state=260727,
            use_raw=False,
        )
        scores["conventional_sasp_scanpy_score"] = score_input.obs["conventional_sasp_scanpy_score"].to_numpy()
        scores["conventional_sasp_n_present"] = len(sasp_genes_present)
    else:
        scores["conventional_sasp_scanpy_score"] = np.nan
        scores["conventional_sasp_n_present"] = 0
    for set_name, genes in gene_sets.items():
        present_genes = [g for g in genes if g in score_input.var_names]
        score_col = f"{set_name}_scanpy_score"
        if present_genes:
            sc.tl.score_genes(
                score_input,
                gene_list=present_genes,
                score_name=score_col,
                random_state=260727,
                use_raw=False,
            )
            scores[score_col] = score_input.obs[score_col].to_numpy()
            scores[f"{set_name}_n_present"] = len(present_genes)
        else:
            scores[score_col] = np.nan
            scores[f"{set_name}_n_present"] = 0
    return scores


def harmonize_scores(scores: pd.DataFrame, dataset: str) -> pd.DataFrame:
    out = scores.copy()
    if dataset == "GTEx_skin":
        out["donor"] = out["Participant ID"].astype(str)
        out["age"] = out["Age_bin"].astype(str)
        out["age_numeric"] = out["Age_bin"].map(age_mid_from_gtex_bin)
        out["site"] = out["Tissue Site Detail"].astype(str)
        out["method"] = out["prep"].astype(str)
        out["sex"] = out["Sex"].astype(str)
        out["major_cell_type"] = out["Broad cell type"].astype(str)
        out["minor_cell_type"] = out["Granular cell type"].astype(str)
    else:
        out["donor"] = out["donor"].astype(str)
        out["age"] = out["age"].astype(str)
        out["age_numeric"] = pd.to_numeric(out["age"], errors="coerce")
        out["site"] = out["anatomical_position"].astype(str)
        out["method"] = out["method"].astype(str)
        out["sex"] = out["sex"].astype(str)
        out["major_cell_type"] = out["broad_cell_class"].astype(str)
        out["minor_cell_type"] = out["cell_ontology_class"].astype(str)

    keep = [
        "dataset",
        "donor",
        "age",
        "age_numeric",
        "site",
        "method",
        "sex",
        "major_cell_type",
        "minor_cell_type",
    ] + [f"{name}_scanpy_score" for name in GENE_SETS]
    keep += [f"{name}_log1p_norm" for name in {**SENESCENCE_MARKERS, **PROLIFERATION_MARKERS}]
    keep += ["conventional_sasp_scanpy_score", "conventional_sasp_n_present"]
    return out[keep].copy()


def summarize_scores(cell_scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    score_cols = [f"{name}_scanpy_score" for name in GENE_SETS]
    marker_cols = [f"{name}_log1p_norm" for name in {**SENESCENCE_MARKERS, **PROLIFERATION_MARKERS}]
    mean_cols = score_cols + marker_cols + ["conventional_sasp_scanpy_score", "conventional_sasp_n_present"]
    donor_summary = (
        cell_scores.groupby(["dataset", "donor", "age", "age_numeric", "site", "method", "sex"], dropna=False)[mean_cols]
        .mean()
        .reset_index()
    )
    celltype_summary = (
        cell_scores.groupby(["dataset", "donor", "age", "age_numeric", "major_cell_type"], dropna=False)
        .agg(n_cells=("dataset", "size"), **{col: (col, "mean") for col in mean_cols})
        .reset_index()
    )
    celltype_summary = celltype_summary.loc[celltype_summary["n_cells"] >= MIN_CELLS_PER_GROUP].copy()

    trend_rows = []
    for dataset, ddf in donor_summary.groupby("dataset"):
        for col in score_cols:
            valid = ddf[["age_numeric", col]].dropna()
            if len(valid) >= 3 and valid["age_numeric"].nunique() >= 3:
                rho, pval = stats.spearmanr(valid["age_numeric"], valid[col])
            else:
                rho, pval = np.nan, np.nan
            trend_rows.append(
                {
                    "dataset": dataset,
                    "score": col,
                    "n_donor_or_age_points": int(len(valid)),
                    "spearman_rho_age": rho,
                    "spearman_p_age": pval,
                    "note": "descriptive only; donor/site/method confounded and low N",
                }
            )
    return donor_summary, celltype_summary, pd.DataFrame(trend_rows)


def zscore_within_groups(values: pd.Series, groups: pd.Series, min_cells: int = MIN_CELLS_PER_QUARTILE_GROUP) -> pd.Series:
    out = pd.Series(np.nan, index=values.index, dtype=float)
    for _, idx in groups.groupby(groups, dropna=False).groups.items():
        idx = list(idx)
        if len(idx) < min_cells:
            continue
        group_values = values.loc[idx].astype(float)
        sd = group_values.std(ddof=0)
        if not np.isfinite(sd) or sd == 0:
            out.loc[idx] = 0.0
        else:
            out.loc[idx] = (group_values - group_values.mean()) / sd
    return out


def top_quartile_within_groups(values: pd.Series, groups: pd.Series, min_cells: int = MIN_CELLS_PER_QUARTILE_GROUP) -> pd.Series:
    out = pd.Series(False, index=values.index, dtype=bool)
    for _, idx in groups.groupby(groups, dropna=False).groups.items():
        idx = list(idx)
        group_values = values.loc[idx].astype(float).dropna()
        if len(group_values) < min_cells or group_values.nunique() < 2:
            continue
        percentile_rank = group_values.rank(method="average", pct=True)
        out.loc[group_values.index] = percentile_rank.gt(0.75)
    return out


def add_snc_and_sssp_flags(cell_scores: pd.DataFrame) -> pd.DataFrame:
    out = cell_scores.copy()
    score_col = f"{PRIMARY_SCORE_SET}_scanpy_score"
    cdkn_cols = ["CDKN1A_p21_log1p_norm", "CDKN2A_p16_log1p_norm"]
    proliferation_cols = ["PCNA_log1p_norm", "MKI67_log1p_norm"]
    required = cdkn_cols + proliferation_cols + [score_col]
    missing = [col for col in required if col not in out.columns]
    if missing:
        raise ValueError(f"Missing columns for SnC/SSSP calls: {missing}")

    out["cdkn_any_positive"] = out[cdkn_cols].fillna(0).gt(0).any(axis=1)
    out["proliferation_any_positive"] = out[proliferation_cols].fillna(0).gt(0).any(axis=1)
    out["snc_original_positive"] = out["cdkn_any_positive"] & ~out["proliferation_any_positive"]

    out["cdkn_mean_log1p_norm"] = out[cdkn_cols].mean(axis=1)
    out["proliferation_mean_log1p_norm"] = out[proliferation_cols].mean(axis=1)
    groups = out["dataset"].astype(str) + " | " + out["major_cell_type"].astype(str)
    out["cdkn_z_within_dataset_celltype"] = zscore_within_groups(out["cdkn_mean_log1p_norm"], groups)
    out["proliferation_z_within_dataset_celltype"] = zscore_within_groups(out["proliferation_mean_log1p_norm"], groups)
    out["sasp_z_within_dataset_celltype"] = zscore_within_groups(out["conventional_sasp_scanpy_score"], groups)
    out["snc_cdkn_cycle_sasp_index"] = (
        out["cdkn_z_within_dataset_celltype"]
        + out["sasp_z_within_dataset_celltype"]
        - out["proliferation_z_within_dataset_celltype"]
    )
    out["snc_cdkn_cycle_sasp_q4"] = top_quartile_within_groups(out["snc_cdkn_cycle_sasp_index"], groups)
    out["sssp_high_q4"] = top_quartile_within_groups(out[score_col], groups)

    out["snc_original_label"] = np.where(out["snc_original_positive"], "SnC+", "SnC-")
    out["snc_q4_label"] = np.where(out["snc_cdkn_cycle_sasp_q4"], "SnC-Q4+", "SnC-Q4-")
    out["sssp_high_label"] = np.where(out["sssp_high_q4"], "SSSP-high", "SSSP-low")
    return out


def summarize_snc_sssp_relations(cell_scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    snc_defs = {
        "original_positive": "snc_original_positive",
        "cdkn_cycle_sasp_q4": "snc_cdkn_cycle_sasp_q4",
    }
    sample_cols = ["dataset", "donor", "age", "age_numeric", "site", "method", "sex"]
    sample_rows = []
    overlap_rows = []
    celltype_rows = []
    for def_name, call_col in snc_defs.items():
        for group_values, sdf in cell_scores.groupby(sample_cols, dropna=False):
            n_cells = len(sdf)
            n_snc = int(sdf[call_col].sum())
            n_sssp = int(sdf["sssp_high_q4"].sum())
            n_overlap = int((sdf[call_col] & sdf["sssp_high_q4"]).sum())
            frac_sssp_in_snc = n_overlap / n_snc if n_snc else np.nan
            frac_sssp_in_non_snc = int((~sdf[call_col] & sdf["sssp_high_q4"]).sum()) / int((~sdf[call_col]).sum()) if int((~sdf[call_col]).sum()) else np.nan
            sample_rows.append(
                {
                    **dict(zip(sample_cols, group_values)),
                    "snc_definition": def_name,
                    "n_cells": n_cells,
                    "n_snc": n_snc,
                    "snc_fraction": n_snc / n_cells if n_cells else np.nan,
                    "n_sssp_high": n_sssp,
                    "sssp_high_fraction": n_sssp / n_cells if n_cells else np.nan,
                    "n_snc_and_sssp_high": n_overlap,
                    "fraction_snc_that_are_sssp_high": frac_sssp_in_snc,
                    "fraction_non_snc_that_are_sssp_high": frac_sssp_in_non_snc,
                    "sssp_high_fraction_snc_minus_non_snc": frac_sssp_in_snc - frac_sssp_in_non_snc if pd.notna(frac_sssp_in_snc) and pd.notna(frac_sssp_in_non_snc) else np.nan,
                }
            )

        for dataset, ddf in cell_scores.groupby("dataset", dropna=False):
            table = pd.crosstab(ddf[call_col], ddf["sssp_high_q4"])
            for value in [False, True]:
                if value not in table.index:
                    table.loc[value, :] = 0
                if value not in table.columns:
                    table.loc[:, value] = 0
            table = table.sort_index().sort_index(axis=1)
            oddsratio, pval = stats.fisher_exact(
                [[table.loc[True, True], table.loc[True, False]], [table.loc[False, True], table.loc[False, False]]]
            )
            overlap_rows.append(
                {
                    "dataset": dataset,
                    "snc_definition": def_name,
                    "n_cells": int(table.to_numpy().sum()),
                    "n_snc": int(table.loc[True].sum()),
                    "n_sssp_high": int(table[True].sum()),
                    "n_snc_and_sssp_high": int(table.loc[True, True]),
                    "fisher_odds_ratio_descriptive": oddsratio,
                    "fisher_p_descriptive": pval,
                    "note": "Cell-level enrichment is descriptive because cells are not independent biological replicates.",
                }
            )

        ct = (
            cell_scores.groupby(["dataset", "major_cell_type"], dropna=False)
            .agg(
                n_cells=(call_col, "size"),
                n_snc=(call_col, "sum"),
                n_sssp_high=("sssp_high_q4", "sum"),
                n_snc_and_sssp_high=(call_col, lambda x, df=cell_scores, col=call_col: int((x.astype(bool) & df.loc[x.index, "sssp_high_q4"].astype(bool)).sum())),
            )
            .reset_index()
        )
        ct["snc_definition"] = def_name
        ct["snc_fraction"] = ct["n_snc"] / ct["n_cells"].replace(0, np.nan)
        ct["sssp_high_fraction"] = ct["n_sssp_high"] / ct["n_cells"].replace(0, np.nan)
        ct["snc_sssp_overlap_fraction"] = ct["n_snc_and_sssp_high"] / ct["n_cells"].replace(0, np.nan)
        celltype_rows.append(ct)

    return pd.DataFrame(sample_rows), pd.DataFrame(overlap_rows), pd.concat(celltype_rows, ignore_index=True)


def calculate_marker_score_correlations(cell_scores: pd.DataFrame, donor_summary: pd.DataFrame) -> pd.DataFrame:
    score_col = f"{PRIMARY_SCORE_SET}_scanpy_score"
    rows = []
    for level_name, table in [("cell", cell_scores), ("donor_site_method", donor_summary)]:
        for dataset, ddf in table.groupby("dataset"):
            for marker_name in SENESCENCE_MARKERS:
                marker_col = f"{marker_name}_log1p_norm"
                valid = ddf[[score_col, marker_col]].replace([np.inf, -np.inf], np.nan).dropna()
                if len(valid) >= 3 and valid[score_col].nunique() >= 2 and valid[marker_col].nunique() >= 2:
                    rho, pval = stats.spearmanr(valid[marker_col], valid[score_col])
                else:
                    rho, pval = np.nan, np.nan
                rows.append(
                    {
                        "level": level_name,
                        "dataset": dataset,
                        "marker": marker_name,
                        "score": score_col,
                        "n_points": int(len(valid)),
                        "spearman_rho": rho,
                        "spearman_p": pval,
                        "note": "cell-level is descriptive because cells are not biological replicates; donor/site/method level has low N",
                    }
                )
    return pd.DataFrame(rows)


def visual_group_col(dataset: str) -> str:
    return "method" if dataset == "GTEx_skin" else "site"


def visual_group_title(dataset: str) -> str:
    return "Prep method" if dataset == "GTEx_skin" else "Site"


def display_label(row: pd.Series, dataset: str) -> str:
    if dataset == "GTEx_skin":
        return f"{row['donor']} | {row['site']} | {row['method']}"
    return str(row["donor"])


def plot_donor_scores(donor_summary: pd.DataFrame) -> None:
    score_col = f"{PRIMARY_SCORE_SET}_scanpy_score"
    if score_col not in donor_summary.columns:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, (dataset, ddf) in zip(axes, donor_summary.groupby("dataset")):
        ddf = ddf.sort_values("age_numeric")
        group_col = visual_group_col(dataset)
        group_values = sorted(ddf[group_col].astype(str).unique())
        colors = dict(zip(group_values, plt.cm.Set2(np.linspace(0, 1, max(len(group_values), 3)))))
        for group_value, sdf in ddf.groupby(group_col, dropna=False):
            ax.scatter(
                sdf["age_numeric"],
                sdf[score_col],
                s=80,
                color=colors[str(group_value)],
                edgecolor="black",
                linewidth=0.4,
                label=str(group_value),
            )
        for _, row in ddf.iterrows():
            ax.text(row["age_numeric"], row[score_col], display_label(row, dataset), fontsize=8, ha="left", va="bottom")
        ax.set_title(dataset)
        ax.set_xlabel("Age / age-bin midpoint")
        ax.grid(alpha=0.25)
        ax.legend(title=visual_group_title(dataset), fontsize=7, title_fontsize=8, frameon=False, loc="best")
    axes[0].set_ylabel("Scanpy SSSP score\n(full detected 164-gene panel)")
    fig.suptitle("External datasets: donor-level SSSP score is descriptive only")
    fig.tight_layout()
    fig.savefig(OUTDIR / "external_gtex_tabula_sssp164_detected_scanpy_score_by_method_or_site.png", dpi=300)
    fig.savefig(OUTDIR / "external_gtex_tabula_sssp164_detected_scanpy_score_by_method_or_site.pdf")
    plt.close(fig)


def plot_marker_score_correlations(donor_summary: pd.DataFrame) -> None:
    score_col = f"{PRIMARY_SCORE_SET}_scanpy_score"
    marker_cols = [f"{name}_log1p_norm" for name in SENESCENCE_MARKERS]
    available_marker_cols = [col for col in marker_cols if col in donor_summary.columns]
    if score_col not in donor_summary.columns or not available_marker_cols:
        return

    datasets = list(donor_summary["dataset"].drop_duplicates())
    fig, axes = plt.subplots(
        len(available_marker_cols),
        len(datasets),
        figsize=(5 * len(datasets), 4 * len(available_marker_cols)),
        squeeze=False,
        sharey=True,
    )
    for col_idx, dataset in enumerate(datasets):
        ddf = donor_summary.loc[donor_summary["dataset"].eq(dataset)].copy()
        group_col = visual_group_col(dataset)
        group_values = sorted(ddf[group_col].astype(str).unique())
        colors = dict(zip(group_values, plt.cm.Set2(np.linspace(0, 1, max(len(group_values), 3)))))
        for row_idx, marker_col in enumerate(available_marker_cols):
            ax = axes[row_idx][col_idx]
            for group_value, sdf in ddf.groupby(group_col, dropna=False):
                ax.scatter(
                    sdf[marker_col],
                    sdf[score_col],
                    s=80,
                    color=colors[str(group_value)],
                    edgecolor="black",
                    linewidth=0.4,
                    label=str(group_value),
                )
            for _, sample_row in ddf.iterrows():
                ax.text(sample_row[marker_col], sample_row[score_col], display_label(sample_row, dataset), fontsize=8, ha="left", va="bottom")
            valid = ddf[[marker_col, score_col]].replace([np.inf, -np.inf], np.nan).dropna()
            if len(valid) >= 3 and valid[marker_col].nunique() >= 2 and valid[score_col].nunique() >= 2:
                rho, pval = stats.spearmanr(valid[marker_col], valid[score_col])
                corr_text = f"Spearman rho={rho:.2f}, p={pval:.3g}"
            else:
                corr_text = "Spearman rho=NA"
            marker_label = marker_col.replace("_log1p_norm", "").replace("_", " / ")
            ax.set_title(f"{dataset}\n{corr_text}")
            ax.set_xlabel(f"{marker_label} log1p normalized expression")
            ax.grid(alpha=0.25)
            if col_idx == 0:
                ax.set_ylabel("Scanpy SSSP score\n(full detected 164-gene panel)")
            if row_idx == 0:
                ax.legend(title=visual_group_title(dataset), fontsize=7, title_fontsize=8, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(OUTDIR / "external_cdkn_marker_vs_sssp164_scanpy_score_by_method_or_site.png", dpi=300)
    fig.savefig(OUTDIR / "external_cdkn_marker_vs_sssp164_scanpy_score_by_method_or_site.pdf")
    plt.close(fig)


def plot_marker_score_correlation_lines(donor_summary: pd.DataFrame) -> None:
    score_col = f"{PRIMARY_SCORE_SET}_scanpy_score"
    marker_cols = [f"{name}_log1p_norm" for name in SENESCENCE_MARKERS]
    available_marker_cols = [col for col in marker_cols if col in donor_summary.columns]
    if score_col not in donor_summary.columns or not available_marker_cols:
        return

    datasets = list(donor_summary["dataset"].drop_duplicates())
    fig, axes = plt.subplots(
        len(available_marker_cols),
        len(datasets),
        figsize=(5 * len(datasets), 4 * len(available_marker_cols)),
        squeeze=False,
        sharey=True,
    )
    for col_idx, dataset in enumerate(datasets):
        ddf = donor_summary.loc[donor_summary["dataset"].eq(dataset)].copy()
        group_col = visual_group_col(dataset)
        group_values = sorted(ddf[group_col].astype(str).unique())
        colors = dict(zip(group_values, plt.cm.Set2(np.linspace(0, 1, max(len(group_values), 3)))))
        for row_idx, marker_col in enumerate(available_marker_cols):
            ax = axes[row_idx][col_idx]
            keep_cols = list(dict.fromkeys([marker_col, score_col, group_col, "donor", "method", "site"]))
            valid = ddf[keep_cols].replace([np.inf, -np.inf], np.nan).dropna()
            if len(valid) >= 2 and valid[marker_col].nunique() >= 2:
                slope, intercept, _, _, _ = stats.linregress(valid[marker_col], valid[score_col])
                x_line = np.linspace(valid[marker_col].min(), valid[marker_col].max(), 100)
                ax.plot(x_line, intercept + slope * x_line, color="black", linewidth=2, alpha=0.8)
            for group_value, sdf in valid.groupby(group_col, dropna=False):
                ax.scatter(
                    sdf[marker_col],
                    sdf[score_col],
                    s=80,
                    color=colors[str(group_value)],
                    edgecolor="black",
                    linewidth=0.4,
                    label=str(group_value),
                    zorder=3,
                )
            if len(valid) >= 3 and valid[marker_col].nunique() >= 2 and valid[score_col].nunique() >= 2:
                rho, pval = stats.spearmanr(valid[marker_col], valid[score_col])
                corr_text = f"Spearman rho={rho:.2f}, p={pval:.3g}"
            else:
                corr_text = "Spearman rho=NA"
            marker_label = marker_col.replace("_log1p_norm", "").replace("_", " / ")
            ax.set_title(f"{dataset}\n{corr_text}")
            ax.set_xlabel(f"{marker_label} log1p normalized expression")
            ax.grid(alpha=0.25)
            if col_idx == 0:
                ax.set_ylabel("Scanpy SSSP score\n(full detected 164-gene panel)")
            if row_idx == 0:
                ax.legend(title=visual_group_title(dataset), fontsize=7, title_fontsize=8, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(OUTDIR / "external_cdkn_marker_vs_sssp164_scanpy_score_correlation_lines.png", dpi=300)
    fig.savefig(OUTDIR / "external_cdkn_marker_vs_sssp164_scanpy_score_correlation_lines.pdf")
    plt.close(fig)


def summarize_cell_binned_marker_trends(
    cell_scores: pd.DataFrame,
    marker_col: str,
    score_col: str,
    min_cells_per_sample: int = 50,
    min_cells_per_bin: int = 10,
    n_positive_bins: int = 4,
) -> pd.DataFrame:
    rows = []
    group_cols = ["dataset", "donor", "site", "method"]
    for group_values, sdf in cell_scores.groupby(group_cols, dropna=False):
        sdf = sdf[[marker_col, score_col]].replace([np.inf, -np.inf], np.nan).dropna().copy()
        if len(sdf) < min_cells_per_sample or sdf[marker_col].nunique() < 2:
            continue

        zero_mask = sdf[marker_col].le(0)
        bin_frames = []
        if zero_mask.sum() >= min_cells_per_bin:
            zero_df = sdf.loc[zero_mask].copy()
            zero_df["marker_bin"] = "zero"
            bin_frames.append(zero_df)

        positive_df = sdf.loc[~zero_mask].copy()
        if len(positive_df) >= min_cells_per_bin:
            n_bins = min(n_positive_bins, positive_df[marker_col].nunique(), len(positive_df) // min_cells_per_bin)
            if n_bins >= 2:
                positive_df["marker_bin"] = pd.qcut(
                    positive_df[marker_col],
                    q=n_bins,
                    labels=False,
                    duplicates="drop",
                ).astype(str)
            else:
                positive_df["marker_bin"] = "positive"
            bin_frames.append(positive_df)

        if not bin_frames:
            continue

        binned = pd.concat(bin_frames, ignore_index=True)
        summary = (
            binned.groupby("marker_bin", dropna=False)
            .agg(
                n_cells=(score_col, "size"),
                mean_marker=(marker_col, "mean"),
                mean_sssp=(score_col, "mean"),
            )
            .reset_index()
        )
        summary = summary.loc[summary["n_cells"].ge(min_cells_per_bin)].sort_values("mean_marker")
        if len(summary) < 2:
            continue

        dataset, donor, site, method = group_values
        sample_label = f"{donor} | {site} | {method}" if dataset == "GTEx_skin" else f"{donor} | {site}"
        for _, row in summary.iterrows():
            rows.append(
                {
                    "dataset": dataset,
                    "donor": donor,
                    "site": site,
                    "method": method,
                    "sample_label": sample_label,
                    "marker": marker_col.replace("_log1p_norm", ""),
                    "marker_bin": row["marker_bin"],
                    "n_cells": int(row["n_cells"]),
                    "mean_marker": row["mean_marker"],
                    "mean_sssp": row["mean_sssp"],
                }
            )
    return pd.DataFrame(rows)


def plot_cell_binned_marker_trend_lines(cell_scores: pd.DataFrame) -> None:
    score_col = f"{PRIMARY_SCORE_SET}_scanpy_score"
    marker_cols = [f"{name}_log1p_norm" for name in SENESCENCE_MARKERS]
    available_marker_cols = [col for col in marker_cols if col in cell_scores.columns]
    if score_col not in cell_scores.columns or not available_marker_cols:
        return

    trend_tables = [
        summarize_cell_binned_marker_trends(cell_scores, marker_col, score_col)
        for marker_col in available_marker_cols
    ]
    trend_tables = [tbl for tbl in trend_tables if not tbl.empty]
    if not trend_tables:
        return
    trends = pd.concat(trend_tables, ignore_index=True)
    trends.to_csv(OUTDIR / "external_cdkn_marker_sssp_cell_binned_sample_trends.csv", index=False)

    datasets = list(trends["dataset"].drop_duplicates())
    markers = list(trends["marker"].drop_duplicates())
    fig, axes = plt.subplots(
        len(markers),
        len(datasets),
        figsize=(5.7 * len(datasets), 4.3 * len(markers)),
        squeeze=False,
        sharey=True,
    )

    for col_idx, dataset in enumerate(datasets):
        ddf = trends.loc[trends["dataset"].eq(dataset)]
        sample_labels = sorted(ddf["sample_label"].unique())
        colors = dict(zip(sample_labels, plt.cm.tab20(np.linspace(0, 1, max(len(sample_labels), 3)))))
        for row_idx, marker in enumerate(markers):
            ax = axes[row_idx][col_idx]
            mdf = ddf.loc[ddf["marker"].eq(marker)]
            for sample_label, sdf in mdf.groupby("sample_label", dropna=False):
                sdf = sdf.sort_values("mean_marker")
                ax.plot(
                    sdf["mean_marker"],
                    sdf["mean_sssp"],
                    marker="o",
                    linewidth=1.8,
                    markersize=4,
                    color=colors[sample_label],
                    label=sample_label,
                )
            marker_label = marker.replace("_", " / ")
            ax.set_title(dataset)
            ax.set_xlabel(f"{marker_label} expression bin mean")
            ax.grid(alpha=0.25)
            if col_idx == 0:
                ax.set_ylabel("Mean Scanpy SSSP score\nwithin expression bin")
            ax.legend(fontsize=6, frameon=False, loc="best")

    fig.tight_layout()
    fig.savefig(OUTDIR / "external_cdkn_marker_vs_sssp164_cell_binned_sample_lines.png", dpi=300)
    fig.savefig(OUTDIR / "external_cdkn_marker_vs_sssp164_cell_binned_sample_lines.pdf")
    plt.close(fig)


def plot_snc_sssp_sample_relations(sample_summary: pd.DataFrame) -> None:
    if sample_summary.empty:
        return
    definitions = list(sample_summary["snc_definition"].drop_duplicates())
    datasets = list(sample_summary["dataset"].drop_duplicates())
    fig, axes = plt.subplots(
        len(definitions),
        len(datasets),
        figsize=(5.4 * len(datasets), 4.2 * len(definitions)),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    for row_idx, def_name in enumerate(definitions):
        for col_idx, dataset in enumerate(datasets):
            ax = axes[row_idx][col_idx]
            ddf = sample_summary.loc[
                sample_summary["snc_definition"].eq(def_name) & sample_summary["dataset"].eq(dataset)
            ].copy()
            group_col = visual_group_col(dataset)
            group_values = sorted(ddf[group_col].astype(str).unique())
            colors = dict(zip(group_values, plt.cm.Set2(np.linspace(0, 1, max(len(group_values), 3)))))
            for group_value, sdf in ddf.groupby(group_col, dropna=False):
                ax.scatter(
                    sdf["snc_fraction"],
                    sdf["sssp_high_fraction"],
                    s=90,
                    color=colors[str(group_value)],
                    edgecolor="black",
                    linewidth=0.4,
                    label=str(group_value),
                )
            for _, sample_row in ddf.iterrows():
                ax.text(
                    sample_row["snc_fraction"],
                    sample_row["sssp_high_fraction"],
                    display_label(sample_row, dataset),
                    fontsize=7,
                    ha="left",
                    va="bottom",
                )
            ax.set_title(f"{dataset}\n{def_name}")
            ax.grid(alpha=0.25)
            if row_idx == len(definitions) - 1:
                ax.set_xlabel("SnC fraction")
            if col_idx == 0:
                ax.set_ylabel("SSSP-high fraction")
            if row_idx == 0:
                ax.legend(title=visual_group_title(dataset), fontsize=7, title_fontsize=8, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(OUTDIR / "external_snc_fraction_vs_sssp_high_fraction_by_sample.png", dpi=300)
    fig.savefig(OUTDIR / "external_snc_fraction_vs_sssp_high_fraction_by_sample.pdf")
    plt.close(fig)


def plot_sssp_high_by_snc_status(cell_scores: pd.DataFrame) -> None:
    rows = []
    definitions = {
        "original_positive": "snc_original_positive",
        "cdkn_cycle_sasp_q4": "snc_cdkn_cycle_sasp_q4",
    }
    sample_cols = ["dataset", "donor", "site", "method", "age", "age_numeric"]
    for def_name, call_col in definitions.items():
        work = cell_scores.copy()
        work["snc_status"] = np.where(work[call_col], "SnC+", "SnC-")
        sample_level = (
            work.groupby(sample_cols + ["snc_status"], dropna=False)
            .agg(n_cells=("sssp_high_q4", "size"), sssp_high_fraction=("sssp_high_q4", "mean"))
            .reset_index()
        )
        sample_level["snc_definition"] = def_name
        rows.append(sample_level)
    plot_df = pd.concat(rows, ignore_index=True)
    plot_df.to_csv(OUTDIR / "external_sssp_high_fraction_by_snc_status_sample_level.csv", index=False)

    definitions = list(plot_df["snc_definition"].drop_duplicates())
    datasets = list(plot_df["dataset"].drop_duplicates())
    fig, axes = plt.subplots(
        len(definitions),
        len(datasets),
        figsize=(5.2 * len(datasets), 4.0 * len(definitions)),
        squeeze=False,
        sharey=True,
    )
    for row_idx, def_name in enumerate(definitions):
        for col_idx, dataset in enumerate(datasets):
            ax = axes[row_idx][col_idx]
            ddf = plot_df.loc[
                plot_df["snc_definition"].eq(def_name) & plot_df["dataset"].eq(dataset)
            ].copy()
            statuses = ["SnC-", "SnC+"]
            x_positions = np.arange(len(statuses))
            for status_idx, status in enumerate(statuses):
                sdf = ddf.loc[ddf["snc_status"].eq(status)]
                jitter = np.linspace(-0.12, 0.12, max(len(sdf), 1)) if len(sdf) else []
                ax.scatter(
                    np.repeat(x_positions[status_idx], len(sdf)) + jitter,
                    sdf["sssp_high_fraction"],
                    s=55,
                    alpha=0.85,
                    edgecolor="black",
                    linewidth=0.3,
                )
                if len(sdf):
                    ax.plot(
                        [x_positions[status_idx] - 0.22, x_positions[status_idx] + 0.22],
                        [sdf["sssp_high_fraction"].mean()] * 2,
                        color="black",
                        linewidth=2,
                    )
            ax.set_xticks(x_positions)
            ax.set_xticklabels(statuses)
            ax.set_title(f"{dataset}\n{def_name}")
            ax.grid(axis="y", alpha=0.25)
            if col_idx == 0:
                ax.set_ylabel("Sample-level fraction SSSP-high")
    fig.tight_layout()
    fig.savefig(OUTDIR / "external_sssp_high_fraction_by_snc_status.png", dpi=300)
    fig.savefig(OUTDIR / "external_sssp_high_fraction_by_snc_status.pdf")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    gene_sets = load_gene_sets()

    print("Reading GTEx skin subset")
    gtex_backed = ad.read_h5ad(GTEX_H5AD, backed="r")
    # Avoid an AnnData backed/raw sparse indexing bug when materializing a
    # subset. We only need X/layers/obs/var for this score audit.
    gtex_backed.raw = None
    gtex_mask_series = gtex_backed.obs["tissue"].astype(str).eq("skin")
    if GTEX_PREP_FILTER:
        gtex_mask_series = gtex_mask_series & gtex_backed.obs["prep"].astype(str).eq(GTEX_PREP_FILTER)
    gtex_mask = gtex_mask_series.to_numpy()
    gtex = gtex_backed[gtex_mask].to_memory()
    gtex_backed.file.close()

    print("Reading Tabula Sapiens skin")
    tabula = ad.read_h5ad(TABULA_H5AD)

    metadata = pd.concat(
        [
            summarize_metadata(gtex.obs, "GTEx_skin"),
            summarize_metadata(tabula.obs, "Tabula_Sapiens_skin"),
        ],
        ignore_index=True,
    )
    metadata.to_csv(OUTDIR / "external_dataset_metadata_summary.csv", index=False)

    gene_overlap = pd.concat(
        [
            summarize_gene_overlap("GTEx_skin", get_var_symbols(gtex, "Approved symbol"), gene_sets),
            summarize_gene_overlap("Tabula_Sapiens_skin", get_var_symbols(tabula, "gene_symbol"), gene_sets),
        ],
        ignore_index=True,
    )
    gene_overlap.to_csv(OUTDIR / "external_sssp_gene_overlap.csv", index=False)

    gtex_scores = harmonize_scores(add_scores(gtex, "GTEx_skin", gene_sets, layer="counts", symbol_col="Approved symbol"), "GTEx_skin")
    tabula_layer = "raw_counts" if "raw_counts" in tabula.layers else None
    tabula_scores = harmonize_scores(add_scores(tabula, "Tabula_Sapiens_skin", gene_sets, layer=tabula_layer, symbol_col="gene_symbol"), "Tabula_Sapiens_skin")
    cell_scores = pd.concat([gtex_scores, tabula_scores], ignore_index=True)
    cell_scores = add_snc_and_sssp_flags(cell_scores)
    cell_scores.to_csv(OUTDIR / "external_sssp_cell_scores.csv.gz", index=False)

    donor_summary, celltype_summary, trends = summarize_scores(cell_scores)
    donor_summary.to_csv(OUTDIR / "external_sssp_score_by_donor.csv", index=False)
    celltype_summary.to_csv(OUTDIR / "external_sssp_score_by_donor_celltype.csv", index=False)
    trends.to_csv(OUTDIR / "external_sssp_descriptive_age_trends.csv", index=False)
    marker_score_correlations = calculate_marker_score_correlations(cell_scores, donor_summary)
    marker_score_correlations.to_csv(OUTDIR / "external_cdkn_marker_sssp_score_correlations.csv", index=False)
    snc_sssp_sample, snc_sssp_overlap, snc_sssp_celltype = summarize_snc_sssp_relations(cell_scores)
    snc_sssp_sample.to_csv(OUTDIR / "external_snc_sssp_relation_by_sample.csv", index=False)
    snc_sssp_overlap.to_csv(OUTDIR / "external_snc_sssp_overlap_enrichment_descriptive.csv", index=False)
    snc_sssp_celltype.to_csv(OUTDIR / "external_snc_sssp_relation_by_major_cell_type.csv", index=False)
    plot_donor_scores(donor_summary)
    plot_marker_score_correlations(donor_summary)
    plot_marker_score_correlation_lines(donor_summary)
    plot_cell_binned_marker_trend_lines(cell_scores)
    plot_snc_sssp_sample_relations(snc_sssp_sample)
    plot_sssp_high_by_snc_status(cell_scores)

    notes = pd.DataFrame(
        [
            {
                "conclusion": "Use as validation/triangulation, not direct integration with Sinai.",
                "reason": "GTEx skin includes multiple nuclei-prep protocols, so GTEx plots separate points by prep method and prep-specific points should not be treated as independent biological donors; Tabula Sapiens skin has 4 donors with age confounded by donor/site/method. Both are useful for checking gene detectability, broad cell-type enrichment, and directionality, but are weak for independent age-trend claims.",
            }
        ]
    )
    notes.to_csv(OUTDIR / "external_validation_recommendation.csv", index=False)
    print(f"Done: {OUTDIR}")


if __name__ == "__main__":
    main()
