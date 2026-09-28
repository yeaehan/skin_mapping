"""
Squidpy ligand-receptor analysis for SSSP ligands by SnC status and minor cell type.

Runs age-stratified squidpy.gr.ligrec analyses using groups such as
"SnC+|Spinous Keratinocytes". Interactions are restricted to OmniPath-derived
LR pairs whose source/ligand gene is in the SSSP gene set and whose ligand and
receptor are both measured by Xenium.
"""

from __future__ import annotations

from pathlib import Path
import sys
import types
import warnings

import anndata as ad
import anndata._core.views as anndata_views
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix, csr_matrix

# Squidpy 1.5 imports spatialdata at package import time, but ligrec only needs
# AnnData expression and cluster labels. This avoids unrelated spatialdata/zarr
# image-stack imports in the current environment.
if "spatialdata" not in sys.modules:
    spatialdata_stub = types.ModuleType("spatialdata")
    spatialdata_stub.SpatialData = type("SpatialData", (), {})
    sys.modules["spatialdata"] = spatialdata_stub

# Squidpy 1.5 also imports sparse view classes that newer AnnData no longer
# exposes under these names.
if not hasattr(anndata_views, "SparseCSRView"):
    anndata_views.SparseCSRView = csr_matrix
if not hasattr(anndata_views, "SparseCSCView"):
    anndata_views.SparseCSCView = csc_matrix

import squidpy as sq


PROJECT_ROOT = Path(__file__).resolve().parents[3]
H5AD_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_sssp_20gene_scanpy_score_genes"
    / "xenium_20gene_sssp_scanpy_scored.h5ad"
)
GENE_CONFIG_PATH = Path(__file__).resolve().with_name("xenium_sssp_gene_set_config.csv")
LR_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "secretory_ligand_receptor_communication"
    / "secretory_ligand_receptor_pairs_from_omnipath.csv"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_squidpy_sssp_ligrec_by_snc_celltype"

GENE_SET_CONFIGS = [
    {
        "gene_set_id": "main_20",
        "membership_col": "in_xenium_main_20",
    },
    {
        "gene_set_id": "diff_filtered_17",
        "membership_col": "in_xenium_diff_filtered",
    },
]
AGE_ORDER = ["Old", "Young"]
MIN_GROUP_CELLS = 5
N_PERMUTATIONS = 1000
EXPR_THRESHOLD = 0.05
RANDOM_SEED = 260727
TOP_N_PLOT = 30
NATIVE_LIGREC_ALPHA = 0.005
TOP_NATIVE_GROUPS = 4
TOP_NATIVE_LR_PAIRS = 12
PRIORITY_PVALUE_THRESHOLD = 0.05


def clean_label(series: pd.Series) -> pd.Series:
    out = series.astype("string").fillna("Unknown").astype(str)
    return out.replace({"": "Unknown", "nan": "Unknown", "None": "Unknown"})


def load_adata() -> ad.AnnData:
    adata = ad.read_h5ad(H5AD_PATH)
    needed = {"age", "snc", "cell_type_major", "cell_type_minor"}
    missing = sorted(needed - set(adata.obs.columns))
    if missing:
        raise KeyError(f"Missing required obs columns in {H5AD_PATH}: {', '.join(missing)}")

    adata.obs["age"] = clean_label(adata.obs["age"])
    adata.obs["snc"] = clean_label(adata.obs["snc"])
    adata.obs["cell_type_major"] = clean_label(adata.obs["cell_type_major"])
    adata.obs["cell_type_minor"] = clean_label(adata.obs["cell_type_minor"])
    adata.obs["snc_minor_group"] = adata.obs["snc"] + "|" + adata.obs["cell_type_minor"]
    adata.var_names = adata.var_names.astype(str).str.upper()
    return adata


def load_gene_sets() -> dict[str, list[str]]:
    config = pd.read_csv(GENE_CONFIG_PATH)
    config["gene"] = config["gene"].astype(str).str.upper()
    gene_sets = {}
    for gene_set_config in GENE_SET_CONFIGS:
        gene_set_id = str(gene_set_config["gene_set_id"])
        membership_col = str(gene_set_config["membership_col"])
        gene_sets[gene_set_id] = config.loc[config[membership_col].astype(bool), "gene"].tolist()
    return gene_sets


def load_interactions(adata: ad.AnnData, sssp_genes: list[str]) -> pd.DataFrame:
    lr = pd.read_csv(LR_PATH)
    required = {"ligand", "receptor"}
    missing = sorted(required - set(lr.columns))
    if missing:
        raise KeyError(f"Missing columns in {LR_PATH}: {', '.join(missing)}")

    lr["source"] = lr["ligand"].astype(str).str.upper()
    lr["target"] = lr["receptor"].astype(str).str.upper()
    var_genes = set(map(str.upper, adata.var_names))
    sssp_gene_set = set(map(str.upper, sssp_genes))
    lr = lr[lr["source"].isin(sssp_gene_set) & lr["source"].isin(var_genes) & lr["target"].isin(var_genes)].copy()
    lr = lr[["source", "target"]].drop_duplicates().sort_values(["source", "target"])
    if lr.empty:
        raise RuntimeError("No SSSP ligand-receptor pairs remain after filtering to Xenium genes.")
    return lr


def subset_age_groups(adata: ad.AnnData, age: str) -> tuple[ad.AnnData, pd.DataFrame]:
    sub = adata[adata.obs["age"].astype(str).eq(age)].copy()
    group_counts = sub.obs["snc_minor_group"].value_counts().rename_axis("snc_minor_group").reset_index(name="n_cells")
    group_counts[["snc", "cell_type_minor"]] = group_counts["snc_minor_group"].str.split("|", n=1, expand=True)
    keep_groups = group_counts.loc[group_counts["n_cells"].ge(MIN_GROUP_CELLS), "snc_minor_group"].tolist()
    sub = sub[sub.obs["snc_minor_group"].isin(keep_groups)].copy()
    sub.obs["snc_minor_group"] = pd.Categorical(
        sub.obs["snc_minor_group"],
        categories=sorted(keep_groups, key=lambda x: (0 if x.startswith("SnC+") else 1, x)),
    )
    group_counts["age"] = age
    group_counts["included_min_cells"] = group_counts["snc_minor_group"].isin(keep_groups)
    return sub, group_counts


def run_ligrec(sub: ad.AnnData, interactions: pd.DataFrame) -> dict[str, pd.DataFrame]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = sq.gr.ligrec(
            sub,
            cluster_key="snc_minor_group",
            interactions=interactions,
            n_perms=N_PERMUTATIONS,
            threshold=EXPR_THRESHOLD,
            corr_method=None,
            use_raw=False,
            copy=True,
            seed=RANDOM_SEED,
            n_jobs=1,
        )
    return res


def snc_plus_groups_by_compartment(sub: ad.AnnData) -> dict[str, list[str]]:
    group_major = (
        sub.obs.groupby("snc_minor_group", observed=True)["cell_type_major"]
        .agg(lambda values: values.value_counts().idxmax())
        .to_dict()
    )
    groups = [str(group) for group in sub.obs["snc_minor_group"].cat.categories if str(group).startswith("SnC+|")]
    skin_groups = [group for group in groups if group_major.get(group) != "Immune Cells"]
    immune_groups = [group for group in groups if group_major.get(group) == "Immune Cells"]
    return {
        "all": groups,
        "skin": skin_groups,
        "immune": immune_groups,
    }


def ligrec_to_long(res: dict[str, pd.DataFrame], gene_set_id: str, age: str) -> pd.DataFrame:
    means = res["means"].copy()
    pvals = res["pvalues"].copy()
    rows = []
    for interaction in means.index:
        ligand, receptor = interaction
        for source_group, target_group in means.columns:
            mean_value = means.loc[interaction, (source_group, target_group)]
            pvalue = pvals.loc[interaction, (source_group, target_group)]
            if pd.isna(mean_value):
                continue
            source_snc, source_cell_type = str(source_group).split("|", 1)
            target_snc, target_cell_type = str(target_group).split("|", 1)
            rows.append(
                {
                    "gene_set_id": gene_set_id,
                    "age": age,
                    "ligand": ligand,
                    "receptor": receptor,
                    "lr_pair": f"{ligand}->{receptor}",
                    "source_group": source_group,
                    "target_group": target_group,
                    "source_snc": source_snc,
                    "source_cell_type_minor": source_cell_type,
                    "target_snc": target_snc,
                    "target_cell_type_minor": target_cell_type,
                    "mean_expression_product": float(mean_value),
                    "pvalue": float(pvalue) if pd.notna(pvalue) else np.nan,
                    "source_is_snc_plus": source_snc == "SnC+",
                    "target_is_snc_plus": target_snc == "SnC+",
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["minus_log10_p"] = -np.log10(out["pvalue"].clip(lower=1e-300))
        out["significant_p05"] = out["pvalue"].le(0.05)
    return out


def summarize_focus(long: pd.DataFrame) -> pd.DataFrame:
    if long.empty:
        return pd.DataFrame()
    focus = long[long["source_is_snc_plus"] | long["target_is_snc_plus"]].copy()
    focus = focus.sort_values(["pvalue", "mean_expression_product"], ascending=[True, False])
    return focus


def plot_focus_dotplot(focus: pd.DataFrame, gene_set_id: str, age: str, mode: str) -> None:
    if focus.empty:
        return
    if mode == "outgoing":
        work = focus[focus["source_is_snc_plus"]].copy()
        title = f"{age} {gene_set_id}: SnC+ cell types sending SSSP ligands"
    elif mode == "incoming":
        work = focus[focus["target_is_snc_plus"]].copy()
        title = f"{age} {gene_set_id}: SnC+ cell types receiving SSSP-ligand signals"
    else:
        raise ValueError(mode)

    if work.empty:
        return
    work = work.sort_values(["pvalue", "mean_expression_product"], ascending=[True, False]).head(TOP_N_PLOT).copy()
    work["interaction_label"] = work["lr_pair"] + "\n" + work["source_group"] + " -> " + work["target_group"]
    work = work.iloc[::-1].copy()

    max_mean = work["mean_expression_product"].max()
    sizes = 35 + 260 * np.sqrt(work["mean_expression_product"] / max(max_mean, 1e-9))
    colors = work["minus_log10_p"].clip(upper=6)

    fig_height = max(5.2, 0.32 * len(work) + 1.2)
    fig, ax = plt.subplots(figsize=(9.5, fig_height))
    scatter = ax.scatter(
        work["mean_expression_product"],
        np.arange(len(work)),
        s=sizes,
        c=colors,
        cmap="magma",
        edgecolor="#333333",
        linewidth=0.4,
        alpha=0.88,
    )
    ax.set_yticks(np.arange(len(work)))
    ax.set_yticklabels(work["interaction_label"], fontsize=8)
    ax.set_xlabel("Squidpy LR mean expression product")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("-log10 permutation p")
    fig.tight_layout()
    fig.savefig(OUTDIR / f"{gene_set_id}_{age.lower()}_snc_plus_{mode}_sssp_ligrec_dotplot.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / f"{gene_set_id}_{age.lower()}_snc_plus_{mode}_sssp_ligrec_dotplot.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_snc_status_summary(long: pd.DataFrame, gene_set_id: str, age: str) -> None:
    if long.empty:
        return
    work = long.copy()
    work["direction_class"] = np.select(
        [
            work["source_is_snc_plus"] & work["target_is_snc_plus"],
            work["source_is_snc_plus"] & ~work["target_is_snc_plus"],
            ~work["source_is_snc_plus"] & work["target_is_snc_plus"],
        ],
        ["SnC+ -> SnC+", "SnC+ -> SnC-", "SnC- -> SnC+"],
        default="SnC- -> SnC-",
    )
    summary = (
        work.groupby(["direction_class", "lr_pair"], observed=True)
        .agg(
            max_mean=("mean_expression_product", "max"),
            min_pvalue=("pvalue", "min"),
            n_group_pairs=("target_group", "size"),
        )
        .reset_index()
        .sort_values(["min_pvalue", "max_mean"], ascending=[True, False])
        .head(40)
    )
    summary["minus_log10_p"] = -np.log10(summary["min_pvalue"].clip(lower=1e-300))
    summary = summary.iloc[::-1].copy()
    colors = {
        "SnC+ -> SnC+": "#6A3D9A",
        "SnC+ -> SnC-": "#D62728",
        "SnC- -> SnC+": "#2CA02C",
        "SnC- -> SnC-": "#9AA0A6",
    }
    fig, ax = plt.subplots(figsize=(9.0, max(5.5, 0.26 * len(summary) + 1.2)))
    ax.barh(
        np.arange(len(summary)),
        summary["max_mean"],
        color=[colors[x] for x in summary["direction_class"]],
        edgecolor="#333333",
        linewidth=0.45,
        alpha=0.86,
    )
    ax.set_yticks(np.arange(len(summary)))
    ax.set_yticklabels(summary["direction_class"] + " | " + summary["lr_pair"], fontsize=8)
    ax.set_xlabel("Max Squidpy LR mean expression product")
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / f"{gene_set_id}_{age.lower()}_snc_status_direction_sssp_ligrec_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / f"{gene_set_id}_{age.lower()}_snc_status_direction_sssp_ligrec_summary.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_native_squidpy_ligrec(
    res: dict[str, pd.DataFrame],
    snc_plus_groups: list[str],
    gene_set_id: str,
    age: str,
    mode: str,
    compartment: str = "all",
) -> None:
    if not snc_plus_groups:
        return

    if mode == "outgoing":
        kwargs = {"source_groups": snc_plus_groups}
    elif mode == "incoming":
        kwargs = {"target_groups": snc_plus_groups}
    else:
        raise ValueError(mode)

    width = max(10.0, min(34.0, 1.4 * len(snc_plus_groups) + 9.0))
    height = 5.8
    sq.pl.ligrec(
        res,
        alpha=NATIVE_LIGREC_ALPHA,
        remove_empty_interactions=True,
        remove_nonsig_interactions=False,
        title="",
        figsize=(width, height),
        **kwargs,
    )
    fig = plt.gcf()
    fig.savefig(
        OUTDIR / f"{gene_set_id}_{age.lower()}_{compartment}_snCplus_{mode}_native_squidpy_ligrec_dotplot.png",
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        OUTDIR / f"{gene_set_id}_{age.lower()}_{compartment}_snCplus_{mode}_native_squidpy_ligrec_dotplot.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)


def rank_snc_plus_groups(
    long: pd.DataFrame,
    snc_plus_groups: list[str],
    mode: str,
) -> pd.DataFrame:
    if long.empty or not snc_plus_groups:
        return pd.DataFrame()

    if mode == "outgoing":
        group_col = "source_group"
        work = long[long[group_col].isin(snc_plus_groups)].copy()
    elif mode == "incoming":
        group_col = "target_group"
        work = long[long[group_col].isin(snc_plus_groups)].copy()
    else:
        raise ValueError(mode)

    if work.empty:
        return pd.DataFrame()

    work["is_priority"] = work["pvalue"].le(PRIORITY_PVALUE_THRESHOLD)
    ranked = (
        work.groupby(group_col, observed=True)
        .agg(
            n_priority_links=("is_priority", "sum"),
            n_priority_lr_pairs=("lr_pair", lambda x: int(work.loc[x.index, "is_priority"].groupby(x).any().sum())),
            min_pvalue=("pvalue", "min"),
            max_mean=("mean_expression_product", "max"),
        )
        .reset_index()
        .rename(columns={group_col: "snc_plus_group"})
        .sort_values(["n_priority_links", "n_priority_lr_pairs", "max_mean", "min_pvalue"], ascending=[False, False, False, True])
    )
    return ranked


def choose_priority_lr_pairs(
    long: pd.DataFrame,
    selected_groups: list[str],
    mode: str,
) -> list[str]:
    if long.empty or not selected_groups:
        return []

    if mode == "outgoing":
        work = long[long["source_group"].isin(selected_groups)].copy()
    elif mode == "incoming":
        work = long[long["target_group"].isin(selected_groups)].copy()
    else:
        raise ValueError(mode)

    if work.empty:
        return []

    work["is_priority"] = work["pvalue"].le(PRIORITY_PVALUE_THRESHOLD)
    ranked = (
        work.groupby("lr_pair", observed=True)
        .agg(
            n_priority_links=("is_priority", "sum"),
            min_pvalue=("pvalue", "min"),
            max_mean=("mean_expression_product", "max"),
        )
        .reset_index()
        .sort_values(["n_priority_links", "max_mean", "min_pvalue"], ascending=[False, False, True])
    )
    return ranked.head(TOP_NATIVE_LR_PAIRS)["lr_pair"].tolist()


def subset_ligrec_result(res: dict[str, pd.DataFrame], selected_lr_pairs: list[str]) -> dict[str, pd.DataFrame]:
    if not selected_lr_pairs:
        return res

    selected = set(selected_lr_pairs)
    out = {}
    for key, table in res.items():
        if not isinstance(table, pd.DataFrame):
            out[key] = table
            continue
        if isinstance(table.index, pd.MultiIndex):
            keep = [f"{idx[0]}->{idx[1]}" in selected for idx in table.index]
            out[key] = table.loc[keep].copy()
        else:
            out[key] = table.copy()
    return out


def remove_obsolete_outputs() -> dict[str, int]:
    patterns = {
        "prioritized": "*prioritized*",
        "incoming": "*incoming*",
        "snc_status_direction": "*snc_status_direction*",
        "legacy_mixed_focus": "*_snc_plus_focus_ligrec_results.csv",
    }
    removed = {}
    for label, pattern in patterns.items():
        removed[label] = 0
        for path in OUTDIR.glob(pattern):
            if path.is_file():
                path.unlink()
                removed[label] += 1
    return removed


def write_plot_interpretation() -> None:
    note = pd.DataFrame(
        [
            {
                "plot_family": "native_squidpy_ligrec_dotplot",
                "included_direction": "SnC+ sender/source groups only",
                "receiver_groups": "All eligible target groups are retained, including SnC- and SnC+ groups.",
                "how_to_read": "Top labels are ligand senders; bottom labels are receptor receivers. Read as ligand from top cell group to receptor on bottom cell group.",
            },
            {
                "plot_family": "snc_plus_outgoing_sssp_ligrec_dotplot",
                "included_direction": "SnC+ sender/source groups only",
                "receiver_groups": "All eligible target groups are retained, including SnC- and SnC+ groups.",
                "how_to_read": "Each row label is ligand->receptor followed by sender group -> receiver group.",
            },
        ]
    )
    note.to_csv(OUTDIR / "xenium_squidpy_sssp_ligrec_plot_interpretation.csv", index=False)


def remove_obsolete_prioritized_outputs() -> int:
    removed = 0
    for path in OUTDIR.glob("*prioritized*"):
        if path.is_file():
            path.unlink()
            removed += 1
    return removed


def plot_prioritized_native_squidpy_ligrec(
    res: dict[str, pd.DataFrame],
    long: pd.DataFrame,
    snc_plus_groups: list[str],
    gene_set_id: str,
    age: str,
    mode: str,
    compartment: str,
) -> pd.DataFrame:
    ranked_groups = rank_snc_plus_groups(long, snc_plus_groups, mode)
    if ranked_groups.empty:
        return ranked_groups

    selected_groups = ranked_groups.head(TOP_NATIVE_GROUPS)["snc_plus_group"].tolist()
    selected_lr_pairs = choose_priority_lr_pairs(long, selected_groups, mode)
    if not selected_lr_pairs:
        return ranked_groups

    ranked_groups["selected_for_prioritized_plot"] = ranked_groups["snc_plus_group"].isin(selected_groups)
    ranked_groups["mode"] = mode
    ranked_groups["compartment"] = compartment
    ranked_groups["age"] = age
    ranked_groups["gene_set_id"] = gene_set_id
    ranked_groups["top_lr_pairs_for_plot"] = "; ".join(selected_lr_pairs)

    prioritized_res = subset_ligrec_result(res, selected_lr_pairs)
    if prioritized_res["means"].empty:
        return ranked_groups

    if mode == "outgoing":
        kwargs = {"source_groups": selected_groups}
        title = f"{age} {gene_set_id}: top {compartment} SnC+ SSSP ligand sources"
    elif mode == "incoming":
        kwargs = {"target_groups": selected_groups}
        title = f"{age} {gene_set_id}: top {compartment} SnC+ SSSP signal receivers"
    else:
        raise ValueError(mode)

    width = max(8.0, min(18.0, 1.6 * len(selected_groups) + 7.0))
    height = max(4.8, min(9.0, 0.45 * len(selected_lr_pairs) + 3.8))
    sq.pl.ligrec(
        prioritized_res,
        alpha=NATIVE_LIGREC_ALPHA,
        remove_empty_interactions=True,
        remove_nonsig_interactions=False,
        title=title,
        figsize=(width, height),
        **kwargs,
    )
    fig = plt.gcf()
    fig.savefig(
        OUTDIR / f"{gene_set_id}_{age.lower()}_{compartment}_snCplus_{mode}_prioritized_native_squidpy_ligrec_dotplot.png",
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        OUTDIR / f"{gene_set_id}_{age.lower()}_{compartment}_snCplus_{mode}_prioritized_native_squidpy_ligrec_dotplot.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)
    return ranked_groups


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    removed_obsolete = remove_obsolete_outputs()
    write_plot_interpretation()
    adata = load_adata()
    gene_sets = load_gene_sets()

    all_long = []
    all_group_counts = []
    all_interactions = []

    for gene_set_id, sssp_genes in gene_sets.items():
        interactions = load_interactions(adata, sssp_genes)
        interactions.assign(gene_set_id=gene_set_id).to_csv(
            OUTDIR / f"{gene_set_id}_sssp_ligand_receptor_pairs_used.csv", index=False
        )
        all_interactions.append(interactions.assign(gene_set_id=gene_set_id))

        for age in AGE_ORDER:
            sub, group_counts = subset_age_groups(adata, age)
            group_counts.insert(0, "gene_set_id", gene_set_id)
            all_group_counts.append(group_counts)
            res = run_ligrec(sub, interactions)
            compartment_groups = snc_plus_groups_by_compartment(sub)
            for compartment, snc_plus_groups in compartment_groups.items():
                plot_native_squidpy_ligrec(res, snc_plus_groups, gene_set_id, age, "outgoing", compartment)
            long = ligrec_to_long(res, gene_set_id, age)
            all_long.append(long)
            focus = summarize_focus(long)
            sender_focus = focus[focus["source_is_snc_plus"]].copy()
            sender_focus.to_csv(OUTDIR / f"{gene_set_id}_{age.lower()}_snc_plus_sender_focus_ligrec_results.csv", index=False)
            plot_focus_dotplot(sender_focus, gene_set_id, age, "outgoing")

    long_all = pd.concat(all_long, ignore_index=True)
    group_counts_all = pd.concat(all_group_counts, ignore_index=True)
    interactions_all = pd.concat(all_interactions, ignore_index=True)

    long_all.to_csv(OUTDIR / "xenium_squidpy_sssp_ligrec_all_group_pairs.csv.gz", index=False)
    group_counts_all.to_csv(OUTDIR / "xenium_squidpy_sssp_ligrec_group_counts.csv", index=False)
    interactions_all.to_csv(OUTDIR / "xenium_squidpy_sssp_ligrec_pairs_used_all_gene_sets.csv", index=False)

    summary = (
        long_all[long_all["source_is_snc_plus"]]
        .groupby(["gene_set_id", "age", "source_is_snc_plus", "target_is_snc_plus"], observed=True)
        .agg(
            n_lr_group_pairs=("lr_pair", "size"),
            n_unique_lr_pairs=("lr_pair", "nunique"),
            min_pvalue=("pvalue", "min"),
            max_mean=("mean_expression_product", "max"),
        )
        .reset_index()
    )
    summary.to_csv(OUTDIR / "xenium_squidpy_sssp_ligrec_snc_focus_summary.csv", index=False)

    run_summary = pd.DataFrame(
        [
            {
                "input_h5ad": str(H5AD_PATH),
                "interaction_source": str(LR_PATH),
                "gene_config": str(GENE_CONFIG_PATH),
                "min_group_cells": MIN_GROUP_CELLS,
                "n_permutations": N_PERMUTATIONS,
                "expression_threshold": EXPR_THRESHOLD,
                "cluster_key": "snc_minor_group",
                "method": "squidpy.gr.ligrec, age-stratified, SnC+ sender/source groups only; target/receiver groups can be SnC- or SnC+",
                "obsolete_prioritized_files_removed": removed_obsolete["prioritized"],
                "obsolete_incoming_files_removed": removed_obsolete["incoming"],
                "obsolete_snc_status_direction_files_removed": removed_obsolete["snc_status_direction"],
                "obsolete_legacy_mixed_focus_files_removed": removed_obsolete["legacy_mixed_focus"],
            }
        ]
    )
    run_summary.to_csv(OUTDIR / "xenium_squidpy_sssp_ligrec_run_summary.csv", index=False)

    print("Squidpy SSSP ligand-receptor by SnC/cell_type_minor complete.")
    print(summary.to_string(index=False))
    print(f"\nWrote outputs to: {OUTDIR}")


if __name__ == "__main__":
    main()
