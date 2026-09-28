"""Squidpy SSSP communication by SnC cell type and CDKN phenotype.

SnC+ source cells are divided into p21+, p16+, and p16+p21+ groups within
minor cell type. All other cells are receiver groups. Analyses are run per
biological sample, using only sssp_xenium ligands and measured OmniPath
receptors. This keeps Squidpy as the single cell-communication package.

Squidpy ligrec is expression/permutation based; it does not itself require
physical adjacency. Sample-level results are therefore the inferential unit.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import anndata as ad
import anndata._core.views as anndata_views
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.sparse import csc_matrix, csr_matrix

if not hasattr(anndata_views, "SparseCSRView"):
    anndata_views.SparseCSRView = csr_matrix
if not hasattr(anndata_views, "SparseCSCView"):
    anndata_views.SparseCSCView = csc_matrix
import squidpy as sq

ROOT = Path(__file__).resolve().parents[3]
H5AD = ROOT / "data" / "Xenium" / "xenium_full.h5ad"
SSSP_JSON = ROOT / "code" / "shared" / "sssp.json"
LR_CSV = ROOT / "results" / "260727" / "secretory_ligand_receptor_communication" / "secretory_ligand_receptor_pairs_from_omnipath.csv"
OUTDIR = ROOT / "results" / "260831" / "xenium_squidpy_sssp_cdkn_celltype_communication"

MIN_GROUP_CELLS = 10
N_PERMUTATIONS = 500
EXPR_THRESHOLD = 0.05
SEED = 260831
EPS = 1e-8
CDKN_ORDER = ["p21+", "p16+", "p16+p21+"]


def clean(x: pd.Series) -> pd.Series:
    return x.astype("string").fillna("Unknown").astype(str).replace({"": "Unknown", "nan": "Unknown"})


def bh(values: pd.Series) -> np.ndarray:
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    ans = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        ids = np.flatnonzero(ok)
        order = ids[np.argsort(p[ok])]
        q = p[order] * len(order) / np.arange(1, len(order) + 1)
        ans[order] = np.clip(np.minimum.accumulate(q[::-1])[::-1], 0, 1)
    return ans


def load_data() -> tuple[ad.AnnData, list[str], pd.DataFrame]:
    a = ad.read_h5ad(H5AD)
    # The exported H5AD uses sample_csv/sample_concat; standardize locally.
    if "sample" not in a.obs.columns:
        for candidate in ["sample_csv", "sample_concat"]:
            if candidate in a.obs.columns:
                a.obs["sample"] = a.obs[candidate].astype(str)
                break
    needed = {"sample", "age", "snc", "cell_type_minor"}
    missing = sorted(needed - set(a.obs.columns))
    if missing:
        raise KeyError(f"Missing Xenium metadata: {missing}")
    a.var_names = a.var_names.astype(str).str.upper()
    if not a.var_names.is_unique:
        raise ValueError("Upper-cased Xenium gene symbols are not unique")
    for col in needed:
        a.obs[col] = clean(a.obs[col])

    # CDKN state is derived from detected counts, matching the R burden plot.
    for gene in ["CDKN1A", "CDKN2A"]:
        if gene not in a.var_names:
            raise KeyError(f"{gene} is absent from Xenium")
    counts = a[:, ["CDKN1A", "CDKN2A"]].X
    counts = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)
    p21, p16 = counts[:, 0] > 0, counts[:, 1] > 0
    a.obs["cdkn_state"] = np.select(
        [p21 & p16, p21, p16], ["p16+p21+", "p21+", "p16+"], default="p16-p21-"
    )
    is_sender = a.obs["snc"].eq("SnC+") & a.obs["cdkn_state"].isin(CDKN_ORDER)
    a.obs["communication_group"] = np.where(
        is_sender,
        "Sender|" + a.obs["cdkn_state"] + "|" + a.obs["cell_type_minor"],
        "Receiver|" + a.obs["cell_type_minor"],
    )

    genes = list(dict.fromkeys(str(g).upper() for g in json.loads(SSSP_JSON.read_text())["sssp_xenium"]))
    lr = pd.read_csv(LR_CSV)
    lr["source"] = lr["ligand"].astype(str).str.upper()
    lr["target"] = lr["receptor"].astype(str).str.upper()
    available = set(a.var_names)
    lr = (lr[lr["source"].isin(genes) & lr["source"].isin(available) & lr["target"].isin(available)]
          [["source", "target"]].drop_duplicates().sort_values(["source", "target"]).reset_index(drop=True))
    if lr.empty:
        raise RuntimeError("No measured sssp_xenium ligand-receptor pairs")

    # Preserve counts for CDKN calling above, then use log-normalized data for ligrec.
    sc.pp.normalize_total(a, target_sum=10000)
    sc.pp.log1p(a)
    return a, genes, lr


def group_audit(a: ad.AnnData) -> pd.DataFrame:
    d = (a.obs.groupby(["sample", "age", "communication_group"], observed=True).size()
         .rename("n_cells").reset_index())
    d["included"] = d["n_cells"].ge(MIN_GROUP_CELLS)
    d["is_sender"] = d["communication_group"].str.startswith("Sender|")
    return d


def run_sample(a: ad.AnnData, lr: pd.DataFrame, sample: str) -> pd.DataFrame:
    sub = a[a.obs["sample"].eq(sample)].copy()
    counts = sub.obs["communication_group"].value_counts()
    groups = sorted(counts[counts >= MIN_GROUP_CELLS].index)
    sender_groups = [g for g in groups if g.startswith("Sender|")]
    receiver_groups = [g for g in groups if g.startswith("Receiver|")]
    if not sender_groups or not receiver_groups:
        return pd.DataFrame()
    sub = sub[sub.obs["communication_group"].isin(groups)].copy()
    sub.obs["communication_group"] = pd.Categorical(sub.obs["communication_group"], categories=groups)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = sq.gr.ligrec(
            sub, cluster_key="communication_group", interactions=lr,
            n_perms=N_PERMUTATIONS, threshold=EXPR_THRESHOLD, corr_method=None,
            use_raw=False, copy=True, seed=SEED, n_jobs=1,
        )
    means, pvals = result["means"], result["pvalues"]
    means.index.names = pvals.index.names = ["ligand", "receptor"]
    means.columns.names = pvals.columns.names = ["source_group", "target_group"]
    m = means.stack(["source_group", "target_group"], future_stack=True).rename("mean_expression_product")
    p = pvals.stack(["source_group", "target_group"], future_stack=True).rename("pvalue")
    d = pd.concat([m, p], axis=1).reset_index().dropna(subset=["mean_expression_product"])
    d = d[d["source_group"].isin(sender_groups) & d["target_group"].isin(receiver_groups)].copy()
    if d.empty:
        return d
    for col in ["mean_expression_product", "pvalue"]:
        if isinstance(d[col].dtype, pd.SparseDtype):
            d[col] = d[col].sparse.to_dense()
        d[col] = pd.to_numeric(d[col], errors="coerce")
    d[["source_role", "cdkn_state", "source_cell_type"]] = d["source_group"].str.split("|", n=2, expand=True)
    d[["target_role", "target_cell_type"]] = d["target_group"].str.split("|", n=1, expand=True)
    d["sample"] = sample
    d["age"] = str(sub.obs["age"].iloc[0])
    d["lr_pair"] = d["ligand"] + "->" + d["receptor"]
    d["fdr_within_sample"] = bh(d["pvalue"])
    d["nominal_p05"] = d["pvalue"].le(0.05)
    d["fdr05"] = d["fdr_within_sample"].le(0.05)
    return d.drop(columns=["source_role", "target_role"])


def compare_p16_p21(edges: pd.DataFrame) -> pd.DataFrame:
    keys = ["sample", "age", "source_cell_type", "target_cell_type", "ligand", "receptor", "lr_pair"]
    value = ["mean_expression_product", "pvalue", "fdr_within_sample"]
    p21 = edges[edges["cdkn_state"].eq("p21+")][keys + value].rename(columns={
        "mean_expression_product": "p21_mean_product", "pvalue": "p21_pvalue",
        "fdr_within_sample": "p21_fdr"})
    p16 = edges[edges["cdkn_state"].eq("p16+")][keys + value].rename(columns={
        "mean_expression_product": "p16_mean_product", "pvalue": "p16_pvalue",
        "fdr_within_sample": "p16_fdr"})
    d = p16.merge(p21, on=keys, how="inner")
    d["p16_minus_p21"] = d["p16_mean_product"] - d["p21_mean_product"]
    d["p16_over_p21"] = (d["p16_mean_product"] + EPS) / (d["p21_mean_product"] + EPS)
    d["log2_p16_over_p21"] = np.log2(d["p16_over_p21"])
    return d


def summarize_profiles(edges: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    ligand = edges.groupby(
        ["sample", "age", "source_cell_type", "cdkn_state", "ligand"], observed=True
    ).agg(
        outgoing_sssp_product=("mean_expression_product", "sum"),
        n_receiver_types=("target_cell_type", "nunique"),
        n_lr_pairs=("lr_pair", "nunique"), n_fdr_edges=("fdr05", "sum")
    ).reset_index()
    receiver = edges.groupby(
        ["sample", "age", "source_cell_type", "cdkn_state", "target_cell_type"], observed=True
    ).agg(
        incoming_product=("mean_expression_product", "sum"),
        n_ligands=("ligand", "nunique"), n_lr_pairs=("lr_pair", "nunique"),
        n_fdr_edges=("fdr05", "sum")
    ).reset_index()
    return ligand, receiver


def plot_ligand_heatmap(profile: pd.DataFrame) -> None:
    if profile.empty:
        return
    d = profile.copy()
    d["sender"] = d["source_cell_type"] + " | " + d["cdkn_state"]
    p = d.pivot_table(index="sender", columns="ligand", values="outgoing_sssp_product", aggfunc="mean", fill_value=0)
    fig, ax = plt.subplots(figsize=(max(9, .65 * len(p.columns)), max(6, .35 * len(p))))
    sns.heatmap(np.log1p(p), cmap="magma", ax=ax)
    ax.set_title("Xenium SnC+ outgoing SSSP profile by cell type and CDKN state")
    ax.set_xlabel("SSSP ligand"); ax.set_ylabel("Sender cell type | CDKN state")
    fig.tight_layout()
    fig.savefig(OUTDIR / "outgoing_sssp_ligand_profile_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "outgoing_sssp_ligand_profile_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_cdkn_contrast(contrast: pd.DataFrame) -> None:
    if contrast.empty:
        return
    d = contrast.groupby(["source_cell_type", "ligand"], observed=True).agg(
        median_log2_p16_over_p21=("log2_p16_over_p21", "median"),
        n_samples=("sample", "nunique")
    ).reset_index()
    d = d[d["n_samples"].ge(2)]
    p = d.pivot(index="source_cell_type", columns="ligand", values="median_log2_p16_over_p21")
    if p.empty:
        return
    lim = np.nanpercentile(np.abs(p.to_numpy()), 95)
    fig, ax = plt.subplots(figsize=(max(9, .65 * len(p.columns)), max(5, .4 * len(p))))
    sns.heatmap(p, cmap="vlag", center=0, vmin=-lim, vmax=lim, ax=ax)
    ax.set_title("Matched p16+ versus p21+ SSSP communication\npositive = higher in p16+")
    ax.set_xlabel("SSSP ligand"); ax.set_ylabel("Sender cell type")
    fig.tight_layout()
    fig.savefig(OUTDIR / "p16_vs_p21_sssp_ligand_contrast_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "p16_vs_p21_sssp_ligand_contrast_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    a, genes, lr = load_data()
    lr.to_csv(OUTDIR / "sssp_xenium_omnipath_lr_pairs_used.csv", index=False)
    pd.DataFrame({"sssp_xenium_gene": genes,
                  "present_in_xenium": [g in set(a.var_names) for g in genes],
                  "used_as_ligand": [g in set(lr.source) for g in genes]}).to_csv(
        OUTDIR / "sssp_xenium_gene_availability.csv", index=False)
    audit = group_audit(a)
    audit.to_csv(OUTDIR / "sample_cdkn_celltype_group_audit.csv", index=False)
    edge_path = OUTDIR / "xenium_squidpy_sssp_cdkn_celltype_edges.csv.gz"
    if edge_path.exists():
        print(f"Reusing completed Squidpy edges: {edge_path}", flush=True)
        edges = pd.read_csv(edge_path)
    else:
        frames = []
        for sample in sorted(a.obs["sample"].unique()):
            print(f"Running Squidpy: {sample}", flush=True)
            d = run_sample(a, lr, sample)
            if not d.empty:
                frames.append(d)
        edges = pd.concat(frames, ignore_index=True)
        edges.to_csv(edge_path, index=False)
    contrast = compare_p16_p21(edges)
    contrast.to_csv(OUTDIR / "matched_p16_vs_p21_sssp_edges.csv.gz", index=False)
    ligand, receiver = summarize_profiles(edges)
    ligand.to_csv(OUTDIR / "outgoing_sssp_ligand_profiles_by_sample.csv", index=False)
    receiver.to_csv(OUTDIR / "receiver_profiles_by_sample.csv", index=False)
    plot_ligand_heatmap(ligand)
    plot_cdkn_contrast(contrast)
    summary = {
        "input": str(H5AD), "samples": int(a.obs["sample"].nunique()),
        "sssp_xenium_genes": len(genes), "sssp_ligands_with_lr_pairs": int(lr.source.nunique()),
        "lr_pairs": len(lr), "min_group_cells": MIN_GROUP_CELLS,
        "permutations_per_sample": N_PERMUTATIONS,
        "communication_package": "Squidpy only",
        "sender_definition": "SnC+ x CDKN state x cell_type_minor",
        "receiver_definition": "SnC- and p16-p21- cells grouped by cell_type_minor",
        "caveat": "Ligrec is expression/permutation based. Use sample-level consistency; cells are not independent replicates."
    }
    (OUTDIR / "run_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
