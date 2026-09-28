"""CellChat-like communication screen for the shared Xenium SSSP.

Cell-cell communication analysis for the shared Xenium SSSP panel.
Scores are computed per biological sample as:
ligand mean log1p(CP10K) * ligand fraction positive *
receptor mean log1p(CP10K) * receptor fraction positive.
This is expression-inferred communication, not a spatial-neighbor test.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parents[3]
H5AD = ROOT / "data" / "Xenium" / "xenium_full.h5ad"
SSSP_JSON = ROOT / "code" / "shared" / "sssp.json"
LR_CACHE = ROOT / "results" / "260727" / "secretory_ligand_receptor_communication" / "omnipath_ligrecextra_human_lr.tsv"
OUTDIR = ROOT / "results" / "260727" / "xenium_shared_sssp_cell_cell_communication"
GROUPINGS = ["cell_type_minor", "snc_cell_type_minor"]
MIN_GROUP_CELLS = 20
MIN_LIGAND_PCT = 0.05
MIN_RECEPTOR_PCT = 0.05
EDGE = ["grouping", "sender_group", "receiver_group", "ligand", "receptor"]


def sem(x):
    x = pd.to_numeric(x, errors="coerce").dropna()
    return float(x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else np.nan


def bh_fdr(values):
    p = pd.to_numeric(values, errors="coerce").to_numpy(float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if valid.any():
        idx = np.flatnonzero(valid)
        order = idx[np.argsort(p[valid])]
        q = p[order] * len(order) / np.arange(1, len(order) + 1)
        out[order] = np.clip(np.minimum.accumulate(q[::-1])[::-1], 0, 1)
    return out


def age_from_sample(sample):
    sample = str(sample).upper()
    return "Young" if sample.endswith("-Y") else "Old" if sample.endswith("-O") else "Unknown"


def load_data():
    adata = sc.read_h5ad(H5AD)
    adata.var_names = adata.var_names.astype(str).str.upper()
    if "sample" not in adata.obs:
        source = "sample_csv" if "sample_csv" in adata.obs else "sample_concat"
        adata.obs["sample"] = adata.obs[source].astype(str)
    adata.obs["sample"] = adata.obs["sample"].astype(str)
    age = adata.obs["age"].astype(str) if "age" in adata.obs else pd.Series("", index=adata.obs_names)
    adata.obs["age"] = [a if a in {"Young", "Old"} else age_from_sample(s) for a, s in zip(age, adata.obs["sample"])]
    minor = adata.obs["cell_type_minor"].astype(str)
    if "cell_type_final" in adata.obs:
        minor = minor.where(~minor.isin(["nan", "none", "Unknown"]), adata.obs["cell_type_final"].astype(str))
    adata.obs["cell_type_minor"] = minor
    adata.obs["snc"] = adata.obs["snc"].astype(str) if "snc" in adata.obs else "Unknown"
    adata.obs["snc_cell_type_minor"] = adata.obs["snc"] + " | " + minor
    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    return adata


def load_pairs(sssp, available):
    lr = pd.read_csv(LR_CACHE, sep="\t").rename(columns={"source_genesymbol": "ligand", "target_genesymbol": "receptor"})
    lr["ligand"] = lr["ligand"].astype(str).str.upper()
    lr["receptor"] = lr["receptor"].astype(str).str.upper()
    keep = lr["ligand"].isin(sssp) & lr["ligand"].isin(available) & lr["receptor"].isin(available) & lr["ligand"].ne(lr["receptor"])
    annotation = [c for c in ["is_directed", "is_stimulation", "is_inhibition", "consensus_direction"] if c in lr]
    return lr.loc[keep, ["ligand", "receptor", *annotation]].drop_duplicates(["ligand", "receptor"])


def group_expression(adata, genes):
    x = adata[:, genes].X
    x = x.toarray() if sparse.issparse(x) else np.asarray(x)
    base = pd.DataFrame(x, columns=genes)
    base["sample"] = adata.obs["sample"].to_numpy()
    base["age"] = adata.obs["age"].to_numpy()
    rows, audit = [], []
    for grouping in GROUPINGS:
        base["group"] = adata.obs[grouping].astype(str).to_numpy()
        for (sample, age, group), sub in base.groupby(["sample", "age", "group"], observed=True):
            included = len(sub) >= MIN_GROUP_CELLS and group not in {"nan", "Unknown"}
            audit.append({"sample": sample, "age": age, "grouping": grouping, "group": group, "n_cells": len(sub), "included": included})
            if not included:
                continue
            values = sub[genes]
            means, pct = values.mean(), values.gt(0).mean()
            rows.extend({"sample": sample, "age": age, "grouping": grouping, "group": group, "gene": gene,
                         "n_cells": len(sub), "mean_log1p_cp10k": means[gene], "pct_expressing": pct[gene]} for gene in genes)
    return pd.DataFrame(rows), pd.DataFrame(audit)


def score_edges(expr, pairs):
    rows = []
    for (sample, age, grouping), sub in expr.groupby(["sample", "age", "grouping"], observed=True):
        lig = sub[sub["gene"].isin(pairs["ligand"]) & sub["pct_expressing"].ge(MIN_LIGAND_PCT)][
            ["group", "gene", "n_cells", "mean_log1p_cp10k", "pct_expressing"]].rename(
            columns={"group": "sender_group", "gene": "ligand", "n_cells": "sender_n_cells", "mean_log1p_cp10k": "ligand_mean", "pct_expressing": "ligand_pct"})
        rec = sub[sub["gene"].isin(pairs["receptor"]) & sub["pct_expressing"].ge(MIN_RECEPTOR_PCT)][
            ["group", "gene", "n_cells", "mean_log1p_cp10k", "pct_expressing"]].rename(
            columns={"group": "receiver_group", "gene": "receptor", "n_cells": "receiver_n_cells", "mean_log1p_cp10k": "receptor_mean", "pct_expressing": "receptor_pct"})
        if lig.empty or rec.empty:
            continue
        scored = lig.merge(pairs, on="ligand").merge(rec, on="receptor")
        scored["communication_score"] = scored["ligand_mean"] * scored["ligand_pct"] * scored["receptor_mean"] * scored["receptor_pct"]
        scored = scored[scored["communication_score"].gt(0)].copy()
        scored.insert(0, "sample", sample); scored.insert(1, "age", age); scored.insert(2, "grouping", grouping)
        rows.append(scored)
    return pd.concat(rows, ignore_index=True)


def summarize(scores):
    by_age = scores.groupby([*EDGE, "age"], observed=True).agg(
        n_samples=("sample", "nunique"), mean_score=("communication_score", "mean"),
        median_score=("communication_score", "median"), sem_score=("communication_score", sem)).reset_index()
    delta = by_age.pivot(index=EDGE, columns="age", values="mean_score").reset_index().rename_axis(None, axis=1)
    if {"Young", "Old"}.issubset(delta):
        delta["old_minus_young"] = delta["Old"] - delta["Young"]
        delta["old_over_young"] = np.divide(delta["Old"], delta["Young"], out=np.full(len(delta), np.nan), where=delta["Young"].to_numpy() > 0)
    tests = []
    for keys, sub in scores.groupby(EDGE, observed=True):
        sample = sub.groupby(["sample", "age"], observed=True)["communication_score"].mean().reset_index()
        old = sample.loc[sample["age"].eq("Old"), "communication_score"]
        young = sample.loc[sample["age"].eq("Young"), "communication_score"]
        if len(old) >= 2 and len(young) >= 2:
            tests.append({**dict(zip(EDGE, keys)), "n_old_samples": len(old), "n_young_samples": len(young),
                          "old_mean": old.mean(), "young_mean": young.mean(), "old_minus_young": old.mean() - young.mean(),
                          "welch_p": stats.ttest_ind(old, young, equal_var=False).pvalue,
                          "mann_whitney_p": stats.mannwhitneyu(old, young, alternative="two-sided").pvalue})
    tests = pd.DataFrame(tests)
    if not tests.empty:
        tests["welch_fdr"] = bh_fdr(tests["welch_p"])
        tests["mann_whitney_fdr"] = bh_fdr(tests["mann_whitney_p"])
        tests = tests.sort_values(["welch_fdr", "old_minus_young"], ascending=[True, False])
    return by_age, delta, tests


def network_summary(scores):
    return scores.groupby(["grouping", "sample", "age", "sender_group", "receiver_group"], observed=True).agg(
        total_score=("communication_score", "sum"), mean_score=("communication_score", "mean"),
        n_lr_edges=("communication_score", "size"), top_ligands=("ligand", lambda x: ";".join(x.value_counts().head(5).index))).reset_index()


def plot_network(network):
    for grouping in GROUPINGS:
        x = network[network["grouping"].eq(grouping)].groupby(
            ["age", "sender_group", "receiver_group"], observed=True)["total_score"].mean().reset_index()
        senders = x.groupby("sender_group")["total_score"].sum().nlargest(12).index
        receivers = x.groupby("receiver_group")["total_score"].sum().nlargest(12).index
        fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
        for ax, age in zip(axes, ["Young", "Old"]):
            p = x[x["age"].eq(age)].pivot(index="sender_group", columns="receiver_group", values="total_score")
            p = p.reindex(index=senders, columns=receivers).fillna(0)
            sns.heatmap(np.log1p(p), cmap="magma", ax=ax, cbar_kws={"label": "log1p mean total score"})
            ax.set_title(age); ax.set_xlabel("Receiver"); ax.set_ylabel("SSSP sender" if age == "Young" else "")
        fig.suptitle(f"Shared SSSP communication: {grouping}", fontweight="bold")
        fig.tight_layout()
        fig.savefig(OUTDIR / f"{grouping}_young_old_network_heatmap.png", dpi=300, bbox_inches="tight")
        fig.savefig(OUTDIR / f"{grouping}_young_old_network_heatmap.pdf", bbox_inches="tight")
        plt.close(fig)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sssp = [str(x).upper() for x in json.loads(SSSP_JSON.read_text())["sssp_xenium"]]
    adata = load_data()
    pairs = load_pairs(sssp, set(adata.var_names))
    if pairs.empty:
        raise RuntimeError("No measured OmniPath pairs for sssp_xenium ligands")
    pairs.to_csv(OUTDIR / "sssp_xenium_omnipath_lr_pairs_used.csv", index=False)
    pd.DataFrame({"gene": sssp, "is_lr_ligand": [g in set(pairs["ligand"]) for g in sssp]}).to_csv(OUTDIR / "sssp_ligand_availability.csv", index=False)
    genes = sorted(set(pairs["ligand"]).union(pairs["receptor"]))
    expr, audit = group_expression(adata, genes)
    expr.to_csv(OUTDIR / "group_gene_expression_by_sample.csv.gz", index=False)
    audit.to_csv(OUTDIR / "group_cell_count_audit.csv", index=False)
    scores = score_edges(expr, pairs)
    scores.to_csv(OUTDIR / "communication_scores_by_sample.csv.gz", index=False)
    by_age, delta, tests = summarize(scores)
    by_age.to_csv(OUTDIR / "communication_edges_by_age.csv.gz", index=False)
    delta.to_csv(OUTDIR / "communication_edges_old_minus_young.csv.gz", index=False)
    tests.to_csv(OUTDIR / "communication_edges_sample_level_tests.csv.gz", index=False)
    network = network_summary(scores)
    network.to_csv(OUTDIR / "sender_receiver_network_by_sample.csv.gz", index=False)
    plot_network(network)
    summary = {"input": str(H5AD), "method": "expression product weighted by fraction positive", "spatial_test": False,
               "n_cells": adata.n_obs, "n_sssp_genes": len(sssp), "n_sssp_ligands": pairs["ligand"].nunique(),
               "n_lr_pairs": len(pairs), "n_sample_edges": len(scores), "min_group_cells": MIN_GROUP_CELLS,
               "min_ligand_pct": MIN_LIGAND_PCT, "min_receptor_pct": MIN_RECEPTOR_PCT,
               "inference_unit": "sample"}
    (OUTDIR / "run_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("\nTop age comparisons:")
    print(tests.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
