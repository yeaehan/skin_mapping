"""SSSP ligand-receptor communication with Sinai SnC+ cells constrained as senders.

For each anatomical site, Squidpy is run on compound age-bin/SnC-state/cell-type
groups.  Results retain within-age communication and compare every SnC+ sender
edge with the matching SnC- sender edge.  The comparison is descriptive: the
Squidpy permutation p-value tests an individual edge, while the SnC+/SnC- ratio
is not itself a donor-level hypothesis test.
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
import seaborn as sns
from scipy.sparse import csc_matrix, csr_matrix

if not hasattr(anndata_views, "SparseCSRView"):
    anndata_views.SparseCSRView = csr_matrix
if not hasattr(anndata_views, "SparseCSCView"):
    anndata_views.SparseCSCView = csc_matrix
import squidpy as sq

ROOT = Path(__file__).resolve().parents[2]
H5AD = ROOT / "data" / "scrna.integrated_annotated_sinai.h5ad"
SNC_OBS = ROOT / "results" / "260727" / "sinai_age_decade_trends_all_sites" / "sinai_all_sites_snc_decade_obs.csv.gz"
SSSP_JSON = ROOT / "code" / "shared" / "sssp.json"
LR_CACHE = ROOT / "results" / "260727" / "secretory_ligand_receptor_communication" / "omnipath_ligrecextra_human_lr.tsv"
OUTDIR = ROOT / "results" / "260727" / "sinai_shared_sssp_snc_sender_ligrec_enrichment"
MIN_GROUP_CELLS = 20
N_PERMUTATIONS = 100
EXPR_THRESHOLD = 0.05
SEED = 260727
EPS = 1e-8


def clean(x):
    return x.astype("string").fillna("Unknown").astype(str).replace({"nan": "Unknown", "": "Unknown"})


def bh_fdr(x):
    p = pd.to_numeric(x, errors="coerce").to_numpy(float)
    ans = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        ids = np.flatnonzero(ok)
        order = ids[np.argsort(p[ok])]
        q = p[order] * len(order) / np.arange(1, len(order) + 1)
        ans[order] = np.clip(np.minimum.accumulate(q[::-1])[::-1], 0, 1)
    return ans


def load_inputs():
    panel = list(dict.fromkeys(str(g).upper() for g in json.loads(SSSP_JSON.read_text())["genes"]))
    full = ad.read_h5ad(H5AD, backed="r")
    lookup = {str(g).upper(): str(g) for g in full.var_names}
    lr = pd.read_csv(LR_CACHE, sep="\t")
    lr["source"] = lr["source_genesymbol"].astype(str).str.upper()
    lr["target"] = lr["target_genesymbol"].astype(str).str.upper()
    lr = (lr[lr["source"].isin(panel) & lr["source"].isin(lookup) & lr["target"].isin(lookup)
             & lr["source"].ne(lr["target"])][["source", "target"]]
          .drop_duplicates().sort_values(["source", "target"]).reset_index(drop=True))
    pd.DataFrame({"gene": panel, "present_in_sinai": [g in lookup for g in panel],
                  "used_as_ligand": [g in set(lr["source"]) for g in panel]}).to_csv(
        OUTDIR / "shared_sssp_gene_availability.csv", index=False)
    genes = sorted(set(lr["source"]).union(lr["target"]))
    x = full[:, [lookup[g] for g in genes]].layers["log1p"]
    obs = full.obs.copy()
    full.file.close()

    calls = pd.read_csv(SNC_OBS, usecols=["cell_id", "snc", "snc_bool"])
    if calls["cell_id"].duplicated().any():
        raise ValueError("Duplicate cell_id values in saved Sinai SnC calls")
    calls = calls.set_index("cell_id")
    keep = obs.index.astype(str).isin(calls.index)
    if not keep.any():
        raise ValueError("No H5AD cell IDs overlap the saved Sinai SnC calls")
    x = x[keep]
    obs = obs.loc[keep].copy()
    aligned = calls.reindex(obs.index.astype(str))
    if aligned["snc_bool"].isna().any():
        raise ValueError("Unexpected missing SnC calls after restricting to overlapping cells")
    obs["snc_state"] = np.where(aligned["snc_bool"].astype(bool).to_numpy(), "SnC+", "SnC-")
    obs["cell_type_group"] = clean(obs["cell_type.broad"])
    obs["site"] = clean(obs["anatomic_site"])
    age = pd.to_numeric(obs["donor_age"], errors="coerce")
    obs["age_bin"] = pd.cut(age, bins=list(range(20, 101, 10)), right=False,
                             labels=[f"{i}-{i+9}" for i in range(20, 100, 10)]).astype(str)
    a = ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=pd.Index(genes)))
    return a, panel, lr, sum(g in lookup for g in panel)


def audit_strata(a):
    rows = []
    for keys, d in a.obs.groupby(["site", "age_bin"], observed=True):
        site, age = keys
        plus_counts = d.loc[d["snc_state"].eq("SnC+"), "cell_type_group"].value_counts()
        minus_counts = d.loc[d["snc_state"].eq("SnC-"), "cell_type_group"].value_counts()
        rows.append({"site": site, "age_bin": age, "n_cells": len(d),
                     "n_snc_plus": int(d["snc_state"].eq("SnC+").sum()),
                     "n_donors": d["donor_id"].nunique(),
                     "n_eligible_snc_plus_types": int((plus_counts >= MIN_GROUP_CELLS).sum()),
                     "n_eligible_snc_minus_types": int((minus_counts >= MIN_GROUP_CELLS).sum()),
                     "donors": ";".join(sorted(set(d["donor_id"].astype(str))))})
    return pd.DataFrame(rows)


def run_site(a, lr, site):
    sub = a[a.obs["site"].eq(site)].copy()
    sub.obs["compound"] = (sub.obs["age_bin"].astype(str) + "|" + sub.obs["snc_state"] + "|"
                           + sub.obs["cell_type_group"])
    counts = sub.obs["compound"].value_counts()
    groups = sorted(counts[counts >= MIN_GROUP_CELLS].index)
    sub = sub[sub.obs["compound"].isin(groups)].copy()
    sub.obs["compound"] = pd.Categorical(sub.obs["compound"], categories=groups)
    if len(groups) < 2:
        return pd.DataFrame()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = sq.gr.ligrec(sub, cluster_key="compound", interactions=lr,
                              n_perms=N_PERMUTATIONS, threshold=EXPR_THRESHOLD,
                              corr_method=None, use_raw=False, copy=True,
                              seed=SEED, n_jobs=1)
    means, pvals = result["means"], result["pvalues"]
    means.index.names = pvals.index.names = ["ligand", "receptor"]
    means.columns.names = pvals.columns.names = ["source_compound", "target_compound"]
    m = means.stack(["source_compound", "target_compound"], future_stack=True).rename("mean_expression_product")
    p = pvals.stack(["source_compound", "target_compound"], future_stack=True).rename("pvalue")
    d = pd.concat([m, p], axis=1).reset_index().dropna(subset=["mean_expression_product"])
    d[["source_age", "source_snc", "source_group"]] = d["source_compound"].str.split("|", n=2, expand=True)
    d[["target_age", "target_snc", "target_group"]] = d["target_compound"].str.split("|", n=2, expand=True)
    d = d[d["source_age"].eq(d["target_age"])].copy()
    d["site"] = site
    d["age_bin"] = d["source_age"]
    d["lr_pair"] = d["ligand"] + "->" + d["receptor"]
    return d.drop(columns=["source_age", "target_age"])


def snc_sender_enrichment(all_edges):
    for col in ["mean_expression_product", "pvalue"]:
        if isinstance(all_edges[col].dtype, pd.SparseDtype):
            all_edges[col] = all_edges[col].sparse.to_dense()
        all_edges[col] = pd.to_numeric(all_edges[col], errors="coerce")
    keys = ["site", "age_bin", "ligand", "receptor", "source_group", "target_snc", "target_group"]
    cols = keys + ["mean_expression_product", "pvalue"]
    plus = all_edges[all_edges["source_snc"].eq("SnC+")][cols].rename(columns={
        "mean_expression_product": "snc_plus_mean_product", "pvalue": "snc_plus_pvalue"})
    minus = all_edges[all_edges["source_snc"].eq("SnC-")][cols].rename(columns={
        "mean_expression_product": "snc_minus_mean_product", "pvalue": "snc_minus_pvalue"})
    d = plus.merge(minus, on=keys, how="left")
    d["matched_snc_minus"] = d["snc_minus_mean_product"].notna()
    d["delta_mean_product"] = d["snc_plus_mean_product"] - d["snc_minus_mean_product"]
    d["fold_enrichment"] = (d["snc_plus_mean_product"] + EPS) / (d["snc_minus_mean_product"] + EPS)
    d["log2_enrichment"] = np.log2(d["fold_enrichment"])
    d["fdr_snc_plus_within_site_age"] = d.groupby(["site", "age_bin"], observed=True)["snc_plus_pvalue"].transform(bh_fdr)
    d["snc_plus_nominal_p05"] = d["snc_plus_pvalue"].le(0.05)
    d["snc_plus_fdr05"] = d["fdr_snc_plus_within_site_age"].le(0.05)
    d["enriched_over_snc_minus"] = d["matched_snc_minus"] & d["log2_enrichment"].gt(0)
    d["lr_pair"] = d["ligand"] + "->" + d["receptor"]
    return d


def involvement_tables(d):
    sender = d.groupby(["site", "age_bin", "source_group"], observed=True).agg(
        total_snc_plus_communication=("snc_plus_mean_product", "sum"),
        positive_enrichment_sum=("delta_mean_product", lambda x: x.clip(lower=0).sum()),
        median_log2_enrichment=("log2_enrichment", "median"),
        n_lr_edges=("lr_pair", "size"), n_unique_ligands=("ligand", "nunique"),
        n_nominal_edges=("snc_plus_nominal_p05", "sum"), n_fdr_edges=("snc_plus_fdr05", "sum"),
        top_ligands=("ligand", lambda x: ";".join(x.value_counts().head(5).index))).reset_index()
    receiver = d.groupby(["site", "age_bin", "target_group"], observed=True).agg(
        total_incoming_snc_plus_communication=("snc_plus_mean_product", "sum"),
        positive_incoming_enrichment_sum=("delta_mean_product", lambda x: x.clip(lower=0).sum()),
        n_lr_edges=("lr_pair", "size"), n_unique_ligands=("ligand", "nunique"),
        n_nominal_edges=("snc_plus_nominal_p05", "sum"), n_fdr_edges=("snc_plus_fdr05", "sum"),
        top_ligands=("ligand", lambda x: ";".join(x.value_counts().head(5).index))).reset_index()
    sender["rank_within_site_age"] = sender.groupby(["site", "age_bin"])["total_snc_plus_communication"].rank(method="dense", ascending=False)
    receiver["rank_within_site_age"] = receiver.groupby(["site", "age_bin"])["total_incoming_snc_plus_communication"].rank(method="dense", ascending=False)
    return sender, receiver


def plot_involvement(table, group_col, metric, stem):
    top = table.groupby(group_col, observed=True)[metric].sum().nlargest(15).index
    d = table[table[group_col].isin(top)].copy()
    d["site_age"] = d["site"] + " | " + d["age_bin"]
    p = d.pivot_table(index="site_age", columns=group_col, values=metric, aggfunc="sum", fill_value=0)
    if p.empty:
        return
    fig, ax = plt.subplots(figsize=(max(10, .65 * len(p.columns)), max(6, .3 * len(p))))
    sns.heatmap(np.log1p(p), cmap="magma", ax=ax)
    ax.set_title(f"SnC+ SSSP {stem} involvement (log1p communication)")
    ax.set_xlabel("Cell type"); ax.set_ylabel("Site | age bin")
    fig.tight_layout()
    fig.savefig(OUTDIR / f"snc_plus_{stem}_cell_type_involvement_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / f"snc_plus_{stem}_cell_type_involvement_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    a, panel, lr, n_present = load_inputs()
    lr.to_csv(OUTDIR / "shared_sssp_omnipath_lr_pairs_used.csv", index=False)
    audit = audit_strata(a)
    audit.to_csv(OUTDIR / "site_age_bin_snc_audit.csv", index=False)
    edge_path = OUTDIR / "sinai_all_snc_state_ligrec_edges.csv.gz"
    if edge_path.exists():
        print(f"Reusing completed Squidpy edges: {edge_path}", flush=True)
        all_edges = pd.read_csv(edge_path)
    else:
        frames = []
        for site in sorted(a.obs["site"].unique()):
            print(f"Running {site}: {int(a.obs['site'].eq(site).sum())} cells", flush=True)
            result = run_site(a, lr, site)
            if not result.empty:
                frames.append(result)
        all_edges = pd.concat(frames, ignore_index=True)
        all_edges.to_csv(edge_path, index=False)
    enriched = snc_sender_enrichment(all_edges)
    enriched.to_csv(OUTDIR / "sinai_snc_plus_sender_sssp_enrichment_edges.csv.gz", index=False)
    sender, receiver = involvement_tables(enriched)
    sender.to_csv(OUTDIR / "sinai_snc_plus_sender_cell_type_involvement_by_site_age.csv", index=False)
    receiver.to_csv(OUTDIR / "sinai_snc_plus_receiver_cell_type_involvement_by_site_age.csv", index=False)
    plot_involvement(sender, "source_group", "total_snc_plus_communication", "sender")
    plot_involvement(receiver, "target_group", "total_incoming_snc_plus_communication", "receiver")
    summary = {
        "input_cells": a.n_obs, "full_sssp_genes": len(panel), "sssp_genes_present": n_present,
        "sssp_ligands_used": int(lr.source.nunique()),
        "lr_pairs_used": len(lr), "permutations": N_PERMUTATIONS,
        "snc_definition_source": str(SNC_OBS), "snc_plus_cells": int(a.obs.snc_state.eq("SnC+").sum()),
        "n_snc_plus_sender_edges": len(enriched), "n_edges_with_matched_snc_minus": int(enriched.matched_snc_minus.sum()),
        "interpretation": "Squidpy p/FDR assess each SnC+ edge; fold enrichment is the matched SnC+/SnC- mean-expression-product ratio and is descriptive, not a direct donor-level test."
    }
    (OUTDIR / "run_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
