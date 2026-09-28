#!/usr/bin/env python
"""AUCell age trends for the SSSP core in Sinai, overall and by study/cell type."""
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numba import njit, prange
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "data" / "scrna.integrated_annotated_sinai.h5ad"
OUT = ROOT / "results" / "public_validation" / "sinai_sssp_core_aucell_age"
GENES = ["CXCL16", "VEGFA", "IL18", "IL1B", "CSF1", "IL4", "CCL19", "CXCL1", "CCL28", "CCL26", "CXCL11", "CCL11", "EGF", "CXCL5", "IL6", "CXCL6", "IL5", "MMP3", "CCL7", "CCL15", "CCL16"]


@njit(parallel=True)
def sparse_aucell(indptr, indices, values, signature_mask, cutoff, maximum):
    scores = np.zeros(len(indptr) - 1, dtype=np.float32)
    for row in prange(len(scores)):
        start, stop = indptr[row], indptr[row + 1]
        order = np.argsort(values[start:stop])[::-1]
        limit = min(cutoff, len(order))
        total = 0.0
        for rank in range(limit):
            if signature_mask[indices[start + order[rank]]]:
                total += cutoff - rank
        scores[row] = total / maximum
    return scores


def bh(pvalues: pd.Series) -> np.ndarray:
    p = pvalues.to_numpy(float)
    out = np.full(len(p), np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid])]
    adjusted = p[order] * len(order) / np.arange(1, len(order) + 1)
    out[order] = np.minimum.accumulate(adjusted[::-1])[::-1].clip(0, 1)
    return out


def correlations(frame: pd.DataFrame, groups: list[str]) -> pd.DataFrame:
    rows = []
    for keys, sub in frame.groupby(groups, observed=True):
        keys = keys if isinstance(keys, tuple) else (keys,)
        row = dict(zip(groups, keys))
        row.update(n_donors=sub["donor_id"].nunique(), n_ages=sub["donor_age"].nunique(),
                   youngest_age=sub["donor_age"].min(), oldest_age=sub["donor_age"].max())
        if row["n_donors"] >= 3 and row["n_ages"] >= 2:
            row.update(pearson_r=stats.pearsonr(sub["donor_age"], sub["aucell_score"]).statistic,
                       pearson_p=stats.pearsonr(sub["donor_age"], sub["aucell_score"]).pvalue,
                       spearman_rho=stats.spearmanr(sub["donor_age"], sub["aucell_score"]).statistic,
                       spearman_p=stats.spearmanr(sub["donor_age"], sub["aucell_score"]).pvalue)
        else:
            row.update(pearson_r=np.nan, pearson_p=np.nan, spearman_rho=np.nan, spearman_p=np.nan)
        rows.append(row)
    result = pd.DataFrame(rows)
    result["pearson_fdr"] = bh(result["pearson_p"])
    result["spearman_fdr"] = bh(result["spearman_p"])
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    atlas = ad.read_h5ad(INPUT, backed="r")
    lookup = {str(g).upper(): i for i, g in enumerate(atlas.var_names)}
    present = [g for g in GENES if g in lookup]
    missing = [g for g in GENES if g not in lookup]
    pd.DataFrame({"gene": GENES, "present": [g in lookup for g in GENES]}).to_csv(OUT / "gene_overlap.csv", index=False)

    signature_mask = np.zeros(atlas.n_vars, dtype=np.bool_)
    signature_mask[[lookup[g] for g in present]] = True
    cutoff = max(1, int(atlas.n_vars * 0.05))
    possible = min(len(present), cutoff)
    maximum = possible * cutoff - possible * (possible - 1) / 2
    scores = np.empty(atlas.n_obs, dtype=np.float32)
    layer = atlas.layers["log1p"]
    for start in range(0, atlas.n_obs, 10_000):
        stop = min(start + 10_000, atlas.n_obs)
        block = sparse.csr_matrix(layer[start:stop])
        scores[start:stop] = sparse_aucell(block.indptr, block.indices, block.data,
                                           signature_mask, cutoff, maximum)

    obs = atlas.obs[["study_id", "donor_id", "donor_age", "cell_type"]].copy()
    site_coverage = atlas.obs[["study_id", "donor_id", "donor_age", "anatomic_site", "anatomic_site.detailed"]].copy()
    atlas.file.close()
    obs["aucell_score"] = scores
    obs["donor_age"] = pd.to_numeric(obs["donor_age"], errors="coerce")
    obs = obs.dropna(subset=["study_id", "donor_id", "donor_age", "cell_type"])
    obs["age_decade"] = (obs["donor_age"] // 10 * 10).astype(int).astype(str) + "s"
    site_coverage["donor_age"] = pd.to_numeric(site_coverage["donor_age"], errors="coerce")
    site_coverage = (site_coverage.groupby(["study_id", "anatomic_site", "anatomic_site.detailed"], observed=True, dropna=False)
                     .agg(n_cells=("donor_id", "size"), n_donors=("donor_id", "nunique"),
                          minimum_age=("donor_age", "min"), maximum_age=("donor_age", "max"))
                     .reset_index())
    site_coverage.to_csv(OUT / "study_anatomic_site_coverage.csv", index=False)

    donor_study = (obs.groupby(["study_id", "donor_id", "donor_age", "age_decade", "cell_type"], observed=True)
                   .agg(aucell_score=("aucell_score", "mean"), n_cells=("aucell_score", "size")).reset_index())
    donor_overall = (donor_study.groupby(["donor_id", "donor_age", "age_decade", "cell_type"], observed=True)
                     .agg(aucell_score=("aucell_score", "mean"), n_cells=("n_cells", "sum")).reset_index())
    overall = correlations(donor_overall, ["cell_type"])
    by_study = correlations(donor_study, ["study_id", "cell_type"])
    age_bins = (donor_overall.groupby(["cell_type", "age_decade"], observed=True)
                .agg(n_donors=("donor_id", "nunique"), mean_score=("aucell_score", "mean"),
                     sem_score=("aucell_score", "sem")).reset_index())

    donor_study.to_csv(OUT / "donor_study_celltype_aucell_scores.csv", index=False)
    donor_overall.to_csv(OUT / "donor_celltype_aucell_scores_overall.csv", index=False)
    overall.to_csv(OUT / "overall_age_correlations_by_celltype.csv", index=False)
    by_study.to_csv(OUT / "study_specific_age_correlations_by_celltype.csv", index=False)
    age_bins.to_csv(OUT / "overall_age_decade_summary_by_celltype.csv", index=False)

    cell_types = sorted(age_bins["cell_type"].unique())
    fig, axes = plt.subplots(3, 4, figsize=(18, 12), squeeze=False)
    for ax, cell_type in zip(axes.flat, cell_types):
        frame = age_bins.loc[age_bins["cell_type"].eq(cell_type)].sort_values("age_decade")
        test = overall.loc[overall["cell_type"].eq(cell_type)].iloc[0]
        ax.errorbar(frame["age_decade"], frame["mean_score"], yerr=frame["sem_score"], marker="o", capsize=2)
        ax.set_title(f"{cell_type}\nρ={test.spearman_rho:.2f}, p={test.spearman_p:.2g}", fontsize=9)
        ax.tick_params(axis="x", rotation=45)
    for ax in axes.flat[len(cell_types):]:
        ax.axis("off")
    fig.suptitle(f"Sinai SSSP core AUCell by age and cell type ({len(present)}/{len(GENES)} genes)")
    fig.tight_layout()
    fig.savefig(OUT / "overall_aucell_age_bins_by_celltype.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(3, 4, figsize=(18, 12), squeeze=False)
    for ax, cell_type in zip(axes.flat, cell_types):
        frame = donor_overall.loc[donor_overall["cell_type"].eq(cell_type)]
        test = overall.loc[overall["cell_type"].eq(cell_type)].iloc[0]
        ax.scatter(frame["donor_age"], frame["aucell_score"], alpha=0.75, s=25)
        if frame["donor_age"].nunique() >= 2:
            slope, intercept = np.polyfit(frame["donor_age"], frame["aucell_score"], 1)
            x = np.array([frame["donor_age"].min(), frame["donor_age"].max()])
            ax.plot(x, slope * x + intercept, color="black", linewidth=1)
        ax.set_title(f"{cell_type}\nPearson r={test.pearson_r:.2f}, p={test.pearson_p:.2g}\nSpearman ρ={test.spearman_rho:.2f}, p={test.spearman_p:.2g}", fontsize=9)
        ax.set(xlabel="Donor age (years)", ylabel="Mean donor AUCell score")
    for ax in axes.flat[len(cell_types):]:
        ax.axis("off")
    fig.suptitle(f"Sinai SSSP core AUCell versus actual donor age ({len(present)}/{len(GENES)} genes)")
    fig.tight_layout()
    fig.savefig(OUT / "overall_aucell_actual_age_by_celltype.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    heat = by_study.pivot(index="cell_type", columns="study_id", values="spearman_rho")
    fig, ax = plt.subplots(figsize=(max(10, heat.shape[1] * 0.8), max(6, heat.shape[0] * 0.5)))
    image = ax.imshow(heat, cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(heat.shape[1]), heat.columns, rotation=45, ha="right")
    ax.set_yticks(range(heat.shape[0]), heat.index)
    ax.set_title("Study-specific Spearman age correlation of SSSP core AUCell")
    fig.colorbar(image, ax=ax, label="Spearman ρ")
    fig.tight_layout()
    fig.savefig(OUT / "study_specific_age_correlation_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    assert len(present) == 17 and len(donor_study) > 0
    print(f"Scored {len(obs):,} cells; present={present}; missing={missing}")


if __name__ == "__main__":
    main()
