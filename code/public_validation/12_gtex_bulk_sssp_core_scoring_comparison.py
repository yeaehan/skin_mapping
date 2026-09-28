#!/usr/bin/env python
"""Compare mean-expression, Scanpy, and AUCell scores for GTEx SSSP core genes."""
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
QC = ROOT / "results" / "260727" / "gtex_bulk_skin_qc"
BASE = ROOT / "results" / "public_validation" / "gtex_bulk_sssp_core_age"
GENES = ["CXCL16", "VEGFA", "IL18", "IL1B", "CSF1", "IL4", "CCL19", "CXCL1", "CCL28", "CCL26", "CXCL11", "CCL11", "EGF", "CXCL5", "IL6", "CXCL6", "IL5", "MMP3", "CCL7", "CCL15", "CCL16"]
AGE_ORDER = ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"]
GCTS = {
    "Skin - Sun Exposed (Lower leg)": "gene_reads_adult_gtex_v11_skin_sun_exposed_lower_leg.gct.gz",
    "Skin - Not Sun Exposed (Suprapubic)": "gene_reads_adult_gtex_v11_skin_not_sun_exposed_suprapubic.gct.gz",
}


def aucell(matrix: np.ndarray, signature: np.ndarray, top_fraction: float = 0.05) -> np.ndarray:
    cutoff = max(1, int(matrix.shape[1] * top_fraction))
    # ponytail: ties at the cutoff are arbitrary; use pySCENIC if exact implementation parity is required.
    top = np.argpartition(-matrix, cutoff - 1, axis=1)[:, :cutoff]
    top = np.take_along_axis(top, np.argsort(-np.take_along_axis(matrix, top, axis=1), axis=1), axis=1)
    hits = np.isin(top, signature)
    n_possible = min(len(signature), cutoff)
    maximum = n_possible * cutoff - n_possible * (n_possible - 1) / 2
    return (hits * np.arange(cutoff, 0, -1)).sum(axis=1) / maximum


def score_site(site: str, filename: str, metadata: pd.DataFrame) -> pd.DataFrame:
    table = pd.read_csv(ROOT / "data" / "external" / filename, sep="\t", skiprows=2)
    sample_ids = table.columns[2:].tolist()
    counts = table.iloc[:, 2:].apply(pd.to_numeric, errors="coerce").fillna(0).astype(np.float32)
    counts.insert(0, "gene", table["Description"].astype(str).str.upper().to_numpy())
    counts = counts.loc[~counts["gene"].isin(["", "NAN", "NONE"])].groupby("gene", observed=True).sum()
    missing = sorted(set(GENES) - set(counts.index))
    if missing:
        raise ValueError(f"{site}: missing genes {missing}")

    size_factors = metadata.set_index("sample_id")["deseq2_size_factor"].reindex(sample_ids)
    normalized = counts.divide(size_factors.to_numpy(), axis=1).T
    work = ad.AnnData(normalized.to_numpy(dtype=np.float32), obs=pd.DataFrame(index=sample_ids),
                     var=pd.DataFrame(index=counts.index))
    sc.pp.log1p(work)
    sc.tl.score_genes(work, GENES, score_name="scanpy_score", random_state=260727)
    signature = counts.index.get_indexer(GENES)
    result = pd.DataFrame({"sample_id": sample_ids,
                           "scanpy_score": work.obs["scanpy_score"].to_numpy(),
                           "aucell_score": aucell(counts.T.to_numpy(), signature),
                           "bulk_skin_site": site})
    return result


def main() -> None:
    metadata = pd.read_csv(QC / "gtex_bulk_skin_hpa_secreted_expression_metadata.csv")
    baseline = pd.read_csv(BASE / "gtex_bulk_sssp_core_sample_scores.csv")
    scored = pd.concat([score_site(site, filename, metadata) for site, filename in GCTS.items()], ignore_index=True)
    scored = scored.merge(baseline[["sample_id", "sssp_core_mean_log1p", "young_old", "age_bin", "age_midpoint"]],
                          on="sample_id", how="left", validate="one_to_one")
    scored["old_indicator"] = scored["young_old"].eq("Old").astype(float).where(scored["young_old"].notna())

    methods = ["sssp_core_mean_log1p", "scanpy_score", "aucell_score"]
    rows = []
    binary = scored.dropna(subset=["old_indicator"])
    for site, frame in list(binary.groupby("bulk_skin_site")) + [("All sites (unadjusted)", binary)]:
        for method in methods:
            pearson_r, pearson_p = stats.pearsonr(frame["old_indicator"], frame[method])
            spearman_r, spearman_p = stats.spearmanr(frame["old_indicator"], frame[method])
            rows.append({"bulk_skin_site": site, "method": method,
                         "n_young": int(frame["young_old"].eq("Young").sum()),
                         "n_old": int(frame["young_old"].eq("Old").sum()),
                         "pearson_r": pearson_r, "pearson_p": pearson_p,
                         "spearman_rho": spearman_r, "spearman_p": spearman_p})

    concordance = []
    for site, frame in list(scored.groupby("bulk_skin_site")) + [("All sites", scored)]:
        for left, right in [("scanpy_score", "aucell_score"), ("sssp_core_mean_log1p", "scanpy_score"),
                            ("sssp_core_mean_log1p", "aucell_score")]:
            concordance.append({"bulk_skin_site": site, "method_1": left, "method_2": right,
                                "pearson_r": stats.pearsonr(frame[left], frame[right]).statistic,
                                "spearman_rho": stats.spearmanr(frame[left], frame[right]).statistic})

    scored["age_bin"] = pd.Categorical(scored["age_bin"], AGE_ORDER, ordered=True)
    age_rows, age_correlations = [], []
    for (site, age_bin), frame in scored.groupby(["bulk_skin_site", "age_bin"], observed=True):
        for method in methods:
            age_rows.append({"bulk_skin_site": site, "age_bin": age_bin, "method": method,
                             "n_samples": len(frame), "mean_score": frame[method].mean(),
                             "sd_score": frame[method].std(), "sem_score": frame[method].sem()})
    age_summary = pd.DataFrame(age_rows).sort_values(["method", "bulk_skin_site", "age_bin"])
    age_summary["change_from_previous_bin"] = age_summary.groupby(["method", "bulk_skin_site"])["mean_score"].diff()
    age_summary["increased_from_previous_bin"] = age_summary["change_from_previous_bin"].gt(0).where(age_summary["change_from_previous_bin"].notna())
    for site, frame in scored.groupby("bulk_skin_site"):
        for method in methods:
            age_correlations.append({"bulk_skin_site": site, "method": method,
                                     "pearson_r": stats.pearsonr(frame["age_midpoint"], frame[method]).statistic,
                                     "pearson_p": stats.pearsonr(frame["age_midpoint"], frame[method]).pvalue,
                                     "spearman_rho": stats.spearmanr(frame["age_midpoint"], frame[method]).statistic,
                                     "spearman_p": stats.spearmanr(frame["age_midpoint"], frame[method]).pvalue})
    age_correlations = pd.DataFrame(age_correlations)

    scored.to_csv(BASE / "gtex_bulk_sssp_core_scoring_method_scores.csv", index=False)
    pd.DataFrame(rows).to_csv(BASE / "gtex_bulk_sssp_core_scoring_method_young_old_correlations.csv", index=False)
    pd.DataFrame(concordance).to_csv(BASE / "gtex_bulk_sssp_core_scoring_method_concordance.csv", index=False)
    age_summary.to_csv(BASE / "gtex_bulk_sssp_core_scoring_method_age_bin_summary.csv", index=False)
    age_correlations.to_csv(BASE / "gtex_bulk_sssp_core_scoring_method_age_correlations.csv", index=False)

    titles = {"sssp_core_mean_log1p": "Mean expression", "scanpy_score": "Scanpy score", "aucell_score": "AUCell score"}
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, method in zip(axes, methods):
        for site, frame in age_summary.loc[age_summary["method"].eq(method)].groupby("bulk_skin_site"):
            corr = age_correlations.loc[(age_correlations["method"].eq(method)) & age_correlations["bulk_skin_site"].eq(site)].iloc[0]
            label = (f"{site}\nPearson r={corr.pearson_r:.3f}, p={corr.pearson_p:.2g}; "
                     f"Spearman ρ={corr.spearman_rho:.3f}, p={corr.spearman_p:.2g}")
            ax.errorbar(frame["age_bin"].astype(str), frame["mean_score"], yerr=frame["sem_score"],
                        marker="o", capsize=3, label=label)
        ax.set(title=titles[method], xlabel="Age bin", ylabel="Mean score")
        ax.tick_params(axis="x", rotation=45)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("GTEx bulk-skin SSSP core scores across consecutive age bins")
    fig.tight_layout()
    fig.savefig(BASE / "gtex_bulk_sssp_core_three_scoring_methods_age_bins.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    assert len(scored) == 1405 and scored[methods].notna().all().all()


if __name__ == "__main__":
    main()
