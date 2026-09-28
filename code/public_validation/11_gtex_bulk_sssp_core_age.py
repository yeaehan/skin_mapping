#!/usr/bin/env python
"""Test the 21-gene SSSP core composite across GTEx skin age bins."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results" / "260727" / "gtex_bulk_skin_qc"
OUT = ROOT / "results" / "public_validation" / "gtex_bulk_sssp_core_age"
GENES = ["CXCL16", "VEGFA", "IL18", "IL1B", "CSF1", "IL4", "CCL19", "CXCL1", "CCL28", "CCL26", "CXCL11", "CCL11", "EGF", "CXCL5", "IL6", "CXCL6", "IL5", "MMP3", "CCL7", "CCL15", "CCL16"]
AGE_ORDER = ["20-29", "30-39", "40-49", "50-59", "60-69", "70-79"]
FILES = {
    "Skin - Sun Exposed (Lower leg)": "gtex_bulk_skin_hpa_secreted_log1p_norm_lower_leg.csv.gz",
    "Skin - Not Sun Exposed (Suprapubic)": "gtex_bulk_skin_hpa_secreted_log1p_norm_suprapubic.csv.gz",
}
GCTS = {
    "Skin - Sun Exposed (Lower leg)": "gene_reads_adult_gtex_v11_skin_sun_exposed_lower_leg.gct.gz",
    "Skin - Not Sun Exposed (Suprapubic)": "gene_reads_adult_gtex_v11_skin_not_sun_exposed_suprapubic.gct.gz",
}


def add_missing_genes(expression: pd.DataFrame, genes: list[str], site: str, metadata: pd.DataFrame) -> None:
    found = []
    gct = ROOT / "data" / "external" / GCTS[site]
    for chunk in pd.read_csv(gct, sep="\t", skiprows=2, chunksize=5000):
        rows = chunk.loc[chunk["Description"].astype(str).str.upper().isin(genes)]
        if not rows.empty:
            found.append(rows)
    if not found:
        raise ValueError(f"{site}: genes absent from raw GTEx counts: {genes}")
    counts = pd.concat(found).assign(Description=lambda x: x["Description"].str.upper()).groupby("Description").sum(numeric_only=True)
    size_factors = metadata.set_index("sample_id")["deseq2_size_factor"].reindex(expression["sample_id"])
    for gene in genes:
        if gene not in counts.index:
            raise ValueError(f"{site}: {gene} absent from raw GTEx counts")
        expression[gene] = np.log1p(counts.loc[gene, expression["sample_id"]].to_numpy() / size_factors.to_numpy())


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(SOURCE / "gtex_bulk_skin_hpa_secreted_expression_metadata.csv")
    parts = []
    for site, filename in FILES.items():
        expression = pd.read_csv(SOURCE / filename, usecols=lambda c: c == "sample_id" or c in GENES)
        missing = sorted(set(GENES) - set(expression.columns))
        if missing:
            add_missing_genes(expression, missing, site, metadata)
        expression["sssp_core_mean_log1p"] = expression[GENES].mean(axis=1)
        expression["bulk_skin_site"] = site
        parts.append(expression)

    scores = pd.concat(parts, ignore_index=True).merge(
        metadata[["sample_id", "donor_id", "bulk_skin_site", "age_bin", "sex"]],
        on=["sample_id", "bulk_skin_site"], how="left", validate="one_to_one",
    )
    scores["age_midpoint"] = scores["age_bin"].str.split("-").map(
        lambda x: np.mean([float(x[0]), float(x[1])]) if isinstance(x, list) and len(x) == 2 else np.nan
    )
    scores["age_bin"] = pd.Categorical(scores["age_bin"], AGE_ORDER, ordered=True)

    summary = (scores.dropna(subset=["age_bin"])
               .groupby(["bulk_skin_site", "age_bin"], observed=True)
               .agg(n_samples=("sssp_core_mean_log1p", "size"), mean_score=("sssp_core_mean_log1p", "mean"),
                    median_score=("sssp_core_mean_log1p", "median"), sd_score=("sssp_core_mean_log1p", "std"))
               .reset_index())
    summary["sem_score"] = summary["sd_score"] / np.sqrt(summary["n_samples"])
    summary["change_from_previous_bin"] = summary.groupby("bulk_skin_site")["mean_score"].diff()
    summary["increased_from_previous_bin"] = summary["change_from_previous_bin"].gt(0).where(summary["change_from_previous_bin"].notna())

    correlations = []
    for site, frame in list(scores.groupby("bulk_skin_site")) + [("All sites (unadjusted)", scores)]:
        valid = frame[["age_midpoint", "sssp_core_mean_log1p"]].dropna()
        pearson_r, pearson_p = stats.pearsonr(valid["age_midpoint"], valid["sssp_core_mean_log1p"])
        spearman_r, spearman_p = stats.spearmanr(valid["age_midpoint"], valid["sssp_core_mean_log1p"])
        correlations.append({"bulk_skin_site": site, "n_samples": len(valid), "pearson_r": pearson_r,
                             "pearson_p": pearson_p, "spearman_rho": spearman_r, "spearman_p": spearman_p})

    scores["young_old"] = scores["age_bin"].astype("string").map(
        {"20-29": "Young", "30-39": "Young", "40-49": "Young", "60-69": "Old", "70-79": "Old"}
    )
    binary_correlations = []
    binary_scores = scores.dropna(subset=["young_old"]).copy()
    binary_scores["old_indicator"] = binary_scores["young_old"].eq("Old").astype(int)
    for site, frame in list(binary_scores.groupby("bulk_skin_site")) + [("All sites (unadjusted)", binary_scores)]:
        pearson_r, pearson_p = stats.pearsonr(frame["old_indicator"], frame["sssp_core_mean_log1p"])
        spearman_r, spearman_p = stats.spearmanr(frame["old_indicator"], frame["sssp_core_mean_log1p"])
        binary_correlations.append({"bulk_skin_site": site, "coding": "Young=0; Old=1; age 50-59 excluded",
                                    "n_young": int(frame["young_old"].eq("Young").sum()),
                                    "n_old": int(frame["young_old"].eq("Old").sum()),
                                    "mean_young": frame.loc[frame["young_old"].eq("Young"), "sssp_core_mean_log1p"].mean(),
                                    "mean_old": frame.loc[frame["young_old"].eq("Old"), "sssp_core_mean_log1p"].mean(),
                                    "pearson_r": pearson_r, "pearson_p": pearson_p,
                                    "spearman_rho": spearman_r, "spearman_p": spearman_p})

    scores.to_csv(OUT / "gtex_bulk_sssp_core_sample_scores.csv", index=False)
    summary.to_csv(OUT / "gtex_bulk_sssp_core_age_bin_summary.csv", index=False)
    correlations = pd.DataFrame(correlations)
    correlations.to_csv(OUT / "gtex_bulk_sssp_core_age_correlations.csv", index=False)
    pd.DataFrame(binary_correlations).to_csv(OUT / "gtex_bulk_sssp_core_young_old_correlations.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 5))
    for site, frame in summary.groupby("bulk_skin_site", observed=True):
        corr = correlations.loc[correlations["bulk_skin_site"].eq(site)].iloc[0]
        label = (f"{site}\nPearson r={corr.pearson_r:.3f}, p={corr.pearson_p:.2g}; "
                 f"Spearman ρ={corr.spearman_rho:.3f}, p={corr.spearman_p:.2g}")
        ax.errorbar(frame["age_bin"].astype(str), frame["mean_score"], yerr=frame["sem_score"], marker="o", capsize=3, label=label)
    ax.set(xlabel="Age bin", ylabel="Mean log1p normalized expression", title="GTEx bulk-skin 21-gene SSSP core composite")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "gtex_bulk_sssp_core_age_bins.png", dpi=300)
    plt.close(fig)
    assert len(scores) == 1405 and scores[GENES].notna().all().all()


if __name__ == "__main__":
    main()
