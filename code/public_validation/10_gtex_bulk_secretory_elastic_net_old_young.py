"""
GTEx bulk skin HPA-secreted-gene elastic-net screen: old vs young.

Preprocessing input is produced by 09_gtex_bulk_skin_qc.Rmd:
  - DESeq2-style size-factor QC
  - HPA secreted-gene log1p-normalized expression matrices
  - sample metadata

Age contrast:
  Young = GTEx AGE 20-29 or 30-39
  Old   = GTEx AGE 60-69 or 70-79
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler



PROJECT_ROOT = Path(__file__).resolve().parents[2]
PREP_DIR = PROJECT_ROOT / "results" / "260727" / "gtex_bulk_skin_qc"
OUTDIR = PROJECT_ROOT / "results" / "260727" / "gtex_bulk_hpa_secreted_elastic_net_old_young_from_16_preprocessed"
HPA_SECRETED_GENES = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "sinai_age_decade_trends_all_sites"
    / "hpa_predicted_secreted_gene_list.csv"
)

EXPR_FILES = {
    "lower_leg": PREP_DIR / "gtex_bulk_skin_hpa_secreted_log1p_norm_lower_leg.csv.gz",
    "suprapubic": PREP_DIR / "gtex_bulk_skin_hpa_secreted_log1p_norm_suprapubic.csv.gz",
}
METADATA_FILE = PREP_DIR / "gtex_bulk_skin_hpa_secreted_expression_metadata.csv"

YOUNG_BINS = {"20-29", "30-39"}
OLD_BINS = {"60-69", "70-79"}
AGE_MIDPOINTS = {
    "20-29": 24.5,
    "30-39": 34.5,
    "40-49": 44.5,
    "50-59": 54.5,
    "60-69": 64.5,
    "70-79": 74.5,
}
PRIMARY_TREND_FILTER_ANALYSIS = "suprapubic_plus_lower_leg"
SPEARMAN_RHO_CUTOFF = 0.1
SPEARMAN_FDR_CUTOFF = 0.05
RANDOM_SEEDS = [260727, 260728]
N_SPLITS = 4
C_GRID = [0.01, 0.1, 1, 10, 100]
L1_RATIO_GRID = [0.1, 0.5, 0.9, 1.0]
MIN_SELECTION_FREQUENCY = 0.75
MIN_DIRECTION_FREQUENCY = 0.75


@dataclass(frozen=True)
class AnalysisData:
    expr: pd.DataFrame
    meta: pd.DataFrame


def load_hpa_annotations() -> pd.DataFrame:
    hpa = pd.read_csv(HPA_SECRETED_GENES)
    hpa["Gene"] = hpa["Gene"].astype(str).str.upper()
    return hpa.drop_duplicates("Gene").rename(
        columns={
            "Gene": "gene",
            "Ensembl": "hpa_ensembl",
            "Gene description": "hpa_gene_description",
            "Protein class": "hpa_protein_class",
            "Secretome location": "hpa_secretome_location",
            "Secretome function": "hpa_secretome_function",
        }
    )


def load_preprocessed_inputs() -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    missing = [path for path in [METADATA_FILE, *EXPR_FILES.values()] if not path.exists()]
    if missing:
        missing_text = "\n".join(str(path) for path in missing)
        raise FileNotFoundError(
            "Missing preprocessed GTEx files. Render/run 09_gtex_bulk_skin_qc.Rmd first:\n"
            f"{missing_text}"
        )

    metadata = pd.read_csv(METADATA_FILE)
    metadata["age_group"] = np.select(
        [metadata["age_bin"].isin(YOUNG_BINS), metadata["age_bin"].isin(OLD_BINS)],
        ["Young", "Old"],
        default="Middle_excluded",
    )

    expr = {}
    for site_key, path in EXPR_FILES.items():
        df = pd.read_csv(path, index_col=0)
        df.columns = df.columns.astype(str).str.upper()
        expr[site_key] = df
    return expr, metadata


def build_analysis(expr_tables: dict[str, pd.DataFrame], metadata: pd.DataFrame, site_keys: list[str]) -> AnalysisData:
    common_genes = sorted(set.intersection(*(set(expr_tables[key].columns) for key in site_keys)))
    expr = pd.concat([expr_tables[key].loc[:, common_genes] for key in site_keys], axis=0)
    meta = (
        pd.DataFrame({"sample_id": expr.index})
        .merge(metadata, on="sample_id", how="left")
        .set_index("sample_id")
        .loc[expr.index]
    )

    keep = meta["age_group"].isin(["Young", "Old"])
    expr = expr.loc[keep]
    meta = meta.loc[keep]
    expr = expr.loc[:, expr.std(axis=0) > 1e-8]
    return AnalysisData(expr=expr, meta=meta)


def build_age_trend_analysis(
    expr_tables: dict[str, pd.DataFrame],
    metadata: pd.DataFrame,
    site_keys: list[str],
) -> AnalysisData:
    common_genes = sorted(set.intersection(*(set(expr_tables[key].columns) for key in site_keys)))
    expr = pd.concat([expr_tables[key].loc[:, common_genes] for key in site_keys], axis=0)
    meta = (
        pd.DataFrame({"sample_id": expr.index})
        .merge(metadata, on="sample_id", how="left")
        .set_index("sample_id")
        .loc[expr.index]
    )

    keep = meta["age_bin"].isin(AGE_MIDPOINTS)
    expr = expr.loc[keep]
    meta = meta.loc[keep].copy()
    meta["age_midpoint"] = meta["age_bin"].map(AGE_MIDPOINTS).astype(float)
    expr = expr.loc[:, expr.std(axis=0) > 1e-8]
    return AnalysisData(expr=expr, meta=meta)


def bh_fdr(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return pd.Series(out, index=pvalues.index)

    valid_idx = np.where(valid)[0]
    order = valid_idx[np.argsort(p[valid])]
    adjusted = p[order] * len(order) / np.arange(1, len(order) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out[order] = np.clip(adjusted, 0, 1)
    return pd.Series(out, index=pvalues.index)


def calculate_spearman_age_trends(data: AnalysisData, analysis: str) -> pd.DataFrame:
    age = data.meta["age_midpoint"].to_numpy(dtype=float)
    age_bins = sorted(data.meta["age_bin"].dropna().unique(), key=lambda x: AGE_MIDPOINTS.get(str(x), np.inf))
    bin_means = data.expr.groupby(data.meta["age_bin"]).mean().reindex(age_bins)

    rows = []
    for gene in data.expr.columns:
        values = data.expr[gene].to_numpy(dtype=float)
        if np.nanstd(values) == 0:
            rho = np.nan
            pvalue = np.nan
        else:
            result = stats.spearmanr(age, values, nan_policy="omit")
            rho = float(result.statistic)
            pvalue = float(result.pvalue)

        row = {
            "analysis": analysis,
            "gene": gene,
            "n_trend_samples": len(data.meta),
            "n_trend_age_bins": len(age_bins),
            "trend_age_bins": ";".join(age_bins),
            "spearman_rho_age": rho,
            "spearman_p_age": pvalue,
            "first_age_bin_mean": float(bin_means.loc[age_bins[0], gene]),
            "last_age_bin_mean": float(bin_means.loc[age_bins[-1], gene]),
            "last_minus_first_age_bin_mean": float(bin_means.loc[age_bins[-1], gene] - bin_means.loc[age_bins[0], gene]),
        }
        for age_bin in age_bins:
            row[f"mean_{age_bin}"] = float(bin_means.loc[age_bin, gene])
        rows.append(row)

    trends = pd.DataFrame(rows)
    trends["spearman_fdr_age"] = bh_fdr(trends["spearman_p_age"])
    trends["positive_significant_age_trend"] = (
        trends["spearman_rho_age"].gt(SPEARMAN_RHO_CUTOFF)
        & trends["spearman_fdr_age"].lt(SPEARMAN_FDR_CUTOFF)
    )
    trends["spearman_filter_rule"] = (
        f"spearman_rho_age > {SPEARMAN_RHO_CUTOFF} and "
        f"spearman_fdr_age < {SPEARMAN_FDR_CUTOFF}"
    )
    return trends.sort_values(["spearman_rho_age", "spearman_fdr_age"], ascending=[False, True])


def make_model(seed: int = 260727):
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(
            penalty="elasticnet",
            solver="saga",
            class_weight="balanced",
            max_iter=5000,
            random_state=seed,
        ),
    )


def run_repeated_cv(data: AnalysisData, analysis: str):
    X = data.expr.to_numpy(dtype=float)
    y = (data.meta["age_group"] == "Old").astype(int).to_numpy()
    groups = data.meta["donor_id"].astype(str).to_numpy()
    genes = np.asarray(data.expr.columns)

    prediction_rows = []
    perf_rows = []
    coef_rows = []

    for seed in RANDOM_SEEDS:
        outer_cv = StratifiedGroupKFold(
            n_splits=N_SPLITS,
            shuffle=True,
            random_state=seed,
        )

        for fold, (train_idx, test_idx) in enumerate(
            outer_cv.split(X, y, groups),
            start=1,
        ):
            inner_cv = StratifiedGroupKFold(
                n_splits=N_SPLITS,
                shuffle=True,
                random_state=seed + 100,
            )

            grid = GridSearchCV(
                make_model(seed),
                {
                    "logisticregression__C": C_GRID,
                    "logisticregression__l1_ratio": L1_RATIO_GRID,
                },
                scoring="roc_auc",
                cv=inner_cv,
                n_jobs=-1,
            )

            grid.fit(
                X[train_idx],
                y[train_idx],
                groups=groups[train_idx],
            )

            model = grid.best_estimator_

            prob = model.predict_proba(X[test_idx])[:, 1]
            pred = model.predict(X[test_idx])
            
            test_meta = data.meta.iloc[test_idx]

            for sample_id, y_true, y_prob, y_pred in zip(
                test_meta.index,
                y[test_idx],
                prob,
                pred,
            ):
                prediction_rows.append({
                    "analysis": analysis,
                    "seed": seed,
                    "fold": fold,
                    "sample_id": sample_id,
                    "donor_id": test_meta.loc[sample_id, "donor_id"],
                    "bulk_skin_site": test_meta.loc[sample_id, "bulk_skin_site"],
                    "age_group": test_meta.loc[sample_id, "age_group"],
                    "y_true": int(y_true),
                    "prob_old": float(y_prob),
                    "predicted_class": "Old" if y_pred == 1 else "Young",
                })

            coef = model.named_steps["logisticregression"].coef_.ravel()
            selected = np.abs(coef) > 1e-9

            perf_rows.append({
                "analysis": analysis,
                "seed": seed,
                "fold": fold,
                "auc": roc_auc_score(y[test_idx], prob),
                "accuracy": accuracy_score(y[test_idx], pred),
                "balanced_accuracy": balanced_accuracy_score(
                    y[test_idx], pred
                ),
                "n_selected": int(selected.sum()),
                "C": grid.best_params_["logisticregression__C"],
                "l1_ratio": grid.best_params_[
                    "logisticregression__l1_ratio"
                ],
            })

            coef_rows.extend(
                {
                    "analysis": analysis,
                    "seed": seed,
                    "fold": fold,
                    "gene": gene,
                    "coef": gene_coef,
                    "abs_coef": abs(gene_coef),
                }
                for gene, gene_coef
                in zip(genes[selected], coef[selected])
            )

    return (
        pd.DataFrame(perf_rows),
        pd.DataFrame(coef_rows),
        pd.DataFrame(prediction_rows),
    )

def fit_final_model(data: AnalysisData, analysis: str):
    X = data.expr.to_numpy(dtype=float)
    y = (data.meta["age_group"] == "Old").astype(int).to_numpy()
    groups = data.meta["donor_id"].astype(str).to_numpy()
    genes = np.asarray(data.expr.columns)

    cv = StratifiedGroupKFold(
        n_splits=N_SPLITS,
        shuffle=True,
        random_state=260727,
    )

    grid = GridSearchCV(
        make_model(),
        {
            "logisticregression__C": C_GRID,
            "logisticregression__l1_ratio": L1_RATIO_GRID,
        },
        scoring="roc_auc",
        cv=cv,
        n_jobs=-1,
    )

    grid.fit(X, y, groups=groups)

    coef = grid.best_estimator_.named_steps[
        "logisticregression"
    ].coef_.ravel()

    final = pd.DataFrame({
        "analysis": analysis,
        "gene": genes,
        "final_coef": coef,
        "final_abs_coef": np.abs(coef),
        "selected_in_final_model": np.abs(coef) > 1e-9,
    })

    final["C"] = grid.best_params_["logisticregression__C"]
    final["l1_ratio"] = grid.best_params_["logisticregression__l1_ratio"]

    return final




def summarize_genes(data, analysis, cv_coefs, final):
    n_models = len(RANDOM_SEEDS) * N_SPLITS

    stability = (
        cv_coefs.groupby(["analysis", "gene"], as_index=False)
        .agg(
            n_selected_cv=("coef", "size"),
            n_positive_cv=("coef", lambda x: (x > 0).sum()),
            mean_cv_coef=("coef", "mean"),
            median_cv_coef=("coef", "median"),
            mean_cv_abs_coef=("abs_coef", "mean"),
        )
    )

    stability["selection_frequency"] = (
        stability["n_selected_cv"] / n_models
    )

    stability["positive_frequency"] = (
        stability["n_positive_cv"] / n_models
    )

    marker = (
        final
        .merge(stability, on=["analysis", "gene"], how="left")
        .fillna({
            "n_selected_cv": 0,
            "n_positive_cv": 0,
            "selection_frequency": 0,
            "positive_frequency": 0,
        })
        .sort_values(
            ["positive_frequency", "final_abs_coef"],
            ascending=False,
        )
    )

    marker["n_negative_cv"] = marker["n_selected_cv"] - marker["n_positive_cv"]
    marker["negative_frequency"] = marker["n_negative_cv"] / n_models
    marker["stable_selected"] = (
        marker["selected_in_final_model"]
        & (marker["selection_frequency"] >= MIN_SELECTION_FREQUENCY)
    )
    marker["stable_old_predictive"] = (
        marker["stable_selected"]
        & (marker["positive_frequency"] >= MIN_DIRECTION_FREQUENCY)
        & (marker["final_coef"] > 0)
    )
    marker["stable_young_predictive"] = (
        marker["stable_selected"]
        & (marker["negative_frequency"] >= MIN_DIRECTION_FREQUENCY)
        & (marker["final_coef"] < 0)
    )
    marker["stable_direction"] = np.select(
        [marker["stable_old_predictive"], marker["stable_young_predictive"]],
        ["Old", "Young"],
        default="Not_stable",
    )
    marker["stability_rule"] = (
        f"selected_in_final_model and selection_frequency >= {MIN_SELECTION_FREQUENCY}; "
        f"Old if positive_frequency >= {MIN_DIRECTION_FREQUENCY}; "
        f"Young if negative_frequency >= {MIN_DIRECTION_FREQUENCY}"
    )
    return marker


def write_plots(performance: pd.DataFrame, markers: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.boxplot(data=performance, x="analysis", y="auc", ax=ax, color="#BFD7EA")
    sns.stripplot(data=performance, x="analysis", y="auc", ax=ax, color="#222222", size=4, alpha=0.55)
    ax.axhline(0.5, color="0.4", linestyle="--", linewidth=1)
    ax.set_xlabel("")
    ax.set_ylabel("Held-out donor-grouped CV AUC")
    ax.set_title("HPA-secreted genes classify old vs young GTEx skin")
    ax.tick_params(axis="x", rotation=20)
    plt.tight_layout()
    fig.savefig(OUTDIR / "elastic_net_old_young_cv_auc.png", dpi=300)
    fig.savefig(OUTDIR / "elastic_net_old_young_cv_auc.pdf")
    plt.close(fig)

    top = (
        markers
        .loc[markers["stable_selected"]]
        .sort_values(["analysis", "selection_frequency", "final_abs_coef"], ascending=[True, False, False])
        .groupby("analysis", group_keys=False)
        .head(20)
    )
    if top.empty:
        return
    fig, axes = plt.subplots(1, top["analysis"].nunique(), figsize=(12, 8), squeeze=False)
    for ax, (analysis, ddf) in zip(axes.ravel(), top.groupby("analysis")):
        ddf = ddf.sort_values("selection_frequency")
        colors = np.where(ddf["final_coef"] > 0, "#B2182B", "#2166AC")
        ax.barh(ddf["gene"], ddf["selection_frequency"], color=colors)
        ax.set_xlim(0, 1)
        ax.set_xlabel("CV selection frequency")
        ax.set_title(analysis)
    axes[0][0].set_ylabel("Gene")
    plt.tight_layout()
    fig.savefig(OUTDIR / "elastic_net_top_secretory_genes_stability.png", dpi=300)
    fig.savefig(OUTDIR / "elastic_net_top_secretory_genes_stability.pdf")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    expr_tables, metadata = load_preprocessed_inputs()
    hpa = load_hpa_annotations()
    analysis_specs = {
        "suprapubic_only": ["suprapubic"],
        "suprapubic_plus_lower_leg": ["suprapubic", "lower_leg"],
    }
    predictions = []

    summaries, performances, cv_coefs, finals, markers, trends = [], [], [], [], [], []
    for analysis, site_keys in analysis_specs.items():
        data = build_analysis(expr_tables, metadata, site_keys)
        perf, coef, pred = run_repeated_cv(data, analysis)
        predictions.append(pred)
        final = fit_final_model(data, analysis)
        marker = summarize_genes(data, analysis, coef, final).merge(hpa, on="gene", how="left")
        trend_data = build_age_trend_analysis(expr_tables, metadata, site_keys)
        trend = calculate_spearman_age_trends(trend_data, analysis)

        summaries.append(
            {
                "analysis": analysis,
                "n_samples": len(data.meta),
                "n_young": int(data.meta["age_group"].eq("Young").sum()),
                "n_old": int(data.meta["age_group"].eq("Old").sum()),
                "n_donors": int(data.meta["donor_id"].nunique()),
                "n_sites": int(data.meta["bulk_skin_site"].nunique()),
                "sites": ";".join(sorted(data.meta["bulk_skin_site"].dropna().unique())),
                "young_bins": ";".join(sorted(YOUNG_BINS)),
                "old_bins": ";".join(sorted(OLD_BINS)),
                "feature_universe": "HPA predicted secreted proteins",
                "n_hpa_secreted_genes_present": data.expr.shape[1],
            }
        )
        performances.append(perf)
        cv_coefs.append(coef)
        finals.append(final)
        markers.append(marker)
        trends.append(trend)

    dataset_summary = pd.DataFrame(summaries)
    performance_all = pd.concat(performances, ignore_index=True)
    cv_coefs_all = pd.concat(cv_coefs, ignore_index=True)
    final_all = pd.concat(finals, ignore_index=True)
    marker_all = pd.concat(markers, ignore_index=True)
    trend_all = pd.concat(trends, ignore_index=True)
    marker_with_trends = marker_all.merge(
        trend_all,
        on=["analysis", "gene"],
        how="left",
    )

    dataset_summary.to_csv(OUTDIR / "elastic_net_dataset_summary.csv", index=False)
    performance_all.to_csv(OUTDIR / "elastic_net_cv_performance.csv", index=False)
    cv_coefs_all.to_csv(OUTDIR / "elastic_net_cv_selected_coefficients.csv", index=False)
    final_all.to_csv(OUTDIR / "elastic_net_final_model_coefficients.csv", index=False)
    marker_all.to_csv(OUTDIR / "elastic_net_candidate_secretory_markers.csv", index=False)
    trend_all.to_csv(OUTDIR / "hpa_secreted_spearman_age_trends.csv", index=False)
    marker_with_trends.to_csv(OUTDIR / "elastic_net_candidate_secretory_markers_with_spearman_age_trends.csv", index=False)
    stable_selected = marker_with_trends.loc[marker_with_trends["stable_selected"]].sort_values(
        ["analysis", "stable_direction", "selection_frequency", "final_abs_coef"],
        ascending=[True, True, False, False],
    )
    stable_old = marker_with_trends.loc[marker_with_trends["stable_old_predictive"]].sort_values(
        ["analysis", "selection_frequency", "final_abs_coef"],
        ascending=[True, False, False],
    )
    stable_young = marker_with_trends.loc[marker_with_trends["stable_young_predictive"]].sort_values(
        ["analysis", "selection_frequency", "final_abs_coef"],
        ascending=[True, False, False],
    )
    final_selected_positive_spearman = marker_with_trends.loc[
        marker_with_trends["selected_in_final_model"]
        & marker_with_trends["final_coef"].gt(0)
        & marker_with_trends["positive_significant_age_trend"].fillna(False)
    ].sort_values(
        ["analysis", "selection_frequency", "final_abs_coef", "spearman_rho_age"],
        ascending=[True, False, False, False],
    )
    stable_old_positive_spearman = stable_old.loc[
        stable_old["positive_significant_age_trend"].fillna(False)
    ].sort_values(
        ["analysis", "selection_frequency", "final_abs_coef", "spearman_rho_age"],
        ascending=[True, False, False, False],
    )
    primary_panel = stable_old_positive_spearman.loc[
        stable_old_positive_spearman["analysis"].eq(PRIMARY_TREND_FILTER_ANALYSIS)
    ].copy()
    primary_panel["panel_rule"] = (
        "HPA-secreted genes selected as stable old-predictive elastic-net features "
        f"in {PRIMARY_TREND_FILTER_ANALYSIS}, then filtered to "
        f"Spearman rho > {SPEARMAN_RHO_CUTOFF} and FDR < {SPEARMAN_FDR_CUTOFF} "
        "using all available GTEx age bins."
    )
    stable_selected.to_csv(OUTDIR / "elastic_net_stable_selected_genes.csv", index=False)
    stable_old.to_csv(OUTDIR / "elastic_net_stable_old_predictive_genes.csv", index=False)
    stable_young.to_csv(OUTDIR / "elastic_net_stable_young_predictive_genes.csv", index=False)
    final_selected_positive_spearman.to_csv(
        OUTDIR / "elastic_net_final_selected_old_positive_spearman_filtered_genes.csv",
        index=False,
    )
    stable_old_positive_spearman.to_csv(
        OUTDIR / "elastic_net_stable_old_predictive_positive_spearman_filtered_genes.csv",
        index=False,
    )
    primary_panel.to_csv(
        OUTDIR / f"{PRIMARY_TREND_FILTER_ANALYSIS}_elastic_net_stable_old_predictive_positive_spearman_panel.csv",
        index=False,
    )
    (
        stable_selected
        .groupby(["analysis", "stable_direction"], as_index=False)
        .agg(n_genes=("gene", "nunique"))
        .to_csv(OUTDIR / "elastic_net_stable_selected_gene_summary.csv", index=False)
    )
    spearman_filter_summary = (
        marker_with_trends
        .assign(
            final_selected_old_positive=lambda x: x["selected_in_final_model"] & x["final_coef"].gt(0),
            stable_old_positive_spearman=lambda x: x["stable_old_predictive"] & x["positive_significant_age_trend"].fillna(False),
            final_selected_old_positive_spearman=lambda x: x["final_selected_old_positive"] & x["positive_significant_age_trend"].fillna(False),
        )
        .groupby("analysis", as_index=False)
        .agg(
            n_hpa_genes=("gene", "nunique"),
            n_final_selected_old_positive=("final_selected_old_positive", "sum"),
            n_final_selected_old_positive_spearman=("final_selected_old_positive_spearman", "sum"),
            n_stable_old_predictive=("stable_old_predictive", "sum"),
            n_stable_old_predictive_positive_spearman=("stable_old_positive_spearman", "sum"),
        )
    )
    spearman_filter_summary.to_csv(OUTDIR / "elastic_net_spearman_filter_gene_summary.csv", index=False)
    write_plots(performance_all, marker_all)

    run_summary = (
        performance_all.groupby("analysis")
        .agg(
            mean_auc=("auc", "mean"),
            sd_auc=("auc", "std"),
            mean_balanced_accuracy=("balanced_accuracy", "mean"),
            median_n_selected=("n_selected", "median"),
        )
        .reset_index()
        .merge(dataset_summary, on="analysis", how="left")
    )

    predictions_all = pd.concat(predictions, ignore_index=True)

    predictions_all.to_csv(
        OUTDIR / "elastic_net_oof_predictions.csv",
        index=False
    )
    
    run_summary = run_summary.merge(spearman_filter_summary, on="analysis", how="left")
    run_summary.to_csv(OUTDIR / "elastic_net_run_summary.csv", index=False)
    print(run_summary.to_string(index=False))


if __name__ == "__main__":
    main()
