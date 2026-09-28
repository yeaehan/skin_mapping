"""
Focused Xenium ligand-receptor scores by SnC status, immune subtype, and minor
cell type.

This script answers: for the selected secretory ligand-receptor pairs, which
SnC/cell-type groups express the ligand and which groups express the receptor?

Scores are computed per sample first:
  sender ligand mean log1p(CP10K) * sender ligand pct positive
  * receiver receptor mean log1p(CP10K) * receiver receptor pct positive
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[3]
XENIUM_DIR = PROJECT_ROOT / "data" / "xenium"
METADATA_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_snc_hotspot_neighborhoods"
    / "xenium_snc_hotspot_cell_neighborhoods.csv"
)
MARKER_EXPR_PATH = (
    PROJECT_ROOT
    / "results"
    / "260727"
    / "xenium_secretory_score_marker_correlations"
    / "xenium_cdkn_proliferation_marker_log_expression.csv.gz"
)
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_secretory_lr_by_snc_immune_minor"

FOCUSED_PAIRS = [
    ("ANXA1", "EGFR", "broad_skin_public_secretory", "SSSP gene as ligand"),
    ("DEFB1", "CCR6", "broad_and_49_gene_skin_secretory", "SSSP gene as ligand"),
    ("CXCL14", "CXCR4", "broad_skin_public_secretory", "SSSP gene as ligand"),
    ("ERBB3", "EGF", "broad_skin_public_secretory", "SSSP/receptor-family gene as LR source; interpret cautiously"),
    ("ERBB3", "AREG", "broad_skin_public_secretory", "SSSP/receptor-family gene as LR source; interpret cautiously"),
    ("EGFR", "EGF", "broad_skin_public_secretory", "SSSP/receptor-family gene as LR source; interpret cautiously"),
    ("EGFR", "AREG", "broad_skin_public_secretory", "SSSP/receptor-family gene as LR source; interpret cautiously"),
    ("CXCL16", "CXCR6", "49_gene_skin_secretory", "SSSP gene as ligand"),
]

GROUP_COLS = [
    "snc",
    "cdkn1a_status",
    "cdkn2a_status",
    "cdkn_any_status",
    "cell_type_minor",
    "immcelltype",
    "snc_cell_type_minor",
    "snc_immcelltype",
    "snc_immcelltype_or_minor",
    "cdkn1a_cell_type_minor",
    "cdkn2a_cell_type_minor",
    "cdkn_any_cell_type_minor",
    "snc_cdkn1a_cell_type_minor",
    "snc_cdkn2a_cell_type_minor",
    "snc_cdkn_any_cell_type_minor",
    "cdkn1a_immcelltype_or_minor",
    "cdkn2a_immcelltype_or_minor",
    "cdkn_any_immcelltype_or_minor",
]

MIN_GROUP_CELLS = 5


def decode_array(values: np.ndarray) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def sem(values: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce").dropna()
    if len(values) <= 1:
        return np.nan
    return float(values.std(ddof=1) / np.sqrt(len(values)))


def write_csv_safely(df: pd.DataFrame, path: Path, **kwargs) -> Path:
    try:
        df.to_csv(path, **kwargs)
        return path
    except PermissionError:
        fallback = path.with_name(path.stem + "_current" + path.suffix)
        df.to_csv(fallback, **kwargs)
        return fallback


def age_from_sample(sample: str) -> str:
    if str(sample).upper().endswith("-O"):
        return "Old"
    if str(sample).upper().endswith("-Y"):
        return "Young"
    return "Unknown"


def load_metadata() -> pd.DataFrame:
    cols = [
        "cell_id",
        "sample",
        "age",
        "snc",
        "immune_cell",
        "cell_type_minor",
        "cell_type_final",
        "immune_subclusters1",
        "zone",
    ]
    meta = pd.read_csv(METADATA_PATH, usecols=cols)
    markers = pd.read_csv(MARKER_EXPR_PATH, usecols=["cell_id", "sample", "CDKN1A", "CDKN2A"])
    markers["cdkn1a_status"] = np.where(markers["CDKN1A"].gt(0), "CDKN1A+", "CDKN1A-")
    markers["cdkn2a_status"] = np.where(markers["CDKN2A"].gt(0), "CDKN2A+", "CDKN2A-")
    markers["cdkn_any_status"] = np.select(
        [
            markers["CDKN1A"].gt(0) & markers["CDKN2A"].gt(0),
            markers["CDKN1A"].gt(0),
            markers["CDKN2A"].gt(0),
        ],
        ["CDKN1A+/CDKN2A+", "CDKN1A+ only", "CDKN2A+ only"],
        default="CDKN-",
    )
    meta = meta.merge(
        markers[["cell_id", "sample", "CDKN1A", "CDKN2A", "cdkn1a_status", "cdkn2a_status", "cdkn_any_status"]],
        on=["cell_id", "sample"],
        how="left",
    )
    meta["age"] = meta["age"].fillna(meta["sample"].map(age_from_sample))
    meta["snc"] = meta["snc"].fillna("SnC-")
    meta["cdkn1a_status"] = meta["cdkn1a_status"].fillna("CDKN1A-")
    meta["cdkn2a_status"] = meta["cdkn2a_status"].fillna("CDKN2A-")
    meta["cdkn_any_status"] = meta["cdkn_any_status"].fillna("CDKN-")
    meta["cell_type_minor"] = meta["cell_type_minor"].fillna(meta["cell_type_final"]).fillna("Unknown")
    meta["immcelltype"] = meta["immune_subclusters1"].replace({"none": np.nan, "nan": np.nan})
    meta["immcelltype"] = meta["immcelltype"].where(meta["immcelltype"].notna(), "non_immune_or_untyped")
    meta["immcelltype_or_minor"] = np.where(
        meta["immcelltype"].eq("non_immune_or_untyped"),
        meta["cell_type_minor"],
        meta["immcelltype"],
    )
    meta["snc_cell_type_minor"] = meta["snc"].astype(str) + " | " + meta["cell_type_minor"].astype(str)
    meta["snc_immcelltype"] = meta["snc"].astype(str) + " | " + meta["immcelltype"].astype(str)
    meta["snc_immcelltype_or_minor"] = meta["snc"].astype(str) + " | " + meta["immcelltype_or_minor"].astype(str)
    meta["cdkn1a_cell_type_minor"] = meta["cdkn1a_status"].astype(str) + " | " + meta["cell_type_minor"].astype(str)
    meta["cdkn2a_cell_type_minor"] = meta["cdkn2a_status"].astype(str) + " | " + meta["cell_type_minor"].astype(str)
    meta["cdkn_any_cell_type_minor"] = meta["cdkn_any_status"].astype(str) + " | " + meta["cell_type_minor"].astype(str)
    meta["snc_cdkn1a_cell_type_minor"] = (
        meta["snc"].astype(str) + " | " + meta["cdkn1a_status"].astype(str) + " | " + meta["cell_type_minor"].astype(str)
    )
    meta["snc_cdkn2a_cell_type_minor"] = (
        meta["snc"].astype(str) + " | " + meta["cdkn2a_status"].astype(str) + " | " + meta["cell_type_minor"].astype(str)
    )
    meta["snc_cdkn_any_cell_type_minor"] = (
        meta["snc"].astype(str) + " | " + meta["cdkn_any_status"].astype(str) + " | " + meta["cell_type_minor"].astype(str)
    )
    meta["cdkn1a_immcelltype_or_minor"] = (
        meta["cdkn1a_status"].astype(str) + " | " + meta["immcelltype_or_minor"].astype(str)
    )
    meta["cdkn2a_immcelltype_or_minor"] = (
        meta["cdkn2a_status"].astype(str) + " | " + meta["immcelltype_or_minor"].astype(str)
    )
    meta["cdkn_any_immcelltype_or_minor"] = (
        meta["cdkn_any_status"].astype(str) + " | " + meta["immcelltype_or_minor"].astype(str)
    )
    return meta


def read_sample_expression(sample: str, genes: list[str], wanted_cells: set[str]) -> pd.DataFrame:
    h5_path = XENIUM_DIR / sample / "cell_feature_matrix.h5"
    with h5py.File(h5_path, "r") as handle:
        matrix_group = handle["matrix"]
        feature_names = [gene.upper() for gene in decode_array(matrix_group["features"]["name"][:])]
        barcodes = decode_array(matrix_group["barcodes"][:])
        shape = tuple(int(x) for x in matrix_group["shape"][:])
        matrix = sparse.csc_matrix(
            (
                matrix_group["data"][:],
                matrix_group["indices"][:],
                matrix_group["indptr"][:],
            ),
            shape=shape,
        )

    lookup = {gene: i for i, gene in enumerate(feature_names)}
    present = [gene for gene in genes if gene in lookup]
    full_cell_ids = np.array([f"{sample}_{barcode}" for barcode in barcodes], dtype=object)
    keep = np.array([cell_id in wanted_cells for cell_id in full_cell_ids])
    if not present or not np.any(keep):
        return pd.DataFrame()

    counts = matrix[[lookup[gene] for gene in present], :][:, keep].astype(np.float32).toarray()
    totals = counts.sum(axis=0)
    norm = np.divide(counts, totals, out=np.zeros_like(counts), where=totals > 0) * 10000
    log_norm = np.log1p(norm).astype(np.float32)

    expr = pd.DataFrame(log_norm.T, columns=present)
    expr.insert(0, "cell_id", full_cell_ids[keep])
    return expr


def build_group_gene_expression(meta: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    for sample, sample_meta in meta.groupby("sample", observed=True):
        expr = read_sample_expression(sample, genes, set(sample_meta["cell_id"]))
        if expr.empty:
            continue
        sample_frame = sample_meta.merge(expr, on="cell_id", how="inner")
        for group_col in GROUP_COLS:
            for group_value, sub in sample_frame.groupby(group_col, observed=True):
                if len(sub) < MIN_GROUP_CELLS:
                    continue
                for gene in genes:
                    if gene not in sub.columns:
                        continue
                    values = sub[gene]
                    rows.append(
                        {
                            "sample": sample,
                            "age": sub["age"].iloc[0],
                            "grouping": group_col,
                            "group": group_value,
                            "gene": gene,
                            "n_cells": len(sub),
                            "mean_log1p_cp10k": float(values.mean()),
                            "pct_expressing": float(values.gt(0).mean()),
                        }
                    )
    return pd.DataFrame(rows)


def score_pairs_by_sample(group_expr: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sample, age, grouping), sub in group_expr.groupby(["sample", "age", "grouping"], observed=True):
        ligand_expr = sub.rename(
            columns={
                "group": "sender_group",
                "gene": "ligand",
                "n_cells": "sender_n_cells",
                "mean_log1p_cp10k": "sender_ligand_mean_log1p_cp10k",
                "pct_expressing": "sender_ligand_pct",
            }
        )
        receptor_expr = sub.rename(
            columns={
                "group": "receiver_group",
                "gene": "receptor",
                "n_cells": "receiver_n_cells",
                "mean_log1p_cp10k": "receiver_receptor_mean_log1p_cp10k",
                "pct_expressing": "receiver_receptor_pct",
            }
        )
        for ligand, receptor, panel_source, lr_role_note in pairs[
            ["ligand", "receptor", "panel_source", "lr_role_note"]
        ].itertuples(index=False):
            lig = ligand_expr[ligand_expr["ligand"].eq(ligand)]
            rec = receptor_expr[receptor_expr["receptor"].eq(receptor)]
            if lig.empty or rec.empty:
                continue
            crossed = lig.merge(rec, how="cross")
            crossed["communication_score"] = (
                crossed["sender_ligand_mean_log1p_cp10k"]
                * crossed["sender_ligand_pct"]
                * crossed["receiver_receptor_mean_log1p_cp10k"]
                * crossed["receiver_receptor_pct"]
            )
            crossed = crossed[crossed["communication_score"].gt(0)].copy()
            if crossed.empty:
                continue
            crossed["sample"] = sample
            crossed["age"] = age
            crossed["grouping"] = grouping
            crossed["panel_source"] = panel_source
            crossed["lr_role_note"] = lr_role_note
            keep_cols = [
                "sample",
                "age",
                "grouping",
                "panel_source",
                "lr_role_note",
                "ligand",
                "receptor",
                "sender_group",
                "receiver_group",
                "communication_score",
                "sender_ligand_mean_log1p_cp10k",
                "sender_ligand_pct",
                "receiver_receptor_mean_log1p_cp10k",
                "receiver_receptor_pct",
                "sender_n_cells",
                "receiver_n_cells",
            ]
            rows.append(crossed[keep_cols])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def summarize_scores(scores: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    group_cols = [
        "grouping",
        "age",
        "panel_source",
        "lr_role_note",
        "ligand",
        "receptor",
        "sender_group",
        "receiver_group",
    ]
    summary = (
        scores.groupby(group_cols, observed=True)
        .agg(
            n_samples=("sample", "nunique"),
            mean_score=("communication_score", "mean"),
            sem_score=("communication_score", sem),
            mean_sender_ligand_pct=("sender_ligand_pct", "mean"),
            mean_receiver_receptor_pct=("receiver_receptor_pct", "mean"),
            mean_sender_n_cells=("sender_n_cells", "mean"),
            mean_receiver_n_cells=("receiver_n_cells", "mean"),
            samples=("sample", lambda x: ", ".join(sorted(x.astype(str).unique()))),
        )
        .reset_index()
    )

    key_cols = [
        "grouping",
        "panel_source",
        "lr_role_note",
        "ligand",
        "receptor",
        "sender_group",
        "receiver_group",
    ]
    delta = (
        summary.pivot_table(
            index=key_cols,
            columns="age",
            values="mean_score",
            aggfunc="first",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    if "Old" in delta.columns and "Young" in delta.columns:
        delta["old_minus_young"] = delta["Old"] - delta["Young"]
        delta["old_over_young"] = np.divide(
            delta["Old"],
            delta["Young"],
            out=np.full(len(delta), np.nan),
            where=delta["Young"].to_numpy() > 0,
        )
    return summary, delta


def write_focus_tables(summary: pd.DataFrame, delta: pd.DataFrame) -> None:
    old = summary[summary["age"].eq("Old")].copy()
    write_csv_safely(
        old.sort_values(["grouping", "ligand", "receptor", "mean_score"], ascending=[True, True, True, False]),
        OUTDIR / "xenium_focused_secretory_lr_old_by_snc_immune_minor_summary.csv",
        index=False,
    )
    if "old_minus_young" in delta.columns:
        write_csv_safely(
            delta.sort_values(
                ["grouping", "ligand", "receptor", "old_minus_young"],
                ascending=[True, True, True, False],
            ),
            OUTDIR / "xenium_focused_secretory_lr_old_minus_young_by_snc_immune_minor.csv",
            index=False,
        )

    for grouping in GROUP_COLS:
        old_group = old[old["grouping"].eq(grouping)].copy()
        if not old_group.empty:
            top = (
                old_group.sort_values("mean_score", ascending=False)
                .groupby(["ligand", "receptor"], observed=True, group_keys=False)
                .head(10)
            )
            write_csv_safely(top, OUTDIR / f"top_old_{grouping}_focused_lr.csv", index=False)
        if "old_minus_young" in delta.columns:
            delta_group = delta[delta["grouping"].eq(grouping)].copy()
            if not delta_group.empty:
                top_delta = (
                    delta_group.sort_values("old_minus_young", ascending=False)
                    .groupby(["ligand", "receptor"], observed=True, group_keys=False)
                    .head(10)
                )
                write_csv_safely(top_delta, OUTDIR / f"top_old_minus_young_{grouping}_focused_lr.csv", index=False)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    pairs = pd.DataFrame(FOCUSED_PAIRS, columns=["ligand", "receptor", "panel_source", "lr_role_note"])
    write_csv_safely(pairs, OUTDIR / "xenium_focused_secretory_lr_pair_role_annotation.csv", index=False)
    genes = sorted(set(pairs["ligand"]).union(pairs["receptor"]))
    meta = load_metadata()

    group_expr = build_group_gene_expression(meta, genes)
    write_csv_safely(group_expr, OUTDIR / "xenium_focused_lr_group_gene_expression.csv", index=False)
    scores = score_pairs_by_sample(group_expr, pairs)
    write_csv_safely(scores, OUTDIR / "xenium_focused_secretory_lr_scores_by_sample.csv.gz", index=False)
    summary, delta = summarize_scores(scores)
    write_csv_safely(summary, OUTDIR / "xenium_focused_secretory_lr_by_age_group_summary.csv", index=False)
    write_csv_safely(delta, OUTDIR / "xenium_focused_secretory_lr_old_minus_young_summary.csv", index=False)
    write_focus_tables(summary, delta)

    print(f"Saved outputs to: {OUTDIR}")
    print("\nTop old SnC + minor cell-type communication scores:")
    top = (
        summary[
            summary["age"].eq("Old")
            & summary["grouping"].eq("snc_cell_type_minor")
            & summary["n_samples"].ge(2)
        ]
        .sort_values("mean_score", ascending=False)
        .head(30)
    )
    print(
        top[
            [
                "ligand",
                "receptor",
                "sender_group",
                "receiver_group",
                "n_samples",
                "mean_score",
                "mean_sender_ligand_pct",
                "mean_receiver_receptor_pct",
                "samples",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
