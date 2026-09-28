"""
Compact prioritized plots for Squidpy SSSP ligand-receptor results.

Reads the full long-form output from
21_xenium_squidpy_sssp_ligrec_by_snc_celltype.py and creates less crowded
dotplots by:
1. ranking SnC+ minor cell types by the number of significant LR links,
2. keeping only the top SnC+ groups,
3. keeping only the strongest LR pairs for those groups,
4. aggregating each LR/group combination to its strongest partner group.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUTDIR = PROJECT_ROOT / "results" / "260727" / "xenium_squidpy_sssp_ligrec_by_snc_celltype"
LONG_RESULTS = OUTDIR / "xenium_squidpy_sssp_ligrec_all_group_pairs.csv.gz"
GROUP_COUNTS = OUTDIR / "xenium_squidpy_sssp_ligrec_group_counts.csv"

AGE_ORDER = ["Old", "Young"]
GENE_SET_ORDER = ["main_20", "diff_filtered_17"]
COMPARTMENTS = ["skin", "immune"]
MODES = ["outgoing", "incoming"]
IMMUNE_MINOR_TYPES = {
    "B Cells",
    "CD4+ T Cells",
    "CD8+ T Cells/NK Cells",
    "Dendritic Cells",
    "Langerhans Cells",
    "Macrophages",
    "Mast Cells",
}
TOP_GROUPS = 4
TOP_LR_PAIRS = 12
PVALUE_THRESHOLD = 0.05


def compartment_of_minor(cell_type: str) -> str:
    return "immune" if str(cell_type) in IMMUNE_MINOR_TYPES else "skin"


def clean_group_label(group: str) -> str:
    return str(group).replace("SnC+|", "").replace("SnC-|", "SnC- ")


def rank_groups(long: pd.DataFrame, candidate_groups: list[str], mode: str) -> pd.DataFrame:
    group_col = "source_group" if mode == "outgoing" else "target_group"
    work = long[long[group_col].isin(candidate_groups)].copy()
    if work.empty:
        return pd.DataFrame()
    work["is_priority"] = work["pvalue"].le(PVALUE_THRESHOLD)
    ranked = (
        work.groupby(group_col, observed=True)
        .agg(
            n_priority_links=("is_priority", "sum"),
            n_priority_lr_pairs=("lr_pair", lambda s: s[work.loc[s.index, "is_priority"]].nunique()),
            n_total_links=("lr_pair", "size"),
            max_mean=("mean_expression_product", "max"),
            min_pvalue=("pvalue", "min"),
        )
        .reset_index()
        .rename(columns={group_col: "snc_plus_group"})
        .sort_values(["n_priority_links", "n_priority_lr_pairs", "max_mean", "min_pvalue"], ascending=[False, False, False, True])
    )
    return ranked


def strongest_partner_rows(long: pd.DataFrame, selected_groups: list[str], mode: str) -> pd.DataFrame:
    group_col = "source_group" if mode == "outgoing" else "target_group"
    partner_col = "target_group" if mode == "outgoing" else "source_group"
    partner_type_col = "target_cell_type_minor" if mode == "outgoing" else "source_cell_type_minor"
    work = long[long[group_col].isin(selected_groups)].copy()
    if work.empty:
        return work
    work["is_priority"] = work["pvalue"].le(PVALUE_THRESHOLD)
    work = work.sort_values(["pvalue", "mean_expression_product"], ascending=[True, False])
    best = work.groupby([group_col, "lr_pair"], observed=True).head(1).copy()
    best = best.rename(
        columns={
            group_col: "snc_plus_group",
            partner_col: "strongest_partner_group",
            partner_type_col: "strongest_partner_cell_type_minor",
        }
    )
    return best


def select_top_lr_pairs(best: pd.DataFrame) -> list[str]:
    if best.empty:
        return []
    ranked = (
        best.groupby("lr_pair", observed=True)
        .agg(
            n_priority_groups=("is_priority", "sum"),
            n_groups=("snc_plus_group", "nunique"),
            max_mean=("mean_expression_product", "max"),
            min_pvalue=("pvalue", "min"),
        )
        .reset_index()
        .sort_values(["n_priority_groups", "n_groups", "max_mean", "min_pvalue"], ascending=[False, False, False, True])
    )
    return ranked.head(TOP_LR_PAIRS)["lr_pair"].tolist()


def plot_prioritized(best: pd.DataFrame, group_order: list[str], lr_order: list[str], title: str, outbase: Path) -> None:
    if best.empty or not group_order or not lr_order:
        return
    work = best[best["snc_plus_group"].isin(group_order) & best["lr_pair"].isin(lr_order)].copy()
    if work.empty:
        return
    group_to_x = {group: i for i, group in enumerate(group_order)}
    lr_to_y = {lr: i for i, lr in enumerate(lr_order[::-1])}
    work["x"] = work["snc_plus_group"].map(group_to_x)
    work["y"] = work["lr_pair"].map(lr_to_y)
    max_mean = max(float(work["mean_expression_product"].max()), 1e-9)
    sizes = 45 + 420 * np.sqrt(work["mean_expression_product"].to_numpy(dtype=float) / max_mean)
    colors = -np.log10(work["pvalue"].clip(lower=1e-300)).clip(upper=8)

    fig_width = max(6.2, 1.45 * len(group_order) + 3.2)
    fig_height = max(5.0, 0.42 * len(lr_order) + 1.7)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    scatter = ax.scatter(
        work["x"],
        work["y"],
        s=sizes,
        c=colors,
        cmap="magma",
        edgecolor="#333333",
        linewidth=0.45,
        alpha=0.88,
    )
    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels([clean_group_label(g) for g in group_order], rotation=35, ha="right")
    ax.set_yticks(range(len(lr_order)))
    ax.set_yticklabels(lr_order[::-1])
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("Prioritized SnC+ minor cell type")
    ax.set_ylabel("Prioritized SSSP ligand-receptor pair")
    ax.grid(color="#E6E6E6", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label("-log10 permutation p, clipped at 8")

    legend_values = np.array([0.25, 0.5, 1.0]) * max_mean
    handles = [
        ax.scatter([], [], s=45 + 420 * np.sqrt(v / max_mean), color="#BBBBBB", edgecolor="#333333", linewidth=0.45)
        for v in legend_values
    ]
    ax.legend(handles, [f"{v:.2f}" for v in legend_values], title="LR mean", frameon=False, bbox_to_anchor=(1.22, 0.5), loc="center left")

    fig.tight_layout()
    fig.savefig(outbase.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    long = pd.read_csv(LONG_RESULTS)
    counts = pd.read_csv(GROUP_COUNTS)
    counts["compartment"] = counts["cell_type_minor"].map(compartment_of_minor)
    counts = counts[counts["included_min_cells"].astype(bool) & counts["snc"].eq("SnC+")]

    all_ranks = []
    all_best_rows = []
    for gene_set_id in GENE_SET_ORDER:
        for age in AGE_ORDER:
            subset = long[long["gene_set_id"].eq(gene_set_id) & long["age"].eq(age)].copy()
            for compartment in COMPARTMENTS:
                candidate_groups = counts.loc[
                    counts["gene_set_id"].eq(gene_set_id)
                    & counts["age"].eq(age)
                    & counts["compartment"].eq(compartment),
                    "snc_minor_group",
                ].tolist()
                if not candidate_groups:
                    continue
                for mode in MODES:
                    ranks = rank_groups(subset, candidate_groups, mode)
                    if ranks.empty:
                        continue
                    ranks.insert(0, "gene_set_id", gene_set_id)
                    ranks.insert(1, "age", age)
                    ranks.insert(2, "compartment", compartment)
                    ranks.insert(3, "mode", mode)
                    ranks["selected_for_plot"] = False
                    selected_groups = ranks.head(TOP_GROUPS)["snc_plus_group"].tolist()
                    ranks.loc[ranks["snc_plus_group"].isin(selected_groups), "selected_for_plot"] = True
                    all_ranks.append(ranks)

                    best = strongest_partner_rows(subset, selected_groups, mode)
                    lr_order = select_top_lr_pairs(best)
                    best = best[best["lr_pair"].isin(lr_order)].copy()
                    if best.empty:
                        continue
                    best.insert(0, "plot_gene_set_id", gene_set_id)
                    best.insert(1, "plot_age", age)
                    best.insert(2, "plot_compartment", compartment)
                    best.insert(3, "plot_mode", mode)
                    all_best_rows.append(best)
                    best.to_csv(
                        OUTDIR / f"{gene_set_id}_{age.lower()}_{compartment}_sncplus_{mode}_prioritized_compact_ligrec_values.csv",
                        index=False,
                    )
                    plot_prioritized(
                        best,
                        selected_groups,
                        lr_order,
                        f"{age} {gene_set_id}: prioritized {compartment} SnC+ {mode} SSSP LR",
                        OUTDIR / f"{gene_set_id}_{age.lower()}_{compartment}_sncplus_{mode}_prioritized_compact_ligrec_dotplot",
                    )

    ranks_all = pd.concat(all_ranks, ignore_index=True) if all_ranks else pd.DataFrame()
    best_all = pd.concat(all_best_rows, ignore_index=True) if all_best_rows else pd.DataFrame()
    ranks_all.to_csv(OUTDIR / "xenium_squidpy_sssp_ligrec_prioritized_group_ranks.csv", index=False)
    best_all.to_csv(OUTDIR / "xenium_squidpy_sssp_ligrec_prioritized_compact_plot_values.csv", index=False)

    print("Prioritized compact SSSP ligrec plots complete.")
    if not ranks_all.empty:
        print(ranks_all[ranks_all["selected_for_plot"]].to_string(index=False))
    print(f"\nWrote outputs to: {OUTDIR}")


if __name__ == "__main__":
    main()
