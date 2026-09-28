#!/usr/bin/env python
"""Run the existing Xenium Squidpy workflow with Core21 ligands and SnC+ senders."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import omnipath as op
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BASE_PATH = Path(__file__).with_name("31_xenium_shared_sssp_squidpy_ligrec.py")
OUT = ROOT / "results" / "public_validation" / "xenium_core21_snc_sender_ligrec"
PAIR_CACHE = OUT / "core21_omnipath_ligand_receptor_pairs.csv"
CORE21 = "CXCL16 VEGFA IL18 IL1B CSF1 IL4 CCL19 CXCL1 CCL28 CCL26 CXCL11 CCL11 EGF CXCL5 IL6 CXCL6 IL5 MMP3 CCL7 CCL15 CCL16".split()


def load_core21_interactions(adata, genes):
    if PAIR_CACHE.exists():
        return pd.read_csv(PAIR_CACHE)[["source", "target"]]
    network = op.interactions.import_intercell_network()
    network["source"] = network["genesymbol_intercell_source"].astype(str).str.upper()
    network["target"] = network["genesymbol_intercell_target"].astype(str).str.upper()
    measured = set(adata.var_names.astype(str).str.upper())
    pairs = network[
        network["source"].isin(genes)
        & network["target"].isin(measured)
        & network["source"].isin(measured)
        & network["transmitter_intercell_source"].fillna(False)
        & network["receiver_intercell_target"].fillna(False)
        & network["source"].ne(network["target"])
    ][["source", "target"]].drop_duplicates().sort_values(["source", "target"])
    if pairs.empty:
        raise RuntimeError("No Xenium-measurable Core21 ligand-receptor pairs found in OmniPath.")
    OUT.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(PAIR_CACHE, index=False)
    return pairs


def summarize_results():
    frames = []
    for age in ("old", "young"):
        frame = pd.read_csv(OUT / f"core21_{age}_snc_plus_sender_focus_ligrec_results.csv")
        frames.append(frame[frame["pvalue"].le(.05)])
    significant = pd.concat(frames, ignore_index=True)
    keys = ["age", "source_cell_type_minor", "target_snc", "target_cell_type_minor"]
    summary = significant.groupby(keys, observed=True).agg(
        n_significant_lr_pairs=("lr_pair", "nunique"),
        max_mean_expression_product=("mean_expression_product", "max"),
        mean_expression_product=("mean_expression_product", "mean"),
        top_ligands=("ligand", lambda x: ";".join(x.value_counts().head(5).index)),
        top_lr_pairs=("lr_pair", lambda x: ";".join(x.value_counts().head(5).index)),
    ).reset_index()
    summary.to_csv(OUT / "core21_significant_snc_sender_receiver_summary.csv", index=False)

    for age, frame in significant.groupby("age"):
        top = frame.sort_values(["pvalue", "mean_expression_product"], ascending=[True, False]).head(20).iloc[::-1]
        labels = (top["ligand"] + "→" + top["receptor"] + "\n" +
                  top["source_cell_type_minor"] + " → " + top["target_cell_type_minor"])
        fig, ax = plt.subplots(figsize=(12, max(7, .65 * len(top) + 2)))
        ax.barh(range(len(top)), top["mean_expression_product"], color="#7A5195")
        ax.set_yticks(range(len(top)), labels, fontsize=8)
        ax.set_xlabel("Squidpy mean-expression product")
        ax.set_title(f"{age}: top significant Core21 signals from SnC+ senders")
        ax.grid(axis="x", color="#E6E6E6")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(OUT / f"core21_{age.lower()}_top20_snc_sender_edges_spaced.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    collapsed = significant.groupby(
        ["age", "source_cell_type_minor", "target_cell_type_minor"], observed=True
    )["mean_expression_product"].max().reset_index()
    receivers = (collapsed.groupby("target_cell_type_minor")["mean_expression_product"]
                 .max().nlargest(20).index)
    tables = {}
    for age in ("Young", "Old"):
        table = (collapsed[collapsed["age"].eq(age) & collapsed["target_cell_type_minor"].isin(receivers)]
                 .pivot(index="source_cell_type_minor", columns="target_cell_type_minor",
                        values="mean_expression_product").reindex(columns=receivers))
        tables[age] = table
    heights = [max(5, .55 * len(tables[age]) + 2) for age in ("Young", "Old")]
    fig, axes = plt.subplots(2, 1, figsize=(15, sum(heights)),
                             gridspec_kw={"height_ratios": heights}, constrained_layout=True)
    for ax, age in zip(axes, ("Young", "Old")):
        table = tables[age]
        image = ax.imshow(table.fillna(0), cmap="magma", aspect="auto", vmin=0)
        ax.set_xticks(range(len(table.columns)), table.columns, rotation=55, ha="right", fontsize=8)
        ax.set_yticks(range(len(table.index)), table.index, fontsize=9)
        ax.set_title(f"{age}: significant Core21 communication from SnC+ senders")
        ax.set_xlabel("Receiver cell type")
        ax.set_ylabel("SnC+ sender cell type")
        fig.colorbar(image, ax=ax, label="Maximum Squidpy mean-expression product", shrink=.75)
    fig.savefig(OUT / "core21_significant_snc_sender_receiver_heatmaps.png", dpi=300)
    plt.close(fig)


def main():
    spec = spec_from_file_location("xenium_squidpy_base", BASE_PATH)
    base = module_from_spec(spec)
    spec.loader.exec_module(base)
    base.OUTDIR = OUT
    base.LR_PATH = PAIR_CACHE
    base.load_gene_sets = lambda: {"core21": CORE21}
    base.load_interactions = load_core21_interactions
    base.main()
    summarize_results()


if __name__ == "__main__":
    main()
