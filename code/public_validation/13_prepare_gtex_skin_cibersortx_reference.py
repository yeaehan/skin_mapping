#!/usr/bin/env python
"""Prepare a donor-balanced GTEx-skin CIBERSORTx single-cell reference.

By default, all observed Broad cell types are retained and cells are sampled
within donor/type strata. Eligibility thresholds are reported for audit but do
not remove scarce skin compartments. A strict mode and an optional MuSiC
manifest remain available for sensitivity analyses.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=root)
    parser.add_argument("--adata", type=Path, default=None)
    parser.add_argument("--music-manifest", type=Path, default=None,
                        help="Optional CSV/TSV restricting donors/types to a MuSiC reference.")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--selection-mode", choices=("all", "strict"), default="all",
                        help="Retain all observed donor/type strata (default), or apply minimum cell/donor thresholds.")
    parser.add_argument("--seed", type=int, default=20260920)
    parser.add_argument("--min-cells-donor-type", type=int, default=20)
    parser.add_argument("--min-donors-type", type=int, default=2)
    parser.add_argument("--max-donors-type", type=int, default=3)
    parser.add_argument("--max-cells-donor-type", type=int, default=50)
    parser.add_argument("--max-total-cells", type=int, default=5000)
    parser.add_argument("--min-detected-cells", type=int, default=3)
    return parser.parse_args()


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t" if path.suffix.lower() in {".tsv", ".txt"} else ",")


def normalize_manifest(table: pd.DataFrame) -> pd.DataFrame:
    aliases = {
        "donor": ["donor", "participant_id", "Participant ID", "individual"],
        "cell_type": ["cell_type", "broad_cell_type", "Broad cell type"],
    }
    rename = {}
    for target, candidates in aliases.items():
        source = next((column for column in candidates if column in table.columns), None)
        if source is None:
            raise ValueError(f"MuSiC manifest needs a {target!r} column; found {table.columns.tolist()}")
        rename[source] = target
    table = table.rename(columns=rename).copy()
    if "retain" in table.columns:
        retain = table["retain"].astype(str).str.lower().isin({"true", "1", "yes", "y"})
        table = table.loc[retain].copy()
    table["donor"] = table["donor"].astype(str).str.strip()
    table["cell_type"] = table["cell_type"].astype(str).str.strip()
    table = table.loc[table.donor.ne("") & table.cell_type.ne("")]
    if table.empty:
        raise ValueError("MuSiC manifest has no retained donor/type rows")
    return table[["donor", "cell_type"]].drop_duplicates().sort_values(["cell_type", "donor"])


def safe_phenotypes(cell_types: list[str]) -> dict[str, str]:
    mapping = {}
    used = set()
    for cell_type in sorted(set(cell_types)):
        # Periods are deliberately excluded: the only period in a CIBERSORTx
        # cell-column label is the final .N suffix.
        safe = re.sub(r"[^A-Za-z0-9_-]+", "_", cell_type).strip("_")
        if not safe:
            raise ValueError(f"Cell type cannot be converted to a phenotype label: {cell_type!r}")
        if safe in used:
            raise ValueError(f"Phenotype-label collision after sanitizing {cell_type!r}: {safe!r}")
        mapping[cell_type] = safe
        used.add(safe)
    return mapping


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_sparse_reference(path: Path, genes: np.ndarray, matrix: sparse.spmatrix,
                           column_names: list[str], chunk_size: int = 250,
                           float_format: str = "%.0f") -> None:
    """Write genes x cells without densifying the complete matrix."""
    genes_by_cells = matrix.T.tocsr()
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("GeneSymbol\t" + "\t".join(column_names) + "\n")
        for start in range(0, len(genes), chunk_size):
            end = min(start + chunk_size, len(genes))
            block = genes_by_cells[start:end].toarray()
            frame = pd.DataFrame(block, index=genes[start:end], columns=column_names)
            frame.index.name = "GeneSymbol"
            frame.to_csv(handle, sep="\t", header=False, float_format=float_format)


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    adata_path = (args.adata or root / "data" / "external" /
                  "GTEx_skin_snRNAseq_atlas_071421.h5ad").resolve()
    output_dir = (args.output_dir or root / "results" / "public_validation" /
                  "gtex_skin_cibersortx" / "reference").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    adata = ad.read_h5ad(adata_path)
    donor_col, type_col = "Participant ID", "Broad cell type"
    for column in [donor_col, type_col]:
        if column not in adata.obs:
            raise KeyError(f"Missing required observation column {column!r}")
    if adata.raw is None:
        raise ValueError("adata.raw is absent; raw.X UMI counts are required")
    obs = adata.obs[[donor_col, type_col]].astype(str).rename(
        columns={donor_col: "donor", type_col: "cell_type"}
    )
    observed = obs.groupby(["donor", "cell_type"], observed=True).size().rename("n_available").reset_index()
    observed["meets_min_cells"] = observed.n_available.ge(args.min_cells_donor_type)
    donors_per_type = observed.loc[observed.meets_min_cells].groupby("cell_type").donor.nunique()
    observed["eligible_donors_for_type"] = observed.cell_type.map(donors_per_type).fillna(0).astype(int)
    observed["type_meets_min_donors"] = observed.eligible_donors_for_type.ge(args.min_donors_type)
    observed.to_csv(output_dir / "all_donor_type_eligibility.csv", index=False)

    template = observed.loc[observed.meets_min_cells & observed.type_meets_min_donors,
                            ["donor", "cell_type", "n_available"]].copy()
    template["retain"] = False
    template.to_csv(output_dir / "music_reference_manifest_template.csv", index=False)
    if args.audit_only:
        logging.info("Audit only: %d donor/type rows; %d eligible template rows", len(observed), len(template))
        logging.info("Wrote %s", output_dir)
        return
    if args.music_manifest is None and args.selection_mode == "all":
        requested = observed.copy()
        logging.info("Selecting every observed GTEx skin Broad cell type; thresholds remain audit-only")
    elif args.music_manifest is None:
        requested = observed.loc[
            observed.meets_min_cells & observed.type_meets_min_donors,
            ["donor", "cell_type", "n_available", "meets_min_cells",
             "eligible_donors_for_type", "type_meets_min_donors"],
        ].copy()
        logging.info("Strict mode: selecting only threshold-eligible donor/type rows")
    else:
        music_manifest = normalize_manifest(read_table(args.music_manifest.resolve()))
        requested = music_manifest.merge(observed, on=["donor", "cell_type"], how="left", validate="one_to_one")
    missing = requested.loc[requested.n_available.isna(), ["donor", "cell_type"]]
    enforce_thresholds = args.selection_mode == "strict" or args.music_manifest is not None
    too_small = requested.loc[
        requested.n_available.fillna(0).lt(args.min_cells_donor_type),
        ["donor", "cell_type", "n_available"],
    ] if enforce_thresholds else requested.iloc[0:0]
    if not missing.empty or not too_small.empty:
        requested.to_csv(output_dir / "music_manifest_validation_failed.csv", index=False)
        raise ValueError(f"MuSiC manifest has {len(missing)} absent and {len(too_small)} <{args.min_cells_donor_type}-cell donor/type rows")
    if enforce_thresholds:
        retained_donors = requested.groupby("cell_type").donor.nunique()
        failing_types = retained_donors[retained_donors.lt(args.min_donors_type)]
        if not failing_types.empty:
            raise ValueError("Retained types failing minimum donors/type: " + failing_types.to_dict().__repr__())

    rng = np.random.default_rng(args.seed)
    selected_parts = []
    selection_rows = []
    for cell_type, group in requested.groupby("cell_type", sort=True):
        donors = sorted(group.donor.unique())
        if len(donors) > args.max_donors_type:
            donors = sorted(rng.choice(donors, size=args.max_donors_type, replace=False).tolist())
        for donor in donors:
            candidates = obs.index[(obs.donor == donor) & (obs.cell_type == cell_type)].to_numpy()
            n_take = min(args.max_cells_donor_type, len(candidates))
            chosen = rng.choice(candidates, size=n_take, replace=False)
            selected_parts.append(chosen)
            selection_rows.append({"donor": donor, "cell_type": cell_type,
                                   "n_available": len(candidates), "n_selected": n_take})
    selected = np.concatenate(selected_parts)
    if len(selected) > args.max_total_cells:
        # Balanced second-stage cap: sample nearly equally from each selected stratum.
        selected_obs = obs.loc[selected]
        strata = list(selected_obs.groupby(["donor", "cell_type"], sort=True).groups.values())
        base = args.max_total_cells // len(strata)
        if enforce_thresholds and base < args.min_cells_donor_type:
            raise ValueError("5,000-cell cap is incompatible with >=20 cells for every selected donor/type stratum")
        balanced = [rng.choice(np.asarray(list(index)), size=min(base, len(index)), replace=False) for index in strata]
        selected = np.concatenate(balanced)
        remaining = args.max_total_cells - len(selected)
        if remaining:
            pool = np.setdiff1d(np.concatenate([np.asarray(list(x)) for x in strata]), selected)
            selected = np.concatenate([selected, rng.choice(pool, size=min(remaining, len(pool)), replace=False)])
    selected = np.asarray(selected)
    selected_obs = obs.loc[selected].copy()

    raw = adata.raw[selected, :]
    raw_var = raw.var.copy()
    if "Approved symbol" not in raw_var:
        raise KeyError("adata.raw.var lacks 'Approved symbol'; unique HGNC symbols are required")
    symbols = raw_var["Approved symbol"].astype("string").str.strip()
    valid = symbols.notna() & symbols.ne("") & ~symbols.str.upper().str.startswith("MT-")
    symbols = symbols.loc[valid].astype(str)
    counts = raw.X[:, valid.to_numpy()]
    counts = sparse.csr_matrix(counts)
    if counts.data.size and (counts.data.min() < 0 or not np.allclose(counts.data, np.rint(counts.data), atol=1e-6)):
        raise ValueError("adata.raw.X is not non-negative integer-like UMI counts")

    unique_symbols = np.asarray(sorted(symbols.unique()))
    symbol_index = pd.Series(np.arange(len(unique_symbols)), index=unique_symbols)
    mapping = sparse.csr_matrix(
        (np.ones(len(symbols), dtype=np.float32),
         (np.arange(len(symbols)), symbol_index.loc[symbols].to_numpy())),
        shape=(len(symbols), len(unique_symbols)),
    )
    summed = counts @ mapping
    detected = np.asarray((summed > 0).sum(axis=0)).ravel()
    keep_genes = detected >= args.min_detected_cells
    summed = summed[:, keep_genes].tocsr()
    genes = unique_symbols[keep_genes]

    phenotype_map = safe_phenotypes(selected_obs.cell_type.tolist())
    selected_obs["cibersortx_phenotype"] = selected_obs.cell_type.map(phenotype_map)
    counters = selected_obs.groupby("cibersortx_phenotype").cumcount().add(1).astype(str)
    column_names = (selected_obs.cibersortx_phenotype + "." + counters).tolist()
    if any(name.count(".") != 1 or not re.search(r"\.\d+$", name) for name in column_names):
        raise AssertionError("CIBERSORTx phenotype columns must contain only a final .N suffix")

    reference_path = output_dir / "gtex_skin_cibersortx_single_cell_reference.tsv"
    write_sparse_reference(reference_path, genes, summed, column_names)
    # CIBERSORTx's single-cell sample format repeats the phenotype label for
    # every cell in row 1. Normalize each cell to CPM so this upload and the
    # GTEx bulk mixture occupy the same non-log CPM space. Keep the raw-count,
    # uniquely named file above for MuSiC and provenance.
    cell_totals = np.asarray(summed.sum(axis=1)).ravel()
    if np.any(cell_totals <= 0):
        raise ValueError("Selected reference contains a zero-library cell after gene filtering")
    cpm = sparse.diags(1_000_000.0 / cell_totals) @ summed.astype(np.float64)
    upload_labels = selected_obs["cibersortx_phenotype"].tolist()
    upload_path = output_dir / "gtex_skin_cibersortx_single_cell_reference_cpm.tsv"
    write_sparse_reference(upload_path, genes, cpm, upload_labels, float_format="%.6f")
    cell_manifest = selected_obs.copy()
    cell_manifest.insert(0, "source_cell_id", selected_obs.index.astype(str))
    cell_manifest["cibersortx_column"] = column_names
    cell_manifest["cibersortx_upload_label"] = upload_labels
    cell_manifest.to_csv(output_dir / "selected_cell_manifest.csv", index=False)
    final_selection = (selected_obs.groupby(["donor", "cell_type"], observed=True).size()
                       .rename("n_selected").reset_index())
    selection_audit = (requested[["donor", "cell_type", "n_available"]]
                       .merge(final_selection, on=["donor", "cell_type"], how="left"))
    selection_audit["n_selected"] = selection_audit["n_selected"].fillna(0).astype(int)
    selection_audit["donor_selected"] = selection_audit.n_selected.gt(0)
    selection_audit.to_csv(output_dir / "donor_type_selection_audit.csv", index=False)
    pd.DataFrame([{"original_cell_type": key, "cibersortx_phenotype": value}
                  for key, value in phenotype_map.items()]).to_csv(
        output_dir / "cell_type_name_map.csv", index=False
    )
    run_manifest = {
        "source_h5ad": str(adata_path),
        "music_manifest": str(args.music_manifest.resolve()) if args.music_manifest else None,
        "selection_source": "MuSiC manifest" if args.music_manifest else f"GTEx skin Broad cell types ({args.selection_mode})",
        "selection_mode": args.selection_mode,
        "thresholds_used_for_exclusion": enforce_thresholds,
        "expression_source": "adata.raw.X",
        "expression_scale": "raw non-log UMI counts",
        "gene_symbol_field": "adata.raw.var['Approved symbol']",
        "duplicate_symbol_policy": "sum",
        "mitochondrial_gene_policy": "drop symbols beginning MT- (case-insensitive)",
        "min_detected_selected_cells": args.min_detected_cells,
        "seed": args.seed,
        "n_selected_cells": int(len(selected)),
        "n_exported_genes": int(len(genes)),
        "reference_file": reference_path.name,
        "reference_sha256": sha256(reference_path),
        "cibersortx_upload_file": upload_path.name,
        "cibersortx_upload_sha256": sha256(upload_path),
        "cibersortx_upload_scale": "per-cell CPM from raw non-log UMI counts",
        "cibersortx_upload_header": "repeated sanitized Broad cell-type phenotype labels",
        "thresholds": {
            "min_cells_per_donor_type": args.min_cells_donor_type,
            "min_donors_per_type": args.min_donors_type,
            "max_donors_per_type": args.max_donors_type,
            "max_cells_per_donor_type": args.max_cells_donor_type,
            "max_total_cells": args.max_total_cells,
        },
    }
    (output_dir / "reference_manifest.json").write_text(json.dumps(run_manifest, indent=2), encoding="utf-8")
    logging.info("Wrote %s (%d genes x %d cells)", reference_path, len(genes), len(selected))
    logging.info("Cell types: %s", sorted(selected_obs.cell_type.unique()))


if __name__ == "__main__":
    main()
