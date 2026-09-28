#!/usr/bin/env python
"""Export GTEx bulk-skin raw read counts as a CIBERSORTx CPM mixture file.

The output is intended for manual upload to CIBERSORTx Fractions after the
matching single-cell file has been used to create a custom signature matrix.
Only genes present in that single-cell reference are exported, but CPM
denominators use every row in each original raw-count GCT.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd


SITES = {
    "Skin - Sun Exposed (Lower leg)": "gene_reads_adult_gtex_v11_skin_sun_exposed_lower_leg.gct.gz",
    "Skin - Not Sun Exposed (Suprapubic)": "gene_reads_adult_gtex_v11_skin_not_sun_exposed_suprapubic.gct.gz",
}


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=root)
    parser.add_argument("--reference", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--chunk-size", type=int, default=1000)
    return parser.parse_args()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reference_genes(path: Path) -> list[str]:
    genes = pd.read_csv(path, sep="\t", usecols=[0])["GeneSymbol"].astype(str)
    if genes.duplicated().any():
        raise ValueError("Reference GeneSymbol column is not unique")
    return genes.tolist()


def gct_to_reference_counts(path: Path, wanted: set[str], chunk_size: int) -> tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    library_sizes: pd.Series | None = None
    retained: list[pd.DataFrame] = []
    n_rows = 0
    n_symbol_rows = 0
    sample_columns: list[str] | None = None

    for chunk in pd.read_csv(path, sep="\t", skiprows=2, chunksize=chunk_size, low_memory=False):
        if chunk.shape[1] < 3:
            raise ValueError(f"Unexpected GCT structure in {path}")
        id_col, symbol_col = chunk.columns[:2]
        current_samples = chunk.columns[2:].astype(str).tolist()
        if sample_columns is None:
            sample_columns = current_samples
            library_sizes = pd.Series(0.0, index=sample_columns)
        elif current_samples != sample_columns:
            raise ValueError(f"Sample columns changed while reading {path}")

        values = chunk.iloc[:, 2:].apply(pd.to_numeric, errors="raise")
        if (values.to_numpy() < 0).any():
            raise ValueError(f"Negative raw counts found in {path}")
        library_sizes = library_sizes.add(values.sum(axis=0), fill_value=0)
        n_rows += len(chunk)

        symbols = chunk[symbol_col].astype("string").str.strip().str.upper()
        valid = symbols.notna() & symbols.ne("") & symbols.ne("NA") & symbols.ne("NAN")
        valid &= ~symbols.str.startswith("MT-")
        n_symbol_rows += int(valid.sum())
        keep = valid & symbols.isin(wanted)
        if keep.any():
            block = values.loc[keep].copy()
            block.index = symbols.loc[keep].astype(str)
            retained.append(block.groupby(level=0, sort=False).sum())

    if sample_columns is None or library_sizes is None or not retained:
        raise ValueError(f"No usable reference-overlap rows found in {path}")
    if (library_sizes <= 0).any():
        raise ValueError(f"Zero-count bulk libraries found in {path}")
    counts = pd.concat(retained).groupby(level=0, sort=False).sum()
    audit = pd.DataFrame({
        "sample_id": sample_columns,
        "source_file": path.name,
        "raw_library_size": library_sizes.loc[sample_columns].to_numpy(),
        "n_gct_rows": n_rows,
        "n_rows_with_symbol": n_symbol_rows,
        "n_reference_genes_exported": len(counts),
    })
    return counts, library_sizes, audit


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    data_dir = root / "data" / "external"
    output_dir = (args.output_dir or root / "results" / "public_validation" /
                  "gtex_skin_cibersortx" / "website_inputs").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    reference = (args.reference or root / "results" / "public_validation" /
                 "gtex_skin_cibersortx" / "reference" /
                 "gtex_skin_cibersortx_single_cell_reference.tsv").resolve()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    genes = reference_genes(reference)
    wanted = set(genes)
    site_matrices: list[pd.DataFrame] = []
    raw_site_matrices: list[pd.DataFrame] = []
    sample_audits: list[pd.DataFrame] = []
    for site, filename in SITES.items():
        path = data_dir / filename
        logging.info("Reading raw counts: %s", path)
        counts, library_sizes, audit = gct_to_reference_counts(path, wanted, args.chunk_size)
        cpm = counts.div(library_sizes, axis=1).mul(1_000_000).astype(np.float32)
        audit.insert(1, "skin_site", site)
        site_matrices.append(cpm)
        raw_site_matrices.append(counts)
        sample_audits.append(audit)

    duplicated_samples = pd.Index(np.concatenate([x.columns for x in site_matrices])).duplicated()
    if duplicated_samples.any():
        raise ValueError("A GTEx sample ID occurs in more than one skin-site GCT")
    mixture = pd.concat(site_matrices, axis=1).reindex(genes).fillna(0.0)
    mixture.index.name = "GeneSymbol"
    output = output_dir / "gtex_bulk_skin_cibersortx_mixture_cpm.tsv"
    mixture.to_csv(output, sep="\t", float_format="%.6f")
    raw_counts = pd.concat(raw_site_matrices, axis=1).reindex(genes).fillna(0)
    raw_counts.index.name = "GeneSymbol"
    raw_output = output_dir / "gtex_bulk_skin_music_raw_counts.tsv"
    raw_counts.to_csv(raw_output, sep="\t", float_format="%.0f")
    sample_audit = pd.concat(sample_audits, ignore_index=True)
    sample_audit.to_csv(output_dir / "bulk_sample_manifest.csv", index=False)

    manifest = {
        "purpose": "manual CIBERSORTx website upload as mixture file",
        "mixture_file": output.name,
        "mixture_sha256": file_sha256(output),
        "music_raw_counts_file": raw_output.name,
        "music_raw_counts_sha256": file_sha256(raw_output),
        "scale": "CPM computed from raw non-log GTEx read counts",
        "cpm_denominator": "all rows in each original GCT",
        "gene_policy": "uppercase symbols; duplicate symbols summed; missing and MT- symbols dropped; reference genes retained",
        "reference_file": str(reference),
        "reference_sha256": file_sha256(reference),
        "n_genes": int(mixture.shape[0]),
        "n_samples": int(mixture.shape[1]),
        "sites": SITES,
    }
    (output_dir / "bulk_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    logging.info("Wrote %s (%d genes x %d samples)", output, *mixture.shape)
    logging.info("Wrote %s for MuSiC", raw_output)


if __name__ == "__main__":
    main()
