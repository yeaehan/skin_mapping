"""
Build data/Xenium/xenium_full.h5ad from the per-sample Xenium outputs
plus the rds-derived cell metadata CSV.

For each sample (folder `data/Xenium/<sample>`):
  - counts: cell_feature_matrix/{barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz}  (genes x cells)
  - coords: cells.parquet  (x_centroid, y_centroid)

Then merge in `results/260727/xenium_rds_metadata/xenium_rds_cell_metadata.csv.gz`
using key `<sample>_<cell_id>` -> adds snc, age, cell_type_final, etc.

Output: data/Xenium/xenium_full.h5ad  (concatenated AnnData, var = gene symbols, obsm['spatial'])
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-7s  %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("build_xenium_full")

ROOT = Path(__file__).resolve().parents[3]                       # skin_mapping_2/
XEN_DIR = ROOT / "data" / "Xenium"
OUT = XEN_DIR / "xenium_full.h5ad"
META_CSV = ROOT / "results" / "260727" / "xenium_rds_metadata" / "xenium_rds_cell_metadata.csv.gz"

SAMPLES = ["1738-Y", "1883-Y", "4890-O", "4891-O", "5138-O", "7503-Y"]


def load_sample(sample: str) -> ad.AnnData:
    sdir = XEN_DIR / sample
    cfm = sdir / "cell_feature_matrix"
    if not cfm.is_dir():
        raise FileNotFoundError(f"Missing cell_feature_matrix for {sample}: {cfm}")

    log.info("[%s] reading barcodes/features/matrix …", sample)
    barcodes = pd.read_csv(cfm / "barcodes.tsv.gz", header=None, sep="\t")[0].astype(str).values
    features = pd.read_csv(cfm / "features.tsv.gz", header=None, sep="\t", names=["ensembl", "symbol", "feature_type"])
    mtx = sio.mmread(str(cfm / "matrix.mtx.gz")).T.tocsr()        # -> cells x genes
    log.info("[%s] mtx shape after transpose: %s", sample, mtx.shape)

    if mtx.shape != (len(barcodes), len(features)):
        raise ValueError(f"[{sample}] shape mismatch mtx={mtx.shape} vs barcodes/features {(len(barcodes), len(features))}")

    coords_df = pd.read_parquet(sdir / "cells.parquet")[["cell_id", "x_centroid", "y_centroid"]]
    if not np.array_equal(barcodes, coords_df["cell_id"].astype(str).values):
        log.warning("[%s] barcode order != cells.parquet order; reordering coords", sample)
        coords_df = coords_df.set_index("cell_id").loc[barcodes].reset_index()

    obs = pd.DataFrame(
        {
            "cell_id": barcodes,
            "sample": sample,
            "x_centroid": coords_df["x_centroid"].to_numpy(),
            "y_centroid": coords_df["y_centroid"].to_numpy(),
        }
    )
    obs.index = [f"{sample}_{b}" for b in barcodes]
    obs["obs_names"] = obs.index

    var = pd.DataFrame(
        {
            "ensembl": features["ensembl"].astype(str).values,
            "symbol": features["symbol"].astype(str).values,
            "feature_type": features["feature_type"].astype(str).values,
        },
    )
    # Use gene symbols as the var index (Xenium convention; we're keeping ensembl too)
    var = var.set_index("symbol", drop=False)
    var.index.name = "symbol_index"
    if var.index.has_duplicates:
        # disambiguate rare dup symbols
        var["symbol"] = _dedup_index(var["symbol"].tolist())
        var.index = pd.Index(var["symbol"], name="symbol_index")

    a = ad.AnnData(
        X=mtx,
        obs=obs,
        var=var,
    )
    a.obsm["spatial"] = obs[["x_centroid", "y_centroid"]].to_numpy()
    return a


def _dedup_index(names: list[str]) -> list[str]:
    seen = {}
    out = []
    for n in names:
        if n not in seen:
            seen[n] = 0
            out.append(n)
        else:
            seen[n] += 1
            out.append(f"{n}.{seen[n]}")
    return out


def main():
    log.info("=" * 60)
    log.info("Building xenium_full.h5ad from %d samples", len(SAMPLES))
    log.info("=" * 60)

    pieces = []
    for s in tqdm(SAMPLES, desc="Loading samples"):
        pieces.append(load_sample(s))

    log.info("Concatenating %d AnnData objects …", len(pieces))
    adata = ad.concat(
        pieces,
        axis=0,
        join="outer",                 # keep all 541 genes (they should match, but safe)
        label="sample_concat",
        keys=SAMPLES,
        index_unique=None,
    )
    log.info("  concat shape: %s", adata.shape)

    # Inner-join to genes present in every sample (Xenium panel is constant, but just in case)
    if adata.var.shape[0] != pieces[0].var.shape[0]:
        n_per = [a.var.shape[0] for a in pieces]
        log.warning("Gene count differs across samples: %s", n_per)

    # Ensure raw counts are stored (no normalization yet); X is the integer matrix
    if sp.issparse(adata.X):
        adata.X = adata.X.tocsr()
        adata.X.eliminate_zeros()

    log.info("Loading rds-derived metadata: %s", META_CSV)
    meta = pd.read_csv(META_CSV)
    log.info("  metadata shape: %s", meta.shape)

    # CSV cell_id is already prefixed with sample (e.g. "1738-Y_aagkcool-1"),
    # while obs.cell_id from the upstream parquet is the bare barcode (e.g. "aagkcool-1").
    # Build a composite key on the obs side so we filter by the same identifier.
    adata.obs_names = [f"{s}_{b}" for s, b in zip(adata.obs["sample"].astype(str).values,
                                                  adata.obs["cell_id"].astype(str).values)]
    keep_cells = set(meta["cell_id"].astype(str).tolist())
    n_before = adata.shape[0]
    mask = adata.obs_names.isin(keep_cells)
    log.info("Filtering to cells present in metadata CSV: %d / %d (drop %d)", int(mask.sum()), n_before, int((~mask).sum()))
    adata = adata[mask].copy()
    log.info("  after filter: %s", adata.shape)

    log.info("Joining metadata onto obs by composite cell_id …")
    meta_idx = meta.set_index("cell_id")
    adata.obs["cell_id"] = adata.obs_names                                    # use composite key for join
    obs_merged = adata.obs.join(meta_idx, how="left", rsuffix="_csv")

    # Drop pre-existing columns we want to overwrite from CSV (x_centroid is fine to keep from upstream)
    overwrite_cols = [c for c in meta_idx.columns if c in adata.obs.columns and c != "x_centroid"]
    if overwrite_cols:
        obs_merged = obs_merged.drop(columns=overwrite_cols)

    n_unmatched = int(obs_merged["snc"].isna().sum()) if "snc" in obs_merged.columns else None
    log.info("  unmatched after join (no metadata row): %d / %d", n_unmatched, len(obs_merged))

    # Convert string-y categorical columns; coerce mixed dtypes to str for safe h5ad write
    for col in obs_merged.columns:
        s = obs_merged[col]
        if s.dtype == "object":
            obs_merged[col] = s.astype(str)

    for col in ["snc", "age", "cell_type_minor", "cell_type_major", "cell_type_final",
                "segmentation_method", "cell_lineage", "majorskincell", "minorskincell",
                "age_major", "age_minor", "celltype_age", "celltype_snc", "celltype_snc1",
                "cd45snc", "immune_cell", "immune_subclusters", "immune_subclusters1"]:
        if col in obs_merged.columns:
            obs_merged[col] = obs_merged[col].astype("category")

    adata.obs = obs_merged
    adata.obs.index.name = "cell_id"

    log.info("Saving: %s", OUT)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT, compression="gzip")
    log.info("Done. shape = %s", adata.shape)


if __name__ == "__main__":
    main()
