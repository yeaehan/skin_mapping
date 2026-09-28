"""
Merge public skin scRNA references at the raw-count level and run Harmony.

Inputs are the local processed/count-containing objects already present in this
workspace:

- Restrepo / NS-Atlas scRNA h5ad from Zenodo
- Tabula Sapiens skin h5ad, excluding cells with undetermined anatomical site
- Optional: GTEx 8-tissue snRNA h5ad, filtered to skin suprapubic by default

The script deliberately does not assign new harmonized cell-type labels. It
keeps original annotations as metadata so annotation can be revisited after
inspecting the integrated clusters.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import h5py


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESTREPO_H5AD = PROJECT_ROOT / "rsc" / "Restrepo_NatGen_Merfish_Skin" / "16795569" / "scrna.integrated_annotated.h5ad"
GTEX_H5AD = PROJECT_ROOT / "data" / "external" / "GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad"
TABULA_H5AD = PROJECT_ROOT / "data" / "external" / "Skin_TSP1_30_version2d_10X_smartseq_scvi_Nov122024_updated.h5ad"
OUTDIR = PROJECT_ROOT / "results" / "260815" / "public_reference_merge_harmony"

UNDETERMINED_SITE_VALUES = {
    "",
    "na",
    "n/a",
    "nan",
    "none",
    "not available",
    "not_applicable",
    "not applicable",
    "unknown",
    "undetermined",
    "unspecified",
}

ad = None
np = None
pd = None
sc = None
sparse = None


def import_analysis_stack() -> None:
    global ad, np, pd, sc, sparse
    if ad is not None:
        return
    import anndata as _ad
    import numpy as _np
    import pandas as _pd
    import scanpy as _sc
    from scipy import sparse as _sparse

    ad = _ad
    np = _np
    pd = _pd
    sc = _sc
    sparse = _sparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge Restrepo and Tabula Sapiens skin raw counts, then run Harmony."
    )
    parser.add_argument("--restrepo-h5ad", type=Path, default=RESTREPO_H5AD)
    parser.add_argument("--gtex-h5ad", type=Path, default=GTEX_H5AD)
    parser.add_argument("--tabula-h5ad", type=Path, default=TABULA_H5AD)
    parser.add_argument("--outdir", type=Path, default=OUTDIR)
    parser.add_argument(
        "--include-gtex",
        action="store_true",
        help="Also merge GTEx skin snRNA-seq. Off by default to avoid mixing scRNA-seq and snRNA-seq technologies.",
    )
    parser.add_argument(
        "--gtex-sites",
        choices=("suprapubic", "lower_leg", "all_skin"),
        default="suprapubic",
        help="GTEx skin subset to merge when --include-gtex is used.",
    )
    parser.add_argument(
        "--exclude-tabula-buccal",
        action="store_true",
        help="Exclude Tabula Sapiens buccal/oral mucosa cells in addition to undetermined sites.",
    )
    parser.add_argument(
        "--join",
        choices=("inner", "outer"),
        default="inner",
        help="Gene merge mode. Inner is more conservative for integration and avoids a very sparse union matrix.",
    )
    parser.add_argument("--harmony-batch-key", default="study_id")
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--neighbors-pcs", type=int, default=40)
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument("--max-cells-per-source", type=int, default=None)
    parser.add_argument("--random-state", type=int, default=260815)
    parser.add_argument(
        "--skip-harmony",
        action="store_true",
        help="Write the merged raw-count object and summary, then stop before Harmony integration.",
    )
    return parser.parse_args()


def require_file(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path}")


def as_csr(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(np.asarray(matrix))


def first_existing(df: pd.DataFrame, columns: Iterable[str]) -> str | None:
    for col in columns:
        if col in df.columns:
            return col
    return None


def get_var_symbols(adata: ad.AnnData, preferred: str | None = None) -> pd.Series:
    candidates = [preferred, "gene_symbol", "gene_name", "Approved symbol", "feature_name"]
    col = first_existing(adata.var, [c for c in candidates if c])
    if col:
        symbols = adata.var[col].astype(str)
    else:
        symbols = pd.Series(adata.var_names.astype(str), index=adata.var_names)
    symbols.index = adata.var_names
    return symbols.str.upper().str.strip()


def get_var_symbols_from_var(var: pd.DataFrame, preferred: str | None = None) -> pd.Series:
    candidates = [preferred, "gene_symbol", "gene_name", "Approved symbol", "feature_name"]
    col = first_existing(var, [c for c in candidates if c])
    if col:
        symbols = var[col].astype(str)
    else:
        symbols = pd.Series(var.index.astype(str), index=var.index)
    symbols.index = var.index
    return symbols.str.upper().str.strip()


def decode_h5_values(values) -> np.ndarray:
    values = np.asarray(values)
    if values.dtype.kind in {"S", "O"}:
        return np.array(
            [x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else x for x in values],
            dtype=object,
        )
    return values


def read_h5ad_column(obj) -> np.ndarray:
    if isinstance(obj, h5py.Group) and {"categories", "codes"}.issubset(set(obj.keys())):
        categories = decode_h5_values(obj["categories"][:])
        codes = np.asarray(obj["codes"][:], dtype=int)
        out = np.full(codes.shape, pd.NA, dtype=object)
        ok = codes >= 0
        out[ok] = categories[codes[ok]]
        return out
    return decode_h5_values(obj[:])


def read_h5ad_dataframe(group: h5py.Group) -> pd.DataFrame:
    if "_index" in group:
        index = read_h5ad_column(group["_index"])
    else:
        first_key = next(iter(group.keys()))
        index = np.arange(len(read_h5ad_column(group[first_key])))
    n_rows = len(index)
    df = pd.DataFrame(index=pd.Index(index.astype(str)))
    for key in group.keys():
        if key == "_index":
            continue
        values = read_h5ad_column(group[key])
        if len(values) == n_rows:
            df[key] = values
    return df


def choose_positions(n_obs: int, max_cells: int | None, random_state: int, mask=None) -> np.ndarray:
    if mask is None:
        positions = np.arange(n_obs)
    else:
        positions = np.flatnonzero(np.asarray(mask, dtype=bool))
    if max_cells is not None and len(positions) > max_cells:
        rng = np.random.default_rng(random_state)
        positions = rng.choice(positions, size=max_cells, replace=False)
    return np.sort(positions.astype(np.int64, copy=False))


def read_matrix_rows(matrix_obj, row_positions: np.ndarray) -> sparse.csr_matrix:
    if isinstance(matrix_obj, h5py.Group) and {"data", "indices", "indptr"}.issubset(set(matrix_obj.keys())):
        shape = tuple(int(x) for x in matrix_obj.attrs["shape"])
        n_vars = shape[1]
        old_indptr = matrix_obj["indptr"][:]
        new_indptr = np.zeros(len(row_positions) + 1, dtype=old_indptr.dtype)
        data_chunks = []
        index_chunks = []
        nnz = 0
        for out_i, row in enumerate(row_positions):
            start = int(old_indptr[row])
            end = int(old_indptr[row + 1])
            if end > start:
                data_chunks.append(matrix_obj["data"][start:end])
                index_chunks.append(matrix_obj["indices"][start:end])
                nnz += end - start
            new_indptr[out_i + 1] = nnz
        data = np.concatenate(data_chunks) if data_chunks else np.array([], dtype=matrix_obj["data"].dtype)
        indices = np.concatenate(index_chunks) if index_chunks else np.array([], dtype=matrix_obj["indices"].dtype)
        return sparse.csr_matrix((data, indices, new_indptr), shape=(len(row_positions), n_vars))

    dense = matrix_obj[row_positions, :]
    return sparse.csr_matrix(dense)


def make_count_adata_from_matrix(
    counts,
    obs: pd.DataFrame,
    var: pd.DataFrame,
    var_symbol_col: str | None,
) -> ad.AnnData:
    symbols = get_var_symbols_from_var(var, var_symbol_col)
    valid = ~symbols.isna() & ~symbols.isin(["", "NAN", "NONE"])
    symbols = symbols.loc[valid]
    counts = as_csr(counts)[:, np.asarray(valid)]

    codes, unique_symbols = pd.factorize(symbols, sort=True)
    gene_map = sparse.csr_matrix(
        (np.ones(len(codes), dtype=np.float32), (np.arange(len(codes)), codes)),
        shape=(len(codes), len(unique_symbols)),
    )
    counts = (counts @ gene_map).tocsr()

    out = ad.AnnData(
        X=counts,
        obs=obs.copy(),
        var=pd.DataFrame(index=pd.Index(unique_symbols, name="gene_symbol")),
    )
    out.var_names_make_unique()
    return out


def make_count_adata(
    adata: ad.AnnData,
    layer: str,
    obs: pd.DataFrame,
    var_symbol_col: str | None,
) -> ad.AnnData:
    if layer and layer in adata.layers:
        counts = adata.layers[layer]
    elif layer:
        raise KeyError(f"Layer {layer!r} not found. Available layers: {list(adata.layers.keys())}")
    else:
        counts = adata.X

    symbols = get_var_symbols(adata, var_symbol_col)
    valid = ~symbols.isna() & ~symbols.isin(["", "NAN", "NONE"])
    symbols = symbols.loc[valid]
    counts = as_csr(counts)[:, np.asarray(valid)]

    codes, unique_symbols = pd.factorize(symbols, sort=True)
    gene_map = sparse.csr_matrix(
        (np.ones(len(codes), dtype=np.float32), (np.arange(len(codes)), codes)),
        shape=(len(codes), len(unique_symbols)),
    )
    counts = counts @ gene_map
    counts = counts.tocsr()

    out = ad.AnnData(
        X=counts,
        obs=obs.copy(),
        var=pd.DataFrame(index=pd.Index(unique_symbols, name="gene_symbol")),
    )
    out.var_names_make_unique()
    return out


def standard_obs(index: pd.Index, dataset_source: str) -> pd.DataFrame:
    obs = pd.DataFrame(index=index.copy())
    obs["dataset_source"] = dataset_source
    obs["source_cell_id"] = index.astype(str)
    return obs


def clean_string(series: pd.Series) -> pd.Series:
    return series.astype(str).replace({"nan": np.nan, "None": np.nan})


def load_restrepo(path: Path, max_cells: int | None, random_state: int) -> ad.AnnData:
    print(f"Reading Restrepo counts: {path}")
    with h5py.File(path, "r") as h5:
        src_obs = read_h5ad_dataframe(h5["obs"])
        src_var = read_h5ad_dataframe(h5["var"])
        positions = choose_positions(src_obs.shape[0], max_cells, random_state)
        print(f"Retaining {len(positions):,} Restrepo cells")
        counts = read_matrix_rows(h5["layers"]["counts"], positions)
        src_obs = src_obs.iloc[positions].copy()

    obs = standard_obs(src_obs.index, "Restrepo_NSAtlas")
    obs["study_id"] = clean_string(src_obs["study_id"])
    obs["donor_id"] = clean_string(src_obs["donor_id"])
    obs["sample_id"] = clean_string(src_obs["sample_barcode"])
    obs["anatomic_site"] = clean_string(src_obs["anatomic_site"])
    obs["anatomic_site_detail"] = clean_string(src_obs["anatomic_site.detailed"])
    obs["donor_age"] = pd.to_numeric(src_obs.get("donor_age"), errors="coerce")
    obs["donor_sex"] = clean_string(src_obs.get("donor_sex", pd.Series(index=src_obs.index, dtype=object)))
    obs["assay"] = "public_scRNA"
    obs["original_broad_cell_type"] = clean_string(src_obs.get("cell_type.broad", src_obs.get("cell_type")))
    obs["original_cell_type"] = clean_string(src_obs.get("cell_type.detailed", src_obs.get("reported.cell_type")))
    return make_count_adata_from_matrix(counts, obs, src_var, var_symbol_col=None)


def gtex_site_mask(obs: pd.DataFrame, gtex_sites: str) -> pd.Series:
    site_col = first_existing(obs, ["Tissue Site Detail", "Tissue", "tissue"])
    if site_col is None:
        raise KeyError("Could not find a GTEx tissue/site column.")
    site = obs[site_col].astype(str).str.lower()
    is_skin = site.str.contains("skin", na=False)
    if gtex_sites == "suprapubic":
        return is_skin & site.str.contains("suprapubic|not sun exposed", regex=True, na=False)
    if gtex_sites == "lower_leg":
        return is_skin & site.str.contains("lower leg|sun exposed", regex=True, na=False)
    return is_skin


def load_gtex(path: Path, gtex_sites: str, max_cells: int | None, random_state: int) -> ad.AnnData:
    print(f"Reading GTEx counts: {path}")
    with h5py.File(path, "r") as h5:
        src_obs = read_h5ad_dataframe(h5["obs"])
        src_var = read_h5ad_dataframe(h5["var"])
        mask = gtex_site_mask(src_obs, gtex_sites)
        positions = choose_positions(src_obs.shape[0], max_cells, random_state, mask=mask)
        if len(positions) == 0:
            available = sorted(src_obs[first_existing(src_obs, ["Tissue Site Detail", "Tissue", "tissue"])].astype(str).unique())
            raise ValueError(
                f"GTEx filter {gtex_sites!r} retained zero cells. Available site values include: {available[:30]}"
            )
        print(f"Retaining {len(positions):,} GTEx cells")
        counts = read_matrix_rows(h5["layers"]["counts"], positions)
        src_obs = src_obs.iloc[positions].copy()

    if int(mask.sum()) == 0:
        available = sorted(src_obs[first_existing(src_obs, ["Tissue Site Detail", "Tissue", "tissue"])].astype(str).unique())
        raise ValueError(
            f"GTEx filter {gtex_sites!r} retained zero cells. Available site values include: {available[:30]}"
        )

    obs = standard_obs(src_obs.index, f"GTEx_skin_{gtex_sites}")
    obs["study_id"] = "GTEx_snRNAseq_atlas"
    obs["donor_id"] = clean_string(src_obs["Participant ID"])
    obs["sample_id"] = clean_string(src_obs.get("Sample ID short", src_obs.get("Sample ID")))
    obs["anatomic_site"] = clean_string(src_obs["Tissue Site Detail"])
    obs["anatomic_site_detail"] = obs["anatomic_site"]
    obs["donor_age"] = clean_string(src_obs.get("Age_bin", pd.Series(index=src_obs.index, dtype=object)))
    obs["donor_sex"] = clean_string(src_obs.get("Sex", pd.Series(index=src_obs.index, dtype=object)))
    obs["assay"] = clean_string(src_obs.get("prep", pd.Series("snRNA", index=src_obs.index)))
    obs["original_broad_cell_type"] = clean_string(src_obs.get("Broad cell type"))
    obs["original_cell_type"] = clean_string(src_obs.get("Granular cell type", src_obs.get("Cell types level 3")))
    return make_count_adata_from_matrix(counts, obs, src_var, var_symbol_col="gene_name")


def tabula_site_mask(obs: pd.DataFrame, exclude_buccal: bool) -> pd.Series:
    if "tissue" not in obs.columns or "anatomical_position" not in obs.columns:
        raise KeyError("Tabula object must contain 'tissue' and 'anatomical_position' columns.")
    tissue = obs["tissue"].astype(str).str.lower()
    site = obs["anatomical_position"].astype(str).str.strip()
    mask = tissue.eq("skin") & ~site.str.lower().isin(UNDETERMINED_SITE_VALUES)
    if exclude_buccal:
        mask &= ~site.str.lower().str.contains("buccal|oral", regex=True, na=False)
    return mask


def load_tabula(path: Path, exclude_buccal: bool, max_cells: int | None, random_state: int) -> ad.AnnData:
    print(f"Reading Tabula Sapiens counts: {path}")
    with h5py.File(path, "r") as h5:
        src_obs = read_h5ad_dataframe(h5["obs"])
        src_var = read_h5ad_dataframe(h5["var"])
        mask = tabula_site_mask(src_obs, exclude_buccal)
        positions = choose_positions(src_obs.shape[0], max_cells, random_state, mask=mask)
        if len(positions) == 0:
            raise ValueError("Tabula site filter retained zero cells.")
        print(f"Retaining {len(positions):,} Tabula Sapiens cells")
        counts = read_matrix_rows(h5["layers"]["raw_counts"], positions)
        src_obs = src_obs.iloc[positions].copy()

    obs = standard_obs(src_obs.index, "TabulaSapiens_skin_site_known")
    obs["study_id"] = "TabulaSapiens"
    obs["donor_id"] = clean_string(src_obs["donor"])
    obs["sample_id"] = clean_string(src_obs["sample_id"])
    obs["anatomic_site"] = clean_string(src_obs["anatomical_position"])
    obs["anatomic_site_detail"] = obs["anatomic_site"]
    obs["donor_age"] = pd.to_numeric(src_obs.get("age"), errors="coerce")
    obs["donor_sex"] = clean_string(src_obs.get("sex", pd.Series(index=src_obs.index, dtype=object)))
    obs["assay"] = clean_string(src_obs.get("method", src_obs.get("assay")))
    obs["original_broad_cell_type"] = clean_string(src_obs.get("broad_cell_class", src_obs.get("compartment")))
    obs["original_cell_type"] = clean_string(src_obs.get("cell_ontology_class", src_obs.get("free_annotation")))
    return make_count_adata_from_matrix(counts, obs, src_var, var_symbol_col="gene_symbol")


def downsample_per_source(adata: ad.AnnData, max_cells: int | None, random_state: int) -> ad.AnnData:
    # Main loaders already sample before reading count layers into memory. This
    # remains as a guard for any future source loader that returns a full object.
    if max_cells is None:
        return adata
    rng = np.random.default_rng(random_state)
    keep = []
    for _, idx in adata.obs.groupby("dataset_source", observed=True).indices.items():
        idx = np.asarray(idx)
        if len(idx) > max_cells:
            idx = rng.choice(idx, size=max_cells, replace=False)
        keep.append(idx)
    keep_idx = np.sort(np.concatenate(keep))
    return adata[keep_idx].copy()


def add_qc(adata: ad.AnnData) -> None:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPL", "RPS"))
    adata.var["hb"] = adata.var_names.str.startswith(("HB",))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, inplace=True)


def run_harmony(adata: ad.AnnData, args: argparse.Namespace) -> ad.AnnData:
    work = adata.copy()
    work.layers["counts"] = work.X.copy()

    sc.pp.normalize_total(work, target_sum=args.target_sum)
    sc.pp.log1p(work)
    sc.pp.highly_variable_genes(
        work,
        n_top_genes=args.n_top_genes,
        flavor="seurat",
        batch_key=args.harmony_batch_key,
    )
    hvg = work[:, work.var["highly_variable"].to_numpy()].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=args.n_pcs, svd_solver="arpack", random_state=args.random_state)

    try:
        import scanpy.external as sce
    except ImportError as exc:
        raise ImportError("scanpy.external is required for Harmony integration.") from exc

    try:
        sce.pp.harmony_integrate(
            hvg,
            key=args.harmony_batch_key,
            basis="X_pca",
            adjusted_basis="X_pca_harmony",
        )
    except ImportError as exc:
        raise ImportError("Harmony integration requires harmonypy. Install it with: pip install harmonypy") from exc

    n_pcs = min(args.neighbors_pcs, hvg.obsm["X_pca_harmony"].shape[1])
    sc.pp.neighbors(hvg, use_rep="X_pca_harmony", n_pcs=n_pcs)
    sc.tl.umap(hvg, random_state=args.random_state)
    hvg.obsm["X_umap_harmony"] = hvg.obsm.pop("X_umap")

    try:
        for res in (0.2, 0.4, 0.6, 0.8, 1.0):
            sc.tl.leiden(hvg, resolution=res, key_added=f"leiden_harmony_{res:g}", random_state=args.random_state)
    except ImportError:
        print("leidenalg/igraph not available; skipping Leiden clustering.")

    return hvg


def write_summary(adata: ad.AnnData, path: Path) -> None:
    rows = []
    for dataset, obs in adata.obs.groupby("dataset_source", observed=True):
        rows.append(
            {
                "dataset_source": dataset,
                "n_cells": int(obs.shape[0]),
                "n_studies": int(obs["study_id"].astype(str).nunique()),
                "n_donors": int(obs["donor_id"].astype(str).nunique()),
                "sites": "; ".join(sorted(obs["anatomic_site"].astype(str).unique())),
            }
        )
    pd.DataFrame(rows).to_csv(path, index=False)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for path in (args.restrepo_h5ad, args.tabula_h5ad):
        require_file(path)
    if args.include_gtex:
        require_file(args.gtex_h5ad)

    import_analysis_stack()

    datasets = [
        load_restrepo(args.restrepo_h5ad, args.max_cells_per_source, args.random_state),
        load_tabula(args.tabula_h5ad, args.exclude_tabula_buccal, args.max_cells_per_source, args.random_state),
    ]
    if args.include_gtex:
        datasets.insert(1, load_gtex(args.gtex_h5ad, args.gtex_sites, args.max_cells_per_source, args.random_state))

    merged = ad.concat(datasets, join=args.join, merge="same", index_unique=None)
    merged.obs_names = merged.obs["dataset_source"].astype(str) + ":" + merged.obs["source_cell_id"].astype(str)
    merged.obs_names_make_unique()
    merged = downsample_per_source(merged, args.max_cells_per_source, args.random_state)
    add_qc(merged)

    prefix = "public_skin_restrepo_tabula"
    if args.include_gtex:
        prefix = f"{prefix}_gtex_{args.gtex_sites}"
    if args.max_cells_per_source is not None:
        prefix = f"{prefix}.max{args.max_cells_per_source}_per_source"

    raw_path = args.outdir / f"{prefix}.raw_counts_merged.h5ad"
    integrated_path = args.outdir / f"{prefix}.harmony_integrated_hvg.h5ad"
    summary_path = args.outdir / f"{prefix}.merge_summary.csv"
    config_path = args.outdir / f"{prefix}.run_config.json"

    merged.write_h5ad(raw_path, compression="gzip")
    write_summary(merged, summary_path)

    run_config = {
        "restrepo_h5ad": str(args.restrepo_h5ad),
        "gtex_h5ad": str(args.gtex_h5ad),
        "tabula_h5ad": str(args.tabula_h5ad),
        "include_gtex": args.include_gtex,
        "gtex_sites": args.gtex_sites,
        "exclude_tabula_buccal": args.exclude_tabula_buccal,
        "join": args.join,
        "harmony_batch_key": args.harmony_batch_key,
        "n_top_genes": args.n_top_genes,
        "n_pcs": args.n_pcs,
        "neighbors_pcs": args.neighbors_pcs,
        "max_cells_per_source": args.max_cells_per_source,
        "skip_harmony": args.skip_harmony,
    }

    if args.skip_harmony:
        with config_path.open("w", encoding="utf-8") as f:
            json.dump(run_config, f, indent=2)
        print(f"Wrote raw merged counts: {raw_path}")
        print(f"Wrote summary: {summary_path}")
        print("Skipped Harmony integration by request.")
        return

    integrated = run_harmony(merged, args)
    integrated.uns["merge_harmony_config"] = run_config
    integrated.write_h5ad(integrated_path, compression="gzip")

    with config_path.open("w", encoding="utf-8") as f:
        json.dump(run_config, f, indent=2)

    print(f"Wrote raw merged counts: {raw_path}")
    print(f"Wrote Harmony integrated object: {integrated_path}")
    print(f"Wrote summary: {summary_path}")


if __name__ == "__main__":
    main()
