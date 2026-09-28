# %%
"""
Combined GEO public-validation SnC quartile scoring.

Datasets:
  - GSE130973: aging human skin scRNA-seq, processed filtered MTX files.
  - GSE202352: human hip/palm/sole skin scRNA-seq, processed Cell Ranger MEX archives.

This is a separate workflow from the Sinai notebook. It downloads processed GEO
files, combines count matrices, normalizes/log-transforms the combined object,
and applies the same continuous-axis senescence index idea:

    SI = Z_growth + Z_SASP - Z_cycle

Then it assigns SI quartiles and top 5% / top 10% calls within configurable
groups. If curated cell-type metadata are not supplied, all cells are scored
within study-level groups to avoid forcing unavailable lineage labels.

Important interpretation:
  These are SI-high cells, not harmonized binary SnC+ cells. Use them as a
  transcriptomic sensitivity analysis alongside rule-based SnC annotations.
"""

from __future__ import annotations

import gzip
import re
import tarfile
import urllib.request
from pathlib import Path
from typing import Iterable

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp

try:
    from IPython.display import display
except ImportError:
    def display(obj):
        print(obj)


# %%
# =========================
# Config
# =========================
PROJECT_ROOT = Path("..").resolve()
DATA_DIR = PROJECT_ROOT / "data" / "public_validation" / "geo_combined"
OUTDIR = PROJECT_ROOT / "results" / "public_validation" / "snc_quartile_scoring_geo_combined"
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUTDIR.mkdir(parents=True, exist_ok=True)

# Optional CSV with one row per cell if you later obtain author annotations.
# Expected columns can include: cell_id, sample_id, cell_type, cell_type.broad,
# cell_type.detailed. Leave as None to score without cell-type labels.
OPTIONAL_CELL_METADATA_CSV: Path | None = None

TOP_PERCENTILES = [0.95, 0.90]  # top 5% and top 10%
MIN_CELLS_PER_GROUP = 50
TARGET_SUM = 1e4

# When no cell-type labels are present, use study-level scoring groups. If you
# add metadata with cell_type/cell_type.broad, these groups become study + cell type.
SCORING_GROUP_CANDIDATES = ["study_id", "cell_type.broad", "cell_type", "cell_type.detailed"]

W_GROWTH = 1.0
W_SASP = 1.0
W_CYCLE = 1.0


GSE130973_URLS = {
    "matrix": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130973&format=file&file=GSE130973_matrix_filtered.mtx.gz",
    "genes": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130973&format=file&file=GSE130973_genes_filtered.tsv.gz",
    "barcodes": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130973&format=file&file=GSE130973_barcodes_filtered.tsv.gz",
}

GSE202352_RAW_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE202352&format=file"


GSE130973_SUFFIX_META = {
    "1": {"sample_id": "y1", "donor_id": "y1", "donor_age": 25, "age_group": "Young", "anatomic_site": "inguinoiliac region"},
    "2": {"sample_id": "y2", "donor_id": "y2", "donor_age": 27, "age_group": "Young", "anatomic_site": "inguinoiliac region"},
    "3": {"sample_id": "o1", "donor_id": "o1", "donor_age": 53, "age_group": "Old", "anatomic_site": "inguinoiliac region"},
    "4": {"sample_id": "o2", "donor_id": "o2", "donor_age": 70, "age_group": "Old", "anatomic_site": "inguinoiliac region"},
    "5": {"sample_id": "o3", "donor_id": "o3", "donor_age": 69, "age_group": "Old", "anatomic_site": "inguinoiliac region"},
}


GSE202352_SAMPLE_META = {
    "GSM6111844": {"sample_id": "Subject1_hip", "donor_id": "Subject 1", "anatomic_site": "hip", "compartment": "combined epidermis + dermis"},
    "GSM6111845": {"sample_id": "Subject1_sole", "donor_id": "Subject 1", "anatomic_site": "sole", "compartment": "combined epidermis + dermis"},
    "GSM6111846": {"sample_id": "Subject2_hip", "donor_id": "Subject 2", "anatomic_site": "hip", "compartment": "combined epidermis + dermis"},
    "GSM6111847": {"sample_id": "Subject2_sole", "donor_id": "Subject 2", "anatomic_site": "sole", "compartment": "combined epidermis + dermis"},
    "GSM6111848": {"sample_id": "Subject3_hip_dermis", "donor_id": "Subject 3", "anatomic_site": "hip", "compartment": "dermis"},
    "GSM6111849": {"sample_id": "Subject3_palm_dermis", "donor_id": "Subject 3", "anatomic_site": "palm", "compartment": "dermis"},
    "GSM6111850": {"sample_id": "Subject3_hip_epidermis", "donor_id": "Subject 3", "anatomic_site": "hip", "compartment": "epidermis"},
    "GSM6111851": {"sample_id": "Subject3_palm_epidermis", "donor_id": "Subject 3", "anatomic_site": "palm", "compartment": "epidermis"},
    "GSM6111852": {"sample_id": "Subject4_hip_dermis", "donor_id": "Subject 4", "anatomic_site": "hip", "compartment": "dermis"},
    "GSM6111853": {"sample_id": "Subject4_palm_dermis", "donor_id": "Subject 4", "anatomic_site": "palm", "compartment": "dermis"},
    "GSM6111854": {"sample_id": "Subject4_hip_epidermis", "donor_id": "Subject 4", "anatomic_site": "hip", "compartment": "epidermis"},
    "GSM6111855": {"sample_id": "Subject4_palm_epidermis", "donor_id": "Subject 4", "anatomic_site": "palm", "compartment": "epidermis"},
}


# %%
# =========================
# Gene sets
# =========================
GROWTH_GENES = ["CDKN1A", "CDKN2A"]

S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
    "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", "RFC2", "RPA2", "NASP",
    "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2",
    "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2",
    "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]

G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", "SMC4", "CCNB2",
    "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
    "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1",
    "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5",
    "CENPA",
]

SASP_GENES = [
    "CXCL1", "CXCL10", "CXCL11", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL5",
    "CXCL6", "CXCL8", "CCL2", "CCL20", "CCL28", "CCL7", "IL11", "IL12A", "IL15",
    "IL16", "IL17D", "IL18", "IL18BP", "IL1A", "IL1B", "IL23A", "IL24", "IL32",
    "IL33", "IL6", "IL6ST", "IL7", "MMP1", "MMP10", "MMP11", "MMP12", "MMP14",
    "MMP15", "MMP16", "MMP17", "MMP19", "MMP22", "MMP24", "MMP24OS", "MMP25-AS1",
    "MMP3", "TNFAIP1", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFAIP8L1",
    "TNFAIP8L3", "TNFSF10", "TNFSF12", "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15",
    "TNFSF18", "TNFSF4", "VEGFA", "VEGFB", "VEGFC", "EGF", "PDGFA", "PDGFB",
    "GDF11", "GDF15", "GDF5", "GDF6", "FGF1", "FGF10", "FGF13", "FGF14", "FGF2",
    "FGF5", "FGF7", "FGF9",
]


# %%
# =========================
# Helpers
# =========================
def download_file(url: str, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 0:
        print(f"Already downloaded: {path}")
        return path
    tmp = path.with_suffix(path.suffix + ".tmp")
    print(f"Downloading {url}\n  -> {path}")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(path)
    return path


def safe_extract_tar(tar_path: Path, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    root = outdir.resolve()
    with tarfile.open(tar_path) as tar:
        for member in tar.getmembers():
            target = (outdir / member.name).resolve()
            if not str(target).startswith(str(root)):
                raise RuntimeError(f"Refusing unsafe tar member: {member.name}")
        tar.extractall(outdir)


def read_tsv_gz(path: Path, header=None) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=header, compression="gzip")


def canonical_gene_name(name: str) -> str:
    return re.sub(r"-\d+$", "", str(name)).upper()


def make_var_names_unique(names: Iterable[str]) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for name in names:
        if name not in seen:
            seen[name] = 0
            out.append(name)
        else:
            seen[name] += 1
            out.append(f"{name}-{seen[name]}")
    return out


def build_gene_lookup(adata: ad.AnnData) -> dict[str, str]:
    lookup = {}
    for gene in adata.var_names:
        lookup.setdefault(canonical_gene_name(gene), gene)
    return lookup


def present_genes(genes: Iterable[str], lookup: dict[str, str]) -> list[str]:
    present = []
    for gene in genes:
        mapped = lookup.get(canonical_gene_name(gene))
        if mapped is not None and mapped not in present:
            present.append(mapped)
    return present


def get_expr_df(adata: ad.AnnData, genes: list[str], layer: str | None = "log1p") -> pd.DataFrame:
    if not genes:
        return pd.DataFrame(index=adata.obs_names)
    matrix = adata[:, genes].layers[layer] if layer else adata[:, genes].X
    if sp.issparse(matrix):
        matrix = matrix.toarray()
    return pd.DataFrame(np.asarray(matrix), index=adata.obs_names, columns=genes)


def module_score(adata: ad.AnnData, genes: list[str], layer: str | None = "log1p") -> pd.Series:
    expr = get_expr_df(adata, genes, layer=layer)
    if expr.shape[1] == 0:
        return pd.Series(np.nan, index=adata.obs_names)
    return expr.mean(axis=1)


def group_key(obs: pd.DataFrame, columns: list[str]) -> pd.Series:
    if not columns:
        return pd.Series("all_cells", index=obs.index)
    parts = [obs[col].astype("string").fillna("Unknown") for col in columns]
    key = parts[0]
    for part in parts[1:]:
        key = key + " | " + part
    return key


def zscore_within_group(values: pd.Series, groups: pd.Series, min_cells: int = 50) -> pd.Series:
    values = pd.Series(values, index=groups.index, dtype="float64")
    groups = groups.astype("string").fillna("Unknown")
    out = pd.Series(np.nan, index=values.index, dtype="float64")

    global_mean = values.mean(skipna=True)
    global_sd = values.std(skipna=True, ddof=0)
    if not np.isfinite(global_sd) or global_sd == 0:
        global_sd = 1.0

    for _, idx in groups.groupby(groups, observed=True).groups.items():
        idx = list(idx)
        subset = values.loc[idx]
        if len(idx) < min_cells:
            out.loc[idx] = (subset - global_mean) / global_sd
            continue
        mean = subset.mean(skipna=True)
        sd = subset.std(skipna=True, ddof=0)
        if not np.isfinite(sd) or sd == 0:
            out.loc[idx] = (subset - global_mean) / global_sd
        else:
            out.loc[idx] = (subset - mean) / sd
    return out


def percentile_call_within_group(
    values: pd.Series,
    groups: pd.Series,
    percentile: float,
    min_cells: int = 50,
) -> tuple[pd.Series, pd.Series]:
    values = pd.Series(values, index=groups.index, dtype="float64")
    groups = groups.astype("string").fillna("Unknown")
    out = pd.Series(False, index=values.index)
    thresholds = {}
    global_threshold = values.quantile(percentile)

    for group, idx in groups.groupby(groups, observed=True).groups.items():
        idx = list(idx)
        subset = values.loc[idx]
        threshold = global_threshold if len(idx) < min_cells else subset.quantile(percentile)
        thresholds[group] = threshold
        out.loc[idx] = subset >= threshold
    return out, pd.Series(thresholds, name=f"SI_q{int(percentile * 100)}")


def quartile_within_group(values: pd.Series, groups: pd.Series, min_cells: int = 50) -> pd.Categorical:
    values = pd.Series(values, index=groups.index, dtype="float64")
    groups = groups.astype("string").fillna("Unknown")
    out = pd.Series(pd.NA, index=values.index, dtype="object")
    global_edges = values.quantile([0.25, 0.50, 0.75]).to_numpy()

    def assign(vals: pd.Series, edges: np.ndarray) -> pd.Series:
        return pd.cut(vals, bins=[-np.inf, *edges, np.inf], labels=["Q1", "Q2", "Q3", "Q4"], include_lowest=True)

    for _, idx in groups.groupby(groups, observed=True).groups.items():
        idx = list(idx)
        subset = values.loc[idx]
        edges = global_edges if len(idx) < min_cells else subset.quantile([0.25, 0.50, 0.75]).to_numpy()
        out.loc[idx] = assign(subset, edges).astype("object")
    return pd.Categorical(out, categories=["Q1", "Q2", "Q3", "Q4"], ordered=True)


def summarize_call(adata: ad.AnnData, call_col: str, group_cols: list[str]) -> pd.DataFrame:
    frame = adata.obs[[*group_cols, call_col]].copy()
    rows = []
    for keys, sub in frame.groupby(group_cols, observed=True, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        n = len(sub)
        n_called = int(sub[call_col].sum())
        row = dict(zip(group_cols, keys))
        row.update({"n_cells": n, "n_called": n_called, "pct_called": 100 * n_called / n if n else np.nan})
        rows.append(row)
    return pd.DataFrame(rows).sort_values(group_cols).reset_index(drop=True)


# %%
# =========================
# Load GSE130973
# =========================
def load_gse130973() -> ad.AnnData:
    target_dir = DATA_DIR / "GSE130973"
    matrix_path = download_file(GSE130973_URLS["matrix"], target_dir / "GSE130973_matrix_filtered.mtx.gz")
    genes_path = download_file(GSE130973_URLS["genes"], target_dir / "GSE130973_genes_filtered.tsv.gz")
    barcodes_path = download_file(GSE130973_URLS["barcodes"], target_dir / "GSE130973_barcodes_filtered.tsv.gz")

    genes = read_tsv_gz(genes_path, header=None)
    barcodes = read_tsv_gz(barcodes_path, header=None)[0].astype(str).tolist()

    gene_symbols = genes.iloc[:, 1].astype(str).tolist() if genes.shape[1] > 1 else genes.iloc[:, 0].astype(str).tolist()
    gene_symbols = make_var_names_unique(gene_symbols)

    with gzip.open(matrix_path, "rb") as handle:
        matrix = scipy.io.mmread(handle).tocsr()

    # GEO MTX is gene x cell; AnnData expects cell x gene.
    adata = ad.AnnData(X=matrix.T.tocsr())
    adata.obs_names = [f"GSE130973_{bc}" for bc in barcodes]
    adata.var_names = gene_symbols
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()

    suffix = pd.Series(barcodes, index=adata.obs_names).str.extract(r"-(\d+)$", expand=False)
    meta_rows = [
        GSE130973_SUFFIX_META.get(str(value), {
            "sample_id": f"unknown_suffix_{value}",
            "donor_id": pd.NA,
            "donor_age": pd.NA,
            "age_group": pd.NA,
            "anatomic_site": "inguinoiliac region",
        })
        for value in suffix
    ]
    meta_frame = pd.DataFrame(meta_rows, index=adata.obs_names)
    adata.obs = adata.obs.join(meta_frame)
    adata.obs["study_id"] = "GSE130973"
    adata.obs["geo_accession"] = "GSE130973"
    adata.obs["compartment"] = "whole skin"
    adata.obs["cell_id"] = adata.obs_names
    return adata


adata_130973 = load_gse130973()
print(adata_130973)
display(adata_130973.obs[["study_id", "sample_id", "donor_id", "donor_age", "age_group"]].drop_duplicates())


# %%
# =========================
# Load GSE202352
# =========================
def find_mex_dirs(root: Path) -> list[Path]:
    matrix_files = list(root.rglob("matrix.mtx.gz")) + list(root.rglob("matrix.mtx"))
    out = []
    for matrix_file in matrix_files:
        parent = matrix_file.parent
        has_barcodes = (parent / "barcodes.tsv.gz").exists() or (parent / "barcodes.tsv").exists()
        has_features = (
            (parent / "features.tsv.gz").exists()
            or (parent / "features.tsv").exists()
            or (parent / "genes.tsv.gz").exists()
            or (parent / "genes.tsv").exists()
        )
        if has_barcodes and has_features:
            out.append(parent)
    return sorted(set(out))


def load_gse202352() -> ad.AnnData:
    target_dir = DATA_DIR / "GSE202352"
    raw_tar = download_file(GSE202352_RAW_URL, target_dir / "GSE202352_RAW.tar")
    outer_extract = target_dir / "raw_extract"
    if not outer_extract.exists() or not any(outer_extract.iterdir()):
        safe_extract_tar(raw_tar, outer_extract)

    inner_archives = sorted(list(outer_extract.rglob("*.tar.gz")) + list(outer_extract.rglob("*.tgz")) + list(outer_extract.rglob("*.tar")))
    extracted_dir = target_dir / "mex_extract"
    extracted_dir.mkdir(parents=True, exist_ok=True)

    for archive in inner_archives:
        sample_extract = extracted_dir / archive.name.replace(".tar.gz", "").replace(".tgz", "").replace(".tar", "")
        if not sample_extract.exists() or not any(sample_extract.iterdir()):
            safe_extract_tar(archive, sample_extract)

    adatas = []
    for mex_dir in find_mex_dirs(extracted_dir):
        archive_context = " ".join(part.name for part in mex_dir.parents)
        accession_match = re.search(r"(GSM\d+)", str(mex_dir))
        if accession_match is None:
            accession_match = re.search(r"(GSM\d+)", archive_context)
        if accession_match is None:
            print(f"Skipping MEX dir without GSM accession: {mex_dir}")
            continue
        accession = accession_match.group(1)
        if accession not in GSE202352_SAMPLE_META:
            print(f"Skipping unrecognized GSE202352 accession: {accession}")
            continue

        one = sc.read_10x_mtx(mex_dir, var_names="gene_symbols", cache=False)
        one.var_names_make_unique()
        one.layers["counts"] = one.X.copy()

        meta = GSE202352_SAMPLE_META[accession].copy()
        one.obs_names = [f"{accession}_{bc}" for bc in one.obs_names.astype(str)]
        for key, value in meta.items():
            one.obs[key] = value
        one.obs["study_id"] = "GSE202352"
        one.obs["geo_accession"] = accession
        one.obs["donor_age"] = pd.NA
        one.obs["age_group"] = pd.NA
        one.obs["cell_id"] = one.obs_names
        adatas.append(one)

    if not adatas:
        raise RuntimeError(f"No GSE202352 MEX directories were found under {extracted_dir}")

    combined = ad.concat(adatas, join="outer", label=None, index_unique=None, merge="same")
    combined.layers["counts"] = combined.X.copy()
    return combined


adata_202352 = load_gse202352()
print(adata_202352)
display(adata_202352.obs[["study_id", "geo_accession", "sample_id", "donor_id", "anatomic_site", "compartment"]].drop_duplicates())


# %%
# =========================
# Combine, optional metadata, normalize
# =========================
adata = ad.concat(
    [adata_130973, adata_202352],
    join="outer",
    label="source_dataset",
    keys=["GSE130973", "GSE202352"],
    index_unique=None,
    merge="same",
)
adata.var_names_make_unique()
adata.layers["counts"] = adata.X.copy()

if OPTIONAL_CELL_METADATA_CSV is not None:
    cell_meta = pd.read_csv(OPTIONAL_CELL_METADATA_CSV)
    if "cell_id" not in cell_meta.columns:
        raise KeyError("OPTIONAL_CELL_METADATA_CSV must contain a cell_id column.")
    cell_meta = cell_meta.drop_duplicates("cell_id").set_index("cell_id")
    adata.obs = adata.obs.join(cell_meta, on="cell_id", rsuffix="_external")

sc.pp.filter_cells(adata, min_counts=1)
sc.pp.filter_genes(adata, min_cells=3)
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
sc.pp.log1p(adata)
adata.layers["log1p"] = adata.X.copy()

print(adata)
display(adata.obs[["study_id", "sample_id", "donor_id", "donor_age", "age_group", "anatomic_site", "compartment"]].drop_duplicates())


# %%
# =========================
# Apply SI / quartile scoring
# =========================
lookup = build_gene_lookup(adata)
growth_genes = present_genes(GROWTH_GENES, lookup)
s_genes = present_genes(S_GENES, lookup)
g2m_genes = present_genes(G2M_GENES, lookup)
sasp_genes = present_genes(SASP_GENES, lookup)

print(f"Growth genes present: {len(growth_genes)}/{len(GROWTH_GENES)} {growth_genes}")
print(f"S genes present: {len(s_genes)}/{len(S_GENES)}")
print(f"G2M genes present: {len(g2m_genes)}/{len(G2M_GENES)}")
print(f"SASP genes present: {len(sasp_genes)}/{len(SASP_GENES)}")

available_group_cols = [col for col in SCORING_GROUP_CANDIDATES if col in adata.obs.columns and adata.obs[col].notna().any()]
if not any(col.startswith("cell_type") for col in available_group_cols):
    available_group_cols = ["study_id"]
    print("No cell-type metadata found. Scoring groups are study_id only.")

adata.obs["quartile_scoring_group"] = pd.Categorical(group_key(adata.obs, available_group_cols))
groups = pd.Series(adata.obs["quartile_scoring_group"], index=adata.obs_names)
print("Scoring groups:", available_group_cols)
display(adata.obs["quartile_scoring_group"].value_counts().to_frame("n_cells"))

adata.obs["growth_raw"] = module_score(adata, growth_genes, layer="log1p")
adata.obs["Z_growth"] = zscore_within_group(adata.obs["growth_raw"], groups, MIN_CELLS_PER_GROUP)

adata.obs["S_score_custom"] = module_score(adata, s_genes, layer="log1p")
adata.obs["G2M_score_custom"] = module_score(adata, g2m_genes, layer="log1p")
adata.obs["cycle_raw"] = adata.obs[["S_score_custom", "G2M_score_custom"]].mean(axis=1)
adata.obs["Z_cycle"] = zscore_within_group(adata.obs["cycle_raw"], groups, MIN_CELLS_PER_GROUP)
adata.obs["Z_cycle_suppression"] = -adata.obs["Z_cycle"]

adata.obs["sasp_raw"] = module_score(adata, sasp_genes, layer="log1p")
adata.obs["Z_SASP"] = zscore_within_group(adata.obs["sasp_raw"], groups, MIN_CELLS_PER_GROUP)

adata.obs["senescence_index"] = (
    W_GROWTH * adata.obs["Z_growth"]
    + W_SASP * adata.obs["Z_SASP"]
    - W_CYCLE * adata.obs["Z_cycle"]
)

adata.obs["SI_quartile_within_group"] = quartile_within_group(
    adata.obs["senescence_index"], groups, MIN_CELLS_PER_GROUP
)

threshold_tables = []
call_cols = []
for percentile in TOP_PERCENTILES:
    pct = int(round((1 - percentile) * 100))
    col = f"SI_top{pct}_pct_within_group"
    adata.obs[col], thresholds = percentile_call_within_group(
        adata.obs["senescence_index"], groups, percentile, MIN_CELLS_PER_GROUP
    )
    call_cols.append(col)
    threshold_tables.append(thresholds.rename(col))

si_thresholds = pd.concat(threshold_tables, axis=1).rename_axis("quartile_scoring_group").reset_index()
display(si_thresholds)
display(adata.obs[["growth_raw", "Z_growth", "cycle_raw", "Z_cycle", "sasp_raw", "Z_SASP", "senescence_index"]].describe().T)


# %%
# =========================
# Summaries
# =========================
summary_group_cols = ["study_id", "sample_id", "donor_id", "anatomic_site", "compartment", "age_group"]
summary_group_cols = [col for col in summary_group_cols if col in adata.obs.columns]

sample_axis_summary = (
    adata.obs.groupby(summary_group_cols, observed=True, dropna=False)
    .agg(
        n_cells=("senescence_index", "size"),
        SI_mean=("senescence_index", "mean"),
        SI_median=("senescence_index", "median"),
        SI_q90=("senescence_index", lambda x: x.quantile(0.90)),
        SI_q95=("senescence_index", lambda x: x.quantile(0.95)),
        growth_Z_mean=("Z_growth", "mean"),
        SASP_Z_mean=("Z_SASP", "mean"),
        cycle_Z_mean=("Z_cycle", "mean"),
    )
    .reset_index()
)
display(sample_axis_summary)

call_summaries = {}
for call_col in call_cols:
    by_sample = summarize_call(adata, call_col, summary_group_cols)
    by_group = summarize_call(adata, call_col, ["quartile_scoring_group"])
    call_summaries[f"{call_col}_by_sample"] = by_sample
    call_summaries[f"{call_col}_by_scoring_group"] = by_group
    print(call_col, "by sample")
    display(by_sample)
    print(call_col, "by scoring group")
    display(by_group)

quartile_summary = (
    adata.obs.groupby(["study_id", "SI_quartile_within_group"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)
quartile_summary["pct_within_study"] = (
    quartile_summary["n_cells"]
    / quartile_summary.groupby("study_id", observed=True)["n_cells"].transform("sum")
    * 100
)
display(quartile_summary)


# %%
# =========================
# Plots
# =========================
sample_axis_summary["sample_label"] = (
    sample_axis_summary["study_id"].astype(str) + " | " + sample_axis_summary["sample_id"].astype(str)
)

plot_order = sample_axis_summary.sort_values(["study_id", "sample_id"])["sample_label"].tolist()
plot_df = sample_axis_summary.set_index("sample_label").loc[plot_order]
plot_cols = ["SI_mean", "growth_Z_mean", "SASP_Z_mean", "cycle_Z_mean"]

fig_height = max(5, 0.35 * len(plot_df) + 2)
fig, ax = plt.subplots(figsize=(8.5, fig_height))
im = ax.imshow(plot_df[plot_cols].to_numpy(dtype=float), aspect="auto", cmap="coolwarm")
ax.set_xticks(range(len(plot_cols)))
ax.set_xticklabels(plot_cols, rotation=45, ha="right")
ax.set_yticks(range(len(plot_df)))
ax.set_yticklabels(plot_df.index)
ax.set_title("Mean SI axes by GEO sample")
fig.colorbar(im, ax=ax, label="Mean score")
fig.tight_layout()
fig.savefig(OUTDIR / "geo_combined_axis_scores_by_sample_heatmap.png", dpi=300, bbox_inches="tight")
plt.show()

for call_col in call_cols:
    pct_col = "pct_called"
    table = call_summaries[f"{call_col}_by_sample"].copy()
    table["sample_label"] = table["study_id"].astype(str) + " | " + table["sample_id"].astype(str)
    table = table.sort_values(["study_id", "sample_id"])

    fig, ax = plt.subplots(figsize=(max(8, 0.45 * len(table)), 4.5))
    ax.bar(table["sample_label"], table[pct_col], color="#4C78A8", edgecolor="black")
    ax.set_ylabel("% SI-high cells")
    ax.set_title(call_col)
    ax.tick_params(axis="x", rotation=90)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / f"{call_col}_by_sample_barplot.png", dpi=300, bbox_inches="tight")
    plt.show()


# %%
# =========================
# Save outputs
# =========================
adata.uns["geo_combined_quartile_snc_scoring"] = {
    "method": "Continuous-axis SI with configurable group-wise z-scoring and top-percentile calls",
    "datasets": ["GSE130973", "GSE202352"],
    "scoring_group_cols": available_group_cols,
    "min_cells_per_group": MIN_CELLS_PER_GROUP,
    "top_percentiles": TOP_PERCENTILES,
    "weights": {"growth": W_GROWTH, "SASP": W_SASP, "cycle_penalty": W_CYCLE},
    "growth_genes_present": growth_genes,
    "sasp_genes_present": sasp_genes,
    "s_genes_present": s_genes,
    "g2m_genes_present": g2m_genes,
}

adata.write_h5ad(OUTDIR / "geo_combined_quartile_snc_scored.h5ad", compression="gzip")
adata.obs.to_csv(OUTDIR / "geo_combined_quartile_snc_obs.csv")
sample_axis_summary.to_csv(OUTDIR / "geo_combined_SI_axis_summary_by_sample.csv", index=False)
si_thresholds.to_csv(OUTDIR / "geo_combined_SI_thresholds_by_group.csv", index=False)
quartile_summary.to_csv(OUTDIR / "geo_combined_SI_quartile_summary_by_study.csv", index=False)
for name, table in call_summaries.items():
    table.to_csv(OUTDIR / f"geo_combined_{name}.csv", index=False)

print(f"Saved outputs to: {OUTDIR}")
