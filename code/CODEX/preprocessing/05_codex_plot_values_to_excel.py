"""
Bundle CODEX plot-ready values into one Excel workbook.

This exporter intentionally avoids the immune-audit outputs. For older CODEX
plots that did not save CSV value tables, it reconstructs the plotted values
from the raw CODEX Excel files using the marker positivity and annotation rules
used by the CODEX R scripts.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import xlsxwriter


PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUTDIR = PROJECT_ROOT / "results" / "260727" / "codex_plot_values_excel"
WORKBOOK = OUTDIR / "codex_plot_values_for_graphs_clean_minor_cells.xlsx"

RAW_DIR = PROJECT_ROOT / "data" / "codex"
LEGACY_RESULT_DIR = PROJECT_ROOT / "results" / "260303"

AGE_ORDER = ["Young", "Old"]
SKIN_LEVELS = [
    "Fibroblast",
    "Melanocyte",
    "Basal Keratinocyte",
    "Suprabasal Keratinocyte",
    "Corneocyte",
    "Smooth Muscle",
]

MARKER_PANEL = [
    "BCL2",
    "SMA",
    "p21",
    "Keratin 14",
    "Vimentin",
    "Keratin 5",
    "p16",
    "CD34",
    "IP10",
    "CD8",
    "hIL-1alpha",
    "CD14",
    "LaminB1",
    "Ki67",
    "Filaggrin",
    "CXCR3",
    "HMGB1",
    "Sirt1",
    "IDO1",
    "FOXP3",
    "IL6",
    "MelanA",
    "IL8",
    "CD4",
    "CD68",
    "HLA-DR",
    "CD45",
    "PCNA",
    "ICOS",
    "Pan-Cytokeratin",
    "CD3e",
    "CD11c",
]

MARKER_ROOT_MAP = {
    "p21": "p21",
    "p16": "p16",
    "p16bc42": "p16",
    "pcna": "PCNA",
    "ki67": "Ki67",
    "hmgb1": "HMGB1",
    "laminb": "LaminB1",
    "laminb1": "LaminB1",
    "ip10": "IP10",
    "cxcr3": "CXCR3",
    "cxcr": "CXCR3",
    "hil1a": "hIL-1alpha",
    "hil1alpha": "hIL-1alpha",
    "il1alpha": "hIL-1alpha",
    "il1a": "hIL-1alpha",
    "il6": "IL6",
    "il8": "IL8",
    "sirt1": "Sirt1",
    "sirtt1": "Sirt1",
    "ido1": "IDO1",
    "cd3e": "CD3e",
    "cd4": "CD4",
    "cd8": "CD8",
    "cd14": "CD14",
    "cd68": "CD68",
    "cd11c": "CD11c",
    "cd34": "CD34",
    "cd45": "CD45",
    "foxp3": "FOXP3",
    "hladr": "HLA-DR",
    "hla": "HLA-DR",
    "icos": "ICOS",
    "sma": "SMA",
    "vimentin": "Vimentin",
    "k5": "Keratin 5",
    "keratin5": "Keratin 5",
    "keratin14": "Keratin 14",
    "k14": "Keratin 14",
    "melana": "MelanA",
    "mela": "MelanA",
    "filaggrin": "Filaggrin",
    "filagrin": "Filaggrin",
    "fillagrin": "Filaggrin",
    "pancyto": "Pan-Cytokeratin",
    "pancytokeratin": "Pan-Cytokeratin",
    "bcl2": "BCL2",
    "bcl2bc42": "BCL2",
    "bcl2bc": "BCL2",
    "bcl02": "BCL2",
}

SAMPLE_EXCLUDE_MARKERS = {
    "202UA": {"CD3e"},
    "207UA": {"CD3e"},
    "7503": {"CD11c"},
    "4614": {"Pan-Cytokeratin"},
}

LEGACY_CSVS = {
    "legacy_skin_celltypes": (
        LEGACY_RESULT_DIR / "codex_skin_celltype_distribution_minorskincell1.csv",
        "Existing CSV exported by 04_codex_raw_values_to_doc.R: age-level minor skin cell-type distribution.",
    ),
    "legacy_snc_counts": (
        LEGACY_RESULT_DIR / "codex_snc2_counts_combined.csv",
        "Existing CSV exported by 04_codex_raw_values_to_doc.R: age/sample SnC2 counts.",
    ),
    "legacy_snc_by_celltype": (
        LEGACY_RESULT_DIR / "codex_snc2_by_celltype_combined_minorskincell1.csv",
        "Existing CSV exported by 04_codex_raw_values_to_doc.R: SnC2 burden and composition by minor skin cell type.",
    ),
    "legacy_cd45_snc": (
        LEGACY_RESULT_DIR / "codex_cd45_snc2_counts_combined.csv",
        "Existing CSV exported by 04_codex_raw_values_to_doc.R: CD45+ SnC2 counts.",
    ),
    "legacy_260310_cor_expr": (
        PROJECT_ROOT / "results" / "260310" / "cor_snc_expr.csv",
        "Older CODEX marker-expression correlation values from results/260310.",
    ),
    "legacy_260310_cor_pos": (
        PROJECT_ROOT / "results" / "260310" / "cor_snc_pos.csv",
        "Older CODEX marker-positivity correlation values from results/260310.",
    ),
    "legacy_260310_long": (
        PROJECT_ROOT / "results" / "260310" / "df_codex_long.csv",
        "Older long-format per-cell CODEX marker table from results/260310.",
    ),
    "legacy_260310_sample": (
        PROJECT_ROOT / "results" / "260310" / "sum_codex_sample.csv",
        "Older sample-level CODEX marker summary from results/260310.",
    ),
}


@dataclass(frozen=True)
class SheetRecord:
    sheet: str
    description: str
    source: str
    rows: int
    columns: int


def normalize_token(value: object) -> str:
    token = str(value)
    token = re.sub(r"^MeanIntensity[_\s-]*", "", token, flags=re.IGNORECASE)
    return re.sub(r"[^a-z0-9]", "", token.lower())


def marker_for_column(column: object) -> str | None:
    raw = str(column).lower()
    if re.search(r"thermal|thermo|bad", raw):
        return None
    return MARKER_ROOT_MAP.get(normalize_token(column))


def infer_age(path: Path) -> str | None:
    lower_parts = {part.lower() for part in path.parts}
    if "young" in lower_parts:
        return "Young"
    if "old" in lower_parts:
        return "Old"
    return None


def infer_sample(path: Path) -> str:
    name = path.name
    match = re.search(r"(CAMS\d+a|Subject\d+|\d{3,4}UA)", name, flags=re.IGNORECASE)
    if match:
        sample = match.group(1)
        sample = re.sub(r"^cams", "CAMS", sample, flags=re.IGNORECASE)
        sample = re.sub(r"^subject", "Subject", sample, flags=re.IGNORECASE)
        sample = re.sub(r"ua$", "UA", sample, flags=re.IGNORECASE)
        return sample
    match = re.search(r"(^|[^0-9])([0-9]{4})(-[0-9]+)?", name)
    if match:
        return match.group(2)
    return path.stem


def pick_sheet(path: Path) -> str:
    excel = pd.ExcelFile(path)
    scored: list[tuple[int, int, int, str]] = []
    for sheet in excel.sheet_names:
        try:
            head = pd.read_excel(path, sheet_name=sheet, nrows=5)
        except Exception:
            continue
        roots = [normalize_token(col) for col in head.columns]
        marker_hits = sum(marker_for_column(col) is not None for col in head.columns)
        coord_hits = int("x" in roots) + int("y" in roots)
        nuclei_bonus = 2 if re.search(r"nuclei.*props|nucleiprops", sheet, re.IGNORECASE) else 0
        sheet1_bonus = 1 if sheet.lower() == "sheet1" else 0
        scored.append((marker_hits, coord_hits, nuclei_bonus + sheet1_bonus, sheet))
    if not scored:
        return excel.sheet_names[0]
    scored.sort(reverse=True)
    return scored[0][3]


def find_column(columns: list[str], accepted_roots: set[str]) -> str | None:
    for col in columns:
        if normalize_token(col) in accepted_roots:
            return col
    return None


def load_codex_sample(path: Path) -> tuple[pd.DataFrame, dict]:
    sample = infer_sample(path)
    age = infer_age(path)
    sheet = pick_sheet(path)
    raw = pd.read_excel(path, sheet_name=sheet)
    raw.columns = [str(col) for col in raw.columns]

    label_col = find_column(raw.columns.tolist(), {"label", "cellid", "objectid", "id", "barcode"})
    if label_col is None:
        labels = pd.Series(np.arange(1, len(raw) + 1), index=raw.index).astype(str)
    else:
        labels = raw[label_col].astype(str)

    x_col = find_column(raw.columns.tolist(), {"x", "umap1"})
    y_col = find_column(raw.columns.tolist(), {"y", "umap2"})

    marker_sources: dict[str, list[str]] = {marker: [] for marker in MARKER_PANEL}
    for col in raw.columns:
        marker = marker_for_column(col)
        if marker is not None:
            marker_sources[marker].append(col)

    excluded = SAMPLE_EXCLUDE_MARKERS.get(sample, set())
    present_markers = [m for m in MARKER_PANEL if marker_sources[m] and m not in excluded]

    marker_values = {}
    numeric_counts = pd.Series(0, index=raw.index)
    for marker in present_markers:
        values = raw[marker_sources[marker]].apply(pd.to_numeric, errors="coerce")
        numeric_counts = numeric_counts + values.notna().any(axis=1).astype(int)
        marker_values[marker] = values.max(axis=1, skipna=True).fillna(0).astype(float)

    if not marker_values:
        raise ValueError(f"No recognized marker columns found in {path.name} / {sheet}")

    cell = pd.DataFrame(marker_values)
    cell.insert(0, "cell_id", [f"{sample}_{idx}" for idx in pd.RangeIndex(1, len(cell) + 1)])
    cell.insert(1, "source_label", labels)
    cell.insert(2, "sample", sample)
    cell.insert(3, "age", age)
    cell.insert(4, "source_file", path.name)
    cell.insert(5, "source_sheet", sheet)
    cell.insert(6, "x", pd.to_numeric(raw[x_col], errors="coerce") if x_col else np.nan)
    cell.insert(7, "y", pd.to_numeric(raw[y_col], errors="coerce") if y_col else np.nan)

    footer_label = labels.str.contains(
        r"threshold|thresh|cutoff|summary|mean|median|sd|note|background",
        case=False,
        regex=True,
        na=False,
    )
    keep = labels.notna() & labels.ne("") & ~footer_label & (numeric_counts >= 2)
    cell = cell.loc[keep].reset_index(drop=True)

    for marker in MARKER_PANEL:
        if marker not in cell.columns:
            cell[marker] = 0.0
        cell[f"{marker}_pos"] = cell[marker] > 0

    p21 = cell["p21_pos"]
    p16 = cell["p16_pos"]
    pcna = cell["PCNA_pos"]
    bcl2 = cell["BCL2_pos"]
    ip10 = cell["IP10_pos"]
    il6 = cell["IL6_pos"]
    il8 = cell["IL8_pos"]

    cell["cdkn_state"] = np.select(
        [p16 & p21, p16 & ~p21, p21 & ~p16],
        ["p16+p21+", "p16+", "p21+"],
        default="p16-p21-",
    )
    cell["snc2"] = np.where((p21 | p16) & ~pcna & (ip10 | il6 | il8), "SnC+", "SnC-")
    cell["snc_bcl2"] = np.where(
        (p21 | p16) & ~pcna & bcl2 & (ip10 | il6 | il8), "SnC+", "SnC-"
    )

    cell["majorskincell"] = np.select(
        [
            cell["MelanA_pos"],
            cell["SMA_pos"],
            cell["Vimentin_pos"],
            cell["Keratin 5_pos"] | cell["Keratin 14_pos"] | cell["Filaggrin_pos"],
        ],
        ["Melanocyte", "Smooth Muscle", "Fibroblast", "Keratinocyte"],
        default="none",
    )
    cell["minorskincell1"] = np.select(
        [
            cell["MelanA_pos"],
            cell["SMA_pos"],
            cell["Vimentin_pos"],
            cell["Keratin 5_pos"] | cell["Keratin 14_pos"],
            cell["Filaggrin_pos"],
            (~cell["Filaggrin_pos"])
            & (~cell["Keratin 5_pos"])
            & (~cell["Keratin 14_pos"])
            & cell["Pan-Cytokeratin_pos"],
        ],
        [
            "Melanocyte",
            "Smooth Muscle",
            "Fibroblast",
            "Basal Keratinocyte",
            "Corneocyte",
            "Suprabasal Keratinocyte",
        ],
        default="none",
    )
    cell["immcelltype"] = np.select(
        [
            cell["CD4_pos"] & cell["FOXP3_pos"],
            cell["CD4_pos"] & cell["CXCR3_pos"],
            cell["CD14_pos"],
            cell["CD68_pos"],
            cell["CD8_pos"],
            cell["CD4_pos"],
            cell["CD45_pos"],
        ],
        ["CD4+ FOXP3+", "CD4+ CXCR3+", "CD14+", "CD68+", "CD8+", "CD4+", "CD45+"],
        default="none",
    )

    meta = {
        "sample": sample,
        "age": age,
        "source_file": path.name,
        "source_path": str(path.relative_to(PROJECT_ROOT)),
        "source_sheet": sheet,
        "n_cells": int(len(cell)),
        "n_markers_present": int(len(present_markers)),
        "markers_present": "; ".join(present_markers),
        "excluded_markers": "; ".join(sorted(excluded)),
    }
    return cell, meta


def ordered(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        if col == "age" and col in out:
            out[col] = pd.Categorical(out[col], categories=AGE_ORDER, ordered=True)
    return out.sort_values([c for c in columns if c in out.columns]).reset_index(drop=True)


def prop_summary(data: pd.DataFrame, group_cols: list[str], flag_col: str, name: str) -> pd.DataFrame:
    grouped = data.groupby(group_cols, dropna=False, observed=False)
    out = grouped.agg(
        n_total_cells=("cell_id", "size"),
        n_positive_cells=(flag_col, lambda x: int((x == "SnC+").sum())),
    ).reset_index()
    out[f"{name}_prop"] = out["n_positive_cells"] / out["n_total_cells"]
    out[f"{name}_pct"] = 100 * out[f"{name}_prop"]
    return ordered(out, group_cols)


def snc_counts(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for definition in ["snc2", "snc_bcl2"]:
        age = prop_summary(cells, ["age"], definition, definition)
        age.insert(0, "definition", definition)
        age.insert(1, "level", "age")
        age["sample"] = pd.NA

        sample = prop_summary(cells, ["age", "sample"], definition, definition)
        sample.insert(0, "definition", definition)
        sample.insert(1, "level", "sample")
        rows.extend([age, sample])
    return pd.concat(rows, ignore_index=True)


def snc_sample_summary(cells: pd.DataFrame) -> pd.DataFrame:
    sample = snc_counts(cells).query("level == 'sample'").copy()
    rows = []
    for definition, part in sample.groupby("definition"):
        prop_col = f"{definition}_prop"
        tmp = part.groupby("age", observed=False)[prop_col].agg(["mean", "std", "count"]).reset_index()
        tmp["sem"] = tmp["std"] / np.sqrt(tmp["count"])
        tmp.insert(0, "definition", definition)
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True)


def celltype_distribution(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cell_col in ["majorskincell", "minorskincell1", "immcelltype"]:
        for level, group_cols in [("age", ["age"]), ("sample", ["age", "sample"])]:
            denom = cells.groupby(group_cols, observed=False)["cell_id"].size().rename("total_cells")
            num = (
                cells.groupby(group_cols + [cell_col], observed=False)["cell_id"]
                .size()
                .rename("n_cells")
                .reset_index()
            )
            out = num.merge(denom.reset_index(), on=group_cols, how="left")
            out["prop_of_total"] = out["n_cells"] / out["total_cells"]
            out["pct_of_total"] = 100 * out["prop_of_total"]
            out.insert(0, "celltype_column", cell_col)
            out.insert(1, "level", level)
            rows.append(out.rename(columns={cell_col: "celltype"}))
    return pd.concat(rows, ignore_index=True)


def snc_by_celltype(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for definition in ["snc2", "snc_bcl2"]:
        for cell_col in ["majorskincell", "minorskincell1", "immcelltype"]:
            for level, group_cols in [("age", ["age"]), ("sample", ["age", "sample"])]:
                grouped = cells.groupby(group_cols + [cell_col], observed=False)
                out = grouped.agg(
                    n_total_cells=("cell_id", "size"),
                    n_snc_cells=(definition, lambda x: int((x == "SnC+").sum())),
                ).reset_index()
                out["snc_prop_within_celltype"] = out["n_snc_cells"] / out["n_total_cells"]
                out["snc_pct_within_celltype"] = 100 * out["snc_prop_within_celltype"]
                total_snc = (
                    out.groupby(group_cols, observed=False)["n_snc_cells"]
                    .sum()
                    .rename("total_snc_in_group")
                    .reset_index()
                )
                out = out.merge(total_snc, on=group_cols, how="left")
                out["snc_composition_prop"] = np.where(
                    out["total_snc_in_group"] > 0,
                    out["n_snc_cells"] / out["total_snc_in_group"],
                    np.nan,
                )
                out["snc_composition_pct"] = 100 * out["snc_composition_prop"]
                out.insert(0, "definition", definition)
                out.insert(1, "celltype_column", cell_col)
                out.insert(2, "level", level)
                rows.append(out.rename(columns={cell_col: "celltype"}))
    return pd.concat(rows, ignore_index=True)


def donut_values(cells: pd.DataFrame) -> dict[str, pd.DataFrame]:
    total_minor = (
        cells[cells["minorskincell1"].isin(SKIN_LEVELS)]
        .groupby(["age", "minorskincell1"], observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    total_minor["total_cells_in_age"] = total_minor.groupby("age", observed=False)["n_cells"].transform("sum")
    total_minor["composition_prop"] = total_minor["n_cells"] / total_minor["total_cells_in_age"]
    total_minor["composition_pct"] = 100 * total_minor["composition_prop"]

    out = {"donut_total_minor": total_minor}
    for definition in ["snc2", "snc_bcl2"]:
        snc_minor = (
            cells[(cells[definition] == "SnC+") & cells["minorskincell1"].isin(SKIN_LEVELS)]
            .groupby(["age", "minorskincell1"], observed=False)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        snc_minor["total_snc_in_age"] = snc_minor.groupby("age", observed=False)["n_cells"].transform("sum")
        snc_minor["composition_prop"] = snc_minor["n_cells"] / snc_minor["total_snc_in_age"]
        snc_minor["composition_pct"] = 100 * snc_minor["composition_prop"]
        snc_minor.insert(0, "definition", definition)
        out[f"donut_snc_minor_{definition}"] = snc_minor

        cdkn = (
            cells[cells[definition] == "SnC+"]
            .groupby(["age", "cdkn_state"], observed=False)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        cdkn["total_snc_in_age"] = cdkn.groupby("age", observed=False)["n_cells"].transform("sum")
        cdkn["composition_prop"] = cdkn["n_cells"] / cdkn["total_snc_in_age"]
        cdkn["composition_pct"] = 100 * cdkn["composition_prop"]
        cdkn.insert(0, "definition", definition)
        out[f"donut_cdkn_{definition}"] = cdkn
    return out


def marker_expression_values(cells: pd.DataFrame) -> dict[str, pd.DataFrame]:
    violin_markers = ["BCL2", "HMGB1", "LaminB1"]
    sen_markers = ["p21", "p16", "PCNA", "BCL2", "IP10", "IL6", "IL8", "CXCR3"]

    violin = cells.melt(
        id_vars=["cell_id", "sample", "age", "snc2", "snc_bcl2", "minorskincell1", "immcelltype"],
        value_vars=violin_markers,
        var_name="marker",
        value_name="expression",
    )

    sample_means = (
        violin.groupby(["marker", "age", "sample", "snc2"], observed=False)["expression"]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
        .rename(columns={"count": "n_cells"})
    )

    rows = []
    for group_col in ["snc2", "snc_bcl2", "cdkn_state", "minorskincell1", "immcelltype"]:
        for marker in sen_markers:
            part = (
                cells.groupby(["age", "sample", group_col], observed=False)
                .agg(
                    avg_expression=(marker, "mean"),
                    pct_positive=(f"{marker}_pos", lambda x: 100 * float(np.mean(x))),
                    n_cells=("cell_id", "size"),
                )
                .reset_index()
                .rename(columns={group_col: "group"})
            )
            part.insert(0, "group_column", group_col)
            part.insert(1, "marker", marker)
            rows.append(part)
    marker_summary = pd.concat(rows, ignore_index=True)

    dotplot = []
    for split_col in ["snc2", "snc_bcl2"]:
        for marker in sen_markers:
            part = (
                cells.groupby(["age", "immcelltype", split_col], observed=False)
                .agg(
                    avg_expression=(marker, "mean"),
                    pct_positive=(f"{marker}_pos", lambda x: 100 * float(np.mean(x))),
                    n_cells=("cell_id", "size"),
                )
                .reset_index()
                .rename(columns={split_col: "snc_status"})
            )
            part.insert(0, "definition", split_col)
            part.insert(1, "marker", marker)
            dotplot.append(part)

    return {
        "marker_violin_cells": violin,
        "marker_violin_sample_means": sample_means,
        "marker_expression_summary": marker_summary,
        "immune_snc_marker_dotplot": pd.concat(dotplot, ignore_index=True),
    }


def cd45_snc_counts(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for definition in ["snc2", "snc_bcl2"]:
        for level, group_cols in [("age", ["age"]), ("sample", ["age", "sample"])]:
            part = cells[cells["CD45_pos"]].copy()
            out = (
                part.groupby(group_cols, observed=False)
                .agg(
                    n_cd45_cells=("cell_id", "size"),
                    n_cd45_snc_cells=(definition, lambda x: int((x == "SnC+").sum())),
                )
                .reset_index()
            )
            out["cd45_snc_prop"] = out["n_cd45_snc_cells"] / out["n_cd45_cells"]
            out["cd45_snc_pct"] = 100 * out["cd45_snc_prop"]
            out.insert(0, "definition", definition)
            out.insert(1, "level", level)
            rows.append(out)
    return pd.concat(rows, ignore_index=True)


def read_legacy_tables() -> dict[str, tuple[pd.DataFrame, str, str]]:
    out = {}
    legacy_minor_col = "minor" + "skincell"
    for sheet, (path, description) in LEGACY_CSVS.items():
        if path.exists():
            table = pd.read_csv(path, low_memory=False)
            table = table.rename(
                columns={
                    legacy_minor_col: "minor_cell_type_legacy",
                    f"{legacy_minor_col}_shared": "minor_cell_type_legacy",
                }
            )
            out[sheet] = (table, description, str(path.relative_to(PROJECT_ROOT)))
    return out


def clean_sheet_name(name: str) -> str:
    cleaned = re.sub(r"[\[\]:*?/\\]", "_", name)
    return cleaned[:31]


def excel_value(value: object) -> object:
    if value is None:
        return None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        if np.isnan(value) or np.isinf(value):
            return None
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if pd.isna(value):
        return None
    return value


def write_sheet(
    workbook: xlsxwriter.Workbook,
    records: list[SheetRecord],
    sheet: str,
    df: pd.DataFrame,
    description: str,
    source: str,
) -> None:
    df = df.copy()
    sheet_name = clean_sheet_name(sheet)
    worksheet = workbook.add_worksheet(sheet_name)
    header_fmt = workbook.add_format({"bold": True, "font_color": "white", "bg_color": "#1F4E79"})

    for col_idx, column in enumerate(df.columns):
        worksheet.write(0, col_idx, str(column), header_fmt)

    for row_idx, row in enumerate(df.itertuples(index=False, name=None), start=1):
        worksheet.write_row(row_idx, 0, [excel_value(value) for value in row])

    if len(df.columns) > 0:
        worksheet.freeze_panes(1, 0)
        worksheet.autofilter(0, 0, max(len(df), 1), len(df.columns) - 1)
        for col_idx, column in enumerate(df.columns):
            sample_values = df[column].head(200).astype(str).tolist()
            width = min(max([len(str(column)), *[len(v) for v in sample_values], 10]) + 2, 45)
            worksheet.set_column(col_idx, col_idx, width)

    records.append(SheetRecord(sheet, description, source, len(df), len(df.columns)))


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    raw_files = sorted(RAW_DIR.glob("*/*.xlsx"))
    if not raw_files:
        raise FileNotFoundError(f"No raw CODEX Excel files found under {RAW_DIR}")

    cell_tables = []
    manifest_rows = []
    for path in raw_files:
        cells, meta = load_codex_sample(path)
        cell_tables.append(cells)
        manifest_rows.append(meta)

    cells = pd.concat(cell_tables, ignore_index=True)
    manifest = pd.DataFrame(manifest_rows).sort_values(["age", "sample"])

    records: list[SheetRecord] = []
    workbook = xlsxwriter.Workbook(
        str(WORKBOOK),
        {"constant_memory": True, "nan_inf_to_errors": True, "tmpdir": str(OUTDIR)},
    )
    try:
        write_sheet(
            workbook,
            records,
            "sample_manifest",
            manifest,
            "Raw CODEX file manifest: sample, age, selected marker sheet, recovered cell count, and markers present.",
            "data/codex/*/*.xlsx",
        )

        per_cell_cols = [
            "cell_id",
            "source_label",
            "sample",
            "age",
            "source_file",
            "source_sheet",
            "x",
            "y",
            *MARKER_PANEL,
            *[f"{m}_pos" for m in MARKER_PANEL],
            "cdkn_state",
            "snc2",
            "snc_bcl2",
            "majorskincell",
            "minorskincell1",
            "immcelltype",
        ]
        write_sheet(
            workbook,
            records,
            "per_cell_annotations",
            cells[per_cell_cols],
            "Recovered per-cell marker values, positivity calls, SnC definitions, and skin/immune annotations.",
            "data/codex/*/*.xlsx",
        )

        write_sheet(
            workbook,
            records,
            "snc_counts",
            snc_counts(cells),
            "Age- and sample-level SnC counts/proportions for snc2 and BCL2-requiring SnC definitions.",
            "recomputed from raw CODEX Excel files",
        )
        write_sheet(
            workbook,
            records,
            "snc_sample_mean_sem",
            snc_sample_summary(cells),
            "Mean/SD/SEM of sample-level SnC proportions by age and SnC definition.",
            "recomputed from raw CODEX Excel files",
        )
        write_sheet(
            workbook,
            records,
            "celltype_distribution",
            celltype_distribution(cells),
            "Age- and sample-level total cell-type composition for major skin labels, minorskincell1 skin labels, and immune labels.",
            "recomputed from raw CODEX Excel files",
        )
        write_sheet(
            workbook,
            records,
            "snc_by_celltype",
            snc_by_celltype(cells),
            "Within-cell-type SnC burden and composition of the SnC pool by age/sample.",
            "recomputed from raw CODEX Excel files",
        )
        write_sheet(
            workbook,
            records,
            "cd45_snc_counts",
            cd45_snc_counts(cells),
            "CD45+ SnC counts/proportions by age/sample and SnC definition.",
            "recomputed from raw CODEX Excel files",
        )

        for sheet, table in donut_values(cells).items():
            write_sheet(
                workbook,
                records,
                sheet,
                table,
                "Donut/stacked composition values for total minor skin cells, SnC+ minor skin cells, or p16/p21 state.",
                "recomputed from raw CODEX Excel files",
            )

        for sheet, table in marker_expression_values(cells).items():
            write_sheet(
                workbook,
                records,
                sheet,
                table,
                "Recovered marker-expression values/summaries used for violin and marker dot-style plots.",
                "recomputed from raw CODEX Excel files",
            )

        for sheet, (table, description, source) in read_legacy_tables().items():
            write_sheet(workbook, records, sheet, table, description, source)

        index = pd.DataFrame([record.__dict__ for record in records])
        index.insert(
            0,
            "note",
            "No immune-audit sheets included. Older non-BANKSY plots are reconstructed from raw CODEX Excel files; BANKSY domain labels are not recoverable from raw marker spreadsheets alone.",
        )
        write_sheet(
            workbook,
            records,
            "index",
            index,
            "Workbook index describing each sheet, source, row count, and column count.",
            "generated by this exporter",
        )
    finally:
        workbook.close()

    print(f"Wrote {WORKBOOK}")
    print(pd.DataFrame([record.__dict__ for record in records]).to_string(index=False))


if __name__ == "__main__":
    main()
