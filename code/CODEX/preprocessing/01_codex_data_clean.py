#!/usr/bin/env python
"""Clean zero-thresholded CODEX Excel exports into an analysis-ready H5AD."""
from __future__ import annotations
import argparse, logging, re
from collections import Counter
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from openpyxl import load_workbook
import matplotlib.pyplot as plt

PANEL = ["BCL2","SMA","p21","Keratin 14","Vimentin","Keratin 5","p16","CD34","IP10","CD8","hIL-1alpha","CD14","LaminB1","Ki67","Filaggrin","CXCR3","HMGB1","Sirt1","IDO1","FOXP3","IL6","MelanA","IL8","CD4","CD68","HLA-DR","CD45","PCNA","ICOS","Pan-Cytokeratin","CD3e","CD11c"]
ROOT_MAP = {"p21":"p21","p16":"p16_non_bc42","p16bc42":"p16","p16layer9":"p16_layer9","pcna":"PCNA","ki67":"Ki67","hmgb1":"HMGB1","laminb":"LaminB1","laminb1":"LaminB1","ip10":"IP10","cxcr3":"CXCR3","hil1a":"hIL-1alpha","hil1alpha":"hIL-1alpha","il6":"IL6","il8":"IL8","sirt1":"Sirt1","sirtt1":"Sirt1","ido1":"IDO1","cd3e":"CD3e","cd4":"CD4","cd8":"CD8","cd14":"CD14","cd68":"CD68","cd11c":"CD11c","cd34":"CD34","cd45":"CD45","foxp3":"FOXP3","hladr":"HLA-DR","icos":"ICOS","sma":"SMA","vimentin":"Vimentin","k5":"Keratin 5","keratin5":"Keratin 5","k14":"Keratin 14","keratin14":"Keratin 14","mela":"MelanA","melana":"MelanA","filagrin":"Filaggrin","fillagrin":"Filaggrin","filaggrin":"Filaggrin","pancyto":"Pan-Cytokeratin","pancytokeratin":"Pan-Cytokeratin","bcl2":"BCL2","bcl2bc42":"BCL2","bcl2bc":"BCL2","bcl02":"BCL2"}
TYPE = {"SMA","Keratin 14","Vimentin","Keratin 5","CD34","CD8","CD14","Filaggrin","FOXP3","MelanA","CD4","CD68","HLA-DR","CD45","Pan-Cytokeratin","CD3e","CD11c"}
# Broad markers are retained for annotation but do not independently establish
# a lineage conflict (Vimentin is broad; HLA-DR is inducible outside immunity).
LINEAGE = {
    "epithelial":["Keratin 14","Keratin 5","Filaggrin","Pan-Cytokeratin"],
    "immune":["CD45","CD3e","CD4","CD8","FOXP3","CD14","CD68","CD11c"],
    "stromal":["SMA","CD34"],
    "melanocytic":["MelanA"],
}
BROAD_TYPE_MARKERS = ["Vimentin","HLA-DR"]
FOOTER = re.compile(r"threshold|thresh|cutoff|summary|mean|median|sd|note|background|^total$|^%$", re.I)

def token(x):
    return re.sub(r"[^a-z0-9]", "", re.sub(r"^MeanIntensity[_ ]*", "", str(x), flags=re.I).lower())

def sample_id(path):
    m = re.search(r"(CAMS\d+a|Subject\d+|\d{3,4}UA)", path.name, re.I)
    if m:
        x = re.sub(r"^cams", "CAMS", m.group(1), flags=re.I)
        x = re.sub(r"^subject", "Subject", x, flags=re.I)
        return re.sub(r"ua$", "UA", x, flags=re.I)
    m = re.search(r"(?:^|[^0-9])([0-9]{4})(?:-[0-9]+)?", path.name)
    return m.group(1) if m else path.stem

def excluded(sample):
    return {"CD3e"} if sample in {"202UA","207UA"} else {"CD11c"} if sample == "7503" else {"Pan-Cytokeratin"} if sample == "4614" else set()

def read_excel(path):
    wb = load_workbook(path, read_only=True, data_only=True)
    sheets = wb.sheetnames
    sheet = next((s for s in sheets if re.search(r"^nuclei\s*props$|nuclei.*props|nucleiprops", s, re.I)), "Sheet1" if "Sheet1" in sheets else sheets[0])
    rows = wb[sheet].iter_rows(values_only=True)
    headers = ["" if x is None else str(x) for x in next(rows)]
    frame = pd.DataFrame(rows, columns=headers)
    wb.close()
    return sheet, frame

def make_unique(values):
    used, next_n, out = set(), Counter(), []
    for value in values.astype(str):
        candidate = value
        if candidate in used:
            n = max(1, next_n[value])
            while f"{value}.{n}" in used: n += 1
            candidate, next_n[value] = f"{value}.{n}", n + 1
        used.add(candidate); out.append(candidate)
    return out

def clean(path):
    sample, (sheet, raw) = sample_id(path), read_excel(path)
    roots = [token(x) for x in raw.columns]
    label_i = next((i for i,x in enumerate(roots) if x == "label"), None)
    if label_i is None: label_i = next((i for i,x in enumerate(roots) if x in {"cellid","objectid","id","barcode"}), None)
    if label_i is None: raise ValueError("missing cell identifier")
    coord_i = (roots.index("umap1"),roots.index("umap2")) if "umap1" in roots and "umap2" in roots else (roots.index("x"),roots.index("y")) if "x" in roots and "y" in roots else None
    names = ["Label" if x == "label" else "X" if x == "x" else "Y" if x == "y" else ROOT_MAP.get(x,str(c)) for c,x in zip(raw.columns,roots)]
    frame = {}
    for i,name in enumerate(names): frame[name] = raw.iloc[:,i] if name not in frame else frame[name].combine_first(raw.iloc[:,i])
    frame = pd.DataFrame(frame, index=raw.index)
    if "Label" not in frame: frame["Label"] = raw.iloc[:,label_i]
    if sample.lower() == "cams18a":
        if "p16_layer9" not in frame: raise ValueError("CAMS18a missing p16 layer 9")
        frame["p16"], p16_source = frame["p16_layer9"], "p16 layer 9"
    else:
        if "p16" not in frame: raise ValueError(f"{sample} missing p16-BC42")
        p16_source = "p16-BC42"
    markers = [x for x in PANEL if x in frame and x not in excluded(sample)]
    labels = frame.Label.astype("string")
    values = frame[markers].apply(pd.to_numeric, errors="coerce")
    keep = labels.notna() & labels.str.strip().ne("") & ~labels.fillna("").str.contains(FOOTER) & values.notna().sum(1).ge(2)
    labels, values = labels[keep].astype(str), values[keep].replace([np.inf,-np.inf],np.nan).fillna(0)
    if (values < 0).any().any(): raise ValueError("negative values found in zero-thresholded data")
    ids = [f"{sample}_{x}" for x in make_unique(labels)]; values.index = ids
    spatial = np.full((len(ids),2),np.nan,np.float32) if coord_i is None else raw.iloc[:,list(coord_i)][keep].apply(pd.to_numeric,errors="coerce").to_numpy(np.float32)
    source_age = next((x.title() for x in path.parts if x.lower() in {"old","young"}), None)
    if source_age is None:
        raise ValueError("age is missing: workbook must be inside an old/ or young/ folder")
    return dict(sample=sample,age=source_age,source_age=source_age,sheet=sheet,file=path.name,p16_source=p16_source,markers=markers,values=values,labels=labels.to_numpy(str),ids=ids,spatial=spatial)

def assemble(items):
    availability = pd.DataFrame(False,index=[x["sample"] for x in items],columns=PANEL)
    matrices, obs, spatial = [], [], []
    for x in items:
        availability.loc[x["sample"],x["markers"]] = True
        matrices.append(x["values"].reindex(columns=PANEL,fill_value=0).to_numpy(np.float32))
        measured = x["values"].to_numpy(float)
        type_markers=[m for m in x["markers"] if m in TYPE]
        type_positive=x["values"][type_markers].gt(0)
        n_type=type_positive.sum(axis=1).to_numpy()
        type_names=type_positive.apply(lambda row: ";".join(row.index[row.to_numpy()]) if row.any() else "none",axis=1).to_numpy()
        obs.append(pd.DataFrame({"sample":x["sample"],"age":x["age"],"source_age":x["source_age"],"p16_source":x["p16_source"],"source_file":x["file"],"source_sheet":x["sheet"],"original_label":x["labels"],"n_panel_markers_measured":len(x["markers"]),"n_panel_markers_positive":(measured>0).sum(1),"total_panel_intensity":measured.sum(1),"n_type_markers_measured":len(type_markers),"n_type_markers_positive":n_type,"type_marker_count_group":np.where(n_type==0,"0",np.where(n_type==1,"1","2+")),"positive_type_markers":type_names},index=pd.Index(x["ids"],name="cell_id")))
        spatial.append(x["spatial"])
    matrix, obs = np.vstack(matrices), pd.concat(obs)
    obs["all_panel_zero"] = obs.n_panel_markers_positive.eq(0)
    for c in ["sample","age","source_age","p16_source"]: obs[c] = obs[c].astype("category")
    def marker_lineage(marker):
        if marker in BROAD_TYPE_MARKERS: return "broad_type"
        return next((lineage for lineage,markers in LINEAGE.items() if marker in markers),"state")
    var = pd.DataFrame({"marker_class":["type" if x in TYPE else "state" for x in PANEL],"marker_lineage":[marker_lineage(x) for x in PANEL],"used_for_lineage_conflict":[any(x in markers for markers in LINEAGE.values()) for x in PANEL]},index=pd.Index(PANEL,name="marker"))
    adata = ad.AnnData(matrix,obs=obs,var=var); adata.layers["log1p"] = np.log1p(matrix).astype(np.float32); adata.obsm["spatial"] = np.vstack(spatial)
    adata.uns["age_definition"] = "Canonical age group read from each workbook's old/ or young/ source folder. Any plot-specific sample restriction must be applied only in that plotting analysis."
    adata.uns["missing_marker_policy"] = "Unavailable markers are zero-filled; consult marker_availability_by_sample.csv."
    return adata, availability

def write_type_marker_qc(adata,availability,result_dir):
    # Type-marker summaries are defined only for cells that passed the
    # all-panel-zero filter. Keep this guard even though main() filters first.
    eligible=~adata.obs["all_panel_zero"].to_numpy()
    adata=adata[eligible].copy()
    obs=adata.obs.copy(); order=["0","1","2+"]
    grouped=obs.groupby(["sample","type_marker_count_group"],observed=True).size().rename("n_cells").reset_index()
    grouped["fraction_of_sample"]=grouped.n_cells/grouped.groupby("sample").n_cells.transform("sum")
    grouped.to_csv(result_dir/"type_marker_count_groups_by_sample.csv",index=False)
    exact=obs.groupby(["sample","n_type_markers_positive"],observed=True).size().rename("n_cells").reset_index()
    exact["fraction_of_sample"]=exact.n_cells/exact.groupby("sample").n_cells.transform("sum")
    exact.to_csv(result_dir/"number_type_markers_positive_by_sample.csv",index=False)
    combos=obs.groupby(["sample","positive_type_markers"],observed=True).size().rename("n_cells").reset_index()
    combos["fraction_of_sample"]=combos.n_cells/combos.groupby("sample").n_cells.transform("sum")
    combos["n_markers_in_combination"]=combos.positive_type_markers.map(lambda x: 0 if x=="none" else len(x.split(";")))
    combos.sort_values(["sample","n_cells"],ascending=[True,False]).to_csv(result_dir/"type_marker_combinations_by_sample.csv",index=False)
    combos[combos.n_markers_in_combination==1].to_csv(result_dir/"single_type_marker_cells_by_sample.csv",index=False)
    combos[combos.n_markers_in_combination>=2].to_csv(result_dir/"multi_type_marker_combinations_by_sample.csv",index=False)
    raw=np.asarray(adata.X); marker_rows=[]
    for sample in obs["sample"].cat.categories:
        cell_mask=obs["sample"].astype(str).eq(str(sample)).to_numpy(); measured={m for m in TYPE if m in availability.columns and bool(availability.loc[str(sample),m])}
        for marker in [m for m in PANEL if m in TYPE]:
            if marker not in measured: continue
            marker_rows.append({"sample":sample,"marker":marker,"n_positive":int((raw[cell_mask,adata.var_names.get_loc(marker)]>0).sum()),"n_cells":int(cell_mask.sum()),"fraction_positive":float((raw[cell_mask,adata.var_names.get_loc(marker)]>0).mean())})
    marker_df=pd.DataFrame(marker_rows); marker_df.to_csv(result_dir/"type_marker_positivity_by_sample.csv",index=False)
    table=grouped.pivot(index="sample",columns="type_marker_count_group",values="fraction_of_sample").reindex(columns=order,fill_value=0).fillna(0)
    ax=table.plot.bar(stacked=True,color=["#D9D9D9","#4C78A8","#E45756"],figsize=(10,5)); ax.set_ylabel("Fraction of retained cells"); ax.set_xlabel("Sample"); ax.legend(title="Positive type markers"); plt.tight_layout(); plt.savefig(result_dir/"type_marker_count_distribution_by_sample.png",dpi=300); plt.close()
    heat=marker_df.pivot(index="marker",columns="sample",values="fraction_positive")
    fig,ax=plt.subplots(figsize=(max(8,0.75*heat.shape[1]),max(5,0.38*heat.shape[0]))); im=ax.imshow(heat.to_numpy(),aspect="auto",vmin=0,vmax=1,cmap="Blues"); ax.set_xticks(range(heat.shape[1]),heat.columns,rotation=45,ha="right"); ax.set_yticks(range(heat.shape[0]),heat.index); ax.set_title("Fraction positive for cell-type markers"); fig.colorbar(im,ax=ax,label="Fraction positive"); plt.tight_layout(); plt.savefig(result_dir/"type_marker_positivity_heatmap_by_sample.png",dpi=300); plt.close()

def write_lineage_conflict_qc(adata,availability,result_dir):
    raw=pd.DataFrame(np.asarray(adata.X),index=adata.obs_names,columns=adata.var_names)
    quartiles=[]; conflict_cells=[]
    for lineage in LINEAGE:
        adata.obs[f"{lineage}_positive"]=False
        adata.obs[f"{lineage}_percentile_score"]=0.0
    for sample in adata.obs["sample"].cat.categories:
        idx=adata.obs_names[adata.obs["sample"].astype(str)==str(sample)]
        percentile=pd.DataFrame(0.0,index=idx,columns=adata.var_names)
        for marker in [m for m in TYPE if m in adata.var_names and bool(availability.loc[str(sample),m])]:
            values=raw.loc[idx,marker]; positive=values>0; pos=values[positive]
            quartiles.append({"sample":sample,"marker":marker,"marker_lineage":next((k for k,v in LINEAGE.items() if marker in v),"broad"),"n_cells":len(values),"n_positive":int(positive.sum()),"fraction_positive":float(positive.mean()),"positive_q25":float(pos.quantile(.25)) if len(pos) else np.nan,"positive_median":float(pos.quantile(.5)) if len(pos) else np.nan,"positive_q75":float(pos.quantile(.75)) if len(pos) else np.nan})
            if positive.any(): percentile.loc[pos.index,marker]=pos.rank(method="average",pct=True)
        for lineage,markers in LINEAGE.items():
            measured=[m for m in markers if m in adata.var_names and bool(availability.loc[str(sample),m])]
            if not measured: continue
            positive=raw.loc[idx,measured].gt(0).any(axis=1)
            # Mean percentile across available lineage markers captures both
            # intensity and breadth while remaining comparable between samples.
            score=percentile.loc[idx,measured].mean(axis=1)
            adata.obs.loc[idx,f"{lineage}_positive"]=positive.to_numpy()
            adata.obs.loc[idx,f"{lineage}_percentile_score"]=score.to_numpy()
    flags=pd.DataFrame({x:adata.obs[f"{x}_positive"].astype(bool) for x in LINEAGE},index=adata.obs_names)
    adata.obs["lineage_pattern"]=flags.apply(lambda r:"+".join(r.index[r.to_numpy()]) if r.any() else "no_specific_lineage",axis=1)
    adata.obs["epithelial_immune_conflict"]=flags.epithelial & flags.immune
    diff=adata.obs["epithelial_percentile_score"]-adata.obs["immune_percentile_score"]
    dominance=np.where(diff>0.10,"epithelial-dominant",np.where(diff<-.10,"immune-dominant","similar-strength"))
    adata.obs["epithelial_immune_dominance"]=np.where(adata.obs["epithelial_immune_conflict"],dominance,"not-a-conflict")
    pd.DataFrame(quartiles).sort_values(["sample","marker_lineage","marker"]).to_csv(result_dir/"type_marker_positive_quartiles_by_sample.csv",index=False)
    summary=adata.obs.groupby(["sample","lineage_pattern"],observed=True).size().rename("n_cells").reset_index(); summary["fraction_of_sample"]=summary.n_cells/summary.groupby("sample",observed=True).n_cells.transform("sum"); summary.to_csv(result_dir/"lineage_patterns_by_sample.csv",index=False)
    conflicts=adata.obs.loc[adata.obs["epithelial_immune_conflict"],["sample","age","positive_type_markers","epithelial_percentile_score","immune_percentile_score","stromal_percentile_score","epithelial_immune_dominance"]].copy(); conflicts.to_csv(result_dir/"epithelial_immune_conflict_cells.csv.gz",compression="gzip")
    conflict_summary=adata.obs.groupby("sample",observed=True).agg(n_cells=("sample","size"),n_epithelial_immune_conflicts=("epithelial_immune_conflict","sum")).reset_index(); conflict_summary["fraction_epithelial_immune_conflicts"]=conflict_summary.n_epithelial_immune_conflicts/conflict_summary.n_cells; conflict_summary.to_csv(result_dir/"epithelial_immune_conflicts_by_sample.csv",index=False)
    ax=conflict_summary.plot.bar(x="sample",y="fraction_epithelial_immune_conflicts",legend=False,color="#B279A2",figsize=(10,5)); ax.set_ylabel("Fraction of retained cells"); ax.set_title("Epithelial + immune marker conflicts"); plt.tight_layout(); plt.savefig(result_dir/"epithelial_immune_conflict_fraction_by_sample.png",dpi=300); plt.close()
    if len(conflicts):
        dom=conflicts.groupby(["sample","epithelial_immune_dominance"],observed=True).size().rename("n").reset_index(); dom["fraction"]=dom.n/dom.groupby("sample",observed=True).n.transform("sum"); tab=dom.pivot(index="sample",columns="epithelial_immune_dominance",values="fraction").fillna(0); ax=tab.plot.bar(stacked=True,figsize=(10,5),color={"epithelial-dominant":"#E45756","immune-dominant":"#4C78A8","similar-strength":"#BDBDBD"}); ax.set_ylabel("Fraction of conflict cells"); ax.set_title("Dominant lineage among epithelial–immune conflicts"); plt.tight_layout(); plt.savefig(result_dir/"epithelial_immune_conflict_dominance_by_sample.png",dpi=300); plt.close()

def main():
    parser = argparse.ArgumentParser(description=__doc__); root = Path(__file__).resolve().parents[3]
    parser.add_argument("--project-root",type=Path,default=root); parser.add_argument("--input-dir",type=Path); parser.add_argument("--output",type=Path); parser.add_argument("--strict-r-parity",action="store_true"); args=parser.parse_args()
    root=args.project_root.resolve(); input_dir=(args.input_dir or root/"data"/"codex").resolve(); output=(args.output or input_dir/"codex_cleaned_nonzero.h5ad").resolve(); sample_dir=input_dir/"cleaned_samples"; sample_dir.mkdir(parents=True,exist_ok=True); result_dir=root/"results"/"codex"/"preprocessing"/"python_cleaning"; result_dir.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.INFO,format="%(levelname)s: %(message)s"); items=[]; rows=[]
    for path in sorted(input_dir.rglob("*.xlsx")):
        try:
            x=clean(path); items.append(x); rows.append(dict(file=path.name,sample_id=x["sample"],status="ok",sheet=x["sheet"],n_cells=len(x["values"]),n_markers=len(x["markers"]),markers=";".join(x["markers"]),p16_source=x["p16_source"],message="")); logging.info("%s: %d cells, %d markers (%s)",x["sample"],len(x["values"]),len(x["markers"]),x["p16_source"])
        except Exception as exc: rows.append(dict(file=path.name,sample_id=sample_id(path),status="error",message=str(exc))); logging.error("%s: %s",path.name,exc)
    log=pd.DataFrame(rows); log.to_csv(result_dir/"codex_preprocessing_log.csv",index=False)
    if log.empty or (log.status=="error").any(): raise RuntimeError(f"Import errors; see {result_dir}")
    # Remove cells with no positive measured marker before writing any H5AD.
    filtered_items=[]; qc_rows=[]
    for x in items:
        keep=x["values"].gt(0).any(axis=1).to_numpy()
        n_original=len(keep); n_removed=int((~keep).sum())
        qc_rows.append({"sample":x["sample"],"n_cells_original":n_original,"n_all_panel_zero_removed":n_removed,"n_cells_retained":int(keep.sum()),"fraction_all_panel_zero_removed":n_removed/n_original,"percent_all_panel_zero_removed":100*n_removed/n_original})
        x=x.copy(); x["values"]=x["values"].iloc[keep].copy(); x["labels"]=x["labels"][keep]; x["ids"]=[v for v,k in zip(x["ids"],keep) if k]; x["spatial"]=x["spatial"][keep]
        filtered_items.append(x)
        sample_adata,_=assemble([x])
        # Per-sample files contain only channels genuinely measured in that sample.
        sample_adata=sample_adata[:,x["markers"]].copy()
        sample_adata.write_h5ad(sample_dir/f"{x['sample']}_cleaned_nonzero.h5ad",compression="gzip")
    adata,availability=assemble(filtered_items); availability.index.name="sample"; availability.to_csv(result_dir/"marker_availability_by_sample.csv")
    pd.DataFrame(qc_rows).to_csv(result_dir/"all_panel_zero_qc_by_sample.csv",index=False); write_type_marker_qc(adata,availability,result_dir); write_lineage_conflict_qc(adata,availability,result_dir); adata.obs.to_csv(result_dir/"cell_qc.csv.gz",compression="gzip")
    # Retained for backward compatibility; step 02 merges the per-sample H5ADs.
    adata.write_h5ad(output,compression="gzip"); primary=adata
    r_path=input_dir/"codex_preprocessing_log.csv"
    if r_path.exists():
        fields=["sample_id","n_cells","n_markers","p16_source"]; parity=pd.read_csv(r_path)[fields].merge(log[fields],on="sample_id",how="outer",suffixes=("_r","_python")); parity["match"]=(parity.n_cells_r==parity.n_cells_python)&(parity.n_markers_r==parity.n_markers_python)&(parity.p16_source_r==parity.p16_source_python); parity.to_csv(result_dir/"r_python_import_parity.csv",index=False)
        if not parity.match.fillna(False).all():
            if args.strict_r_parity: raise RuntimeError("R/Python import mismatch")
            logging.warning("R/Python counts differ; inspect parity CSV")
    logging.info("Wrote %s: %d nonzero cells x %d markers",output,primary.n_obs,primary.n_vars); logging.info("Wrote %d filtered per-sample H5ADs to %s",len(filtered_items),sample_dir)
if __name__ == "__main__": main()
