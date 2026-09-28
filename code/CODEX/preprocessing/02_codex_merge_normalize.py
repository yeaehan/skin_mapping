#!/usr/bin/env python
"""Merge cleaned CODEX samples and prepare log1p/scaled protein intensities."""
from __future__ import annotations
import argparse, logging
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    root=Path(__file__).resolve().parents[3]
    parser.add_argument("--project-root",type=Path,default=root)
    parser.add_argument("--input-dir",type=Path)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--clip-z",type=float,default=10.0)
    args=parser.parse_args(); root=args.project_root.resolve()
    input_dir=(args.input_dir or root/"data"/"codex"/"cleaned_samples").resolve()
    output=(args.output or root/"data"/"codex"/"codex_merged_log1p.h5ad").resolve()
    results=root/"results"/"codex"/"preprocessing"/"merge_normalize"; results.mkdir(parents=True,exist_ok=True); output.parent.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=logging.INFO,format="%(levelname)s: %(message)s")
    paths=sorted(input_dir.glob("*_cleaned_nonzero.h5ad"))
    if not paths: raise FileNotFoundError(f"No cleaned per-sample H5ADs in {input_dir}; run step 01 first")
    objects={}
    for path in paths:
        obj=ad.read_h5ad(path)
        samples=obj.obs["sample"].astype(str).unique()
        if len(samples)!=1: raise ValueError(f"{path.name} contains {len(samples)} sample labels")
        if obj.obs["all_panel_zero"].any(): raise ValueError(f"{path.name} still contains all-panel-zero cells")
        objects[samples[0]]=obj
    common=set.intersection(*(set(x.var_names) for x in objects.values()))
    common_markers=[x for x in next(iter(objects.values())).var_names if x in common]
    if not common_markers: raise ValueError("No markers are shared by every sample")
    input_qc=pd.DataFrame([{"sample":s,"n_cells_input":x.n_obs,"n_markers_input":x.n_vars} for s,x in objects.items()])
    input_qc.to_csv(results/"input_samples_qc.csv",index=False)
    merged=ad.concat({s:x[:,common_markers].copy() for s,x in objects.items()},join="inner",merge="same",index_unique=None)
    raw=np.asarray(merged.X,dtype=np.float32)
    finite=np.isfinite(raw).all(axis=1); nonnegative=(raw>=0).all(axis=1); shared_nonzero=(raw>0).any(axis=1)
    if not finite.all(): raise ValueError(f"Step-01 inputs contain {int((~finite).sum())} cells with non-finite values")
    if not nonnegative.all(): raise ValueError(f"Step-01 inputs contain {int((~nonnegative).sum())} cells with negative values")
    # Step 01 is authoritative for cell filtering. Cells positive only for an
    # optional marker can be zero over the shared-marker subset; flag them but
    # do not silently remove them during merge.
    merged.obs["all_shared_markers_zero"] = ~shared_nonzero
    filter_qc=pd.DataFrame({"criterion":["input_cells_from_step01","nonfinite_cells","negative_cells","all_shared_markers_zero_flagged","cells_removed_in_step02","output_cells"],"n_cells":[len(shared_nonzero),int((~finite).sum()),int((~nonnegative).sum()),int((~shared_nonzero).sum()),0,len(shared_nonzero)]})
    filter_qc.to_csv(results/"postmerge_cell_qc_no_additional_filtering.csv",index=False)
    merged.layers["raw_intensity"]=raw.copy()
    log1p=np.log1p(raw).astype(np.float32); merged.X=log1p; merged.layers["log1p"]=log1p.copy()
    means=log1p.mean(axis=0); std=log1p.std(axis=0); std[std==0]=1
    merged.layers["scaled"] = np.clip((log1p-means)/std,-args.clip_z,args.clip_z).astype(np.float32)
    merged.var["log1p_mean"]=means; merged.var["log1p_sd"]=std
    merged.uns["preprocessing"]={"cell_filter":"No cells removed in step 02; step 01 removed all-panel-zero cells","normalization":"log1p raw zero-thresholded protein intensity; no library-size normalization","scaled_layer":f"per-marker z-score of log1p, clipped to +/-{args.clip_z}","shared_markers_only":True,"all_shared_markers_zero_policy":"flagged in obs, not removed","n_input_samples":len(objects)}
    merged.write_h5ad(output,compression="gzip")
    pd.DataFrame({"marker":common_markers}).to_csv(results/"shared_markers_used_no_additional_filtering.csv",index=False)
    logging.info("Merged %d samples: %d cells x %d shared markers",len(objects),merged.n_obs,merged.n_vars)
    logging.info("X/log1p=log1p intensity; layers: raw_intensity, log1p, scaled")
    logging.info("Wrote %s",output)

if __name__=="__main__": main()
