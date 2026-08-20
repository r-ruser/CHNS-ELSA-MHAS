from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import shap, shapiq, sklearn


def ess(w):
    w = np.asarray(w, float)
    return float(w.sum() ** 2 / np.sum(w ** 2))


def smd_table(z, features, weight_name, imp):
    x = pd.get_dummies(z[features], drop_first=False, dtype=float).fillna(0)
    a = z.A.to_numpy(); w = z[weight_name].to_numpy(float)
    rows = []
    for c in x:
        v = x[c].to_numpy(float)
        m1 = np.average(v[a == 1], weights=w[a == 1]); m0 = np.average(v[a == 0], weights=w[a == 0])
        s = np.sqrt((np.average((v[a == 1]-m1)**2, weights=w[a == 1]) + np.average((v[a == 0]-m0)**2, weights=w[a == 0]))/2)
        rows.append({"imputation": imp, "weight": weight_name, "covariate": c, "SMD": (m1-m0)/s if s > 1e-12 else 0.0})
    return rows


def main(project):
    p = Path(project) / "causal_dml"; dat = p/"data"; res=p/"results"; fig=p/"figures"; logs=p/"logs"
    features = (dat/"frozen_feature_list.txt").read_text(encoding="utf-8").splitlines()
    balances=[]; essrows=[]
    for j in range(1,6):
        d=pd.read_parquet(dat/f"dml_imputation_{j}.parquet")
        q=pd.read_parquet(res/f"dml_crossfit_predictions_imp{j}.parquet")
        z=d.merge(q[["dml_id","overlap_weight","observed_combined_weight"]],on="dml_id",validate="one_to_one")
        balances += smd_table(z,features,"overlap_weight",j)
        balances += smd_table(z,features,"observed_combined_weight",j)
        for wn in ["overlap_weight","observed_combined_weight"]:
            for a in [0,1]: essrows.append({"imputation":j,"weight":wn,"A":a,"ESS":ess(z.loc[z.A==a,wn])})
    bal=pd.DataFrame(balances); bal.to_csv(res/"DML_balance_diagnostics.csv",index=False)
    summary=(bal.assign(abs_SMD=bal.SMD.abs()).groupby(["imputation","weight"],as_index=False)
             .agg(max_abs_SMD=("abs_SMD","max"),n_ge_0_10=("abs_SMD",lambda x:int((x>=.10).sum()))))
    summary.to_csv(res/"DML_balance_summary.csv",index=False)
    pd.DataFrame(essrows).to_csv(res/"DML_ESS_summary.csv",index=False)

    primary=pd.read_csv(res/"Table_DML_primary_effects.csv")
    assert set(primary.outcome)=={"HD_z","KDM_BAA"} and (primary.m==5).all()
    for y in ["HD_z","KDM_BAA"]:
        assert len(pd.read_parquet(res/f"SHAP_values_{y}.parquet"))>0
        assert len(pd.read_csv(res/f"SHAPIQ_interactions_{y}.csv"))>0
        for stem in [f"Figure_SHAP_{y}",f"Figure_SHAPIQ_{y}",f"Figure_SHAP_beeswarm_{y}",
                     f"Figure_SHAP_heatmap_{y}",f"Figure_SHAP_force_{y}",f"Figure_SHAP_waterfall_{y}"]:
            for ext in ["svg","pdf","png","tiff"]: assert (fig/f"{stem}.{ext}").stat().st_size>1000
        local=pd.read_csv(res/f"SHAP_local_representative_{y}.csv")
        assert len(local)==1 and local.selection_rule.iloc[0]=="closest predicted conditional contrast to sample median"
    manifest=[]
    for f in sorted([x for x in p.rglob("*") if x.is_file() and "internal_restricted" not in x.name and "crosswalk" not in x.name]):
        manifest.append({"path":str(f.relative_to(p)).replace("\\","/"),"bytes":f.stat().st_size,"md5":hashlib.md5(f.read_bytes()).hexdigest()})
    pd.DataFrame(manifest).to_csv(p/"public_code_output_manifest_md5.csv",index=False)
    note={"validation":"PASS","python":sys.version,"sklearn":sklearn.__version__,"shap":shap.__version__,"shapiq":shapiq.__version__,
          "balance_interpretation":"DML orthogonality does not impose exact weighted balance; diagnostics are reported, not used as a post-outcome tuning gate.",
          "max_treatment_overlap_SMD":float(summary.loc[summary.weight=="overlap_weight","max_abs_SMD"].max()),
          "max_observed_combined_SMD":float(summary.loc[summary.weight=="observed_combined_weight","max_abs_SMD"].max())}
    (logs/"23_diagnostics_validate.json").write_text(json.dumps(note,indent=2),encoding="utf-8")
    print(json.dumps(note,indent=2))

if __name__=="__main__": main(sys.argv[1] if len(sys.argv)>1 else ".")
