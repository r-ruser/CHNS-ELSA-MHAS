from __future__ import annotations

import json
import math
import pickle
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from sklearn.base import clone
from sklearn.ensemble import (
    ExtraTreesClassifier, ExtraTreesRegressor,
    HistGradientBoostingClassifier, HistGradientBoostingRegressor,
    RandomForestClassifier, RandomForestRegressor,
)
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import log_loss, mean_squared_error, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("error", category=RuntimeWarning)
SEED = 20260817
EPS = 0.01


@dataclass
class SuperLearner:
    task: str
    learners: dict
    weights: np.ndarray | None = None
    fitted: list | None = None

    def _pred(self, model, x):
        return model.predict_proba(x)[:, 1] if self.task == "binary" else model.predict(x)

    def fit(self, x, y, groups, inner_splits=5):
        unique_groups = np.unique(groups)
        nsplit = min(inner_splits, len(unique_groups))
        if nsplit < 2:
            raise ValueError("Fewer than two groups for Super Learner cross-validation")
        oof = np.zeros((len(y), len(self.learners)))
        cv = GroupKFold(n_splits=nsplit)
        for train, valid in cv.split(x, y, groups):
            for k, model in enumerate(self.learners.values()):
                fit = clone(model).fit(x[train], y[train])
                oof[valid, k] = self._pred(fit, x[valid])
        if self.task == "binary":
            # NNLS on probabilities is the classical convex Super Learner metalearner.
            target = y.astype(float)
        else:
            target = y.astype(float)
        coef, _ = nnls(oof, target)
        if not np.isfinite(coef).all() or coef.sum() <= 0:
            coef = np.repeat(1.0 / oof.shape[1], oof.shape[1])
        self.weights = coef / coef.sum()
        self.fitted = [clone(m).fit(x, y) for m in self.learners.values()]
        pred = np.clip(oof @ self.weights, EPS, 1 - EPS) if self.task == "binary" else oof @ self.weights
        return pred

    def predict(self, x):
        p = np.column_stack([self._pred(m, x) for m in self.fitted]) @ self.weights
        return np.clip(p, EPS, 1 - EPS) if self.task == "binary" else p


def learner_library(task: str):
    if task == "binary":
        return {
            "logistic": make_pipeline(StandardScaler(), LogisticRegression(C=1.0, max_iter=3000, random_state=SEED)),
            "random_forest": RandomForestClassifier(n_estimators=350, min_samples_leaf=20, max_features="sqrt", n_jobs=-1, random_state=SEED),
            "extra_trees": ExtraTreesClassifier(n_estimators=350, min_samples_leaf=20, max_features="sqrt", n_jobs=-1, random_state=SEED + 1),
            "hist_gradient_boosting": HistGradientBoostingClassifier(max_iter=180, learning_rate=.04, max_leaf_nodes=15, l2_regularization=1.0, random_state=SEED),
        }
    return {
        "ridge": make_pipeline(StandardScaler(), Ridge(alpha=1.0)),
        "random_forest": RandomForestRegressor(n_estimators=350, min_samples_leaf=15, max_features=.7, n_jobs=-1, random_state=SEED),
        "extra_trees": ExtraTreesRegressor(n_estimators=350, min_samples_leaf=15, max_features=.7, n_jobs=-1, random_state=SEED + 1),
        "hist_gradient_boosting": HistGradientBoostingRegressor(max_iter=180, learning_rate=.04, max_leaf_nodes=15, l2_regularization=1.0, random_state=SEED),
    }


def design_matrix(df, frozen_columns=None):
    x = df.copy()
    for c in x.columns:
        if x[c].dtype == "object" or str(x[c].dtype).startswith("category"):
            x[c] = x[c].astype(str).fillna("Missing")
    x = pd.get_dummies(x, drop_first=False, dtype=float)
    if frozen_columns is not None:
        x = x.reindex(columns=frozen_columns, fill_value=0.0)
    x = x.replace([np.inf, -np.inf], np.nan)
    x = x.fillna(x.median(numeric_only=True)).fillna(0.0)
    arr = x.to_numpy(dtype=float)
    return x.columns.tolist(), arr


def cluster_variance(num, h, psi, groups):
    den = h.sum()
    ifv = num - psi * h
    cs = pd.Series(ifv).groupby(pd.Series(groups).astype(str)).sum().to_numpy()
    g = len(cs)
    return (g / (g - 1.0)) * np.sum((cs - cs.mean()) ** 2) / (den ** 2)


def ess(w):
    w = np.asarray(w, dtype=float)
    return float(w.sum() ** 2 / np.sum(w ** 2))


def run_imputation(df, feature_names, imp, outdir, modeldir):
    rawx = df[feature_names].copy()
    columns, x = design_matrix(rawx)
    a = df.A.to_numpy(dtype=int)
    r = df.R.to_numpy(dtype=int)
    groups = df.commid2004.astype(str).to_numpy()
    outer = GroupKFold(n_splits=5)
    e = np.zeros(len(df)); rp = np.zeros(len(df))
    outcome_pred = {y: {"m0": np.zeros(len(df)), "m1": np.zeros(len(df))} for y in ["HD_z", "KDM_BAA"]}
    fold_id = np.zeros(len(df), dtype=int)
    diag = []
    weights_records = []

    for fold, (train, valid) in enumerate(outer.split(x, a, groups), start=1):
        fold_id[valid] = fold
        sl_e = SuperLearner("binary", learner_library("binary"))
        pe_train = sl_e.fit(x[train], a[train], groups[train])
        e[valid] = sl_e.predict(x[valid])
        sl_r = SuperLearner("binary", learner_library("binary"))
        pr_train = sl_r.fit(np.column_stack([x[train], a[train]]), r[train], groups[train])
        rp[valid] = sl_r.predict(np.column_stack([x[valid], a[valid]]))
        diag.append({"imputation": imp, "fold": fold, "n_train": len(train), "n_valid": len(valid),
                     "treatment_auc_train_oof": roc_auc_score(a[train], pe_train),
                     "treatment_logloss_train_oof": log_loss(a[train], pe_train),
                     "observation_auc_train_oof": roc_auc_score(r[train], pr_train),
                     "observation_logloss_train_oof": log_loss(r[train], pr_train)})
        for name, sl in [("treatment", sl_e), ("observation", sl_r)]:
            for learner, weight in zip(sl.learners.keys(), sl.weights):
                weights_records.append({"imputation": imp, "fold": fold, "nuisance": name, "outcome": "", "arm": "", "learner": learner, "weight": weight})
        for yname in outcome_pred:
            y = df[yname].to_numpy(dtype=float)
            for arm in [0, 1]:
                tr = train[(a[train] == arm) & (r[train] == 1) & np.isfinite(y[train])]
                if len(tr) < 100:
                    raise ValueError(f"Insufficient observed outcome rows: imp={imp}, fold={fold}, y={yname}, arm={arm}, n={len(tr)}")
                sl_y = SuperLearner("continuous", learner_library("continuous"))
                py_train = sl_y.fit(x[tr], y[tr], groups[tr])
                outcome_pred[yname][f"m{arm}"][valid] = sl_y.predict(x[valid])
                diag[-1][f"{yname}_arm{arm}_rmse_train_oof"] = math.sqrt(mean_squared_error(y[tr], py_train))
                for learner, weight in zip(sl_y.learners.keys(), sl_y.weights):
                    weights_records.append({"imputation": imp, "fold": fold, "nuisance": "outcome", "outcome": yname, "arm": arm, "learner": learner, "weight": weight})

    e = np.clip(e, EPS, 1 - EPS); rp = np.clip(rp, EPS, 1 - EPS)
    h = e * (1 - e)
    ow = np.where(a == 1, 1 - e, e)
    results = []
    predout = pd.DataFrame({"dml_id": df.dml_id, "imputation": imp, "fold": fold_id, "A": a, "R": r,
                            "commid2004": groups, "e_hat": e, "r_hat": rp, "overlap_weight": ow,
                            "observed_combined_weight": r * ow / rp})
    for yname in outcome_pred:
        y = df[yname].to_numpy(dtype=float)
        m0 = outcome_pred[yname]["m0"]; m1 = outcome_pred[yname]["m1"]
        resid = np.where(a == 1, (1 - e) * (y - m1), -e * (y - m0))
        resid = np.where((r == 1) & np.isfinite(y), resid / rp, 0.0)
        num = h * (m1 - m0) + resid
        psi = num.sum() / h.sum()
        var = cluster_variance(num, h, psi, groups)
        se = math.sqrt(var)
        results.append({"imputation": imp, "outcome": yname, "estimate": psi, "variance": var, "SE": se,
                        "lower95": psi - 1.95996398454 * se, "upper95": psi + 1.95996398454 * se,
                        "n": len(df), "observed_n": int(np.sum((r == 1) & np.isfinite(y))),
                        "communities": len(np.unique(groups)), "estimand": "ATO"})
        predout[f"{yname}_m0"] = m0; predout[f"{yname}_m1"] = m1
        predout[f"{yname}_conditional_contrast"] = m1 - m0
        predout[f"{yname}_orthogonal_numerator"] = num

    diagnostics = pd.DataFrame(diag)
    diagnostics["e_min"] = e.min(); diagnostics["e_p01"] = np.quantile(e, .01); diagnostics["e_p99"] = np.quantile(e, .99); diagnostics["e_max"] = e.max()
    diagnostics["r_min"] = rp.min(); diagnostics["r_p01"] = np.quantile(rp, .01); diagnostics["r_p99"] = np.quantile(rp, .99); diagnostics["r_max"] = rp.max()
    diagnostics["ess_clean_overlap"] = ess(ow[a == 1]); diagnostics["ess_anysolid_overlap"] = ess(ow[a == 0])
    diagnostics["ess_clean_observed_combined"] = ess((ow / rp)[(a == 1) & (r == 1)])
    diagnostics["ess_anysolid_observed_combined"] = ess((ow / rp)[(a == 0) & (r == 1)])

    predout.to_parquet(outdir / f"dml_crossfit_predictions_imp{imp}.parquet", index=False)
    pd.DataFrame(results).to_csv(outdir / f"dml_effects_imp{imp}.csv", index=False)
    diagnostics.to_csv(outdir / f"dml_diagnostics_imp{imp}.csv", index=False)
    pd.DataFrame(weights_records).to_csv(outdir / f"superlearner_weights_imp{imp}.csv", index=False)

    # Full-data nuisance fits on imputation 1 are interpretation models only.
    if imp == 1:
        full_models = {"feature_names": feature_names, "encoded_columns": columns}
        sl_e = SuperLearner("binary", learner_library("binary")); sl_e.fit(x, a, groups); full_models["treatment"] = sl_e
        sl_r = SuperLearner("binary", learner_library("binary")); sl_r.fit(np.column_stack([x, a]), r, groups); full_models["observation"] = sl_r
        for yname in outcome_pred:
            y = df[yname].to_numpy(dtype=float)
            for arm in [0, 1]:
                ii = (a == arm) & (r == 1) & np.isfinite(y)
                sl = SuperLearner("continuous", learner_library("continuous")); sl.fit(x[ii], y[ii], groups[ii])
                full_models[f"{yname}_m{arm}"] = sl
        with open(modeldir / "superlearner_full_models_imp1.pkl", "wb") as f:
            pickle.dump(full_models, f)
        pd.DataFrame(x, columns=columns).assign(dml_id=df.dml_id.values).to_parquet(modeldir / "encoded_features_imp1.parquet", index=False)

    return results


def rubin_pool(frame):
    out = []
    for outcome, z in frame.groupby("outcome"):
        q = z.estimate.to_numpy(); u = z.variance.to_numpy(); m = len(z)
        qb = q.mean(); within = u.mean(); between = q.var(ddof=1); total = within + (1 + 1/m) * between
        se = math.sqrt(total); p = math.erfc(abs(qb / se) / math.sqrt(2))
        out.append({"outcome": outcome, "method": "Cross-fitted Super Learner DML",
                    "contrast": "clean-only 2006 versus continued any-solid 2006",
                    "estimand": "ATO in the 2004 any-solid overlap population",
                    "estimate": qb, "SE": se, "lower95": qb - 1.95996398454*se, "upper95": qb + 1.95996398454*se,
                    "p_value": p, "within_variance": within, "between_variance": between, "m": m})
    return pd.DataFrame(out)


def main(project):
    project = Path(project)
    datadir = project / "causal_dml" / "data"; outdir = project / "causal_dml" / "results"; modeldir = project / "causal_dml" / "models"; logdir = project / "causal_dml" / "logs"
    for p in [outdir, modeldir, logdir]: p.mkdir(parents=True, exist_ok=True)
    feature_names = [x.strip() for x in (datadir / "frozen_feature_list.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
    all_results = []
    for imp in range(1, 6):
        df = pd.read_parquet(datadir / f"dml_imputation_{imp}.parquet")
        all_results.extend(run_imputation(df, feature_names, imp, outdir, modeldir))
        print(f"DML_IMPUTATION_DONE={imp}/5", flush=True)
    perimp = pd.DataFrame(all_results); perimp.to_csv(outdir / "dml_effects_all_imputations.csv", index=False)
    pooled = rubin_pool(perimp); pooled.to_csv(outdir / "Table_DML_primary_effects.csv", index=False)
    versions = {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__}
    (logdir / "software_versions.json").write_text(json.dumps(versions, indent=2), encoding="utf-8")
    (logdir / "21_superlearner_dml.log").write_text("SUPER_LEARNER_DML=PASS\nMICE=5; outer community GroupKFold=5; epsilon=0.01; cluster-robust orthogonal score.\n", encoding="utf-8")
    print(pooled.to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
