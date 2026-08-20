from __future__ import annotations

import json
import importlib.util
import pickle
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
import shapiq

SEED = 20260817
rng = np.random.default_rng(SEED)


def save_figure(fig, stem: Path, width=7.2, height=5.2):
    fig.set_size_inches(width, height)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")


def feature_group(encoded_name: str, originals: list[str]):
    matches = [x for x in originals if encoded_name == x or encoded_name.startswith(x + "_")]
    return max(matches, key=len) if matches else encoded_name


def explain_outcome(outcome, models, x, columns, dml_ids, outdir, figdir):
    m0 = models[f"{outcome}_m0"]
    m1 = models[f"{outcome}_m1"]
    display_outcome = "HD" if outcome == "HD_z" else "KDM-BAA"
    predict_tau = lambda a: m1.predict(np.asarray(a)) - m0.predict(np.asarray(a))

    bg_idx = rng.choice(len(x), size=min(100, len(x)), replace=False)
    ex_idx = rng.choice(len(x), size=min(300, len(x)), replace=False)
    background = x[bg_idx]
    explain_x = x[ex_idx]
    shap_path = outdir / f"SHAP_values_{outcome}.parquet"
    if shap_path.exists():
        shap_long = pd.read_parquet(shap_path)
    else:
        explainer = shap.Explainer(predict_tau, shap.maskers.Independent(background, max_samples=100), algorithm="permutation")
        sv = explainer(explain_x, max_evals=2 * x.shape[1] + 1, batch_size=20)
        vals = np.asarray(sv.values)
        shap_long = pd.DataFrame(vals, columns=columns).assign(dml_id=np.asarray(dml_ids)[ex_idx]).melt(id_vars="dml_id", var_name="encoded_feature", value_name="shap_value")
        shap_long["outcome"] = outcome
        shap_long["feature"] = shap_long.encoded_feature.map(lambda z: feature_group(z, models["feature_names"]))
        shap_long.to_parquet(shap_path, index=False)

    # Aggregate dummy columns to their frozen original feature for publication-level importance.
    importance = (shap_long.groupby("feature", as_index=False).shap_value
                  .apply(lambda z: np.mean(np.abs(z))).rename(columns={"shap_value": "mean_abs_SHAP"})
                  .sort_values("mean_abs_SHAP", ascending=False))
    importance["outcome"] = outcome
    importance.to_csv(outdir / f"SHAP_importance_{outcome}.csv", index=False)

    # Rebuild a SHAP Explanation from the cached long-form values. Base values are
    # recovered by the additivity identity for the actual Super Learner contrast.
    wide = shap_long.pivot(index="dml_id", columns="encoded_feature", values="shap_value").reindex(columns=columns, fill_value=0.0)
    id_to_row = {str(v): i for i, v in enumerate(dml_ids)}
    row_idx = np.asarray([id_to_row[str(v)] for v in wide.index])
    x_shap = x[row_idx]
    vals_shap = wide.to_numpy(dtype=float)
    pred_shap = predict_tau(x_shap)
    base_shap = pred_shap - vals_shap.sum(axis=1)
    explanation = shap.Explanation(values=vals_shap, base_values=base_shap, data=x_shap, feature_names=columns)

    plt.figure()
    shap.summary_plot(explanation.values, features=explanation.data, feature_names=columns, max_display=15,
                      color_bar_label="Encoded feature value", show=False)
    ax = plt.gca(); ax.set_title(f"{display_outcome}: SHAP beeswarm", loc="left", fontweight="bold")
    ax.set_xlabel("SHAP value for predicted conditional contrast")
    plt.figtext(.01, .005, "Prediction explanation; not an individual causal effect", fontsize=8, color="#555555")
    save_figure(plt.gcf(), figdir / f"Figure_SHAP_beeswarm_{outcome}", width=8.0, height=6.0); plt.close()

    shap.plots.heatmap(explanation, max_display=15, instance_order=explanation.sum(1), show=False)
    ax = plt.gca(); ax.set_title(f"{display_outcome}: SHAP heatmap", loc="left", fontweight="bold")
    plt.figtext(.01, .005, "Rows are the prespecified 300-person explanation sample; prediction explanation only", fontsize=8, color="#555555")
    save_figure(plt.gcf(), figdir / f"Figure_SHAP_heatmap_{outcome}", width=8.2, height=6.2); plt.close()

    # Representative local explanation: closest predicted contrast to the sample median.
    local_pos = int(np.argmin(np.abs(pred_shap - np.median(pred_shap))))
    local_abs = np.abs(vals_shap[local_pos]); keep = np.argsort(local_abs)[::-1][:12]
    other = np.setdiff1d(np.arange(len(columns)), keep)
    local_values = np.r_[vals_shap[local_pos, keep], vals_shap[local_pos, other].sum()]
    local_data = np.r_[x_shap[local_pos, keep], np.nan]
    local_names = [columns[i] for i in keep] + ["Other encoded features"]
    local_exp = shap.Explanation(values=local_values, base_values=base_shap[local_pos], data=local_data, feature_names=local_names)
    pd.DataFrame({"outcome":[outcome], "dml_id":[wide.index[local_pos]], "selection_rule":["closest predicted conditional contrast to sample median"],
                  "prediction":[pred_shap[local_pos]], "base_value":[base_shap[local_pos]]}).to_csv(outdir/f"SHAP_local_representative_{outcome}.csv",index=False)

    shap.plots.waterfall(local_exp, max_display=13, show=False)
    ax=plt.gca(); ax.set_title(f"{display_outcome}: representative SHAP waterfall", loc="left", fontweight="bold")
    plt.figtext(.01,.005,"Representative model prediction; not an individual treatment effect",fontsize=8,color="#555555")
    save_figure(plt.gcf(),figdir/f"Figure_SHAP_waterfall_{outcome}",width=8.0,height=6.2);plt.close()

    # Journal-readable signed force layout. It uses the same additive SHAP values
    # as the official force plot but avoids label collisions for long CHNS names.
    force_keep = np.argsort(np.abs(local_values))[::-1][:8]
    force_other = np.setdiff1d(np.arange(len(local_values)), force_keep)
    force_values = np.r_[local_values[force_keep], local_values[force_other].sum()]
    def force_short(s):
        return (s.replace("proportion_prior", "prop.prior").replace("number_", "n.")
                .replace("last_two_wave_pattern", "last-2-wave").replace("_", " "))[:25]
    force_names = [force_short(local_names[i]) for i in force_keep] + ["Other"]
    cumulative = [float(base_shap[local_pos])]
    for v in force_values: cumulative.append(cumulative[-1] + float(v))
    fig, ax = plt.subplots(figsize=(10.5, 3.6))
    for i, (name, value) in enumerate(zip(force_names, force_values), start=1):
        start, end = cumulative[i-1], cumulative[i]
        color = "#D55E00" if value >= 0 else "#0072B2"
        ax.barh(0, width=value, left=start, height=.32, color=color, edgecolor="white", linewidth=.8)
        ax.text((start+end)/2, 0, str(i), ha="center", va="center", fontsize=7, color="white", fontweight="bold")
    lo,hi=min(cumulative),max(cumulative);pad=max((hi-lo)*.12,.01);ax.set_xlim(lo-pad,hi+pad)
    ax.axvline(cumulative[0], color="#777777", ls="--", lw=.8); ax.axvline(cumulative[-1], color="black", ls="--", lw=.8)
    ax.text(cumulative[0], .30, f"base={cumulative[0]:.3f}", ha="center", fontsize=8, color="#666666")
    ax.text(cumulative[-1], .30, f"prediction={cumulative[-1]:.3f}", ha="center", fontsize=8, fontweight="bold")
    ax.set_ylim(-.35,.42); ax.set_yticks([]); ax.set_xlabel("Predicted conditional contrast")
    ax.set_title(f"{display_outcome}: representative SHAP force plot",loc="left",fontweight="bold")
    ax.spines[["left","right","top"]].set_visible(False)
    legend_lines=[f"{i}. {n}: {v:+.3f}" for i,(n,v) in enumerate(zip(force_names,force_values),start=1)]
    for i,line in enumerate(legend_lines):
        col=i%3;row=i//3;fig.text(.03+col*.33,.19-row*.055,line,fontsize=7,color="#D55E00" if force_values[i]>=0 else "#0072B2")
    fig.text(.01,.01,"Representative model prediction; red raises and blue lowers the predicted contrast",fontsize=8,color="#555555")
    fig.subplots_adjust(left=.08,right=.98,top=.82,bottom=.34); save_figure(fig,figdir/f"Figure_SHAP_force_{outcome}",width=10.5,height=3.6);plt.close(fig)

    top = importance.head(15).sort_values("mean_abs_SHAP")
    fig, ax = plt.subplots()
    ax.barh(top.feature, top.mean_abs_SHAP, color="#0072B2")
    ax.set_xlabel("Mean |SHAP value| for predicted conditional contrast")
    ax.set_ylabel("")
    ax.set_title(f"{display_outcome}: Super Learner contrast model")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0, -0.16, "Prediction explanation; not a causal-effect decomposition", transform=ax.transAxes, fontsize=8, color="#555555")
    fig.tight_layout()
    save_figure(fig, figdir / f"Figure_SHAP_{outcome}")
    plt.close(fig)

    # SHAP-IQ is evaluated on an 8-feature conditional slice of the actual SL contrast.
    # Non-selected encoded features are fixed at the background median; no surrogate is fitted.
    encoded_mean = shap_long.assign(abs_value=shap_long.shap_value.abs()).groupby("encoded_feature").abs_value.mean()
    encoded_imp = np.asarray([encoded_mean.get(c, 0.0) for c in columns])
    top_idx = np.argsort(encoded_imp)[::-1][:8]
    base = np.median(background, axis=0)
    def sliced_model(z):
        z = np.asarray(z)
        full = np.tile(base, (len(z), 1))
        full[:, top_idx] = z
        return predict_tau(full)
    sx = background[:, top_idx]
    iq_explainer = shapiq.TabularExplainer(model=sliced_model, data=sx, index="k-SII", max_order=2, random_state=SEED)
    iq_path = outdir / f"SHAPIQ_interactions_{outcome}.csv"
    if iq_path.exists():
        iq = pd.read_csv(iq_path)
    else:
        selected = explain_x[:5, top_idx]
        records = []
        for row_no, row in enumerate(selected):
            iv = iq_explainer.explain(row, budget=512)
            for interaction, pos in iv.interaction_lookup.items():
                if len(interaction) != 2:
                    continue
                records.append({"outcome": outcome, "row": row_no, "feature_1": columns[top_idx[interaction[0]]],
                                "feature_2": columns[top_idx[interaction[1]]], "kSII": float(iv.values[pos])})
        iq = pd.DataFrame(records)
        iq.to_csv(iq_path, index=False)
    agg = (iq.assign(abs_kSII=iq.kSII.abs()).groupby(["feature_1", "feature_2"], as_index=False)
           .agg(mean_abs_kSII=("abs_kSII", "mean"), signed_mean_kSII=("kSII", "mean"))
           .sort_values("mean_abs_kSII", ascending=False).head(12).sort_values("mean_abs_kSII"))
    def short_label(s):
        s = s.replace("_", " ").replace("proportion prior", "prop. prior").replace("province", "prov.")
        s = s.replace("last two wave pattern", "last-2-wave").replace("ever clean before 2006", "prior clean").replace("ever mixed before 2006", "prior mixed")
        return s[:34]
    labels = agg.feature_1.map(short_label) + " x " + agg.feature_2.map(short_label)
    colors = np.where(agg.signed_mean_kSII >= 0, "#D55E00", "#0072B2")
    fig, ax = plt.subplots()
    ax.barh(labels, agg.mean_abs_kSII, color=colors)
    ax.set_xlabel("Mean |k-SII| for predicted conditional contrast")
    ax.set_ylabel("")
    ax.set_title(f"{display_outcome}: SHAP-IQ pairwise interactions")
    ax.tick_params(axis="y", labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0, -0.16, "Top-8 encoded-feature conditional slice; prediction interaction, not effect modification", transform=ax.transAxes, fontsize=8, color="#555555")
    fig.tight_layout()
    save_figure(fig, figdir / f"Figure_SHAPIQ_{outcome}", width=8.0, height=5.4)
    plt.close(fig)
    return importance


def main(project):
    project = Path(project)
    modeldir = project / "causal_dml" / "models"; outdir = project / "causal_dml" / "results"; figdir = project / "causal_dml" / "figures"; logdir = project / "causal_dml" / "logs"
    figdir.mkdir(parents=True, exist_ok=True)
    # Models were serialized when the training script was executed as __main__.
    # Bind the identical class into this interpreter before unpickling.
    train_script = project / "causal_dml" / "scripts" / "21_superlearner_dml.py"
    spec = importlib.util.spec_from_file_location("dml_training", train_script)
    training = importlib.util.module_from_spec(spec)
    sys.modules["dml_training"] = training
    spec.loader.exec_module(training)
    sys.modules["__main__"].SuperLearner = training.SuperLearner
    with open(modeldir / "superlearner_full_models_imp1.pkl", "rb") as f: models = pickle.load(f)
    xf = pd.read_parquet(modeldir / "encoded_features_imp1.parquet")
    ids = xf.pop("dml_id").to_numpy(); columns = xf.columns.tolist(); x = xf.to_numpy(dtype=float)
    all_imp = []
    for outcome in ["HD_z", "KDM_BAA"]:
        all_imp.append(explain_outcome(outcome, models, x, columns, ids, outdir, figdir))
        print(f"EXPLANATION_DONE={outcome}", flush=True)
    pd.concat(all_imp, ignore_index=True).to_csv(outdir / "SHAP_importance_all_outcomes.csv", index=False)
    versions = {"shap": shap.__version__, "shapiq": shapiq.__version__}
    (logdir / "22_shap_shapiq.log").write_text("SHAP_SHAPIQ=PASS\n" + json.dumps(versions) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
