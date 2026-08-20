#!/usr/bin/env python
"""Read-only feasibility audit for CHNS fuel trajectories and 2009 biomarkers.

No row-level derived data are written. Outputs contain aggregate counts and a
variable dictionary only, under C2/results and C2/logs.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
WORKSPACE = PROJECT.parent
RESULTS = PROJECT / "results"
LOGS = PROJECT / "logs"
RESULTS.mkdir(exist_ok=True)
LOGS.mkdir(exist_ok=True)

logging.basicConfig(
    filename=LOGS / "01_audit_chns.log",
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    encoding="utf-8",
    filemode="w",
)


def locate_one(pattern: str) -> Path:
    roots = list(WORKSPACE.glob("CHNS*/A023*/CHNS-STATA*"))
    if len(roots) != 1:
        raise RuntimeError(f"Expected one CHNS-STATA root, found {len(roots)}")
    hits = list(roots[0].glob(pattern))
    if len(hits) != 1:
        raise RuntimeError(f"Expected one file for {pattern}, found {len(hits)}")
    return hits[0]


FILES = {
    "asset": locate_one("Master_Asset*STATA/asset_12.dta"),
    "biomarker": locate_one("biomarker_09-STATA/biomarker_09.dta"),
    "master_id": locate_one("Master_ID*STATA/mast_pub_12.dta"),
    "surveys": locate_one("Master_ID*STATA/surveys_pub_12.dta"),
    "physical_exam": locate_one("Master_PE_PA*STATA/pexam_pub_12.dta"),
}


def read_stata(path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    reader = pd.io.stata.StataReader(str(path), convert_categoricals=False)
    labels = reader.variable_labels()
    return reader.read(), labels


logging.info("Starting CHNS read-only audit")
frames: dict[str, pd.DataFrame] = {}
labels: dict[str, dict[str, str]] = {}
inventory = []
for key, path in FILES.items():
    df, lab = read_stata(path)
    frames[key], labels[key] = df, lab
    inventory.append(
        {
            "dataset": key,
            "path": str(path),
            "bytes": path.stat().st_size,
            "rows": len(df),
            "columns": len(df.columns),
        }
    )
    logging.info("Read %s rows=%d cols=%d", key, len(df), len(df.columns))
pd.DataFrame(inventory).to_csv(RESULTS / "file_inventory.csv", index=False, encoding="utf-8-sig")

asset = frames["asset"]
bio = frames["biomarker"]
master = frames["master_id"]
surveys = frames["surveys"]
pexam = frames["physical_exam"]

# Uniqueness checks at documented keys.
key_checks = pd.DataFrame(
    [
        ["asset", "hhid+wave", "all waves", len(asset), asset.duplicated(["hhid", "wave"]).sum()],
        ["asset", "hhid+wave", "through 2009", int((asset.wave <= 2009).sum()), asset.loc[asset.wave <= 2009].duplicated(["hhid", "wave"]).sum()],
        ["biomarker", "IDind", "2009", len(bio), bio.duplicated(["IDind"]).sum()],
        ["master_id", "Idind", "master", len(master), master.duplicated(["Idind"]).sum()],
        ["surveys", "Idind+wave", "all waves", len(surveys), surveys.duplicated(["Idind", "wave"]).sum()],
        ["physical_exam", "IDind+wave", "all waves", len(pexam), pexam.duplicated(["IDind", "wave"]).sum()],
    ],
    columns=["dataset", "key", "scope", "rows", "duplicate_key_rows"],
)
key_checks.to_csv(RESULTS / "key_uniqueness.csv", index=False, encoding="utf-8-sig")

# Fuel codebook is taken from the local asset_12.pdf frequency table.
fuel_labels = {
    1: "Coal",
    2: "Electricity",
    3: "Kerosene",
    4: "Liquefied petroleum gas",
    5: "Natural gas",
    6: "Wood/sticks/straw",
    7: "Charcoal",
    8: "Other",
}
fuel_class = {1: "solid", 2: "clean", 3: "other", 4: "clean", 5: "clean", 6: "solid", 7: "solid", 8: "other"}

fuel_rows = []
for wave, group in asset[asset["wave"] <= 2009].groupby("wave"):
    counts = group["L8_1"].value_counts(dropna=False)
    for code, n in counts.items():
        if pd.isna(code):
            label, cls = "Missing", "missing"
        else:
            label = fuel_labels.get(int(code), "Unexpected code")
            cls = fuel_class.get(int(code), "unexpected")
        fuel_rows.append([int(wave), code, label, cls, int(n)])
fuel_distribution = pd.DataFrame(fuel_rows, columns=["wave", "code", "label", "class", "n_households"])
fuel_distribution.to_csv(RESULTS / "fuel_distribution_by_wave.csv", index=False, encoding="utf-8-sig")

# Link individual wave records to household fuel; this accommodates household changes.
s = surveys.loc[surveys["wave"] <= 2009, ["Idind", "hhid", "wave", "age", "urban", "stratum"]].copy()
a = asset.loc[asset["wave"] <= 2009, ["hhid", "wave", "L8_1", "L8_2"]].copy()
long = s.merge(a, on=["hhid", "wave"], how="left", validate="many_to_one")
long["fuel_state"] = long["L8_1"].map(fuel_class)

# Main exposure accounts for fuel stacking. L8_2 code 0/NaN is treated as no
# reported secondary fuel; 3 (kerosene) and 8 (other) remain unclassified.
long["secondary_state"] = long["L8_2"].map(fuel_class)
long.loc[long["L8_2"].isna() | long["L8_2"].eq(0), "secondary_state"] = "none"
long["fuel_stacking"] = (
    long["fuel_state"].isin(["solid", "clean"])
    & long["secondary_state"].isin(["solid", "clean"])
    & long["fuel_state"].ne(long["secondary_state"])
)
long["effective_fuel_state"] = "other_or_unknown"
long.loc[
    long["fuel_state"].eq("clean") & long["secondary_state"].isin(["clean", "none"]),
    "effective_fuel_state",
] = "clean_only"
long.loc[
    long["fuel_state"].eq("solid") | long["secondary_state"].eq("solid"),
    "effective_fuel_state",
] = "any_solid"

bio_ids = set(bio["IDind"].dropna().astype("int64"))
long = long[long["Idind"].astype("int64").isin(bio_ids)].copy()
valid = long[long["fuel_state"].isin(["solid", "clean"])].sort_values(["Idind", "wave"])


def classify_strict(group: pd.DataFrame) -> str:
    states = group["fuel_state"].tolist()
    if len(states) < 2:
        return "fewer_than_2_valid_waves"
    if all(x == "solid" for x in states):
        return "persistent_solid"
    if all(x == "clean" for x in states):
        return "persistent_clean"
    transitions = sum(a != b for a, b in zip(states[:-1], states[1:]))
    if states[0] == "solid" and states[-1] == "clean" and transitions == 1:
        return "solid_to_clean_monotonic"
    if states[0] == "clean" and states[-1] == "solid" and transitions == 1:
        return "clean_to_solid_monotonic"
    return "complex_switching"


trajectory = valid.groupby("Idind", sort=False).apply(classify_strict, include_groups=False).rename("strict_trajectory")
valid_summary = valid.groupby("Idind").agg(
    n_valid_fuel_waves=("wave", "size"), first_valid_wave=("wave", "min"), last_valid_wave=("wave", "max")
)
valid_summary["has_valid_2009_fuel"] = valid.groupby("Idind")["wave"].apply(lambda x: bool((x == 2009).any()))
person_traj = valid_summary.join(trajectory)

# Add all biomarker participants so absent histories are counted explicitly.
person_traj = pd.DataFrame(index=pd.Index(sorted(bio_ids), name="Idind")).join(person_traj)
person_traj["strict_trajectory"] = person_traj["strict_trajectory"].fillna("no_valid_fuel_wave")
person_traj["n_valid_fuel_waves"] = person_traj["n_valid_fuel_waves"].fillna(0).astype(int)
person_traj["has_valid_2009_fuel"] = person_traj["has_valid_2009_fuel"].eq(True)
traj_counts = person_traj.groupby(["strict_trajectory", "has_valid_2009_fuel"], dropna=False).size().rename("n").reset_index()
traj_counts.to_csv(RESULTS / "trajectory_counts.csv", index=False, encoding="utf-8-sig")

wave_history = person_traj["n_valid_fuel_waves"].value_counts().sort_index().rename_axis("n_valid_fuel_waves").reset_index(name="n")
wave_history.to_csv(RESULTS / "fuel_history_length.csv", index=False, encoding="utf-8-sig")

# Effective trajectories using both primary and secondary cooking fuel.
effective_valid = long[long["effective_fuel_state"].isin(["any_solid", "clean_only"])].sort_values(["Idind", "wave"])


def classify_effective(group: pd.DataFrame) -> str:
    states = group["effective_fuel_state"].tolist()
    if len(states) < 2:
        return "fewer_than_2_valid_waves"
    if all(x == "any_solid" for x in states):
        return "persistent_any_solid"
    if all(x == "clean_only" for x in states):
        return "persistent_clean_only"
    transitions = sum(a != b for a, b in zip(states[:-1], states[1:]))
    if states[0] == "any_solid" and states[-1] == "clean_only" and transitions == 1:
        return "solid_to_clean_only_monotonic"
    if states[0] == "clean_only" and states[-1] == "any_solid" and transitions == 1:
        return "clean_only_to_solid_monotonic"
    return "complex_switching"


effective_summary = effective_valid.groupby("Idind").agg(
    n_valid_effective_waves=("wave", "size"),
    n_valid_prior_waves=("wave", lambda x: int((x < 2009).sum())),
    first_valid_effective_wave=("wave", "min"),
    last_valid_effective_wave=("wave", "max"),
)
effective_summary["has_known_effective_2009"] = effective_valid.groupby("Idind")["wave"].apply(lambda x: bool((x == 2009).any()))
effective_summary["stacking_ever"] = long.groupby("Idind")["fuel_stacking"].any()
effective_summary = effective_summary.join(
    effective_valid.groupby("Idind", sort=False).apply(classify_effective, include_groups=False).rename("effective_trajectory")
)
effective_summary = pd.DataFrame(index=pd.Index(sorted(bio_ids), name="Idind")).join(effective_summary)
effective_summary["effective_trajectory"] = effective_summary["effective_trajectory"].fillna("no_valid_effective_wave")
effective_summary["n_valid_effective_waves"] = effective_summary["n_valid_effective_waves"].fillna(0).astype(int)
effective_summary["n_valid_prior_waves"] = effective_summary["n_valid_prior_waves"].fillna(0).astype(int)
effective_summary["has_known_effective_2009"] = effective_summary["has_known_effective_2009"].eq(True)
effective_summary["stacking_ever"] = effective_summary["stacking_ever"].eq(True)
effective_summary.groupby(["effective_trajectory", "has_known_effective_2009"], dropna=False).size().rename("n").reset_index().to_csv(
    RESULTS / "effective_trajectory_counts.csv", index=False, encoding="utf-8-sig"
)

stacking_rows = []
for wave, group in long.groupby("wave"):
    stacking_rows.append([
        int(wave), len(group), int(group["fuel_stacking"].sum()),
        int(group["effective_fuel_state"].eq("any_solid").sum()),
        int(group["effective_fuel_state"].eq("clean_only").sum()),
        int(group["effective_fuel_state"].eq("other_or_unknown").sum()),
    ])
pd.DataFrame(
    stacking_rows,
    columns=["wave", "n_biomarker_participant_records", "n_mixed_clean_solid_stacking", "n_any_solid", "n_clean_only", "n_other_or_unknown"],
).to_csv(RESULTS / "fuel_stacking_by_wave.csv", index=False, encoding="utf-8-sig")

# 2009 demographic and physical-exam links.
s09 = surveys.loc[surveys["wave"] == 2009, ["Idind", "hhid", "age", "urban", "stratum"]].copy()
m = master[["Idind", "gender", "WEST_DOB_Y"]].copy()
p09 = pexam.loc[pexam["wave"] == 2009, [
    "IDind", "SYSTOL1", "SYSTOL2", "SYSTOL3", "DIASTOL1", "DIASTOL2", "DIASTOL3", "height", "weight", "U10"
]].copy()
p09["sbp_mean"] = p09[["SYSTOL1", "SYSTOL2", "SYSTOL3"]].mean(axis=1)
p09["dbp_mean"] = p09[["DIASTOL1", "DIASTOL2", "DIASTOL3"]].mean(axis=1)
p09["bmi"] = p09["weight"] / (p09["height"] / 100) ** 2

linked = bio.rename(columns={"IDind": "Idind"}).merge(s09, on="Idind", how="left", validate="one_to_one")
linked = linked.merge(m, on="Idind", how="left", validate="one_to_one")
linked = linked.merge(p09.rename(columns={"IDind": "Idind"})[["Idind", "sbp_mean", "dbp_mean", "bmi", "U10"]], on="Idind", how="left", validate="one_to_one")
linked = linked.merge(person_traj.reset_index(), on="Idind", how="left", validate="one_to_one")
linked = linked.merge(effective_summary.reset_index(), on="Idind", how="left", validate="one_to_one")

labs = ["HS_CRP", "alb", "cre", "glucose", "hgb"]
expanded = ["HS_CRP", "alb", "cre", "glucose", "HbA1c", "tc", "HDL_C", "tg", "hgb", "wbc", "plt", "ua"]
for c in labs + expanded + ["age", "gender", "sbp_mean", "dbp_mean", "bmi", "U10"]:
    if c not in linked:
        raise KeyError(f"Expected variable missing after link: {c}")

lab_rows = []
for c in expanded + ["sbp_mean", "dbp_mean", "bmi", "U10"]:
    source = "biomarker_09.dta" if c in bio.columns else "pexam_pub_12.dta"
    label = labels["biomarker"].get(c, labels["physical_exam"].get(c, c))
    lab_rows.append([c, source, label, int(linked[c].notna().sum()), float(linked[c].notna().mean())])
pd.DataFrame(lab_rows, columns=["variable", "source", "label", "n_nonmissing", "proportion_nonmissing"]).to_csv(
    RESULTS / "candidate_biomarker_completeness.csv", index=False, encoding="utf-8-sig"
)

complete5 = linked[labs].notna().all(axis=1)
complete12 = linked[expanded].notna().all(axis=1)
adult = linked["age"].ge(18)
age_sex = linked[["age", "gender"]].notna().all(axis=1)
valid2 = linked["n_valid_fuel_waves"].fillna(0).ge(2)
valid2009 = linked["has_valid_2009_fuel"].fillna(False)
known_effective_2009 = linked["has_known_effective_2009"].eq(True)
two_prior_effective = linked["n_valid_prior_waves"].fillna(0).ge(2)

# Descriptive age associations inform (but do not automatically select) KDM candidates.
age_assoc_rows = []
adult_rows = linked["age"].ge(18) & linked["age"].notna()
for c in expanded + ["sbp_mean", "dbp_mean", "bmi", "U10"]:
    pair = linked.loc[adult_rows, ["age", c]].dropna()
    age_assoc_rows.append([c, len(pair), pair["age"].corr(pair[c], method="spearman")])
pd.DataFrame(age_assoc_rows, columns=["variable", "n_pairwise", "spearman_rho_with_age"]).to_csv(
    RESULTS / "candidate_age_associations.csv", index=False, encoding="utf-8-sig"
)

flow = [
    ["Biomarker file", len(linked)],
    ["Linked to 2009 survey age/household", int(linked[["age", "hhid"]].notna().all(axis=1).sum())],
    ["Linked age and sex", int(age_sex.sum())],
    ["Adults age >=18 with age and sex", int((adult & age_sex).sum())],
    ["Complete core 5 biomarkers", int(complete5.sum())],
    ["Complete expanded 12 laboratory biomarkers", int(complete12.sum())],
    ["Adults + age/sex + core 5", int((adult & age_sex & complete5).sum())],
    ["Adults + age/sex + core 5 + >=2 valid fuel waves", int((adult & age_sex & complete5 & valid2).sum())],
    ["Adults + age/sex + core 5 + >=2 valid waves + valid 2009 fuel", int((adult & age_sex & complete5 & valid2 & valid2009).sum())],
    ["Primary analysis: adults + age/sex + core 5 + >=2 prior effective waves + known effective 2009 fuel", int((adult & age_sex & complete5 & two_prior_effective & known_effective_2009).sum())],
]
pd.DataFrame(flow, columns=["stage", "n"]).to_csv(RESULTS / "sample_flow.csv", index=False, encoding="utf-8-sig")

eligible = adult & age_sex & complete5 & valid2 & valid2009
eligible_traj = linked.loc[eligible, "strict_trajectory"].value_counts(dropna=False).rename_axis("strict_trajectory").reset_index(name="n")
eligible_traj.to_csv(RESULTS / "trajectory_counts_analysis_eligible.csv", index=False, encoding="utf-8-sig")

main_eligible = adult & age_sex & complete5 & two_prior_effective & known_effective_2009
linked.loc[main_eligible, "effective_trajectory"].value_counts(dropna=False).rename_axis("effective_trajectory").reset_index(name="n").to_csv(
    RESULTS / "effective_trajectory_counts_main_eligible.csv", index=False, encoding="utf-8-sig"
)

linkage = pd.DataFrame(
    [
        ["biomarker unique IDs", len(linked)],
        ["2009 survey link", int(linked["age"].notna().sum())],
        ["master sex link", int(linked["gender"].notna().sum())],
        ["2009 physical exam link", int(linked[["sbp_mean", "dbp_mean", "bmi", "U10"]].notna().any(axis=1).sum())],
        ["at least one valid fuel wave through 2009", int(linked["n_valid_fuel_waves"].fillna(0).ge(1).sum())],
        ["at least two valid fuel waves through 2009", int(valid2.sum())],
        ["valid primary fuel in 2009", int(valid2009.sum())],
        ["known effective fuel state in 2009 (primary + secondary)", int(known_effective_2009.sum())],
        ["at least two valid effective pre-2009 waves", int(two_prior_effective.sum())],
    ],
    columns=["linkage_metric", "n"],
)
linkage.to_csv(RESULTS / "linkage_summary.csv", index=False, encoding="utf-8-sig")

# Candidate dictionary used to lock variable roles before modeling.
dictionary = []
for source_key, variables in {
    "asset": ["hhid", "wave", "L8_1", "L8_2"],
    "biomarker": ["IDind"] + expanded,
    "master_id": ["Idind", "gender", "WEST_DOB_Y"],
    "surveys": ["Idind", "hhid", "wave", "age", "urban", "stratum"],
    "physical_exam": ["IDind", "wave", "SYSTOL1", "SYSTOL2", "SYSTOL3", "DIASTOL1", "DIASTOL2", "DIASTOL3", "height", "weight", "U10"],
}.items():
    for var in variables:
        dictionary.append(
            {
                "domain": "fuel" if var.startswith("L8") else "identifier/demographic" if var.lower() in {"idind", "hhid", "wave", "age", "gender", "west_dob_y", "urban", "stratum"} else "biomarker/physical",
                "file": FILES[source_key].name,
                "variable": var,
                "label": labels[source_key].get(var, "Derived during audit"),
                "wave": "1989-2009" if source_key in {"asset", "surveys"} else "2009" if source_key in {"biomarker", "physical_exam"} else "master",
                "level": "household-wave" if source_key == "asset" else "individual-wave" if source_key in {"surveys", "physical_exam"} else "individual",
                "role": "exposure" if var.startswith("L8") else "outcome candidate" if var in expanded or var in {"SYSTOL1", "SYSTOL2", "SYSTOL3", "DIASTOL1", "DIASTOL2", "DIASTOL3", "height", "weight", "U10"} else "link/covariate",
                "status": "audited",
                "notes": "See protocol and analysis plan",
            }
        )
pd.DataFrame(dictionary).to_csv(PROJECT / "variable_dictionary.csv", index=False, encoding="utf-8-sig")

summary = {
    "n_biomarker": len(linked),
    "n_complete_core5": int(complete5.sum()),
    "n_age_sex": int(age_sex.sum()),
    "n_adult_core5_two_fuel_waves": int((adult & age_sex & complete5 & valid2).sum()),
    "n_adult_core5_two_fuel_waves_valid2009": int((adult & age_sex & complete5 & valid2 & valid2009).sum()),
    "n_main_eligible_stacking_aware": int(main_eligible.sum()),
    "strict_trajectory_counts_all_biomarker": person_traj["strict_trajectory"].value_counts().to_dict(),
}
(RESULTS / "audit_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
logging.info("Audit complete: %s", summary)
print(json.dumps(summary, ensure_ascii=False, indent=2))
