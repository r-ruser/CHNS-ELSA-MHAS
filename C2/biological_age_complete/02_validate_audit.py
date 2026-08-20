"""Validate aggregate audit outputs and write a reproducibility receipt."""
import json
from pathlib import Path
import pandas as pd

project = Path(__file__).resolve().parents[1]
r = project / "results"
summary = json.loads((r / "audit_summary.json").read_text(encoding="utf-8"))
keys = pd.read_csv(r / "key_uniqueness.csv")
flow = pd.read_csv(r / "sample_flow.csv")
traj = pd.read_csv(r / "effective_trajectory_counts_main_eligible.csv")

assert summary["n_biomarker"] == 9549
assert summary["n_complete_core5"] == 9456
assert int(keys.loc[(keys.dataset == "asset") & (keys.scope == "through 2009"), "duplicate_key_rows"].iloc[0]) == 0
assert int(traj.n.sum()) == summary["n_main_eligible_stacking_aware"]
assert summary["n_main_eligible_stacking_aware"] == int(flow.iloc[-1].n)
assert int(traj.loc[traj.effective_trajectory == "solid_to_clean_only_monotonic", "n"].iloc[0]) >= 1000

receipt = "PASS\nAll aggregate linkage, uniqueness, sample-flow and trajectory invariants passed.\n"
(project / "logs" / "02_validate_audit.txt").write_text(receipt, encoding="utf-8")
print(receipt, end="")

