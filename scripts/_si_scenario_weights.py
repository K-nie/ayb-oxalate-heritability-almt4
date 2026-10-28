#!/usr/bin/env python
"""
Build Supplementary Table S12: breeder-priority scenario weights for the two
weighted selection tools, side by side.

Weights are transcribed verbatim from the analysis scripts (verified):
  - OCS scenarios: scripts/68_ocs.py   SCENARIOS (baseline 1 for every trait,
    with the priority trait(s) up-weighted).
  - Smith-Hazel scenarios: scripts/70_selection_index.py  SCENARIOS (only the
    listed traits carry weight; unlisted traits are 0).

Output: results/70_selection_index/tables/selection_scenario_weights.csv
(long format: tool, scenario, trait, weight), the deposit file behind Table S12.

Author: Benjamin Narh-Madey
"""
from pathlib import Path
import pandas as pd

ROOT = Path("/Users/black_einstein/Desktop/Other_Projects/Prof Adewale_s Data")
OUT = ROOT / "results" / "70_selection_index" / "tables" / "selection_scenario_weights.csv"

INDEX_TRAITS = ["Tannin", "Phenol", "Flavonoid", "Antioxidant", "Mass_of_Seeds",
                "Seed_Coat_Tannin", "Crude_Protein", "Total_Oxalate",
                "Soluble_Oxalate", "Insoluble_Oxalate"]

# --- OCS (scripts/68_ocs.py): baseline 1, priority trait(s) up-weighted ------
ocs = {
    "balanced": {t: 1.0 for t in INDEX_TRAITS},
    "soluble_oxalate_priority": {**{t: 1.0 for t in INDEX_TRAITS},
                                 "Soluble_Oxalate": 3.0, "Insoluble_Oxalate": 2.0},
    "crude_protein_priority": {**{t: 1.0 for t in INDEX_TRAITS},
                               "Crude_Protein": 3.0},
}
# --- Smith-Hazel (scripts/70_selection_index.py): unlisted traits = 0 --------
sh = {
    "nutritional": {"Crude_Protein": 3, "Antioxidant": 1, "Flavonoid": 1,
                    "Phenol": 1, "Mass_of_Seeds": 1, "Soluble_Oxalate": 2, "Tannin": 1},
    "anti_nutritional": {"Soluble_Oxalate": 3, "Total_Oxalate": 1, "Seed_Coat_Tannin": 1,
                         "Tannin": 1, "Crude_Protein": 1, "Insoluble_Oxalate": 1},
    "balanced": {t: 1 for t in INDEX_TRAITS},
}

rows = []
for tool, scens in [("Optimal contribution selection", ocs),
                    ("Smith-Hazel index", sh)]:
    for sc, w in scens.items():
        for t in INDEX_TRAITS:
            rows.append({"tool": tool, "scenario": sc, "trait": t,
                         "weight": float(w.get(t, 0))})
df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)
print(f"wrote {OUT}  ({len(df)} rows: {df.tool.nunique()} tools x "
      f"{df.scenario.nunique()} scenario-rows x {len(INDEX_TRAITS)} traits)")
# wide preview for the record
print(df.pivot_table(index=["tool", "scenario"], columns="trait", values="weight").to_string())
