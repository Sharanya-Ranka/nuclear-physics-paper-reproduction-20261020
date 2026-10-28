import json
import os
import re
from typing import Any, Dict, List, Optional
from matplotlib.pyplot import ion
import numpy as np
import pandas as pd

# breakpoint()

# Usage reference at:
# https://github.com/afedynitch/x4i3/blob/master/doc/x4i/x4i.pdf
from x4i3 import exfor_manager, exfor_entry

OUTPUT_DIR = "exfor_g_n_curated"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Selection criteria specified in Section 3.1
MIN_A = 30
MAX_A = 250
MIN_VALID_POINTS = 20
MAX_ENERGY_MEV = 30.0


db = exfor_manager.X4DBManagerDefault()

# breakpoint()
# # Initialize the EXFOR search interface
# exfor = ExforMaster()


def parse_target_a(target_str: str) -> Optional[int]:
    """Extracts mass number A from EXFOR target string (e.g., '60-NI-60' -> 60, '55-MN-55' -> 55)."""
    match = re.search(r"-(\d+)$", target_str.strip())
    if match:
        return int(match.group(1))
    return None


def extract_data_points(dataset) -> Optional[pd.DataFrame]:
    """Extracts photon energy (MeV) and cross-section (millibarns) from an EXFOR dataset."""
    try:
        # Load measurement points into an x4i3 structured table
        table = dataset.data
        if table is None or table.empty:
            return None

        # Identify photon energy and cross-section columns
        energy_col = None
        cs_col = None
        unc_col = None

        breakpoint()

        for col in table.columns:
            col_upper = str(col).upper()
            if "EN" in col_upper and energy_col is None:
                energy_col = col
            elif ("DATA" in col_upper or "SIG" in col_upper) and cs_col is None:
                cs_col = col
            elif "ERR" in col_upper and unc_col is None:
                unc_col = col

        if not energy_col or not cs_col:
            return None

        energies = table[energy_col].to_numpy(dtype=float)
        cross_sections = table[cs_col].to_numpy(dtype=float)

        if unc_col and unc_col in table.columns:
            uncertainties = table[unc_col].to_numpy(dtype=float)
            uncertainties = np.nan_to_num(uncertainties, nan=0.0)
        else:
            uncertainties = np.zeros_like(cross_sections)

        # Standardize units: EXFOR cross-sections in barns -> millibarns (mb)
        # Note: Check dataset units header if variable; default EXFOR cross section is barns (B)
        units_cs = getattr(dataset, "units", {}).get(cs_col, "B").upper()
        if units_cs in ["B", "BARNS"]:
            cross_sections = cross_sections * 1000.0
            uncertainties = uncertainties * 1000.0
        elif units_cs in ["MICRO-B", "MUB"]:
            cross_sections = cross_sections / 1000.0
            uncertainties = uncertainties / 1000.0

        # Standardize photon energy units to MeV
        units_en = getattr(dataset, "units", {}).get(energy_col, "MEV").upper()
        if units_en in ["KEV"]:
            energies = energies / 1000.0
        elif units_en in ["GEV"]:
            energies = energies * 1000.0

        df = pd.DataFrame(
            {
                "energy_MeV": energies,
                "cross_section_mb": cross_sections,
                "uncertainty_mb": uncertainties,
            }
        )

        # Apply filtering:
        # 1. Zero or negative cross-section values removed
        # 2. Maximum energy range bounded up to 30 MeV
        df = df[(df["cross_section_mb"] > 0) & (df["energy_MeV"] <= MAX_ENERGY_MEV)]
        df = df.dropna().sort_values(by="energy_MeV").reset_index(drop=True)
        return df

    except Exception as exc:
        return None


def make_df(dataset):
    headings = pd.MultiIndex.from_tuples(
        map(lambda x: (x[0], x[1]), zip(dataset.labels, dataset.units))
    )
    df = pd.DataFrame(dataset.data, columns=headings)
    return df


def make_json(df, reaction, dataset_key):
    obj = dict(
        metadata=dict(reaction=str(reaction), dataset_key=dataset_key),
        columns=df.columns.tolist(),
        values=df.values.tolist(),
    )

    return obj


def main():
    print(
        f"Searching EXFOR for (g,n) reactions within mass range {MIN_A} <= A <= {MAX_A}..."
    )

    # Search EXFOR records for incident gamma and emitted neutron reactions: '(G,N)'
    search_query: dict[str, Any] = {
        "reaction": "G,N",
        "quantity": "SIG",  # Cross section
    }

    results = db.query(**search_query)
    collect = []
    not_collect = []
    all_reactions = []
    all_data = []
    breakpoint()

    for ent, subentries in list(results.items()):
        print(f"On entry: {ent}")
        for subent in subentries:
            print(f"On subentry: {subent}")
            res = db.retrieve(ENTRY=ent, SUBENT=subent)

            try:
                entry = res[ent]
                # subentry = entry[subent]
                # reaction = subentry["BIB"]["REACTION"]
                # print(f"Enry: {ent}, Subentry: {subent}")
                # print(f"Reaction: {reaction}")
                # all_reactions.append(reaction)

                datasets = entry.getSimplifiedDataSets()
                for dataset_key, dataset in datasets.items():
                    reaction = dataset.reaction[0]
                    if not re.match(
                        pattern=r"Cross section for .*", string=str(reaction)
                    ) or re.search(
                        pattern=r"electric|magnetic",
                        string=str(reaction),
                        flags=re.IGNORECASE,
                    ):
                        continue

                    df = make_df(dataset)
                    print(f"{str(dataset.reaction[0])}")
                    print(df)
                    json_obj = make_json(df, reaction, dataset_key)
                    all_data.append(json_obj)
                    # breakpoint()

                # if re.match(pattern=r"Cross section for .*", string=reaction):
                #     collect.append(reaction)
                # else:
                #     not_collect.append(reaction)
            except Exception as e:
                pass

        # breakpoint()
        # for subent in subentries:
        #     res = db.retrieve(SUBENT=subent)
        #     if len(res) > 0:
        #         breakpoint()

    with open("data_processed/exfor_data_extracted.json", "w") as fp:
        json.dump(all_data, fp, indent=2)

    # breakpoint()
    # y = exfor_entry.X4Entry(results["40950"])
    # # y.
    # # results = exfor.search(**search_query)
    # print(f"Total matching EXFOR reaction records retrieved: {len(results)}")

    # breakpoint()

    # # Group datasets by target isotope to combine multi-lab EXFOR studies per isotope
    # isotopes_data: Dict[str, List[pd.DataFrame]] = {}

    # for entry in results:
    #     target = entry.target
    #     mass_a = parse_target_a(target)

    #     # Filter: 30 <= A <= 250
    #     if mass_a is None or not (MIN_A <= mass_a <= MAX_A):
    #         continue

    #     df_pts = extract_data_points(entry)
    #     if df_pts is None or df_pts.empty:
    #         continue

    #     if target not in isotopes_data:
    #         isotopes_data[target] = []
    #     isotopes_data[target].append(df_pts)

    # retained_isotopes = 0
    # manifest = []

    # # Filter and export isotopes satisfying dataset threshold (> 20 valid points)
    # for target, dfs in isotopes_data.items():
    #     combined_df = pd.concat(dfs, ignore_index=True)
    #     # Drop duplicates across energy coordinate
    #     combined_df = combined_df.drop_duplicates(
    #         subset=["energy_MeV", "cross_section_mb"]
    #     )
    #     combined_df = combined_df.sort_values(by="energy_MeV").reset_index(drop=True)

    #     if len(combined_df) > MIN_VALID_POINTS:
    #         csv_path = os.path.join(OUTPUT_DIR, f"{target}_g_n.csv")
    #         combined_df.to_csv(csv_path, index=False)
    #         retained_isotopes += 1

    #         manifest.append(
    #             {
    #                 "target": target,
    #                 "mass_A": parse_target_a(target),
    #                 "valid_points": len(combined_df),
    #                 "min_energy": combined_df["energy_MeV"].min(),
    #                 "max_energy": combined_df["energy_MeV"].max(),
    #                 "file_path": csv_path,
    #             }
    #         )

    # manifest_df = pd.DataFrame(manifest).sort_values(by="mass_A").reset_index(drop=True)
    # manifest_df.to_csv(os.path.join(OUTPUT_DIR, "dataset_manifest.csv"), index=False)

    # print("\nExtraction Complete:")
    # print(f"- Total qualified isotopes retained: {retained_isotopes}")
    # print(f"- Processed files saved in: {os.path.abspath(OUTPUT_DIR)}")
    # print(
    #     f"- Manifest summary written to: {os.path.join(OUTPUT_DIR, 'dataset_manifest.csv')}"
    # )


if __name__ == "__main__":
    main()
