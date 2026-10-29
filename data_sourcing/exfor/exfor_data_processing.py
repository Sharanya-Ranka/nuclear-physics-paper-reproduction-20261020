import json
import pandas as pd

# In this script, we choose to process tables with the following headers
# [['Energy', 'MeV'], ['Data', 'barns'], ['d(Energy)', 'MeV'], ['d(Data)', 'barns']] --- 36
# [['Energy', 'MeV'], ['Data', 'barns'], ['d(Data)', 'barns']] --- 140
# [['Energy', 'MeV'], ['Data', 'barns']] --- 10

# PENDING
# [['EN-RES', 'MEV'], ['DATA', 'MB'], ['FLAG', 'NO-DIM']] --- 1
# [['E-MIN', 'MEV'], ['E-MAX', 'MEV'], ['EN', 'MEV'], ['DATA', 'MB'], ['DATA-ERR', 'MB']] --- 1
# [['E-MIN', 'MEV'], ['EN-MAX', 'MEV'], ['DATA', 'ARB-UNITS']] --- 1
# [['E-LVL', 'MEV'], ['EN', 'MEV'], ['DATA', 'MB'], ['DATA-ERR', 'MB']] --- 3
# [['EN-MAX', 'MEV'], ['DATA', 'ARB-UNITS'], ['DATA-ERR', 'ARB-UNITS']] --- 2
# [['POLAR', 'NO-DIM'], ['EN-MIN', 'MEV'], ['EN-MAX', 'MEV'], ['DATA', 'MB'], ['DATA-ERR', 'MB']] --- 1
# [['E-MIN', 'MEV'], ['EN-MAX', 'MEV'], ['DATA', 'ARB-UNITS'], ['ERR-S', 'ARB-UNITS']] --- 2
# [['EN-MAX', 'MEV'], ['DATA', 'ARB-UNITS'], ['ERR-S', 'ARB-UNITS']] --- 5
# [['EN-RES', 'MEV'], ['DATA', 'MB']] --- 1
# [['Energy', 'MeV'], ['Data', 'barns'], ['d(Energy)', 'MeV']] --- 1


import re


def get_standardised_data(json_data):
    headers_signature = str(json_data["columns"])
    if headers_signature not in [
        "[['Energy', 'MeV'], ['Data', 'barns'], ['d(Energy)', 'MeV'], ['d(Data)', 'barns']]",
        "[['Energy', 'MeV'], ['Data', 'barns'], ['d(Data)', 'barns']]",
        "[['Energy', 'MeV'], ['Data', 'barns']]",
    ]:
        return None

    # Example
    # (93-NP-237(G,N)93-NP-236,SIG)
    mo = re.search(
        pattern=r"(?P<target_z>\d+)-(?P<target_el>.+)-(?P<target_a>\d+)\(G,N\)(?P<daughter_z>\d+)-(?P<daughter_el>.+)-(?P<daughter_a>\d+)",
        string=json_data["metadata"]["reaction_repr"],
    )
    if not mo:
        return None
        # breakpoint()

    columns = list(map(lambda mic: tuple(mic), json_data["columns"]))

    if mo:
        df = pd.DataFrame(json_data["values"], columns=columns)
        df[["TargetA", "TargetZ", "Target", "DaughterA", "DaughterZ", "Daughter"]] = (
            mo.group("target_a"),
            mo.group("target_z"),
            mo.group("target_el"),
            mo.group("daughter_a"),
            mo.group("daughter_z"),
            mo.group("daughter_el"),
        )
    else:
        df = pd.DataFrame()

    # breakpoint()

    return df


def main():
    filepath = "data_processed/exfor_data_extracted.json"

    with open(filepath, "r") as fp:
        json_data = json.load(fp)

    all_processed_data = []

    for jd in json_data:
        processed_df = get_standardised_data(jd)
        if processed_df is not None:
            all_processed_data.append(processed_df)

    consolidated_df = pd.concat(all_processed_data).reset_index()[
        [
            "TargetA",
            "TargetZ",
            "Target",
            "DaughterA",
            "DaughterZ",
            "Daughter",
            ("Energy", "MeV"),
            ("Data", "barns"),
            ("d(Data)", "barns"),
            ("d(Energy)", "MeV"),
        ]
    ]

    consolidated_df.to_csv("data_processed/exfor_data_processed.csv")

    breakpoint()
    print("hi")
