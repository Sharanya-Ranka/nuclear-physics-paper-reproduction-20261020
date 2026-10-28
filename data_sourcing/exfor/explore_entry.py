import os
import re
from typing import Any, Dict, List, Optional
from matplotlib.pyplot import ion
import numpy as np
import pandas as pd

from x4i3 import exfor_manager, exfor_entry

OUTPUT_DIR = "exfor_g_n_curated"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Selection criteria specified in Section 3.1
MIN_A = 30
MAX_A = 250
MIN_VALID_POINTS = 20
MAX_ENERGY_MEV = 30.0


db = exfor_manager.X4DBManagerDefault()


def make_df(dataset):
    headings = pd.MultiIndex.from_tuples(
        map(lambda x: (x[0], x[1]), zip(dataset.labels, dataset.units))
    )
    df = pd.DataFrame(dataset.data, columns=headings)
    return df


def main():
    print(
        f"Searching EXFOR for (g,n) reactions within mass range {MIN_A} <= A <= {MAX_A}..."
    )

    # Search EXFOR records for incident gamma and emitted neutron reactions: '(G,N)'
    search_query: dict[str, Any] = {
        "reaction": "G,N",
        "quantity": "SIG",  # Cross section
    }

    # results = db.query(**search_query)
    collect = []
    not_collect = []
    all_reactions = []

    res = db.retrieve(ENTRY="G0010", SUBENT="G0010002")
    entry = res["G0010"]
    subentry = entry["G0010002"]

    datasets = entry.getSimplifiedDataSets()
    for dataset_key, dataset in datasets.items():
        df = make_df(dataset)
        print(f"{str(dataset.reaction[0])}")
        print(df)

        # breakpoint()

    print("hi")
