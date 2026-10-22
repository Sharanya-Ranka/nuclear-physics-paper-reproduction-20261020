import re
import pandas as pd

df = pd.read_csv("data_processed/reactions_view.csv")

INCLUDE_REGEX = r"Cross section for .*|Partial cross section for .*"
IGNORE_CHECKED = r".*(B|b)remsstrahlung.*|\(.*\)/\(.*\)|\(.*\)\+\(.*\)"

included_reactions = df[df["0"].str.match(INCLUDE_REGEX)]
ignored_reactions = df[~df["0"].str.match(INCLUDE_REGEX)]
ignored_reactions_checked = ignored_reactions[
    ignored_reactions["0"].str.match(IGNORE_CHECKED)
]
ignored_reactions_unchecked = ignored_reactions[
    ~ignored_reactions["0"].str.match(IGNORE_CHECKED)
]
breakpoint()
