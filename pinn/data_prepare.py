import pandas as pd

EXFOR_DATA_FILEPATH = "data_sourcing/data_processed/exfor_data_processed.csv"
AME2020_MASS_DATA_FILEPATH = "data_sourcing/data_processed/ame2020_mass_features.csv"
RIPL_GDR_DATA_FILEPATH = "data_sourcing/data_processed/ripl3_gdr_exp.csv"
OUTPUT_FILEPATH = "pinn/data/prepared_data.csv"


def main():
    exfor_df = pd.read_csv(EXFOR_DATA_FILEPATH)
    ame2020_mass_df = pd.read_csv(AME2020_MASS_DATA_FILEPATH)
    ripl_gdr_df = pd.read_csv(RIPL_GDR_DATA_FILEPATH)

    exfor_df.rename(columns={"TargetA": "A", "TargetZ": "Z"}, inplace=True)

    df_merged = exfor_df.merge(ame2020_mass_df, on=["A", "Z"], how="left").merge(
        ripl_gdr_df, on=["A", "Z"], how="left"
    )
    required_columns = [
        "A",
        "Z",
        "('Energy', 'MeV')",
        "('Data', 'barns')",
        "Element",
        "Sn",
        "E1",
        "CS1",
        "W1",
    ]
    df_merged = df_merged.drop_duplicates(
        subset=["A", "Z", "('Energy', 'MeV')", "('Data', 'barns')"]
    )[required_columns]
    df_merged.rename(columns={"A": "TargetA", "Z": "TargetZ"}, inplace=True)

    breakpoint()
    df_merged.to_csv(OUTPUT_FILEPATH, index=False)
