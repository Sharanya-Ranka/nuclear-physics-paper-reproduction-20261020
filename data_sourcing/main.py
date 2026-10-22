from other_nuclear_features.download_files import main
from other_nuclear_features.parse_ame2020 import parse_mass_table
from other_nuclear_features.parse_nubase_features import parse_nubase
from exfor.extract_data import run_pipeline

# import exfor.check_reactions

# main(target_dir="./data")
# parse_nubase(
#     input_path="data/nubase_4.mas20.txt",
#     output_csv_path="data_processed/nubase_features.csv",
# )

# parse_mass_table("data/mass_1.mas20.txt", "data_processed/ame2020_mass_features.csv")
run_pipeline()
