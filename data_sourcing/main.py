from other_nuclear_features.download_files import main
from other_nuclear_features.parse_ame2020 import parse_mass_table
from other_nuclear_features.parse_nubase_features import parse_nubase

from exfor.extract_data import main
from exfor.check_exfor_data_variety import main

# import exfor.check_reactions

# from exfor.explore_entry import main

# main(target_dir="./data")
# parse_nubase(
#     input_path="data/nubase_4.mas20.txt",
#     output_csv_path="data_processed/nubase_features.csv",
# )

# parse_mass_table("data/mass_1.mas20.txt", "data_processed/ame2020_mass_features.csv")
# run_pipeline()
main()
