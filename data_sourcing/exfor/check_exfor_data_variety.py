import json


def main():
    filepath = "data_processed/exfor_data_extracted.json"

    with open(filepath, "r") as fp:
        json_data = json.load(fp)

    all_columns = [str(jd["columns"]) for jd in json_data]
    all_columns_unique = list(set(all_columns))
    print("\n".join(all_columns_unique))

    breakpoint()
    print("hi")
