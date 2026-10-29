import json


def main():
    filepath = "data_processed/exfor_data_extracted2.json"

    with open(filepath, "r") as fp:
        json_data = json.load(fp)

    all_columns = [str(jd["columns"]) for jd in json_data]
    all_columns_unique = list(set(all_columns))

    stats = dict()
    for ds in json_data:
        col_sig = str(ds["columns"])
        if col_sig not in stats:
            stats[col_sig] = dict(datasets_count=0, datapoints_count=0)

        stats[col_sig]["datasets_count"] += 1
        stats[col_sig]["datapoints_count"] += len(ds["values"])

    for col in all_columns_unique:
        print(f"{col} --- {stats[col]}")

    breakpoint()
    print("hi")
