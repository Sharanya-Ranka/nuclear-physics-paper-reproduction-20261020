import pandas as pd
from plotly import graph_objects as go

data_filepath = "data_processed/exfor_data_processed.csv"
A = 60
Z = 28


def main():
    df = pd.read_csv(data_filepath)
    df_filtered = df[(df["TargetA"] == A) & (df["TargetZ"] == Z)]
    print(df_filtered)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=df_filtered["('Energy', 'MeV')"],
            y=df_filtered["('Data', 'barns')"],
            mode="markers",
            name=f"A={A}, Z={Z}",
        )
    )
    fig.update_layout(
        title=f"Exfor Data for A={A}, Z={Z}",
        xaxis_title="Energy",
        yaxis_title="Cross Section",
    )
    fig.write_html(f"plots/exfor_data_A{A}_Z{Z}.html")
