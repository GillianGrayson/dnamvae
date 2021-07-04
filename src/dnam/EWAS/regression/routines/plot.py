import pandas as pd
import plotly.graph_objects as go
from src.dnam.routines.plot.routines import save_figure
from src.dnam.routines.plot.scatter import add_scatter_trace
from src.dnam.routines.plot.layout import add_layout
import os


def plot_regression_scatter(
        df: pd.DataFrame,
        continuous_column: tuple,
        categorical_column: str,
        categorical_values: list,
        result: pd.DataFrame,
        num_cpgs: int,
        path: str,
):
    d = {}
    for (real, show) in categorical_values:
        d[show] = df.loc[df[categorical_column] == real, :]

    head = result.head(num_cpgs)
    for cpg_id, (cpg, row) in enumerate(head.iterrows()):

        fig = go.Figure()
        for key in d:
            add_scatter_trace(fig, d[key][continuous_column[0]].values, d[key][cpg].values, key)

        if isinstance(row['Gene'], str):
            gene = row['Gene']
        else:
            gene = 'non-genic'

        add_layout(fig, continuous_column[1], 'Methylation Level', f"{cpg} ({gene})")

        save_path = f"{path}/figs"
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        save_figure(fig, f"{save_path}/{cpg_id}_{cpg}")
