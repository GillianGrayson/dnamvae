import pandas as pd
import plotly.graph_objects as go
from src.dnam.routines.plot.save import save_figure
from src.dnam.routines.plot.layout import add_layout, get_axis
from src.dnam.routines.plot.box import add_box_trace
import os


def plot_mann_whitney_u_test(
        df_1: pd.DataFrame,
        df_2: pd.DataFrame,
        result: pd.DataFrame,
        path: str,
        num_cpgs: int,
        legend: list
):
    result = result.head(num_cpgs)
    for cpg_id, (cpg, row) in enumerate(result.iterrows()):

        fig = go.Figure()
        add_box_trace(fig, df_1[cpg].values, legend[0])
        add_box_trace(fig, df_2[cpg].values, legend[1])

        if isinstance(row['Gene'], str):
            gene = row['Gene']
        else:
            gene = 'non-genic'

        add_layout(fig, '', "Methylation Level", f"{cpg} ({gene}): {row['pvalue']:0.4e}")

        save_path = f"{path}/EWAS/mann_whitney_u_test/figs"
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        save_figure(fig, f"{save_path}/{cpg_id}_{cpg}")
