import pandas as pd
import plotly.graph_objects as go
from src.dnam.EWAS.routines.plot import get_axis, save_figure
import os


def plot_mann_whitney_u_test(
        df_1: pd.DataFrame,
        df_2: pd.DataFrame,
        result: pd.DataFrame,
        path: str,
        num_cpgs: int,
        legend: list
):
    colors = ['rgba(93, 164, 214, 0.5)', 'rgba(255, 65, 54, 0.5)']

    result = result.head(num_cpgs)
    for cpg_id, (cpg, row) in enumerate(result.iterrows()):

        data = [
            go.Box(
                y=df_1[cpg].values,
                name=legend[0],
                boxpoints='all',
                jitter=0.75,
                pointpos=0.0,
                fillcolor=colors[0],
                marker=dict(
                    size=5
                ),
                line=dict(
                    width=3
                )
            ),

            go.Box(
                y=df_2[cpg].values,
                name=legend[1],
                boxpoints='all',
                jitter=0.75,
                pointpos=0.0,
                fillcolor=colors[1],
                marker=dict(
                    size=5
                ),
                line=dict(
                    width=3
                )
            )
        ]

        layout = go.Layout(
            template="none",
            title=dict(
                text=f"{cpg} ({row['Gene']}): {row['pvalue']:0.4e}",
                font=dict(
                    size=25
                )
            ),
            autosize=True,
            margin=go.layout.Margin(
                l=120,
                r=20,
                b=80,
                t=100,
                pad=0
            ),
            showlegend=False,
            xaxis=get_axis('', 20, 20),
            yaxis=get_axis(r'$\Huge\beta$', 15, 20),
        )

        save_path = f"{path}/EWAS/MannWhitneyUTest/figs"
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        fig = go.Figure(data=data, layout=layout)
        save_figure(fig, f"{save_path}/{cpg_id}_{cpg}")
