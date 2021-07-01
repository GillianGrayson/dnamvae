import pandas as pd
import plotly.graph_objects as go
from src.dnam.EWAS.routines.plot import get_axis, save_figure
import os


def plot_regression_scatter(
        df: pd.DataFrame,
        continuous_column: tuple,
        categorical_column: str,
        categorical_values: dict,
        result: pd.DataFrame,
        num_cpgs: int,
        path: str,
):
    d = {}
    for key, value in categorical_values.items():
        d[key] = df.loc[df[categorical_column] == value, :]

    head = result.head(num_cpgs)
    for cpg_id, (cpg, row) in enumerate(head.iterrows()):

        fig = go.Figure()
        for key in d:
            fig.add_trace(
                go.Scatter(
                    x=d[key][continuous_column[0]].values,
                    y=d[key][cpg].values,
                    name=key,
                    mode='markers',
                    marker=dict(
                        size=8,
                        opacity=0.7,
                        line=dict(
                            width=1
                        )
                    )
                )
            )


        if isinstance(row['Gene'], str):
            gene = row['Gene']
        else:
            gene = 'non-genic'

        fig.update_layout(
            template="none",
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="center",
                x=0.5
            ),
            title=dict(
                text=f"{cpg} ({gene})",
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
            showlegend=True,
            xaxis=get_axis(continuous_column[1], 20, 20),
            yaxis=get_axis('Methylation Level', 20, 20),
        )

        save_path = f"{path}/figs"
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        save_figure(fig, f"{save_path}/{cpg_id}_{cpg}")
