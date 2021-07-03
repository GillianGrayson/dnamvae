from plotly import graph_objects as go
from src.dnam.routines.plot.routines import get_axis


def add_layout(fig, x_label, y_label, title):
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
            text=title,
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
        xaxis=get_axis(x_label, 20, 20),
        yaxis=get_axis(y_label, 20, 20),
    )