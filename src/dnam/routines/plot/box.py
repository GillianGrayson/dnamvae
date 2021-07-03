import plotly.graph_objects as go


def add_box_trace(fig, y, name):
    showlegend = False if name == "" else True
    fig.add_trace(
        go.Box(
            y=y,
            name=name,
            showlegend=showlegend,
            boxpoints='all',
            jitter=0.75,
            pointpos=0.0,
            marker=dict(
                size=5
            ),
            line=dict(
                width=3
            )
        )
    )
