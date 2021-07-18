import plotly.graph_objects as go


def add_violin_trace(fig, y, name):
    showlegend = False if name == "" else True
    fig.add_trace(
        go.Violin(
            y=y,
            name=name,
            box_visible=True,
            showlegend=showlegend,
        )
    )
