import plotly.graph_objects as go


def add_histogram_trace(fig, x, name):
    showlegend = False if name == "" else True
    fig.add_trace(
        go.Histogram(
            x=x,
            name=name,
            showlegend=showlegend,
            marker=dict(
                opacity=0.7,
                line=dict(
                    width=1
                ),
            ),
        )
    )
