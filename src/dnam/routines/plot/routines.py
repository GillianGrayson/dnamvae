import plotly

def get_axis(title, title_size, tick_size):
    axis = dict(
        title=title,
        autorange=True,
        showgrid=True,
        zeroline=False,
        linecolor='black',
        showline=True,
        gridcolor='gray',
        gridwidth=0.1,
        mirror=True,
        ticks='outside',
        titlefont=dict(
            color='black',
            size=title_size
        ),
        showticklabels=True,
        tickangle=0,
        tickfont=dict(
            color='black',
            size=tick_size
        ),
        exponentformat='e',
        showexponent='all'
    )
    return axis


def save_figure(fig, fn, width=800, height=600, scale=2):
    plotly.io.write_image(fig, f"{fn}.png", width=width, height=height, scale=scale)
    plotly.io.write_image(fig, f"{fn}.pdf", width=width, height=height, scale=scale)
