import plotly


def save_figure(fig, fn, width=800, height=600, scale=2):
    plotly.io.write_image(fig, f"{fn}.png", width=width, height=height, scale=scale)
    plotly.io.write_image(fig, f"{fn}.pdf", width=width, height=height, scale=scale)
