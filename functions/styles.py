from matplotlib import rc
from matplotlib import rcParams
import warnings

colors = ["#173F5F", "#20639B", "#3CAEA3", "#F6D55C", "#BF6B63", "#D9A384"]

def lighten_color(color, amount):
    """
       Lightens the given color by multiplying (1-luminosity) by the given amount.
       Input can be matplotlib color string, hex string, or RGB tuple.

       Examples:
       >> lighten_color('g', 0.3)
       >> lighten_color('#F034A3', 0.6)
       >> lighten_color((.3,.55,.1), 0.5)
       """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


class Colors:
    palette = ["#173F5F", "#19647E", "#3CAEA3", "#95C623" , "#E55812", "#D81E5B"]
    green =  ["#2b6636", "#409951", "#56cc6c", "#6bff87"]
    red = ["#660f0f", "#991717", "#cc1f1f", "#ff2626"]
    pink = ["#663055", "#99487f", "#cc60a9", "#ff78d4"]
    blue = ["#0f4366", "#176599", "#1f87cc", "#26a8ff"]
    cyan = ["#236660", "#349990", "#45ccc0", "#57fff0"]
    orange = ["#ffc100", "#ff9a00", "#ff7400", "#ff4d00"]


def figsize(shape, size):
    """
    :param shape: List with number of columns and rows. I.e. [3, 2] would create a plot with 3 columns and 2 rows
    :param size: String that can be "small" or "large"
    :return: width, height the resulting plot size
    """
    if size == "large":
        height = 200
        width = 300
    elif size == "small":
        height = 150
        width = 200
    else:
        return warnings.warn("No valid plot size specified. Options are \"small\" and \"large\". ")
    x = shape[0]
    y = shape[1]
    height = y*height
    width = x*width
    return width, height


def set_default_plot_style():
    rcParams['text.latex.preamble'] = r'\usepackage{lmodern}'
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = 'Latin Modern'
    rcParams.update({'font.size': 11})
    rc('text', usetex=True)


def remove_top_right_axis(axis):
    for ax in axis:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)


def remove_all_axis(axis):
    for ax in axis:
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)