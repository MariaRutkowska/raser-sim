"""Girlypop plot style: pinks and purples. - prepared by Claude"""
import matplotlib.pyplot as plt
from cycler import cycler

PINK = "#FF4FA3"
PURPLE = "#9B5DE5"
ROSE = "#F15BB5"
LILAC = "#C77DFF"
PLUM = "#5B1A5E"     # dark, for reference lines and text
BLUSH = "#FFF4FA"    # very light pink background


def use_girlypop():
    """Apply the pink-and-purple style to every plot that follows."""
    plt.rcParams.update({
        "axes.prop_cycle": cycler(color=[PINK, PURPLE, ROSE, LILAC]),
        "figure.facecolor": "white",
        "axes.facecolor": BLUSH,
        "axes.edgecolor": PLUM,
        "axes.labelcolor": PLUM,
        "axes.titlecolor": PLUM,
        "axes.titleweight": "bold",
        "axes.grid": True,
        "grid.color": "#F7C6E0",
        "grid.linewidth": 0.8,
        "xtick.color": PLUM,
        "ytick.color": PLUM,
        "legend.frameon": False,
        "legend.labelcolor": PLUM,
        "lines.linewidth": 2.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })