"""Plotting utilities"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

from typing import Any
import matplotlib.pyplot as plt
from typing import Callable
import boost_histogram as bh
import hist
from numpy.typing import ArrayLike
import numpy as np
from pathlib import Path

plt.style.use("science")

# fig, ax factory functions
# -------------------------
def simple_ax(
    title: str = "LHCb Unofficial",
    ylabel: str = "Candidates",
    normalised: bool = False,
) -> tuple[object, object]:
    """Book simple ax

    Parameters
    ----------
    title: str
        Title of the plot (default: 'LHCb Unofficial')

    ylabel: str
        Y-axis label (default: 'Candidates')

    normalised: bool
        If true, normalise histograms to unity (default: False)

    Returns
    -------
    tuple[Callable, Callable]
        Fig, Ax plt.Axes objects
    """
    fig, ax = plt.subplots()

    ax.set_title(title)
    ax.set_ylabel(ylabel)

    return fig, ax


# save plots in multiple formats
def save_to(
    outdir: str,
    name: str,
) -> None:
    """Save the current figure to a path in multiple formats

    Generate directory path if unexeistent

    Parameters
    ----------
    outdir: str
        Directory path to save the figure
    name: str
        Name of the plot

    Returns
    -------
    None
        Saves the figure to the path in pdf and png formats
    """
    Path(outdir).mkdir(parents=True, exist_ok=True)
    [plt.savefig(f"{outdir}/{name}.{ext}") for ext in ["pdf", "png"]]


# plot data and pdfs
# ------------------
def plot_data(
    data: ArrayLike,
    ax: plt.Axes,
    range: tuple[float, float],
    bins: int = 50,
    weights: ArrayLike | None = None,
    label: str | None = None,
    color: str | None = "black",
) -> None:
    """Plot the data, accounting for weights if provided

    Parameters
    ----------
    data: ArrayLike
        Data to be plotted

    ax: plt.Axes
        Axes to plot on

    bins: int
        Number of bins (default: 50)

    range: tuple[float, float]
        Range of the data

    weights: ArrayLike | None
        Weights for the data (default: None)

    label: str | None
        Legend label for the data (default: None)

    color: str | None
        Color for the datapoints (default: black)

    Returns
    -------
    None
        Plots the data on the axes
    """
    nh, xe = np.histogram(data, bins=bins, range=range)
    cx = 0.5 * (xe[1:] + xe[:-1])
    err = nh**0.5
    if weights is not None:
        whist = bh.Histogram(bh.axis.Regular(bins, *range), storage=bh.storage.Weight())
        whist.fill(data, weight=weights)
        cx = whist.axes[0].centers
        nh = whist.view().value
        err = whist.view().variance ** 0.5

    ax.errorbar(cx, nh, err, label=label, color=f"{color}", fmt=".")
