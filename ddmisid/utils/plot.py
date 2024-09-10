"""Plotting utilities"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

from typing import Any, Optional, Union, List, Tuple
import matplotlib.pyplot as plt
from typing import Callable
import boost_histogram as bh
import hist
from numpy.typing import ArrayLike
import numpy as np
from pathlib import Path
from typing import Any
import mplhep as hep
import matplotlib as mpl
import scienceplots

plt.style.use(["science", "no-latex"])

# custom color map
mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
    color=[
        "#d6604d",
        "#4393c3",
        "#b2182b",
        "#2166ac",
        "#f4a582",
        "#053061",
    ]
)


# software version
VERSION = "0.1"


# book the axes for postfit plot 
# ------------------------------
def data_pull_plot(
    title: str | None = f"DDmisID v{VERSION}",
    ylabel: str = "Candidates",
    scale: str | None = None,
    axp_xlabel: str | None = None,
    axp_ylabel: str = r"Pulls $[\sigma]$",
    annotation: str | None = None,
    is_pull: bool = True,
) -> tuple[Any, plt.Axes, plt.Axes]:
    """
    Generate 2-ax fig with upper ax for data, lower ax for pulls

    Parameters
    ----------
    title: str | None
        Title of the plot (default: 'LHCb Prelminary')

    ylabel: str
        Y-axis label (default: 'Candidates')

    normalised: bool
        If true, normalise histograms to unity (default: False)

    scale: str | None
        If not None, scale the y-axis by the given scale (default: None)

    axp_xlabel: str
        x-axis label for the pull plot (default: r'$m_{\textrm{corr}}(D^0\mu^+)$ [MeV/$c^2$]')

    axp_ylabel: str
        y-axis label for the pull plot (default: r'Pulls $[\sigma]$')

    annotataion: str | None
        If not None, annotate the plot with the given string (default: None)

    is_pull: bool 
        plot fit pulls; if False, data/model (default: True)

    Returns
    -------
    tuple[Any, Callable, Callable]
        fig, ax plt.Axes objects, ax_p plt.Axes
    """
    fig, (ax, ax_p) = plt.subplots(
        2,
        1,
        gridspec_kw={"height_ratios": [5, 1.2]},
        sharex=True,
    )

    # reduce the white space between the two subplots
    fig.subplots_adjust(hspace=0.05)  # Adjust this value as needed

    # main plt config
    ax.set_title(title, loc="right", color="tab:grey")
    ax.set_ylabel(ylabel)

    # additional annotation
    if annotation is not None:
        ax.text(
            0.03,
            0.93,
            annotation,
            ha="left",
            va="top",
            transform=ax.transAxes,
            color="tab:grey",
        )

    if scale is not None:
        ax.set_yscale(scale)

    # pull plt config
    ax_p.set_xlabel(axp_xlabel)
    ax_p.set_ylabel(axp_ylabel)

    if is_pull == True:
        ax_p.set_ylim(-5.5, 5.5)

        ax_p.set_yticks([-3, 0, 3])  # Set y-axis ticks at -3, 0, 3
        ax_p.tick_params(
            axis="y", labelsize="small"
        )  # Set y-axis tick label font size to small

        ax_p.axhline(0, color="tab:grey", lw=0.5, ls="-")
        ax_p.axhline(3, color="firebrick", lw=0.5, ls="--")
        ax_p.axhline(-3, color="firebrick", lw=0.5, ls="--")

    return fig, ax, ax_p


def make_legend(
    ax: plt.Axes,
    on_plot: bool = True,
    ycoord: float = -0.6,
    ncols: None | int = None,
) -> None:
    """
    Place the legend below the plot, adjusting number of columns

    Parameters
    ----------
    ax: plt.Axes
        Axes object to place the legend on

    on_plot: bool
        If true, place the legend on the plot (default: True)

    ycoord: float
        Y-coordinate of the legend (default: -0.6)

    Returns
    -------
    None
        Places the legend on the axes
    """
    # count entries
    handles, labels = ax.get_legend_handles_labels()

    # decide the number of columns accordingly
    if ncols is None:
        match len(labels):
            case 2:
                ncols = 2
            case _ if math.fmod(float(len(labels)), 3.0) == 0:
                ncols = 3
            case _:
                ncols = 1

    # place the legend
    ax.legend(loc="best")
    if on_plot is False:
        ax.legend(
            bbox_to_anchor=(0.5, ycoord),
            loc="lower center",
            ncol=ncols,
            frameon=False,
        )


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