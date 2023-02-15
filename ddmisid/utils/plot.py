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
from .fit import SimpleModelFactory
from typing import Any

plt.style.use("science")

# plot a fitted model
# -------------------
def viz_signal(
    x: ArrayLike,
    y: ArrayLike,
    ax: plt.Axes,
) -> Callable:
    """Partially set the ax and datapoints of the visualiser, leaving freedom for color and label"""

    def inner(
        color: str = "tab:red",
        label: str = "Signal",
    ) -> None:
        """Viz the signal"""
        ax.plot(x, y, label=label, color=color, lw=1.0)

    return inner


def viz_bkg(
    x: ArrayLike,
    ax: plt.Axes,
    y_hi: Any,
) -> Callable:
    """Set the common cosmetics for bkg component(s) in the plot"""

    def inner(
        color: str,
        y_lo: Any = 0,
        label: str = "Background",
    ) -> None:
        """Viz the bkg"""
        ax.fill_between(x, y_lo, y_hi, label=label, color=color, alpha=0.33)

    return inner


class VisualizerFactory:
    """
    A factory class to visualise the component(s) for a fit model.
    """

    def __init__(self, mi: Callable, model_config: Callable) -> None:
        """
        Initialise the ViewFitRes class
        Note: mi and model are private attributes
        """
        self._mi = mi  # minimiser
        self._model_config = (
            model_config  # sources a generic model fed to the minimiser
        )

    @property
    def mi(self) -> object:
        return self._mi

    @property
    def model(self) -> Callable:
        return self._model_config

    def plot(
        self,
        components: str | list[str],
        mrange: tuple[float, float],
        ax: plt.Axes,
        bins: int = 100,
        **kwargs: Any,
    ) -> Any:

        x = np.linspace(*mrange, bins)
        N = (mrange[1] - mrange[0]) / bins

        # set some cosmetics for the plot
        ax.set_prop_cycle(plt.cycler("color", plt.cm.viridis(np.linspace(0, 1, 10))))

        match components:
            case "signal":
                pdf = self._model_config(mrange=mrange, components="signal")(
                    x, *self._mi.values
                )[
                    1
                ]  # original function returns yield and pdf
                breakpoint()
                # return viz_sig(x, pdf, ax=ax)
                ax.plot(x, pdf)
            case "total":
                pdf = self._model_config(mrange=mrange, components="total")(
                    x, *self._mi.values
                )[
                    1
                ]  # original function returns yield and pdf]
                # return viz_sig(
                #     x=x,
                #     y=100.0 * pdf,
                #     ax=ax,
                # )(color="tab:blue", label="Total pdf")
                ax.plot(x, pdf)


# fig, ax factory functions
# -------------------------
def simple_ax(
    title: str = "LHCb Unofficial",
    ylabel: str = "Candidates",
    normalised: bool = False,
) -> tuple[plt.Axes, plt.Axes]:
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

    ax.errorbar(
        cx,
        nh,
        yerr=err,
        xerr=(xe[1] - xe[0]) / 2,
        label=label,
        color=f"{color}",
        fmt=".",
        markersize=3,
    )
