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
from typing import Any
import mplhep as hep
import matplotlib as mpl
import scienceplots

plt.style.use("science")

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

# plot a fitted model
# -------------------
def viz_signal(
    x: ArrayLike,
    y: ArrayLike,
    ax: plt.Axes,
) -> Callable:
    """Partially set the ax and datapoints of the visualiser, leaving freedom for color and label"""

    def inner(
        color: str = None,
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
    title: str | None = "LHCb Unofficial",
    ylabel: str = "Candidates",
    normalised: bool = False,
    scale: str | None = None,
) -> tuple[Any, plt.Axes]:
    """Book simple ax

    Parameters
    ----------
    title: str | None
        Title of the plot (default: 'LHCb Unofficial')

    ylabel: str
        Y-axis label (default: 'Candidates')

    normalised: bool
        If true, normalise histograms to unity (default: False)

    Returns
    -------
    tuple[Any, Callable]
        Fig, ax plt.Axes objects
    """
    fig, ax = plt.subplots()

    ax.set_title(title, loc="right")
    ax.set_ylabel(ylabel)
    # logo
    ax.text(
        0.0,
        1.07,
        r"\textbf{DD}\textit{mis}ID \texttt{v" + VERSION + "}",
        ha="left",
        va="top",
        transform=ax.transAxes,
        color="grey",
    )
    if scale is not None:
        ax.set_yscale(scale)

    return fig, ax


def make_legend(
    ax: plt.Axes,
    on_plot: bool = True,
    ycoord: float = -0.6,
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
    match len(labels):
        case 2:
            ncols = 2
        case other:
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


# plot data and pdfs
# ------------------
def fill_hist_w(
    data: ArrayLike, bins: int, range: tuple, weights: ArrayLike | None = None
) -> tuple[Any, ...]:
    """Fill histogram accounting weights

    Parameters
    ----------
    data: ArrayLike
        Data to be histogrammed
    bins: int
        Number of bins
    range: tuple
        Range of the histogram
    weights: ArrayLike | None
        Weights to be applied to the data (default: None)

    Returns
    -------
    nh: ArrayLike
        Bin contents
    xe: ArrayLike
        Bin edges
    xc: ArrayLike
        Bin centers
    nh_err: ArrayLike
        Bin errors
    """
    # bin contents, bin edges, bin centers
    nh, xe = np.histogram(data, bins=bins, range=range)
    xc = 0.5 * (xe[1:] + xe[:-1])

    # if no weights, Poisson errors
    nh_err = nh**0.5

    # if weights, careful treatment of bin contents and errors via boost-histogram
    # https://www.zeuthen.desy.de/~wischnew/amanda/discussion/wgterror/working.html
    if weights is not None:
        whist = bh.Histogram(bh.axis.Regular(bins, *range), storage=bh.storage.Weight())
        whist.fill(data, weight=weights)
        cx = whist.axes[0].centers
        nh = whist.view().value
        nh_err = whist.view().variance ** 0.5

    return nh, xe, xc, nh_err


def plot_data(
    data: ArrayLike,
    ax: plt.Axes,
    range: tuple[float, float],
    bins: int = 50,
    weights: ArrayLike | None = None,
    label: str | None = None,
    norm: bool = False,
    color: str = "black",
    **kwargs: Any,
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

    kwargs: Any
        Keyword arguments to be passed to the errorbar plot

    norm: bool
        If true, normalise the data to unit area (default: False)

    Returns
    -------
    None
        Plots the data on the axes
    """
    nh, xe, cx, err = fill_hist_w(data, bins, range, weights)

    # normalise to unity
    _normalisation = np.sum(nh)
    if norm:
        nh = nh / _normalisation
        err = err / _normalisation

    ax.errorbar(
        cx,
        nh,
        yerr=err,
        xerr=(xe[1] - xe[0]) / 2,
        label=label,
        fmt=".",
        markersize=3,
        color=color,
        **kwargs,
    )


def hist_err(
    data: ArrayLike,
    ax: plt.Axes,
    range: tuple[float, float],
    bins: int = 50,
    weights: ArrayLike | None = None,
    label: str | None = None,
    norm: bool = False,
    **kwargs,
) -> None:
    """Wrapper for mplhep histplot method

    Paramaters
    ----------
    nh: ArrayLike
        Histogram bin contents
    xe: ArrayLike
        Histogram bin edges
    ax: plt.Axes
        Axes to plot on
    label: str | None
        Legend label for the data (default: None)
    kwargs: Any
        Keyword arguments to be passed to the errorbar plot
    yerr: ArrayLike | None
        Error for the histogram bin contents (default: None)
    norm: bool
        If true, normalise the data to unit area (default: False)

    Returns
    -------
    None
        Plots the histogram on the axes with errorbars, if not None
    """
    nh, xe, xc, err = fill_hist_w(data, bins, range, weights)

    # normalise to unity
    _normalisation = np.sum(nh)
    if norm is True:
        nh = nh / _normalisation
        err = err / _normalisation

    hep.histplot(nh, xe, yerr=err, ax=ax, label=label, **kwargs)


def hist_step_fill(
    data: ArrayLike,
    range: tuple,
    ax: plt.Axes,
    bins: int = 50,
    weights: ArrayLike | None = None,
    label: str | None = None,
    norm: bool = False,
    **kwargs,
) -> None:
    """Histogram with shaded area

    Paramaters
    ----------
    x: ArrayLike
        Bin centers
    y: ArrayLike
        Bin population
    ax: plt.Axes
        Axes to plot on
    label: str | None
        Legend label for the data (default: None)
    kwargs: Any
        Keyword arguments to be passed to the errorbar plot
    yerr: ArrayLike | None
        Error for the histogram bin contents (default: None)
    norm: bool
        If true, normalise the data to unit area (default: False)

    Returns
    -------
    None
        Plots the histogram on the axes with errorbars, if not None
    """
    nh, xe, xc, err = fill_hist_w(data, bins, range, weights)

    # normalise to unity
    _normalisation = np.sum(nh)
    if norm is True:
        nh = nh / _normalisation
        err = err / _normalisation

    hep.histplot(nh, xe, ax=ax, **kwargs)
    lo_y = nh - err
    ax.bar(
        xc,
        bottom=lo_y,
        height=err * 2,
        alpha=0.33,
        width=xe[1] - xe[0],
        label=label,
        **kwargs,
    )


def eff_plot(
    x: ArrayLike,
    eff: Any,
    eff_err: Any,
    xerr: Any,
    label: str | None,
    ax: plt.Axes,
    fmt: str = ".",
    markersize: int = 5,
    draw_band: bool = True,
    **kwargs,
) -> None:
    """Efficiency plot with shaded error bands

    Parameters
    ----------
    x: ArrayLike
        Bin centers
    eff: Any
        Efficiency values
    eff_err: Any
        Efficiency errors
    label: str | None
        Legend label for the data (default: None)
    ax: plt.Axes
        Axes to plot on
    fmt: str
        Format of the plot (default: '.')
    markersize: int
        Size of the markers (default: 5)
    draw_band: bool
        If true, draw shaded error band (default: True)

    Returns
    -------
    None
        Plots the efficiency on the axes, with error band if draw_band is True
    """
    # quote in percentage
    eff = eff * 100.0
    eff_err = eff_err * 100.0

    if draw_band is True:
        _xerr = 0
        ax.fill_between(
            x, eff - eff_err, eff + eff_err, alpha=0.25, linewidth=0, **kwargs
        )
    else:
        _xerr = xerr

    ax.errorbar(
        x=x,
        y=eff,
        yerr=eff_err,
        xerr=_xerr,
        label=label,
        fmt=fmt,
        markersize=markersize,
    )
