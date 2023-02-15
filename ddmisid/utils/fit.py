"""Template models for fitting."""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at mit.edu"

from typing import Any
from numpy.typing import ArrayLike
from scipy.stats import expon
from collections.abc import Callable
import warnings
from termcolor2 import c as tc
import numpy as np
from .models import dcbwg, dcbwg_cdf


# commonly used fit models
# ------------------------
def make_expon(
    mrange: ArrayLike, lb: float, x: ArrayLike
) -> tuple[Callable, Callable, Callable]:
    """Generate the exponential pdf and cdf for the combinatorial background

    Parameters
    ----------
    mrange : ArrayLike
        The observation interval considered in the fit

    lb: float
        The lambda parameter of the exponential

    x: ArrayLike
        The datapoints at which the pdf and cdf are evaluated

    Returns
    -------
    expon_pdf: ArrayLike
        The exponential pdf evaluated at x

    expon_cdf: ArrayLike
        The exponential cdf evaluated at mrange

    norm_expon_pdf: ArrayLike
        The normalised exponential pdf evaluated across the observational interval; serves to appropriately combined pdf components in a mixture model.
    """
    _expon = lambda low_m_edge, lb: expon(low_m_edge, lb)

    expon_pdf = _expon.pdf(x)
    expon_cdf = _expon.cdf(mrange)
    norm_expon_pdf = expon_pdf / np.diff(expon_cdf)

    return expon_pdf, expon_cdf, norm_expon_pdf


def expon_factory(key: str) -> Callable:
    """Factory for the exponential pdf, cdf and normalised pdf"""

    match key:
        case "pdf":
            return lambda mrange, lb, x: make_expon(mrange, lb, x)[0]
        case "cdf":
            return lambda mrange, lb, x: make_expon(mrange, lb, x)[1]
        case "norm_pdf":
            return lambda mrange, lb, x: make_expon(mrange, lb, x)[2]
        case other:
            raise ValueError(f"Invalid key: {key}")


def pdf_factory(
    mrange: tuple[float, float],
    key: tuple[str, ...] | str,
) -> Callable:
    """Factory method to select the desired simple model in the fit.

    Parameters
    ----------
    mrange: tuple[float, float]
        Range of the fitted variable (typically invariant mass)

    components: tuple[str, ...] | str
        Names of the components in the model

    Returns
    -------
    Callable
        Model for unbinned maximum-likelihood fits
    """
    match key:
        case None:
            raise TypeError("Please specify at least one component identifier [str]")
        case "signal":  # mixture of two one-sided crystal ball functions and a gaussian
            return lambda x, f1, f2, mug, sgg, sgl, sgr, al, ar, nl, nr: dcbwg(
                x, f1, f2, mug, mug, mug, sgg, sgl, sgr, al, ar, nl, nr, mrange
            )
        case "combinatorial":
            return lambda expon: expon_factory(key="norm_pdf")
        case _:
            raise ValueError("Invalid component identifier(s)")


def cdf_factory(
    mrange: tuple[float, float],
    key: tuple[str, ...] | str,
) -> Callable:
    """Factory method to select the desired model cdf.

    Parameters
    ----------
    mrange: tuple[float, float]
        Range of the fitted variable (typically invariant mass)

    components: tuple[str, ...] | str
        Names of the components in the model

    Returns
    -------
    Callable
        Model cdf evaluated across the mrange
    """
    match key:
        case None:
            raise TypeError("Please specify at least one component identifier [str]")
        case "signal":  # mixture of two one-sided crystal ball functions and a gaussian
            return lambda x, f1, f2, mug, sgg, sgl, sgr, al, ar, nl, nr: dcbwg_cdf(
                x, f1, f2, mug, mug, mug, sgg, sgl, sgr, al, ar, nl, nr, mrange
            )
        case "combinatorial":
            return lambda expon: expon_factory(key="cdf")
        case _:
            raise ValueError("Invalid component identifier(s)")


# verify correctness of the fit
# -----------------------------
class SanityChecks:
    def __init__(self, mi):
        self.mi = mi  # Minuit object

    def __call__(self):
        """Perform the sanity checks with the Minuit object"""

        # fmin is valid
        try:
            assert self.mi.fmin.is_valid
            print(tc("Minuit fmin is valid", "green"))
        except:
            print(tc("Minuit fmin is NOT valid", "red"))

        # covariance matrix is accurate
        try:
            assert self.mi.fmin.has_accurate_covar
            print(tc("Minuit fmin has accurate covariance matrix").green)
        except:
            print(tc("Minuit fmin does not have accurate covariance matrix").red)

        # warn if parameters are at limit
        if self.mi.fmin.has_parameters_at_limit is True:
            print(tc("Warning: Minuit fmin has parameters at limit").yellow)
