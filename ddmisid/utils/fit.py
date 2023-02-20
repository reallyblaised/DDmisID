"""Template models for fitting."""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at mit.edu"

from typing import Any
from numpy.typing import ArrayLike
from collections.abc import Callable
import warnings
from termcolor2 import c as tc
import numpy as np
from .models import dcbwg, dcbwg_cdf, expon_pdf, expon_cdf
from functools import partial

# commonly used fit models
# ------------------------
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
            return lambda x, lb: expon_pdf(mrange=mrange, lb=lb, x=x)
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
            return lambda x, lb: expon_cdf(x, lb, mrange)
        case _:
            raise ValueError("Invalid component identifier(s)")


def twoclass_pdf(
    comps: tuple[str, ...] | list | None,
    x: ArrayLike | list | tuple,
    f1: float,
    f2: float,
    mug: float,
    mul: float,
    mur: float,
    sgg: float,
    sgl: float,
    sgr: float,
    al: float,
    ar: float,
    nl: float,
    nr: float,
    lb: float,
    mrange: tuple[float, float] | ArrayLike,
    sig_yield: int,
    comb_yield: int,
) -> Any:
    """Generate a mixture of pdfs, suitably normalised, to fit a realistic mass spectrum:
    - a signal modelled by DCB + gaussian
    - combinatorial bkg, modelled by a exp
    """
    if comps is None:
        comps = ["signal", "combinatorial"]

    model = 0
    # successively generate the linear composition of pdfs
    if "signal" in comps:
        model += sig_yield * dcbwg(
            x, f1, f2, mug, mul, mur, sgg, sgl, sgr, al, ar, nl, nr, mrange
        )

    if "combinatorial" in comps:
        model += comb_yield * expon_pdf(x, lb, mrange)

    return model


def twoclass_cdf(
    comps: tuple[str, ...] | list | None,
    x: ArrayLike | list | tuple,
    f1: float,
    f2: float,
    mug: float,  # shared
    sgg: float,
    sgl: float,
    sgr: float,
    al: float,
    ar: float,
    nl: float,
    nr: float,
    lb: float,
    mrange: tuple[float, float] | ArrayLike,
    sig_yield: int,
    comb_yield: int,
) -> Any:
    """Generate a mixture of pdfs, suitably normalised, to fit a realistic mass spectrum:
    - a signal modelled by DCB + gaussian
    - combinatorial bkg, modelled by a exp
    """
    if comps is None:
        comps = ["signal", "combinatorial"]

    model = 0
    # successively generate the linear composition of pdfs
    if "signal" in comps:
        model += sig_yield * dcbwg_cdf(
            x, f1, f2, mug, mug, mug, sgg, sgl, sgr, al, ar, nl, nr, mrange
        )

    if "combinatorial" in comps:
        model += comb_yield * expon_cdf(x, lb, mrange)

    return model


def composite_pdf_factory(
    key: str,
    mrange: tuple[float, float] | ArrayLike,
) -> Any:
    """factory method to select the desired model"""
    match key:
        case "twoclass":
            return partial(
                twoclass_pdf, mrange=mrange, comps=["signal", "combinatorial"]
            )


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
            print(self.mi.params)
