"""Template models for fitting."""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at mit.edu"

from typing import Any
from numpy.typing import ArrayLike
from numba_stats import norm, truncexpon
from collections.abc import Callable
import warnings


# commonly used fit models
# ------------------------
def simple_twoclass_model(
    mrange: tuple[float, float],
) -> Callable:
    """Closure to specify the fitting range for the model.

    Parameters
    ----------
    mrange: tuple[float, float]
        Range of the fitted variable (typically invariant mass)

    Returns
    -------
    Callable
        Two-class model for unbinned maximum-likelihood fits
    """

    def _simple_twoclass_model(
        x: ArrayLike,
        ns: int,
        nb: int,
        mu: float,
        sg: float,
        lb: float,
    ) -> tuple[int, Any]:
        """Simple model for fitting.
        The total pdf is given the mixture of a gaussian (signal)
        and an exponential (background).

        This is a template model for extended unbinned maximum likelihood fits.

        Parameters
        ----------
        ns: int
            Number of signal candidates
        nb: int
            Number of background candidates
        mu: float
            Mean of the signal gaussian
        sg: float
            Standard deviation of the signal gaussian
        lb: float
            Lambda of the background exponential

        Returns
        -------
        tuple[int, Any]
            The abundance of sig and bkg, and joint pdf
        """
        return ns + nb, ns * norm.pdf(x, mu, sg) + nb * truncexpon.pdf(
            x, *mrange, 0, lb
        )

    return _simple_twoclass_model


def simple_comb_model(
    mrange: tuple[float, float],
) -> Callable:
    """Closure to specify the fitting range for the model.

    Parameters
    ----------
    mrange: tuple[float, float]
        Range of the fitted variable (typically invariant mass)

    Returns
    -------
    Callable
        Model for the combinatorial bkg component in unbinned maximum-likelihood fits
    """

    def _simple_comb_model(
        x: ArrayLike,
        nb: int,
        lb: float,
    ) -> tuple[int, Any]:
        """Simple model for fitting, amounting to a (truncated) exponential pdf.

        This is a template model for extended unbinned maximum likelihood fits.

        Parameters
        ----------
        nb: int
            Number of background candidates
        lb: float
            Lambda of the background exponential

        Returns
        -------
        tuple[int, Any]
            The abundance of bkg, and pdf
        """
        return nb, nb * truncexpon.pdf(x, *mrange, 0, lb)

    return _simple_comb_model


def simple_signal_model() -> Callable:
    """Dummy closure to return the gaussian pdf for the signal; the closure keeps syntax consistent with other models.

    Parameters
    ----------

    Returns
    -------
    Callable
        Model for the signal component in unbinned maximum-likelihood fits
    """

    def _simple_signal_model(
        x: ArrayLike,
        ns: int,
        mu: float,
        sg: float,
    ) -> tuple[int, Any]:
        """Simple model for signal fitting, amounting to a gaussian pdf.

        This is a template model for extended unbinned maximum likelihood fits.

        Parameters
        ----------
        ns: int
            Number of signal candidates
        mu: float
            Mean of the signal gaussian
        sg: float
            Standard deviation of the signal gaussian

        Returns
        -------
        tuple[int, Any]
            The abundance of sig, and datal
        """
        return ns, ns * norm.pdf(x, mu, sg)

    return _simple_signal_model


def SimpleModelFactory(
    mrange: tuple[float, float],
    components: tuple[str, ...] | str = ("total"),
) -> Callable:
    """Factory method to select the desired simple model in the fit.

    Parameters
    ----------
    mrange: tuple[float, float]
        Range of the fitted variable (typically invariant mass)

    components: tuple[str, ...] | str
        Names of the components in the model (default: ('sig', 'bkg'))

    Returns
    -------
    Callable[P, R]
        Simple model for unbinned maximum-likelihood fits
    """
    match components:
        case None:
            raise TypeError("Please specify at least one component identifier [str]")
        case "signal":
            return simple_signal_model()  # no need to specify the mrange here
        case "comb":
            return simple_comb_model(mrange)
        case "total" | ("signal", "comb"):
            return simple_twoclass_model(mrange)
        case _:
            raise ValueError("Invalid component identifier(s)")


# verify correctness of the fit
# -----------------------------
class SanityChecks:
    def __init__(self, mi):
        self.mi = mi  # Minuit object

    def __call__(self, args):
        """Perform the sanity checks with the Minuit object"""

        assert self.mi.fmin.is_valid, "Minuit fmin is not valid"
        assert (
            self.mi.fmin.has_accurate_covar
        ), "Minuit fmin does not have accurate covariance matrix"
        if self.mi.fmin.has_parameters_at_limit is True:
            warnings.warn("Minuit fmin has parameters at limit")
