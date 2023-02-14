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
def config_simple_model(
    mrange: tuple[float, float],
) -> Callable:
    """Closure to configure simple (sig+bkg) model for fitting.

    Parameters
    ----------

    mrange: tuple[float, float]
        Range of the fitted variable (typically invariant mass)

    Returns
    -------
    simplemodel: Callable[P, R]
        Two-class modle for unbinned maximum-likelihood fits
    """

    def simple_model(
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
            The abundance of sig and bkg, and datal
        """
        return ns + nb, ns * norm.pdf(x, mu, sg) + nb * truncexpon.pdf(
            x, *mrange, 0, lb
        )

    return simple_model


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
