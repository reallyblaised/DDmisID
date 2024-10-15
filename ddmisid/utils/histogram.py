"""Histogram-related operations."""

import pickle
import warnings
import numpy as np
from hist import Hist
import boost_histogram as bh
from abc import ABC, abstractmethod
from typing import Union
from uncertainties import ufloat

# Constants
_pid_eff_tolerance = 0.1
_epsilon = 1e-6


def load_hist(f: str) -> Union[Hist, bh.Histogram]:
    """Load boost_histogram/Hist from a pickle file."""
    with open(f, "rb") as f_in:
        return pickle.load(f_in)


class BaseEffHistProcessor(ABC):
    """Abstract base class for all efficiency histogram processors."""

    @abstractmethod
    def process(self, hist: Union[Hist, bh.Histogram]) -> Union[Hist, bh.Histogram]:
        """Process the histogram."""
        pass


class NullEffProcessor(BaseEffHistProcessor):
    """Processor to handle efficiency entries with null central value by setting them to epsilon."""

    def process(self, hist: Union[Hist, bh.Histogram]) -> Union[Hist, bh.Histogram]:
        view = hist.view()

        # Iterate over multi-dimensional indices in the histogram
        for index in np.ndindex(view.shape):
            cval = view[index].value
            variance = view[index].variance

            # shift central value and variance accordingly to avoid null values
            if cval == 0.0:
                # If variance is zero, set it to a small value (but don't add to an existing non-zero variance)
                if variance == 0.0:
                    variance = _epsilon**2
                view[index] = (_epsilon, variance)

        return hist


class NegativeEffProcessor(BaseEffHistProcessor):
    """Processor to handle negative efficiency values based on certain rules."""

    def process(self, hist: Union[Hist, bh.Histogram]) -> Union[Hist, bh.Histogram]:
        view = hist.view()

        # Iterate over multi-dimensional indices in the histogram
        for index in np.ndindex(view.shape):
            cval = view[index].value
            variance = view[index].variance

            # Case by case handling of negative PID efficiency values
            if cval < 0.0:
                if abs(cval) < variance**0.5:  # <0 and within 1 sigma from 0
                    view[index] = (_epsilon, variance)  # proxy for 0.0; store variance
                elif abs(cval) > variance**0.5 and abs(cval) < _pid_eff_tolerance:
                    # <0 but within tolerance
                    view[index] = (_epsilon, abs(cval))  # proxy for 0.0
                elif abs(cval) > 3 * (variance**0.5) and abs(cval) > _pid_eff_tolerance:
                    # Kill process if below 0 and outside tolerance
                    raise ValueError(
                        f"PID efficiency value below 0.0 detected at index {index} outside 3 sigma of 0.0 and above tolerance: {view[index]}. Aborting."
                    )

        return hist


class EfficiencyHistogramProcessor(BaseEffHistProcessor):
    """Processor that combines both zero-value and negative-value handling."""

    def __init__(self):
        self.null_eff_processor = NullEffProcessor()
        self.negative_eff_processor = NegativeEffProcessor()

    def process(self, hist: Union[Hist, bh.Histogram]) -> Union[Hist, bh.Histogram]:
        """First handle zero values, then handle negative values."""
        hist = self.null_eff_processor.process(hist)
        hist = self.negative_eff_processor.process(hist)
        return hist
