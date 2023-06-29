"""
Bespoke classes to discretize pandas dataframes into user-specified binnings.
"""

from abc import ABC, abstractmethod
import pandas as pd
import numpy as np
import hist
from collections import OrderedDict

import matplotlib.pyplot as plt # might not have a need for this import 
from utils import load_ntuple
from pathlib import Path


class Discretizer(ABC):
    """
    Abstract base class for discretizing some dataset into bins desired by the user.
    """
    @abstractmethod
    def discretize(self, data: pd.DataFrame, sel_cut: str) -> None:
        """
        Abstract method enforcing that all subclasses of Discretizer must have the discretize method.
        """
        pass

# TODO: consider whether OrderedDicts would be desirable here
class DataDiscretizer(Discretizer):
    """
    User interface for discretizing some dataset into bins as specified by the user with a config file.
    """
    def __init__(self, binning: dict[str: list] | OrderedDict[str: list], data_cuts: dict[str: str] | OrderedDict[str: str], data: pd.DataFrame) -> None:
        """
        Initializes a DataDiscretizer object equipped with parameters for discretizing data.
        """
        self.binning = binning
        self.data_cuts = data_cuts
        self.data = data
        self.preprocess = pd.DataFrame()
        self._df_preprocessing(self.data_cuts)
        self.hist = None

    def _df_preprocessing(self, data_cuts: dict[str: str] | OrderedDict[str: str]) -> None:
        """
        Internal method for preprocessing of data according to user-specified data cuts.
        Explicit calling by user is never necessary.
        """
        none_apply_expr = ""
        for cut, expr in data_cuts.items():
            self.preprocess[cut] = self.data.eval(expr)
            if none_apply_expr == "":
                none_apply_expr += f"not (({expr})" # beginning of none apply expression
            else:
                none_apply_expr += f"|({expr})" 
        none_apply_expr += ")" # close outer expression
        self.preprocess["none"] = self.data.eval(none_apply_expr)

    def add_cuts(self, data_cuts: dict[str: str] | OrderedDict[str: str]) -> None:
        """
        Add additional cuts to consider in Discretizer.
        """
        self.data_cuts.update(data_cuts)
        self._df_preprocessing(data_cuts)

    def discretize(self, verbose: bool = False) -> hist.Hist:
        """
        Discretizes data into user-specified bins. Returns a Hist histogram.
        """
        hist_axes = []
        # binning axes
        for dim in self.binning:
            hist_axes.append(hist.axis.Variable(
                    edges=self.binning[dim], name=dim, label=dim
            ))
        # data cuts axis
        hist_axes.append(
            hist.axis.StrCategory(self.preprocess.columns, name="data cuts", label="data cuts")
        )
        histogram = hist.Hist(*hist_axes)
        # create array with strings as values
        data_cuts_labels = np.where(self.preprocess, self.preprocess.columns, "").sum(axis=1) # NOTE: assumes cuts are orthogonal
        # fill histogram
        histogram.fill(
            *[self.data[col] for col in self.binning.keys()], data_cuts_labels
        )
        self.hist = histogram
        if verbose:
            print(histogram)
        return histogram
    
    def get_hist(self) -> hist.Hist:
        """
        Returns histogram created after discretizing data using the discretize method.
        """
        assert self.hist is not None, f"must run discretize method on {self}"
        return self.hist

    def __str__(self) -> str:
        return f"DataDiscretizer({self.hist})"


def plot_hist_1d(h: hist.Hist, axis: str, path: str) -> None:
    _, ax = plt.subplots(figsize=(6, 4))
    h.plot1d(ax=ax, ls="--", color="teal", lw=3)

    Path(f"{path}").mkdir(parents=True, exist_ok=True)
    plt.savefig(f"{path}/{axis}.png")


def plot_hist_2d(h: hist.Hist, axis_1: str, axis_2: str, path: str) -> None:
    _, ax = plt.subplots(figsize=(6, 6))
    h.plot2d(ax=ax, cmap="plasma")

    Path(f"{path}").mkdir(parents=True, exist_ok=True)
    plt.savefig(f"{path}/{axis_1}_{axis_2}.png")


if __name__ == "__main__":
    pass