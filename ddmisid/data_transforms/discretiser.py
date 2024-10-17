"""
Discretiser objects to partition the data into reco partitions, in bins of kinematics, topology and occupancy.
"""

from ddmisid.utils import PIDCalibAliasFactory
import pandas as pd
import awkward as ak
from typing import Union, Tuple, List, Dict
import pickle
import hist
from hist import Hist
from pathlib import Path
import boost_histogram as bh


class Discretiser:
    """
    Base class for discretizing data into user-specified binnings.
    """

    def __init__(
        self,
        binning: Dict[str, List[float]],
        reco_cuts: Dict[str, str],
        data: Union[pd.DataFrame, ak.Array],
        reco_label: str = "reco",
    ):
        self._binning = binning  # kinematics and occupancy binning
        self._reco_cuts = (
            reco_cuts  # ideally mutually exclusive reco categories pure in one species
        )
        self._data = data  # control dataset
        self.reco_label = reco_label
        self._hist = self.discretize()

    @property
    def binning(self):
        return self._binning

    @binning.setter
    def binning(self, binning: dict):
        self._binning = binning

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data: pd.DataFrame):
        self._data = data

    @property
    def reco_cuts(self):
        return self._reco_cuts

    @reco_cuts.setter
    def reco_cuts(self, reco_cuts: dict):
        self._reco_cuts = reco_cuts

    @property
    def hist(self):
        return self._hist

    def discretize(self):
        """
        Fill the histogram with the data, in bins of occupancy, kinematics and reco-category
        """
        self.assign_reco_id()
        hist = self.book_histogram()

        # book binning keys + reco axis
        binning_keys = list(self.binning.keys())
        reco_binning = [self.data[key] for key in binning_keys]
        reco_binning.append(self.data["reco"])

        # book empty histogram with n_{pidcalib axes} + reco axis
        hist.fill(*reco_binning)

        return hist

    def assign_reco_id(self):
        """
        Add a new column 'reco' to the data, categorizing each row based on reco_cuts.
        Each row is evaluated against the conditions, and the first match assigns the label.
        If no conditions match, the label remains 'ghost' (which acts as the default category).

        NOTE: we can safely apply this label regardless of explicit ghost category inclusion in misID modeling.
        """
        # Initialize the 'reco' column with 'ghost' -> all evading the reco partitions are defined as ghosts
        self.data["reco"] = "ghost"

        # Loop through each particle type and its associated condition
        for particle, condition in self.reco_cuts.items():
            # Create a mask where the condition is true
            mask = self.data.eval(condition)
            # Assign the particle label where the mask is true and the current 'reco' is still 'ghost'
            self.data.loc[mask & (self.data["reco"] == "ghost"), "reco"] = particle

    def book_histogram(self):
        """
        Book a histogram with the specified kinematics, occupancy, and reco-category bins/partitions.
        This version supports dynamic binning with a flexible number of variables.
        """
        # Initialize a list to store the axes for the histogram
        axes = []

        # Dynamically create axes for each binning variable
        for bin_name, bin_values in self.binning.items():
            axes.append(hist.axis.Variable(bin_values, name=bin_name))

        # Add the reco-category axis, including the 'ghost' category
        axes.append(
            hist.axis.StrCategory(
                list(self.reco_cuts.keys()) + ["ghost"], name=self.reco_label
            )
        )

        # Book the histogram with the dynamically created axes
        return Hist(*axes)

    @staticmethod
    def save_histogram(hist: Union[Hist, bh.Histogram], path: str):
        """
        Save the discretizer to a pickle file.
        """
        parent_dir = Path(path).parent
        parent_dir.mkdir(parents=True, exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump(hist, f)

    @staticmethod
    def load_hist(path: str):
        """
        Load a discretizer from a pickle file.
        """
        with open(path, "rb") as f:
            return pickle.load(f)
