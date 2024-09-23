"""
Bespoke classes to discretize pandas dataframes into user-specified binnings.
"""

from abc import ABC, abstractmethod
import pandas as pd
import numpy as np
import hist
import matplotlib.pyplot as plt  # might not have a need for this import
from ddmisid.utils import read_config, load_ntuple, extract_sel_dict_branches
from pathlib import Path
import argparse
import pickle
import pandas as pd
import hist
from hist import Hist
from functools import partial
from pathlib import Path
import scienceplots
import matplotlib.pyplot as plt

plt.style.use(["science", "no-latex"])


def binning_key_matcher(data_key: str, year: str) -> str:
    """Match data binning key to the PIDCalib binning scheme"""
    assert year in ("2011", "2012", "2015", "2016", "2017", "2018"), "ERROR: unrecognised `year` identifier"

    if year in ("2011", "2012"):
        return data_key # no Brunel prefix for Run 1
    else:
        match data_key: 
            case "P": 
                return "Brunel_P"
            case "ETA":
                return "Brunel_ETA"
            case "nTracks":
                return "nTracks_Brunel"
    

class Discretizer:
    """
    Base class for discretizing data into user-specified binnings.
    """

    def __init__(
        self,
        binning: dict,
        reco_cuts: dict,
        data: pd.DataFrame,
        reco_label: str = "reco",
    ):
        self._binning = binning  # kinematics and occupancy binning
        self._reco_cuts = (
            reco_cuts  # ideally mutually exclusive reco categories pure in one species
        )
        self._data = data  # hadron-enriched dataset
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

    def assign_reco_id(self):
        """
        Adds a new column 'reco' to the data, categorizing each row based on reco_cuts.
        Each row is evaluated against the conditions, and the first match assigns the label.
        If no conditions match, the label is 'Unclassified'.
        """
        # Initialize the 'reco' column with 'ghost' -> all evading the reco partitions are defined as ghosts
        self.data["reco"] = "ghost"

        # Loop through each particle type and its associated condition
        for particle, condition in self.reco_cuts.items():
            # Create a mask where the condition is true
            mask = self.data.eval(condition)
            # Where the mask is true and 'reco' is still 'Unclassified', assign the particle label
            self.data.loc[mask & (self.data["reco"] == "ghost"), "reco"] = particle

    def book_histogram(self):
        """
        Book a 4D hist with the specified kinematics, occupancy and reco-category bins/partitions
        """
        return Hist(
            hist.axis.Variable(
                list(self.binning.values())[0], name=list(self.binning.keys())[0]
            ),
            hist.axis.Variable(
                list(self.binning.values())[1], name=list(self.binning.keys())[1]
            ),
            hist.axis.Variable(
                list(self.binning.values())[2], name=list(self.binning.keys())[2]
            ),
            hist.axis.StrCategory(
                list(self.reco_cuts.keys()) + ["ghost"],
                name=self.reco_label,
            ),
        )

    def discretize(self):
        """
        Fill the histogram with the data, in bins of occupancy, kinematics and reco-category
        """
        self.assign_reco_id()
        hist = self.book_histogram()
        hist.fill(
            self.data[list(self.binning.keys())[0]],
            self.data[list(self.binning.keys())[1]],
            self.data[list(self.binning.keys())[2]],
            self.data["reco"],
        )
        return hist

    @staticmethod
    def save_histogram(hist, path: str):
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


if __name__ == "__main__":
    # get path to data file
    parser = argparse.ArgumentParser(description="data discretizer")
    parser.add_argument(
        "data_file_path", type=str, help="provide path to data for discretizer"
    )
    data_path_arg = parser.parse_args().data_file_path
    data_path = r"{}".format(data_path_arg)  # create raw string from arg

    # get binning, data cuts, selection cuts from config file
    pid_config = read_config("config/main.yml", key="pid")
    binning_no_prefix = pid_config["binning"]

    # data-specific info
    obs_config = read_config("config/main.yml", key="data")
    data_cuts = obs_config[
        "data_cuts"
    ]  # binning (no prefix); full set of cuts for hadron-enriched partiions (comprises common sel between HE and signal)
    root_key, root_tree_name = (
        obs_config["root_config"]["root_key"],
        obs_config["root_config"]["root_tree_name"],
    )
    data_prefixes = obs_config[
        "data_prefixes"
    ]  # prepended to each binning variable -> get data-compatible binning selections

    # rename binning variables to include prefixes used in the hadron-enriched data file
    binning = {}
    for key, prefix in data_prefixes.items():
        if prefix:
            binning[f"{prefix}_{key}"] = binning_no_prefix[binning_key_matcher(key, *pid_config["years"])]
        else:
            binning[f"{key}"] = binning_no_prefix[binning_key_matcher(key, *pid_config["years"])]

    # load data into awkward array for increased speeds - read in only BOI, as will fill obs hist only
    hadron_enriched_dataset = load_ntuple(
        file_path=data_path,
        key=root_key,
        tree_name=root_tree_name,
        max_entries=None,
        library="pd",
        branches=list(
            set(  # avoid duplication of branches
                extract_sel_dict_branches(data_cuts) + list(binning.keys())
            )
        ),  # extract the inclusive set of branches used in any selection definig the hadron-enriched partitions + binning variables
    )

    # -------------------------------------------
    # Discretise the dataset into reco categories
    # -------------------------------------------
    d = Discretizer(
        binning=binning,
        reco_cuts=data_cuts,  # partition into categorical splits of hadrons, electrongs, ghosts
        data=hadron_enriched_dataset,
        reco_label="reco",
    )
    reco_h = d.discretize()

    # save the histograms, looping through bins, projecting out the reco category in each bin
    ax_i, ax_j, ax_k = list(binning.keys())
    for i in range(len(reco_h.axes[ax_i])):
        for j in range(len(reco_h.axes[ax_j])):
            for k in range(len(reco_h.axes[ax_k])):

                p_bin = f"{int(reco_h.axes[ax_i][i][0])}-{int(reco_h.axes[ax_i][i][1])}"
                eta_bin = f"{reco_h.axes[ax_j][j][0]}-{reco_h.axes[ax_j][j][1]}"
                ntracks_bin = f"{int(reco_h.axes[ax_k][k][0])}-{int(reco_h.axes[ax_k][k][1])}"

                # save the histogram
                d.save_histogram(
                    hist=reco_h[i, j, k, ...],
                    path=f"obs/{p_bin}/{eta_bin}/{ntracks_bin}/obs.pkl",
                )
