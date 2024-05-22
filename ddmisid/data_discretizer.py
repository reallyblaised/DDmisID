"""
Bespoke classes to discretize pandas dataframes into user-specified binnings.
"""

from abc import ABC, abstractmethod
import pandas as pd
import numpy as np
import hist

import matplotlib.pyplot as plt  # might not have a need for this import
from ddmisid.utils import read_config, load_ntuple
from pathlib import Path
import argparse
import pickle


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

    def __init__(
        self,
        binning: "dict[str: list]",
        data_cuts: "dict[str: str]",
        data: pd.DataFrame,
    ) -> None:
        """
        Initializes a DataDiscretizer object equipped with parameters for discretizing data.
        """
        self.binning = binning
        self.data_cuts = data_cuts
        self.data = data
        self.preprocess = pd.DataFrame()
        self._df_preprocessing(self.data_cuts)
        self.hist = None

    def _df_preprocessing(self, data_cuts: "dict[str: str]") -> None:
        """
        Internal method for preprocessing of data according to user-specified data cuts.
        Explicit calling by user is never necessary.
        """
        none_apply_expr = ""
        for cut, expr in data_cuts.items():
            self.preprocess[cut] = self.data.eval(expr)
            if none_apply_expr == "":
                none_apply_expr += (
                    f"not (({expr})"  # beginning of none apply expression
                )
            else:
                none_apply_expr += f"|({expr})"
        none_apply_expr += ")"  # close outer expression
        self.preprocess["none"] = self.data.eval(none_apply_expr)

    def add_cuts(self, data_cuts: "dict[str: str]") -> None:
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
            hist_axes.append(
                hist.axis.Variable(edges=self.binning[dim], name=dim, label=dim)
            )
        # data cuts axis
        hist_axes.append(
            hist.axis.StrCategory(
                self.preprocess.columns, name="data cuts", label="data cuts"
            )
        )
        histogram = hist.Hist(*hist_axes)
        # create array with strings as values
        data_cuts_labels = np.where(self.preprocess, self.preprocess.columns, "").sum(
            axis=1
        )  # NOTE: assumes cuts are orthogonal
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


def pkl_2d(binning: "dict[str: list]", hist: hist.Hist) -> None:
    keys = list(binning.keys())
    num_p_bins, num_eta_bins = len(binning[keys[0]]), len(binning[keys[1]])
    for i, p_bin in enumerate(binning[keys[0]]):
        if i < num_p_bins - 1:
            file_path_p = f"obs/{p_bin}-{binning[keys[0]][i+1]}"
        else:
            continue  # no more bins
        for j, eta_bin in enumerate(binning[keys[1]]):
            if j < num_eta_bins - 1:
                file_path = file_path_p + f"/{eta_bin}-{binning[keys[1]][j+1]}"
            else:
                continue  # no more bins
            obs = hist[i, j, :]  # include all reco categories

            Path(file_path).mkdir(parents=True, exist_ok=True)
            with open(f"{file_path}/obs.pkl", "wb") as file:
                pickle.dump(obs, file)


def pkl_3d(binning: "dict[str: list]", hist: hist.Hist) -> None:
    keys = list(binning.keys())
    num_p_bins, num_eta_bins, num_ntracks_bins = (
        len(binning[keys[0]]),
        len(binning[keys[1]]),
        len(binning[keys[2]]),
    )
    for i, p_bin in enumerate(binning[keys[0]]):
        if i < num_p_bins - 1:
            file_path_p = f"obs/{p_bin}-{binning[keys[0]][i+1]}"
        else:
            continue  # no more bins
        for j, eta_bin in enumerate(binning[keys[1]]):
            if j < num_eta_bins - 1:
                file_path_eta = file_path_p + f"/{eta_bin}-{binning[keys[1]][j+1]}"
            else:
                continue  # no more bins
            for k, ntracks_bin in enumerate(binning[keys[2]]):
                if k < num_ntracks_bins - 1:
                    file_path = (
                        file_path_eta + f"/{ntracks_bin}-{binning[keys[2]][k+1]}"
                    )
                else:
                    continue  # no more bins
                obs = hist[i, j, k, :]  # include all reco categories

                Path(file_path).mkdir(parents=True, exist_ok=True)
                with open(f"{file_path}/obs.pkl", "wb") as file:
                    pickle.dump(obs, file)


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
    antimu_id, common_sel = pid_config["antimu_id"], pid_config["common_sel"]
    binning_no_prefix, data_cuts = pid_config["binning"], pid_config["data_cuts"]
    root_key, root_tree_name = (
        pid_config["root_config"]["root_key"],
        pid_config["root_config"]["root_tree_name"],
    )
    data_prefixes = pid_config["data_prefixes"]

    # rename binning variables to include prefixes used in root file
    binning = {}
    for key, prefix in data_prefixes.items():
        if prefix:
            binning[f"{prefix}_{key}"] = binning_no_prefix[key]
        else:
            binning[f"{key}"] = binning_no_prefix[key]

    # load data, apply cuts, and discretize
    data = load_ntuple(
        file_path=data_path,
        key=root_key,
        tree_name=root_tree_name,
        library="pd",
        batch_size="200 MB",
        max_entries=None,
    )

    # ---------------------------------
    # define the hadron-enriched region 
    # ---------------------------------
    # if common_sel != "": # if common selection tags both the hadron-enriched and signal region, include it
    #     data_sel = data.query(
    #         f"{antimu_id} & {common_sel}"
    #     )  
    # else:
    #     data_sel = data.query(f"{antimu_id}")
    # FIXME
    data_sel = data # HACK

    # discretize into 
    discretizer = DataDiscretizer(binning, data_cuts, data_sel)
    h = discretizer.discretize()

    # save obs for every P, ETA, nTracks to pkl file
    if len(binning) == 2:
        pkl_2d(binning, h)
    elif len(binning) == 3:
        pkl_3d(binning, h)
    else:
        assert len(binning) not in [2, 3], "inappropriate binning dimensions"
