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
    # test_binning = {
    #     "P": [10_000, 15_000, 20_000, 25_000, 30_000, 35_000, 40_000, 100_000],
    #     "ETA": [1.5, 3.0, 5.0],
    #     "nTracks": [0, 250, 1000],
    # }
    # test_cut = {
    #     "antimu_id": "DLLmu<0 & IsMuon==0",
    #     "mu_id": "DLLmu>3 & IsMuon==1",
    # }
    # test_data = pd.DataFrame({
    #     "P": [5000, 11000, 12000, 40000],
    #     "ETA": [1.5, 2.5, 5.5, 3.4],
    #     "nTracks": [0, 400, 600, 400],
    #     "DLLmu": [-1, 0.5, -0.4, 4],
    #     "IsMuon": [0, 0, 1, 1],
    # })

    # test_discretizer = DataDiscretizer(test_binning, test_cut, test_data)
    # test_discretizer.discretize()
    # test_hist = test_discretizer.get_hist()
    # test_proj = test_hist.project("data cuts")
    # print(test_proj)
    # fig, ax = plt.subplots(figsize=(6, 6))

    # test_proj.plot2d(ax=ax, cmap="plasma")
    # plt.show()

    file_path = r"C:\Users\kkura\Desktop\UROP_SU23\99_2021_01_443970_443970746_Bc2D0MuNuXSlim.root"
    test_data = load_ntuple(
        file_path = file_path,
        key = "B2DMuNuX_D02KPi_FakeMuonTuple",
        tree_name="DecayTree",
        library="pd",
        batch_size="200 MB",
        max_entries=None
    )

    # first_five = test_data.head(5)
    # first_five.to_csv("first_five.csv")

    antimu_id = "Mu_plus_CombDLLMu<0 & Mu_plus_isMuon==0"
    common_sel = "Mu_plus_InMuonAcc==1.0 & Mu_plus_NShared==0 & Mu_plus_P>10000 & Mu_plus_P<100000 & Mu_plus_PT>1500"
    test_data = test_data.query(f"({antimu_id}) & ({common_sel})")

    test_binning = {
        "Mu_plus_P": [10_000, 15_000, 20_000, 25_000, 30_000, 35_000, 40_000, 100_000],
        "Mu_plus_LK_ETA": [1.5, 3.0, 5.0],
        "nTracks": [0, 250, 1000],
    }
    test_cuts = {
        "kaon-like": "Mu_plus_PIDmu<0 & Mu_plus_isMuon==0 & Mu_plus_ProbNNghost<0.1 & Mu_plus_PIDK>0.0 & (Mu_plus_PIDK-Mu_plus_PIDp)>0.0 & (Mu_plus_PIDK-Mu_plus_PIDe)>0.0",
        "pion-like": "Mu_plus_PIDmu<0 & Mu_plus_isMuon==0 & Mu_plus_ProbNNghost<0.1 & Mu_plus_PIDK<0.0 & Mu_plus_PIDp<0.0 & Mu_plus_PIDe<0.0",
        "proton-like": "Mu_plus_PIDmu<0 & Mu_plus_isMuon==0 & Mu_plus_ProbNNghost<0.1 & Mu_plus_PIDp>0.0 & (Mu_plus_PIDp-Mu_plus_PIDK)>0.0 & (Mu_plus_PIDp-Mu_plus_PIDe)>0.0",
        "electron-like": "Mu_plus_PIDmu<0 & Mu_plus_isMuon==0 & Mu_plus_ProbNNghost<0.1 & Mu_plus_PIDe>0.0 & (Mu_plus_PIDe-Mu_plus_PIDp)>0.0 & (Mu_plus_PIDe-Mu_plus_PIDK)>0.0",
    }

    test_discretizer = DataDiscretizer(test_binning, test_cuts, test_data)
    test_hist = test_discretizer.discretize()

    axes = ("Mu_plus_P", "Mu_plus_LK_ETA", "nTracks", "data cuts")
    for key in axes:
        print(test_hist.project(key))
        plot_hist_1d(test_hist.project(key), key, "test_hists")
        for key_2 in axes:
            if key != key_2:
                plot_hist_2d(test_hist.project(key, key_2), key, key_2, "test_hists_2d")

    # sanity checks:
    print(f"-----CHECKING TEST CUTS-----")
    for index, (key, val) in enumerate(test_cuts.items()):
        cut_df = test_data.query(val)
        print(f"{key} count (test data): {cut_df.shape[0]}")


    print(f"\n-----CHECKING SELECTED BINS-----")

    print(f"bin test 1: P in (10000, 15000), ETA in (1.5, 3.0), nTracks in (0, 250), pion-like:")
    sel_cut_1 = f"Mu_plus_P>10000 & Mu_plus_P<15000 & Mu_plus_LK_ETA>1.5 & Mu_plus_LK_ETA<3.0 & nTracks>0 & nTracks<250 & {test_cuts['pion-like']}"
    cut_df = test_data.query(sel_cut_1)
    print(f" count (test data): {cut_df.shape[0]}")
    print(f" count (hist bin): {test_hist[(0, 0, 0, 1)]}")

    print(f"bin test 1: P in (15000, 20000), ETA in (1.5, 3.0), nTracks in (0, 250), pion-like:")
    sel_cut_1 = f"Mu_plus_P>15000 & Mu_plus_P<20000 & Mu_plus_LK_ETA>1.5 & Mu_plus_LK_ETA<3.0 & nTracks>0 & nTracks<250 & {test_cuts['pion-like']}"
    cut_df = test_data.query(sel_cut_1)
    print(f" count (test data): {cut_df.shape[0]}")
    print(f" count (hist bin): {test_hist[(1, 0, 0, 1)]}")

    print(f"bin test 1: P in (10000, 15000), ETA in (3.0, 5.0), nTracks in (0, 250), pion-like:")
    sel_cut_1 = f"Mu_plus_P>10000 & Mu_plus_P<15000 & Mu_plus_LK_ETA>3.0 & Mu_plus_LK_ETA<5.0 & nTracks>0 & nTracks<250 & {test_cuts['pion-like']}"
    cut_df = test_data.query(sel_cut_1)
    print(f" count (test data): {cut_df.shape[0]}")
    print(f" count (hist bin): {test_hist[(0, 1, 0, 1)]}")

    print(f"bin test 1: P in (10000, 15000), ETA in (3.0, 5.0), nTracks in (0, 250), kaon-like:")
    sel_cut_1 = f"Mu_plus_P>10000 & Mu_plus_P<15000 & Mu_plus_LK_ETA>1.5 & Mu_plus_LK_ETA<3.0 & nTracks>0 & nTracks<250 & {test_cuts['kaon-like']}"
    cut_df = test_data.query(sel_cut_1)
    print(f" count (test data): {cut_df.shape[0]}")
    print(f" count (hist bin): {test_hist[(0, 0, 0, 0)]}")

    print(f"bin test 1: P in (35000, 40000), ETA in (1.5, 3.0), nTracks in (0, 250), pion-like:")
    sel_cut_1 = f"Mu_plus_P>35000 & Mu_plus_P<40000 & Mu_plus_LK_ETA>1.5 & Mu_plus_LK_ETA<3.0 & nTracks>0 & nTracks<250 & {test_cuts['pion-like']}"
    cut_df = test_data.query(sel_cut_1)
    print(f" count (test data): {cut_df.shape[0]}")
    print(f" count (hist bin): {test_hist[(5, 0, 0, 1)]}")

    print(f"bin test 1: P in (35000, 40000), ETA in (1.5, 3.0), nTracks in (0, 250), proton-like:")
    sel_cut_1 = f"Mu_plus_P>35000 & Mu_plus_P<40000 & Mu_plus_LK_ETA>1.5 & Mu_plus_LK_ETA<3.0 & nTracks>0 & nTracks<250 & {test_cuts['proton-like']}"
    cut_df = test_data.query(sel_cut_1)
    print(f" count (test data): {cut_df.shape[0]}")
    print(f" count (hist bin): {test_hist[(5, 0, 0, 2)]}")

    print(f"bin test 1: P in (30000, 35000), ETA in (1.5, 3.0), nTracks in (0, 250), electron-like:")
    sel_cut_1 = f"Mu_plus_P>30000 & Mu_plus_P<35000 & Mu_plus_LK_ETA>1.5 & Mu_plus_LK_ETA<3.0 & nTracks>0 & nTracks<250 & {test_cuts['electron-like']}"
    cut_df = test_data.query(sel_cut_1)
    print(f" count (test data): {cut_df.shape[0]}")
    print(f" count (hist bin): {test_hist[(4, 0, 0, 3)]}")