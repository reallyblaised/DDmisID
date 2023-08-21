"""
Checks hadron to hadron-enriched sample efficiencies against sum of efficiencies of reco as a validity
check of species considered in analysis.
"""

from abc import ABC, abstractmethod
import hist
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scienceplots

from ddmisid.utils import read_config
import pickle
import argparse
from pathlib import Path
from typing import Iterable, Union


plt.style.use(["science", "no-latex"])

# custom color map
mpl.rcParams["axes.prop_cycle"] = mpl.cycler(
    color=[
        "#d6604d",
        "#4393c3",
        "#b2182b",
        "#2166ac",
        "#f4a582",
        "#053061",
    ]
)



class ConsistencyChecker(ABC):
    """
    Abstract base class for checking consistency between different measurements of efficiency.
    """
    @abstractmethod
    def perform_check(self) -> None:
        """
        Abstract method enforcing that all subclasses of ConsistencyChecker must include the perform_check method.
        """
    

class HadronConsistencyChecker(ConsistencyChecker):
    """
    Checks whether the selection of species being considered in the analysis is valid for a given species.
    """
    def __init__(self, h_eff_path_prefix: str, template_path_prefix: str, species: str, binning: "dict[str: list]") -> None:
        self.species = species
        self.binning = binning
        self.bin_indices = None
        self.bin_mapping = self._make_bin_mapping()
        self.h_eff = self._open_hist(f"{h_eff_path_prefix}/{species}/all/perf.pkl")
        self.template_path_prefix = template_path_prefix
        self.viz = self._make_viz()

    def _make_bin_mapping(self) -> "list[dict[int:str]]":
        mapping = []
        for key in self.binning.keys():
            inner = {}
            for i in range(len(self.binning[key])-1):
                inner[i] = f"{self.binning[key][i]}-{self.binning[key][i+1]}"
            mapping.append(inner)
        return mapping

    def _make_viz(self) -> plt.Axes:
        _, ax = plt.subplots(figsize=(5.5,3)) # TODO: figure out how to make figsize dynamic to number of bins
        ax.set_title("LHCb Unofficial", loc="right")
        keys = list(self.binning.keys())
        if len(self.binning) == 2:
            bin_indices = [ (x, y) for x in range(len(self.binning[keys[0]])-1) for y in range(len(self.binning[keys[1]])-1) ]
            bin_labels = [ r"$p_{i=" + f"{x}" + r"}, \eta_{j=" + f"{y}" + r"}$" for x, y in bin_indices ] # TODO: pull vars from config
        elif len(self.binning) == 3:
            bin_indices = [ (x, y, z) for x in range(len(self.binning[keys[0]])-1)
                for y in range(len(self.binning[keys[1]])-1)
                    for z in range(len(self.binning[keys[2]])-1) ]
            bin_labels = [ r"$p_{i=" + f"{x}" + r"}, \eta_{j=" + f"{y}" + r"}, n^{tr}_{k=" + f"{z}" + r"}$"
                for x, y, z in bin_indices ]
        else:
            assert len(self.binning) not in [2, 3], "invalid binning dimensions"
        self.bin_indices = bin_indices
        ax.set_xticks(np.arange(len(bin_indices)))
        ax.set_xticklabels(
            bin_labels, rotation="vertical"
            )
        ax.set_ylabel("efficiencies")
        ax.set_xlabel("kinematic bin indices")
        return ax
    
    def _fill_viz(self, h_eff_data: Iterable[float], h_eff_var: Iterable[float], sum_data: Iterable[float], sum_var: Iterable[float]) -> None:
        x_vals = np.arange(len(self.bin_indices))
        self.viz.legend(
            [self.viz.errorbar(x_vals, h_eff_data, h_eff_var, fmt='.', capsize=3), 
                self.viz.errorbar(x_vals, sum_data, sum_var, fmt='.', capsize=3)],
            ["efficiency histogram", "templates"]
        )

    def _open_hist(self, hist_path: str) -> hist.Hist:
        with open(f"{hist_path}", "rb") as file:
            histogram = pickle.load(file)
        return histogram

    def _add_template_effs(self, template_path: str) -> "tuple[float]":
        template = self._open_hist(template_path)
        eff, err = 0, 0
        for val, var in zip(template.values(), template.variances()):
            eff += val
            err += var # variances add in quadrature
        return (eff, err**0.5)

    def perform_check(self) -> None: # TODO: decide on whether or not to include persist data option
        h_eff_vals, h_eff_errs = [], []
        sum_vals, sum_errs = [], []
        if len(self.binning) == 2:
            for x, inner_h in enumerate(self.h_eff.view()):
                for y, bin in enumerate(inner_h):
                    h_eff_vals.append(bin.value)
                    h_eff_errs.append(bin.variance**0.5) # append error, not variance
                    template_path = f"""{self.template_path_prefix}/{self.bin_mapping[0][x]}/{self.bin_mapping[1][y]}/{self.species}.pkl"""
                    sum_val, sum_err = self._add_template_effs(template_path)
                    sum_vals.append(sum_val)
                    sum_errs.append(sum_err)
        elif len(self.binning) == 3:
            for x, inner_h in enumerate(self.h_eff.view()):
                for y, inner_inner_h in enumerate(inner_h):
                    for z, bin in enumerate(inner_inner_h):
                        h_eff_vals.append(bin.value)
                        h_eff_errs.append(bin.variance**0.5) # append error, not variance
                        template_path = f"""{self.template_path_prefix}/{self.bin_mapping[0][x]}/{self.bin_mapping[1][y]}/{self.bin_mapping[2][z]}/{self.species}.pkl"""
                        sum_val, sum_err = self._add_template_effs(template_path)
                        sum_vals.append(sum_val)
                        sum_errs.append(sum_err)
        else:
            assert len(self.binning) not in [2, 3], "invalid binning dimensions"
        self._fill_viz(h_eff_vals, h_eff_errs, sum_vals, sum_errs)

    def save_viz(self, path_prefix: Union[str, None] = None) -> None:
        if path_prefix is None:
            file_path = f"consistency_checks"
        else:
            file_path = f"{path_prefix}"
        Path(file_path).mkdir(parents=True, exist_ok=True)
        [ self.viz.figure.savefig(f"{path_prefix}/{species}.{ext}") for ext in ("pdf", "png") ]


class KaonConsistencyChecker(HadronConsistencyChecker):
    pass


class PionConsistencyChecker(HadronConsistencyChecker):
    pass


class ProtonConsistencyChecker(HadronConsistencyChecker):
    pass


class ElectronConsistencyChecker(HadronConsistencyChecker):
    pass


class ConsistencyCheckerFactory:
    """
    Factory class for ConsistencyChecker selection.
    """
    def __init__(self):
        """
        Initialize with a dictionary of selector creators
        """
        self._creators = {
            "kaon": KaonConsistencyChecker,
            "pion": PionConsistencyChecker,
            "proton": ProtonConsistencyChecker,
            "electron": ElectronConsistencyChecker,
        }
    
    def create_maker(self, mode: str) -> ConsistencyChecker:
        """
        Create a template maker based on the mode given.
        """
        creator = self._creators.get(mode)
        if not creator:
            raise ValueError(f"Invalid mode: {mode}")
        return creator


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="consistency checker")
    parser.add_argument("h_eff_path_prefix", type=str,
        help="provide path prefix for h to HE sample efficiency histograms")
    parser.add_argument("template_path_prefix", type=str,
        help="provide path prefix for template histograms")
    h_eff_path_prefix_arg = parser.parse_args().h_eff_path_prefix
    template_path_prefix_arg = parser.parse_args().template_path_prefix
    h_eff_path_prefix = r"{}".format(h_eff_path_prefix_arg)
    template_path_prefix = r"{}".format(template_path_prefix_arg)

    # h_eff_path_prefix = r"bin/2018/down/antimu_id"
    # template_path_prefix = r"templates"

    # get binning, species from config file
    pid_config = read_config("config/main.yml", key="pid")
    binning, reco_species = pid_config["binning"], list(pid_config["species"].keys())

    # perform checks for every species
    for species in reco_species:
        checker = ConsistencyCheckerFactory().create_maker(species)(
            h_eff_path_prefix, template_path_prefix, species, binning
        )
        checker.perform_check()
        checker.save_viz(path_prefix=r"checks")