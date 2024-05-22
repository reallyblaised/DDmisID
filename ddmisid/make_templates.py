"""
Functionality to generate templates from PIDCalib jobs
"""

from abc import ABC, abstractmethod
import hist

from ddmisid.utils import read_config
import argparse
import pickle
from pathlib import Path
from typing import Union


class TemplateMaker(ABC):
    """
    Abstract base class for creating templates to be fitted from efficiency histograms.
    """
    @abstractmethod
    def make_hist(self) -> None:
        """
        Abstract method enforcing that all subclasses of TemplateMaker must include the make_hist method.
        """
        pass


class HadronTemplateMaker(TemplateMaker):
    """
    Interface for creating fit templates from PIDCalib efficiency histograms.
    """
    def __init__(self, path_prefix: str, species: str, reco_categories: "list[str]") -> None:
        self.species = species
        self.template_bins = [ f"{species}_to_{reco}_like" for reco in reco_categories ]
        self.path_prefix = path_prefix
        self.eff_hists = self._get_eff_hists(f"{self.path_prefix}/{species}")
        self.hist = None

    
    def _get_eff_hists(self, species_path: str) -> "dict[str: hist.Hist]":
        """Source the PIDCalib2 efficiency histograms for the given species."""
        eff_hists_dict = {}
        for reco in self.template_bins:
            eff_hists_dict[reco] = self._open_eff_hist(f"{species_path}/{reco}/perf.pkl")
        return eff_hists_dict


    def _open_eff_hist(self, eff_hist_path: str) -> hist.Hist:
        """Load the efficiency histogram from the given path."""
        with open(f"{eff_hist_path}", "rb") as file:
            eff_hist = pickle.load(file)
        return eff_hist
    

    def make_hist(self) -> hist.Hist:
        """
        Book the histogram axes
        """
        hist_axes = [ binning_axis for binning_axis in self.eff_hists[self.template_bins[0]].axes ]
        hist_axes.append(
            hist.axis.StrCategory(self.template_bins, name = "reco cuts", label = "reco cuts")
        )
        self.hist = hist.Hist(*hist_axes, storage=hist.storage.Weight())

        for reco, eff_hist in self.eff_hists.items():
            if reco not in self.template_bins:
                raise ValueError(f"{reco} is not a valid category value")
            elif len(hist_axes)-1 == 2:
                self.hist[:, :, reco] = eff_hist.view()
            elif len(hist_axes)-1 == 3:
                self.hist[:, :, :, reco] = eff_hist.view()
            else:
                assert len(hist_axes)-1 == 3, "inappropriate binning dimensions"
        return self.hist

    
    def get_hist(self) -> hist.Hist:
        """Source hist"""
        assert self.hist is not None, f"must run get_hist method on {self}"
        return self.hist
    

    def save_templates(self, binning: "dict[str: list]", path: Union[str, None] = None) -> None:
        """Save the templates to a directory structure based on the binning."""
        assert self.hist is not None, f"must run get_hist method on {self}"
        
        hist_axes = [ binning_axis for binning_axis in self.eff_hists[self.template_bins[0]].axes ]
        if len(hist_axes) == 2:
            pkl_2d(binning, self.hist, self.species, path)
        elif len(hist_axes) == 3:
            pkl_3d(binning, self.hist, self.species, path)
        else:
            assert len(hist_axes) not in [2, 3], "inappropriate binning dimensions"


    def __str__(self) -> str:
        return f"HadronTemplateMaker({self.hist})"


def pkl_2d(binning: "dict[str: list]", hist: hist.Hist, species: str, path: Union[str, None] = None) -> None:
    keys = list(binning.keys())
    num_p_bins, num_eta_bins = len(binning[keys[0]]), len(binning[keys[1]])
    for i, p_bin in enumerate(binning[keys[0]]):
        if i < num_p_bins-1:
            if path is None:
                file_path_p = f"templates/{p_bin}-{binning[keys[0]][i+1]}"
            else:
                file_path_p = f"{path}/{p_bin}-{binning[keys[0]][i+1]}"
        else:
            continue # no more bins
        for j, eta_bin in enumerate(binning[keys[1]]):
            if j < num_eta_bins-1:
                file_path = file_path_p + f"/{eta_bin}-{binning[keys[1]][j+1]}"
            else:
                continue # no more bins
            obs = hist[i, j, :] # include all reco categories

            Path(file_path).mkdir(parents=True, exist_ok=True)
            with open(f"{file_path}/{species}.pkl", "wb") as file:
                pickle.dump(obs, file)


def pkl_3d(binning: "dict[str: list]", hist: hist.Hist, species: str, path: Union[str, None] = None) -> None:
    keys = list(binning.keys())
    num_p_bins, num_eta_bins, num_ntracks_bins = len(binning[keys[0]]), len(binning[keys[1]]), len(binning[keys[2]])
    for i, p_bin in enumerate(binning[keys[0]]):
        if i < num_p_bins-1:
            if path is None:
                file_path_p = f"templates/{p_bin}-{binning[keys[0]][i+1]}"
            else:
                file_path_p = f"{path}/{p_bin}-{binning[keys[0]][i+1]}"
        else:
            continue # no more bins
        for j, eta_bin in enumerate(binning[keys[1]]):
            if j < num_eta_bins-1:
                file_path_eta = file_path_p + f"/{eta_bin}-{binning[keys[1]][j+1]}"
            else:
                continue # no more bins
            for k, ntracks_bin in enumerate(binning[keys[2]]):
                if k < num_ntracks_bins-1:
                    file_path = file_path_eta + f"/{ntracks_bin}-{binning[keys[2]][k+1]}"
                else:
                    continue # no more bins
                obs = hist[i, j, k, :] # include all reco categories

                Path(file_path).mkdir(parents=True, exist_ok=True)
                with open(f"{file_path}/{species}.pkl", "wb") as file:
                    pickle.dump(obs, file)


class KaonTemplateMaker(HadronTemplateMaker):
    pass


class PionTemplateMaker(HadronTemplateMaker):
    pass


class ProtonTemplateMaker(HadronTemplateMaker):
    pass


class ElectronTemplateMaker(HadronTemplateMaker):
    pass


class TemplateMakerFactory:
    """
    Factory class for TemplateMaker selection.
    """
    def __init__(self):
        """
        Initialize with a dictionary of selector creators
        """
        self._creators = {
            "kaon": KaonTemplateMaker,
            "pion": PionTemplateMaker,
            "proton": ProtonTemplateMaker,
            "electron": ElectronTemplateMaker,
        }
    

    def create_maker(self, mode: str) -> TemplateMaker:
        """
        Create a template maker based on the mode given.
        """
        creator = self._creators.get(mode)
        if not creator:
            raise ValueError(f"Invalid mode: {mode}")
        return creator


if __name__ == "__main__":
    # get path prefix for efficiency histograms
    parser = argparse.ArgumentParser(description="template maker")
    parser.add_argument("eff_hist_path_prefix", type=str, 
        help="provide path prefix to efficiency histograms for template maker")
    path_prefix_arg = parser.parse_args().eff_hist_path_prefix
    path_prefix = r"{}".format(path_prefix_arg)

    # get binning, species from config file
    pid_config = read_config("config/main.yml", key="pid")
    binning, reco_categories = pid_config["binning"], list(pid_config["species"].keys())

    # make templates
    for species in reco_categories:
        templates = TemplateMakerFactory().create_maker(species)(
            path_prefix, species, reco_categories
        ) # each reco species in the hadron-enriched sample -> each reco category PID eff, eg K->{K, π, p, e} @ !mu
        templates.make_hist()
        templates.save_templates(binning)
        
