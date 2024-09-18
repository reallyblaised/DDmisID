"""
Adopt a suitably truth-matched ghost MC to produce ghost-specific PID efficiency histograms.
"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

import argparse
import numpy as np
import boost_histogram as bh
import pickle
from ddmisid.utils import read_config, simple_load
import re

PID_CONFIG = read_config("config/main.yml", key="pid")

def convert_pid_alias(pidcalib_alias: str, data_prefix: str) -> str:
    """Convert the PIDCalib2 aliases to the branches present in the ghost MC sample"""
    match pidcalib_alias:
        case "DLLmu": 
            return f"{data_prefix}_PIDmu_corr" # ideally we'd like to add _corr to all PID branches, but let's work with this only for the moment
        case "DLLK": 
            return f"{data_prefix}_PIDK"
        case "DLLp": 
            return f"{data_prefix}_PIDp"
        case "DLLe": 
            return f"{data_prefix}_PIDe"
        case _ if "IsMuon" in pidcalib_alias:
            return f"{data_prefix}_isMuon"
        case _ if "hasMuon" in pidcalib_alias:
            return f"{data_prefix}_hasMuon"
        case _ if "InMuonAcc" in pidcalib_alias:
            return f"{data_prefix}_InMuonAcc"
        case _ if "ProbNNghost" in pidcalib_alias:
            return f"{data_prefix}_ProbNNghost"
        case _ if pidcalib_alias.endswith("_P"):
            return f"{data_prefix}_P"
        case _ if pidcalib_alias.endswith("_PT"):
            return f"{data_prefix}_PT"
        case _ if pidcalib_alias.endswith("_ETA"):
            return f"{data_prefix}_LK_ETA"
        case _ if "nTracks" in pidcalib_alias:
            return f"nTracks"
        case _ if "nSPDHits" in pidcalib_alias:
            return f"nSPDHits"
        case _:
            raise ValueError(f"PID branch alias {pidcalib_alias} not recognised")


def convert_pidcalib_selection(selection: str, data_prefix: str) -> str:
    """Translate a selection string by replacing PID branches with data-specific prefixes"""
    
    # Define a function to handle the replacement
    def replace_branch(match):
        pid_branch = match.group(0)  # Extract matched PID branch
        try:
            # Use the pid_branch_convert function to translate it
            return convert_pid_alias(pid_branch, data_prefix)
        except ValueError:
            # If pid_branch is not recognized, keep it unchanged
            return pid_branch
    
    # Regex to find potential pid_branch names (this assumes valid PID names are alphanumeric/underscores)
    pid_branch_pattern = r'[A-Za-z_]+'

    # Use re.sub to replace all matches of pid_branch with their translated form
    translated_selection = re.sub(pid_branch_pattern, replace_branch, selection)

    return translated_selection


def book_pid_hist(binning_dict: dict = PID_CONFIG["binning"], weighted: bool = False) -> bh.Histogram:
    axes = []
    
    # Loop through the binning dictionary and create axes with metadata
    for axis_name, bins in binning_dict.items():
        # Create a regular axis with the bin edges, setting the label as metadata
        axis = bh.axis.Variable(bins, metadata={"name": axis_name})
        axes.append(axis)
    
    # Create the histogram with the axes with weighted storage to save variance values
    if weighted:
        histogram = bh.Histogram(*axes, storage=bh.storage.Weight())
    else:
        histogram = bh.Histogram(*axes)
    
    return histogram


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Produce ghost-specific PID efficiency histograms"
    )
    # parser.add_argument("-o", "--output", help="Path to efficiency histogram")
    # parser.add_argument(
    #     "-k", "--key", help="Directory key in ghost MC file", default=None
    # )
    # parser.add_argument("-t", "--tree", help="Decay tree name", default="DecayTree")
    # parser.add_argument(
    #     "--cut",
    #     help="Cut to apply to ghost MC (denominator & numerator alike in the efficiency)",
    # )
    # parser.add_argument(
    #     "--pid-cut",
    #     help="PID varibale cuts, of which the efficiency histogram is produced",
    # )
    args = parser.parse_args()

    ghost_input_config = PID_CONFIG["ghost_config"]

    # load the ghost MC with the requisite PID-less cuts
    mc = simple_load(
        path=ghost_input_config["path"],
        cut=ghost_input_config["hadron_enriched_def_sel"], # kinematics and geometry cuts (PID-less) bringing the ghost sample in line with the hadron-enriched data
        max_events=100,
        library="pd",
    )

    # ===============================
    # partition-less pid efficiencies 
    # ===============================
    mc_binning_branches = [convert_pid_alias(alias, ghost_input_config["branch_prefix"]) for alias in PID_CONFIG["binning"].keys()]
    breakpoint()

    DENOM_HIST = book_pid_hist().fill(
        mc[mc_binning_branches]        
    )
    breakpoint()

    # h->mu 
    # -----
    tomu_hist_pass = book_pid_hist().fill(
        mc.query(
            convert_pidcalib_selection(PID_CONFIG["mu_id"]["pid_cut"], ghost_input_config["branch_prefix"])
        )
    )
    # efficiency
    #mu_eff = bh.divide(tomu_hist_pass, DENOM_HIST)

    # h->!mu 
    # ------
    toantimu_hist_pass = book_pid_hist().fill(
        mc.query(
            convert_pidcalib_selection(PID_CONFIG["antimu_id"]["pid_cut"], ghost_input_config["branch_prefix"])
        )
    )
    # efficiency
    #antimu_eff = bh.divide(toantimu_hist_pass, DENOM_HIST)


