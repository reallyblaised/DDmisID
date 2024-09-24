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
from pathlib import Path 

PID_CONFIG = read_config("config/main.yml", key="pid")


def compute_efficiency_hist(
    pass_hist: bh.Histogram, 
    denom_hist: bh.Histogram
) -> bh.Histogram:
    """
    Compute the efficiency histogram (pass / denom) with binomial uncertainties.

    Args:
    pass_hist (bh.Histogram): Histogram for the number of events that pass the selection.
    denom_hist (bh.Histogram): Histogram for the total number of events (denominator).
    
    Returns:
    bh.Histogram: Histogram where each bin contains the efficiency and its variance.
    """
    
    # Ensure the histograms are of the same shape
    assert pass_hist.axes.size == denom_hist.axes.size, "Histograms must have the same number of dimensions"
    for i in range(len(pass_hist.axes)):
        assert pass_hist.axes[i].size == denom_hist.axes[i].size, f"Axis {i} dimensions do not match"

    # Create a new histogram to store the efficiency
    efficiency_hist: bh.Histogram = bh.Histogram(*pass_hist.axes, storage=bh.storage.Weight())
    
    # Get the views (arrays of bin values) for pass and denom histograms
    pass_values: bh.accumulators.Weight = pass_hist.view()
    denom_values: bh.accumulators.Weight = denom_hist.view()

    # Iterate through each bin and calculate the efficiency and variance
    for idx in np.ndindex(pass_values.shape):
        pass_count: float = pass_values[idx]
        denom_count: float = denom_values[idx]

        if denom_count > 0:
            # Efficiency = Pass / Denom
            efficiency: float = pass_count / denom_count
            
            # Variance = Pass * (1 - Efficiency) / Denom
            variance: float = (pass_count * (1 - efficiency)) / denom_count
            
            # Store the result in the new histogram (central value and variance)
            efficiency_hist[idx] = (efficiency, variance)
        else:
            # If the denominator is zero, set efficiency to zero and variance to zero
            efficiency_hist[idx] = (1e-6, 1e-6)

    return efficiency_hist
    

def convert_pid_alias(pidcalib_alias: str, data_prefix: str) -> str:
    """Convert the PIDCalib2 aliases to the branches present in the ghost MC sample"""
    match pidcalib_alias:
        case "DLLmu": 
            #return f"{data_prefix}_PIDmu_corr" # ideally we'd like to add _corr to all PID branches, but let's work with this only for the moment
            return f"{data_prefix}_PIDmu" # ideally we'd like to add _corr to all PID branches, but let's work with this only for the moment
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
            #return f"{data_prefix}_PIDmu_corr" # ideally we'd like to add _corr to all PID branches, but let's work with this only for the moment
            return f"{data_prefix}_LOKI_PP_InAccMuon"
        case _ if "ProbNNghost" in pidcalib_alias:
            return f"{data_prefix}_ProbNNghost"
        case _ if pidcalib_alias.endswith("_P"):
            return f"{data_prefix}_P"
        case _ if pidcalib_alias.endswith("_PT"):
            return f"{data_prefix}_PT"
        case _ if pidcalib_alias.endswith("_ETA"):
            #return f"{data_prefix}_LK_ETA"
            return f"{data_prefix}_LOKI_ETA"
        case _ if "nTracks" in pidcalib_alias:
            return f"nTracks"
        case _ if "nSPDHits" in pidcalib_alias:
            return f"nSPDHits"
        case _ if "NShared" in pidcalib_alias:
            return f"{data_prefix}_NShared"
        case _:
            raise ValueError(f"PID branch alias {pidcalib_alias} not recognised")


def convert_pidcalib_selection(selection: str, data_prefix: str = PID_CONFIG['ghost_config']["branch_prefix"]) -> str:
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


def wrap_selection_criteria(selection_string: str) -> str:
    # Split the string by '&' and remove any surrounding spaces for each condition
    conditions = [condition.strip() for condition in selection_string.split("&")]
    
    # wrap each condition in parentheses
    wrapped_conditions = [f"({condition})" for condition in conditions]
    
    # join them back together with '&'
    return " & ".join(wrapped_conditions)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Produce ghost-specific PID efficiency histograms"
    )
    parser.add_argument(
        "-d", "--outdir", type=str, default="bin",
    )
    parser.add_argument(
        "-y", "--year", type=str, help="Year of data taking",
    )
    parser.add_argument(
        "-m", "--magpol", type=str, help="Magnet polarity", choices=["up", "down"],
    )
    opts = parser.parse_args()

    ghost_input_config = PID_CONFIG["ghost_config"]
    MC_BINNING_BRANCHES = [convert_pid_alias(alias, ghost_input_config["branch_prefix"]) for alias in PID_CONFIG["binning"].keys()]

    # ===============================
    # partition-less pid efficiencies 
    # ===============================
    DENOM_HIST = book_pid_hist().fill(
        *simple_load(
            path=ghost_input_config["path"],
            key=ghost_input_config["key"],
            tree=ghost_input_config["tree"],
            cut=ghost_input_config['hadron_enriched_def_sel'], # kinematics and geometry cuts (PID-less) bringing the ghost sample in line with the hadron-enriched data
            library="np",
            branches=MC_BINNING_BRANCHES,
        ).values()        
    )

    # h->mu 
    # -----
    tomu_hist_pass = book_pid_hist().fill(
        *simple_load(
            path=ghost_input_config["path"],
            tree=ghost_input_config["tree"],
            cut=f"{ghost_input_config['hadron_enriched_def_sel']} & ( {wrap_selection_criteria(convert_pidcalib_selection(PID_CONFIG['mu_id']['pid_cut']))} )", 
            library="np",
            branches=MC_BINNING_BRANCHES,
        ).values()       
    )
    # efficiency
    ghost_to_signal_eff = compute_efficiency_hist(
        pass_hist=tomu_hist_pass, 
        denom_hist=DENOM_HIST
    )

    # write to file 
    outpath = Path(f"{opts.outdir}/{opts.year}/{opts.magpol}/mu_id/ghost/ghost_to_muon_like/perf_postprocessed.pkl")
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "wb") as f:
        pickle.dump(ghost_to_signal_eff, f)

    # h->!mu 
    # ------
    toantimu_hist_pass = book_pid_hist().fill(
        *simple_load(
            path=ghost_input_config["path"],
            tree=ghost_input_config["tree"],
            cut=f"{ghost_input_config['hadron_enriched_def_sel']} & ( {wrap_selection_criteria(convert_pidcalib_selection(PID_CONFIG['antimu_id']['pid_cut']))} )", 
            library="np",
            branches=MC_BINNING_BRANCHES,
        ).values()       
    )    
    # efficiency
    ghost_to_hadron_enriched_eff = compute_efficiency_hist(
        pass_hist=toantimu_hist_pass, 
        denom_hist=DENOM_HIST
    )
    # write to file 
    outpath = Path(f"{opts.outdir}/{opts.year}/{opts.magpol}/antimu_id/ghost/all/perf_postprocessed.pkl")
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "wb") as f:
        pickle.dump(ghost_to_hadron_enriched_eff, f)

    # =======================
    # partition effs | h->!mu 
    # =======================
    for target, id_cut in PID_CONFIG['reco_cuts'].items(): 
        # pass
        target_hist_pass = book_pid_hist().fill(
            *simple_load(
                path=ghost_input_config["path"],
                tree=ghost_input_config["tree"],
                # apply the same cuts as the h->!mu efficiency *and* add the partition cut
                cut=f"{ghost_input_config['hadron_enriched_def_sel']} & ( {wrap_selection_criteria(convert_pidcalib_selection(PID_CONFIG['antimu_id']['pid_cut']))} ) & \
                    ( {wrap_selection_criteria(convert_pidcalib_selection(id_cut))} )", 
                library="np",
                branches=MC_BINNING_BRANCHES,
            ).values()       
        )    
        # the conditional probability is expressed by dividing by the h->!mu histogram
        ghost_to_target_eff = compute_efficiency_hist(
            pass_hist=target_hist_pass, 
            denom_hist=toantimu_hist_pass
        )

        # write to file 
        outpath = Path(f"{opts.outdir}/{opts.year}/{opts.magpol}/antimu_id/ghost/ghost_to_{target}/perf_postprocessed.pkl")
        outpath.parent.mkdir(parents=True, exist_ok=True)
        with open(outpath, "wb") as f:
            pickle.dump(ghost_to_target_eff, f)