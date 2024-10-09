"""
Collect all species-specific yields across all bins of kinematics and occupancy into a weighted histogram
"""

import itertools
import hist
from ddmisid.engine import config
from ddmisid.utils import load_hist, load_root, extract_sel_dict_branches
from pathlib import Path
import argparse
import pickle
import hist
from hist import Hist
import numpy as np


def pairwise(iterable):
    """Generate pairwise combinations from an iterable."""
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)


# Generate file paths
def generate_path(p_bin, r_bin, n_bin, base_dir="scratch/postfit"):
    return f"{base_dir}/{p_bin[0]}-{p_bin[1]}/{r_bin[0]}-{r_bin[1]}/{n_bin[0]}-{n_bin[1]}/yields.pkl"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Collected fitted yields across kinematics and occupancy"
    )
    parser.add_argument(
        "-i",
        "--input",
        help="Yield .pkl files obtained from BML fits to the hadron-enriched data",
        required=True,
        nargs="+",  # store into list
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output path for the collected yields",
        required=True,
    )
    opts = parser.parse_args()

    # book the container histogram with weighted storage; required structure should match pidcalib, eg:
    # Histogram(
    #     Variable([10000, 100000], metadata=...), -> momentum-like binning var
    #     Variable([1.5, 3, 5], metadata=...), -> eta-like binning var
    #     Variable([0, 250, 1000], metadata=...), -> occupancy-like binning var
    #     storage=Weight()
    # )
    # NOTE: as this matches the config file spec, follow this for consistency
    data = load_hist(
        opts.input[0]
    )  # open first yields file to fetch the keys in the right order

    pidcalib_keys, pidcalib_values = zip(*config.pid.sweight_binning.items())
    hist4D = (
        Hist.new.Variable(pidcalib_values[0], name=pidcalib_keys[0])
        .Variable(pidcalib_values[1], name=pidcalib_keys[1])
        .Variable(pidcalib_values[2], name=pidcalib_keys[2])
        .StrCat([k for k in data.keys()], name="species")
        .Weight()
    )  # NOTE: convention: store variance, ie std**2

    # look over ascending bins
    _p_bins = list(pairwise(pidcalib_values[0]))
    _eta_bins = list(pairwise(pidcalib_values[1]))
    _n_bins = list(pairwise(pidcalib_values[2]))

    # fill one histogram containing the cval and std^2
    for i, _ in enumerate(_p_bins):  # momentum
        for j, _ in enumerate(_eta_bins):  # eta
            for k, _ in enumerate(_n_bins):  # multiplicity
                for species in data.keys():
                    yields = load_hist(
                        generate_path(_p_bins[i], _eta_bins[j], _n_bins[k])
                    )
                    hist4D[i, j, k, species] = [
                        yields[species].n,
                        yields[species].s ** 2,
                    ]  # NOTE: store variance by convention

    # write out the collected yields
    with open(f"{opts.output}", "wb") as f:
        pickle.dump(hist4D, f)
