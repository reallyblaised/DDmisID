"""
Assign the w_misid weight following the main equation
   
    w_misid = sum_i Ni/Nref eff_i(h->mu)/eff_i(h->!mu)

following SLB WG methodology; see LHCb-ANA-2021-052
"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

import polars as pl
from ddmisid.utils import read_config, simple_load, load_hist
import argparse
import numpy as np
import boost_histogram as bh
from typing import Dict
import warnings
import pprint


def validate_binning(
    histogram: bh.Histogram, binning: Dict[str, np.ndarray | list]
) -> None:
    """Axis validation against the nominal user-defined binning scheme"""
    for axis in histogram.axes:
        name = axis.metadata.get("name")

        # verify the axis in the nominal binning scheme
        if name not in binning:
            raise ValueError(
                f"Value '{name}' from metadata is not a valid key in the binning dictionary"
            )

        # verify edges are as nominal
        if not np.array_equal(axis.edges, binning[name]):
            raise ValueError(f"Binning mismatch in '{name}' projection")


def validate_n_axes(
    histogram: bh.Histogram, binning: Dict[str, np.ndarray | list]
) -> None:
    """Validate axes multiplicity"""
    if len(histogram.axes) != len(binning):
        raise ValueError(
            f"Histogram has {len(histogram.axes)} axes, but binning specifies {len(binning)} - axes multiplicity mismatch"
        )


def process_pid_hists(
    opts: argparse.Namespace,
    binning: Dict[str, np.ndarray | list],
    verbose: bool = False,
) -> Dict[str, bh.Histogram]:
    """Validate PID histogram structure against user-defined nominal binning"""
    pid_container = {}

    for key, value in vars(opts).items():
        if value is None:
            warnings.warn("Warning: None file passed - skipping", UserWarning)
            continue  # skip [accounting for ghost exclusion]

        if "mu" in key:

            pid_hist = load_hist(value)  # Load muon/antimuon PID efficiencies

            # Validate the binning before adding it to the container
            validate_n_axes(pid_hist, binning)
            validate_binning(pid_hist, binning)

            # If test passed, load histograms into container
            pid_container[key] = pid_hist

    # verbose print of dict
    if verbose:
        pprint.pprint(pid_container, width=40, indent=4)

    return pid_container


if __name__ == "__main__":
    # get path to data file
    parser = argparse.ArgumentParser(
        description="tally all effs, yields and observations for the assignment of w_misid"
    )
    parser.add_argument("--obs", help="hadron-enriched data .root file")
    parser.add_argument("--rel_abundances", help="histogram of Ni/Nref")
    parser.add_argument("--p_to_mu", help="path to proton->mu eff", default=None)
    parser.add_argument("--pi_to_mu", help="path to pion->mu eff", default=None)
    parser.add_argument("--k_to_mu", help="path to kaon->mu eff", default=None)
    parser.add_argument("--e_to_mu", help="path to electron->mu eff", default=None)
    parser.add_argument("--g_to_mu", help="path to ghost->mu eff", default=None)
    parser.add_argument(
        "--p_to_antimu",
        help="path to proton->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    parser.add_argument(
        "--pi_to_antimu",
        help="path to pion->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    parser.add_argument(
        "--k_to_antimu",
        help="path to kaon->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    parser.add_argument(
        "--e_to_antimu",
        help="path to electron->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    parser.add_argument(
        "--g_to_antimu",
        help="path to ghost->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    opts = parser.parse_args()

    # user-defined binning
    binning = read_config("config/main.yml", key="pid")["binning"]

    # process input in dedicated containers
    pid_effs = process_pid_hists(opts, binning, False)

    breakpoint()

    # hadron-enriched observations as lazyframe
    data = pl.from_pandas(
        simple_load(
            path=read_config("config/main.yml", key="data")["path"],
            key=read_config("config/main.yml", key="data")["root_config"]["root_key"],
            library="pd",
        )
    ).lazy()

    # Define the custom function
    def custom_function(value, bin_idx=0):
        # Example: use the bin index to modify the value
        # Here, we divide the value by 1000 and multiply by the bin index + 1
        if bin_idx is not None:
            return (value / 1000) * (bin_idx + 1)
        else:
            return None  # Handle out-of-range values

    # Apply binning and pass both the value and bin index to the custom function
    data = data.with_columns(
        [
            pl.col("Mu_plus_P")
            .cut(
                breaks=binning["Brunel_P"],
                labels=[str(i) for i in range(len(binning["Brunel_P"]) + 1)],
            )
            .alias("P_idx"),
            pl.col("Mu_plus_LK_ETA")
            .cut(
                breaks=binning["Brunel_ETA"],
                labels=[str(i) for i in range(len(binning["Brunel_ETA"]) + 1)],
            )
            .alias("ETA_idx"),
            pl.col("nTracks")
            .cut(
                breaks=binning["nTracks_Brunel"],
                labels=[str(i) for i in range(len(binning["nTracks_Brunel"]) + 1)],
            )
            .alias("nTracks_idx"),
        ]
    ).with_columns(
        pl.struct(["P_idx", "ETA_idx", "nTracks_idx"])
        .map_elements(
            lambda x: int(x["P_idx"]) + int(x["ETA_idx"]) + int(x["nTracks_idx"]),
            return_dtype=pl.Float64,
        )
        .alias("dummy_var")
    )
