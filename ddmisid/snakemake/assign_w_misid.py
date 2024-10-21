"""
Assign the w_misid weight following the main equation
   
    w_misid = sum_i Ni/Nref eff_i(h->mu)/eff_i(h->!mu)

following SLB WG methodology; see LHCb-ANA-2021-052
"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

import polars as pl
from ddmisid.utils import load_root, load_hist, write_key_to_df
from ddmisid.engine import config
import argparse
import numpy as np
import warnings
import pprint
from uncertainties import ufloat
from functools import partial
from itertools import product
import boost_histogram as bh
import hist
from typing import Union, List, Dict
from loguru import logger


SPECIES = [
    recocat.replace("_like", "") for recocat in config.pid.reco_partitions.keys()
]


def extract_axis_name(axis: bh.axis) -> str:
    """Extract the name of histogram axis"""
    try:
        return axis.metadata.get("name")
    except:
        return axis.name


def validate_binning(
    histogram: bh.Histogram, binning: Dict[str, Union[List[float], np.ndarray]]
) -> None:
    """Axis validation against the nominal user-defined binning scheme"""
    for axis in histogram.axes:
        name = extract_axis_name(axis)

        # verify the axis in the nominal binning scheme
        if name not in binning:
            raise ValueError(
                f"Value '{name}' from metadata is not a valid key in the binning dictionary"
            )

        # verify edges are as nominal
        if not np.array_equal(axis.edges, binning[name]):
            raise ValueError(f"Binning mismatch in '{name}' projection")


def validate_n_axes(
    histogram: bh.Histogram, binning: Dict[str, Union[List[float], np.ndarray]]
) -> None:
    """Validate axes multiplicity"""
    if len(histogram.axes) != len(binning):
        raise ValueError(
            f"Histogram has {len(histogram.axes)} axes, but binning specifies {len(binning)} - axes multiplicity mismatch"
        )


def process_rel_n_hist(
    hist: bh.Histogram,
    binning: Dict[str, List[float]] = config.pid.sweight_binning,
) -> Dict[str, bh.Histogram]:
    """Validate the relative abundances, binwise"""
    rel_h = load_hist(hist)  # load the relative abundances histogram

    # Validate the binning before adding it to the container, excluding the species axes
    validate_n_axes(rel_h[..., -1], binning)
    validate_binning(rel_h[..., -1], binning)

    return rel_h


def register_pid_maps(
    opts: argparse.Namespace,
    binning: Dict[str, List[float]] = config.pid.pid_extrap_binning,
) -> Dict[str, bh.Histogram]:
    """Validate PID histogram structure against user-defined nominal binning"""
    pid_container = {}

    for key, value in vars(opts).items():
        if value is None:
            warnings.warn(
                f"Warning: {key}:{value} match passed - skipping histograms evaluating to 'None'",
                UserWarning,
            )
            continue  # skip [accounting for ghost exclusion]

        if "mu" in key:
            pid_hist = load_hist(value)  # Load muon/antimuon PID efficiencies

            # # Validate the binning before adding it to the container
            # validate_n_axes(pid_hist, binning)
            # validate_binning(pid_hist, binning)

            # If test passed, load histograms into container
            pid_container[key] = pid_hist

    return pid_container


def compute_misid_weight(
    p_val: float,
    eta_val: float,
    ntracks_val: int,
    pid_effs: Dict[str, bh.Histogram],
    relative_yields_hist: Union[bh.Histogram, hist.Hist],
    species: list = SPECIES,
) -> float:
    """Element-wise weight assignment, accounting for the fact that the per-species abundance and control- and signal-0 efficiecy maps may adopt different binning"""
    # validate binning of each
    misid_w = ufloat(0.0, 0.0)

    for spc in species:
        # relative abundance
        sweight_prefactor = relative_yields_hist[
            bh.loc(p_val), bh.loc(eta_val), bh.loc(ntracks_val), f"{spc}_yield"
        ]
        sw_spc = ufloat(sweight_prefactor.value, sweight_prefactor.variance**0.5)

        # control-channel efficiency
        ctrl_pideff = pid_effs[f"{spc}_to_antimu"][
            bh.loc(p_val), bh.loc(eta_val), bh.loc(ntracks_val)
        ]
        ctrl_pideff_v = ufloat(ctrl_pideff.value, ctrl_pideff.variance**0.5)

        # signal-channel efficiency
        sig_pideff = pid_effs[f"{spc}_to_mu"][
            bh.loc(p_val), bh.loc(eta_val), bh.loc(ntracks_val)
        ]
        sig_pideff_v = ufloat(sig_pideff.value, sig_pideff.variance**0.5)

        # linear combination: species-abundance prefactor (+) unfold the control PID eff (+) fold in the signal PID eff
        try:
            misid_w += sw_spc * (1.0 / ctrl_pideff_v) * sig_pideff_v
        except ZeroDivisionError:
            logger.warning(
                f"Division by zero encountered. ctrl_pideff_v: {ctrl_pideff_v}. Skipping this calculation."
            )

    return misid_w.n


if __name__ == "__main__":
    # get path to data file
    parser = argparse.ArgumentParser(
        description="tally all effs, yields and observations for the assignment of w_misid"
    )
    parser.add_argument("--obs", help="hadron-enriched data .root file")
    parser.add_argument("--rel_abundances", help="histogram of Ni/Nref")
    parser.add_argument("--proton_to_mu", help="path to proton->mu eff", default=None)
    parser.add_argument("--pion_to_mu", help="path to pion->mu eff", default=None)
    parser.add_argument("--kaon_to_mu", help="path to kaon->mu eff", default=None)
    parser.add_argument("--muon_to_mu", help="path to kaon->mu eff", default=None)
    parser.add_argument(
        "--electron_to_mu", help="path to electron->mu eff", default=None
    )
    parser.add_argument("--ghost_to_mu", help="path to ghost->mu eff", default=None)
    parser.add_argument(
        "--proton_to_antimu",
        help="path to proton->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    parser.add_argument(
        "--pion_to_antimu",
        help="path to pion->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    parser.add_argument(
        "--kaon_to_antimu",
        help="path to kaon->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    parser.add_argument(
        "--muon_to_antimu",
        help="path to kaon->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    parser.add_argument(
        "--electron_to_antimu",
        help="path to electron->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    parser.add_argument(
        "--ghost_to_antimu",
        help="path to ghost->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    parser.add_argument(
        "--output",
        help="path to output .root file",
    )
    parser.add_argument(
        "--key",
        help="key where trees are stored within a .root file",
    )
    parser.add_argument(
        "--tree", help="tree name within the I/O .root files", default="DecayTree"
    )
    opts = parser.parse_args()

    # pid maps for (un)folding into the misID weight
    pid_effs = register_pid_maps(opts)

    # relative abundance prefactors
    rel_n_sp = process_rel_n_hist(opts.rel_abundances)

    # partial for compliance with the polars API
    assign_misid_w = partial(
        compute_misid_weight, pid_effs=pid_effs, relative_yields_hist=rel_n_sp
    )

    # hadron-enriched observations as lazyframe for fast API
    data = pl.from_pandas(
        load_root(
            path=config.data.input_path,
            key=config.data.data_key,
            tree=config.data.data_tree,
            library="pd",
        )
    ).lazy()

    # compute misID weights and store into new branch `misid_w`
    data_update = data.with_columns(
        pl.struct(["Mu_plus_P", "Mu_plus_LK_ETA", "nTracks"])
        .map_elements(
            lambda x: assign_misid_w(
                p_val=x["Mu_plus_P"],
                eta_val=x["Mu_plus_LK_ETA"],
                ntracks_val=x["nTracks"],
            ),
            return_dtype=pl.Float64,
        )
        .alias("misid_w")
    ).collect()  # materialise lazyframe to dataframe

    # write outfile anew to avoid issues with mis-matched updates (plays better with Bc2DMuNu analysis pipeline)ary port to pandas to exploit uproot ability to write out files
    write_key_to_df(
        data=data_update.to_pandas(),
        path=config.data.output_path,
        key=config.data.data_key,
        treename=config.data.data_tree,
    )
