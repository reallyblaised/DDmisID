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
from uncertainties import ufloat
from functools import partial


def extract_axis_name(axis: bh.axis) -> str:
    """Extract the name of histogram axis"""
    try:
        return axis.metadata.get("name")
    except:
        return axis.name


def validate_binning(
    histogram: bh.Histogram,
    binning: Dict[str, np.ndarray | list],
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


def process_rel_n_hist(
    hist: bh.Histogram,
    binning: Dict[str, np.ndarray | list],
    verbose: bool = False,
) -> Dict[str, bh.Histogram]:
    """Validate the relative abundances, binwise"""
    rel_h = load_hist(hist)  # load the relative abundances histogram

    # Validate the binning before adding it to the container, excluding the species axes
    validate_n_axes(rel_h[..., -1], binning)
    validate_binning(rel_h[..., -1], binning)

    # verbose print of dict
    if verbose:
        pprint.pprint(rel_h, width=40, indent=4)

    return rel_h


def compute_misid_w(
    i: int,
    j: int,
    k: int,
    pid_effs: Dict[str, bh.Histogram],
    relative_yields: bh.Histogram,
    species: tuple = ("electron", "kaon", "pion", "proton"),
) -> float:
    """Exract the misID weight, assuming index order p, eta, ntracks [checks in place]"""

    # sanity
    assert (
        abs(relative_yields[i, j, k, ...].sum().value - 1.0) < 1e-5
    ), f"Normalisation check fail: relative yields in bin ({i}, {j}, {k}) are incompatible with normalisation to unity and the 1e-5 level"

    # FIXME: consistency of binning order and edges among the pid hs
    # FIXME: make sure i,j,k are correctly mapped on the pid axes

    # consistency check
    misid_w = 0.0
    for s in species:
        # fetch the factors in the species-specific contributions to the total misID weight, incorporating the respective uncertainties
        ni_nref = ufloat(
            relative_yields[i, j, k, f"{s}_yield"].value,
            relative_yields[i, j, k, f"{s}_yield"].variance ** 0.5,
        )
        num_pid_eff = ufloat(
            pid_effs[f"{s}_to_mu"][i, j, k].value,
            pid_effs[f"{s}_to_mu"][i, j, k].variance ** 0.5,
        )
        denom_pid_eff = ufloat(
            pid_effs[f"{s}_to_antimu"][i, j, k].value,
            pid_effs[f"{s}_to_antimu"][i, j, k].variance ** 0.5,
        )

        # species weight factor
        species_w = ni_nref * (1 / denom_pid_eff) * num_pid_eff

        # linear combination of species factors
        misid_w += species_w

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
        "--electron_to_antimu",
        help="path to electron->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    parser.add_argument(
        "--ghost_to_antimu",
        help="path to ghost->!mu eff [ie containted in hadron-enriched data]",
        default=None,
    )
    opts = parser.parse_args()

    # user-defined binning
    binning = read_config("config/main.yml", key="pid")["binning"]

    # process input in dedicated containers
    pid_effs = process_pid_hists(opts, binning, False)

    # relative abundance prefactors
    rel_n_sp = process_rel_n_hist(opts.rel_abundances, binning)

    # hadron-enriched observations as lazyframe
    data = pl.from_pandas(
        simple_load(
            path=read_config("config/main.yml", key="data")["path"],
            key=read_config("config/main.yml", key="data")["root_config"]["root_key"],
            library="pd",
        )
    ).lazy()

    # facilitate the lazy lambda syntax by assignining the kwargs independnent of bin coordinates
    compute_misid_w_binwise = partial(
        compute_misid_w,
        pid_effs=pid_effs,
        relative_yields=rel_n_sp,
    )

    # Apply binning and pass both the value and bin index to the custom function
    data_update = (
        data.with_columns(
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
        )
        .with_columns(
            pl.struct(["P_idx", "ETA_idx", "nTracks_idx"])
            .map_elements(
                lambda x: compute_misid_w_binwise(
                    int(x["P_idx"]), int(x["ETA_idx"]), int(x["nTracks_idx"])
                ),
                return_dtype=pl.Float64,
            )
            .alias("misid_w")
        )
        .collect()  # materialise lazyframe to dataframe
        # .to_pandas()  # covert back to pandas to be able to write out via uproot
    )

    breakpoint()
