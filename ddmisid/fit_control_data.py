"""
Fit control data to extract trigger and PID efficiencies.
"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

import argparse
from iminuit import cost, Minuit
from sweights import SWeight  # for classic sweights
from sweights import Cow  # for custom orthogonal weight functions
from sweights import cov_correct, approx_cov_correct  # for covariance corrections
from ddmisid.utils import SimpleModelFactory, SanityChecks, dcbwg, dcbwg_cdf  # fitting
from ddmisid.utils import load_ntuple  # I/O
from ddmisid.utils import (
    simple_ax,
    plot_data,
    save_to,
    VisualizerFactory,
    viz_signal,
)  # plotting
import numpy as np
from typing import Callable

# number of mins considered in the NLL minimisation in the fit
FIT_BINS = 40
PLOT_BINS = 100

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fit control data")
    parser.add_argument(
        "-d", "--data", help="Path to data considered in the fit", required=True
    )
    args = parser.parse_args()

    # load the data
    # pep8 convention to signify immutable variables
    MC = load_ntuple(
        file_path=f"root://eoslhcb.cern.ch/{args.data}",
        key="B2JpsiMuNuXTuple",
        tree_name="DecayTree",
        max_entries=-1,
        branches=["B_plus_M"],
        cut="((Jpsi_L0MuonDecision_TOS==1) | (Jpsi_L0DiMuonDecision_TOS==1)) \
            & ( (B_plus_Hlt1TrackMVADecision_TOS==1) | (B_plus_Hlt1TwoTrackMVADecision_TOS==1) ) \
            & ( (B_plus_Hlt2Topo2BodyDecision_TOS==1) | (B_plus_Hlt2Topo3BodyDecision_TOS==1) )",  # companion muon trigger-unbiased
    )

    # mass rage considered - remove low-mass bkg for simplicity, taking advantage of the high-mass Bc peak
    MRANGE = (6000, 6600)

    # execute the MC unbinned ML fit to extract the parameters of interest, later fed to data models
    # --------------------------------------------------------------------------------------------
    # use DCB to fit the signal MC & perform a unbinned ML fit for speadup
    mc_model_pdf = lambda x, f1, f2, mug, sgg, sgl, sgr, al, ar, nl, nr: dcbwg(
        x, f1, f2, mug, mug, mug, sgg, sgl, sgr, al, ar, nl, nr, MRANGE
    )
    mc_model_cdf = lambda x, f1, f2, mug, sgg, sgl, sgr, al, ar, nl, nr: dcbwg_cdf(
        x, f1, f2, mug, mug, mug, sgg, sgl, sgr, al, ar, nl, nr, MRANGE
    )

    # minimise the cost [https://iminuit.readthedocs.io/en/stable/notebooks/cost_functions.html#Binned-Fit]
    nh, xe = np.histogram(MC["B_plus_M"].to_numpy(), bins=FIT_BINS, range=MRANGE)
    cost_obj = cost.BinnedNLL(nh, xe, mc_model_cdf)
    mi = Minuit(
        cost_obj,
        f1=0.3,
        f2=0.5,
        mug=6275,
        sgg=20,
        sgl=20,
        sgr=20,
        al=1,
        ar=1,
        nl=10,
        nr=10,
    )

    # define the parameter ranges
    mi.limits["f1"] = (0, 1)
    mi.limits["f2"] = (0, 1)
    mi.limits["mug"] = (6200, 6340)
    mi.limits["sgg"] = (0, 50)
    mi.limits["sgl"] = (0, 50)
    mi.limits["sgr"] = (0, 50)
    mi.limits["al"] = (0, 10)
    mi.limits["ar"] = (0, 10)
    mi.limits["nl"] = (0, 10)
    mi.limits["nr"] = (0, 10)

    # minimise and error estimation
    mi.migrad()
    mi.hesse()

    print(mi.params)
    print(mi.values)

    # sanity checks
    print("FCN min valid? ", mi.fmin.is_valid)

    # viz results
    # -----------
    fig, ax = simple_ax()

    # fit model; handle the norm of the pdf
    nhp, _ = np.histogram(MC.B_plus_M, range=MRANGE, bins=PLOT_BINS)
    _norm = np.sum(nhp) * (MRANGE[1] - MRANGE[0]) / PLOT_BINS
    _x = np.linspace(*MRANGE, PLOT_BINS)
    viz_signal(x=_x, y=_norm * mc_model_pdf(_x, *mi.values), ax=ax)(label="Fit")

    # observation
    plot_data(
        data=MC["B_plus_M"],
        range=MRANGE,
        bins=PLOT_BINS,
        ax=ax,
        label="Simulation 2016",
    )
    # cosmetics & save
    ax.legend()
    ax.set_yscale("log")
    ax.set_xlabel(r"$m(J/\psi \, \pi^+)$ [MeV$/c^2$]")
    save_to(outdir="test_plots", name="test")
