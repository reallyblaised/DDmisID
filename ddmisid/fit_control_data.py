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
from ddmisid.utils import (
    pdf_factory,
    cdf_factory,
    SanityChecks,
    pdf_factory,
    composite_pdf_factory,
)  # fitting
from ddmisid.utils import load_ntuple, load_root  # I/O
from ddmisid.utils import (
    simple_ax,
    plot_data,
    save_to,
    viz_signal,
    hist_err,
    hist_step_fill,
    fill_hist_w,
    make_legend,
)  # plotting
import numpy as np
from typing import Callable, Any
from uncertainties import ufloat
from scipy.stats import expon
from numpy.typing import ArrayLike
from functools import partial


# number of mins considered in the NLL minimisation in the fit
FIT_BINS = 40
PLOT_BINS = 100

# event section
# -------------
# @L0: trigger on jpsi, leave muon unbiased
# @HLT1: trigger on B inclusively without using the singleµ lines, leave muon unbiased
# @HLT2: trigger on B inclusively without using the singleµ lines, leave muon unbiased
# @Muon ID: enforce anti-muon criteria to retain a pion-like track
# @Kinematics and tracking: cuts aligned with event selection
# @B mass: keep in (6.0, 6.6) GeV window to retain Bc+ signal but cut out misID and other low-mass bkgs

# & (Mu_plus_isMuon==False) & (Mu_plus_PIDmu<0) & (Mu_plus_PIDK<0) \
# & (Mu_plus_IPCHI2_OWNPV>4.8) \
# & (Mu_plus_PT>750) & (Mu_plus_P>1e3) & (Jpsi_PT>2000) & (Jpsi_M>3040) & (Jpsi_M<3150) & ((Jpsi_ENDVERTEX_CHI2/Jpsi_ENDVERTEX_NDOF)<9) & (Mu_1_PT>900) & (Mu_2_PT>900) \
# & (B_plus_DOCA_Jpsi_Mu_plus<0.15) & (B_plus_ISOLATION_BDT<0.2) & ((B_plus_ENDVERTEX_Z-B_plus_PV_Z)>0) & (B_plus_ENDVERTEX_CHI2<25) \

EVT_SELECTION = "((muplus_L0MuonDecision_TOS==1) | (muminus_L0MuonDecision_TOS==1)) \
            & (K_PIDK>3) \
            & (B_M>5200) \
            & (B_M<5400) \
            & (abs(B_TRUEID)==521) & (Jpsi_TRUEID==443) & (abs(Jpsi_MC_MOTHER_ID)==521) & (abs(muplus_TRUEID)==13) & (muplus_MC_MOTHER_ID==443) & (abs(muminus_TRUEID)==13) & (muminus_MC_MOTHER_ID==443) & (abs(K_TRUEID)==321) & (abs(K_MC_MOTHER_ID)==521)"

# branches of interest
BRANCHES = [
    "B_M",
    "K_P",
    "K_PT",
    "K_PIDmu",
    "K_PIDK",
    "K_TRACK_MatchCHI2",
    "K_TRACK_GhostProb",
    "K_ProbNNghost",
]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fit control data")
    parser.add_argument(
        "-s", "--sim", help="Path to mc simulation considered in the fit", required=True
    )
    parser.add_argument(
        "-d", "--data", help="Path to data considered in the fit", required=True
    )
    args = parser.parse_args()

    # load the data
    # pep8 convention to signify immutable variables
    MC = load_ntuple(
        file_path=f"root://eoslhcb.cern.ch/{args.sim}",
        key="B2DsMuNuTuple",
        tree_name="DecayTree",
        max_entries=-1,
        branches=BRANCHES
        + [
            "K_TRUEID",
            "B_TRUEID",
            "K_MC_MOTHER_ID",
        ],
        cut=EVT_SELECTION,
    )
    # mass rage considered - remove low-mass bkg for simplicity, taking advantage of the high-mass Bc peak
    MRANGE = (5200, 5400)

    # execute the MC unbinned ML fit to extract the parameters of interest, later fed to data models
    # --------------------------------------------------------------------------------------------
    # use DCB to fit the signal MC & perform a unbinned ML fit for speadup
    dcbgauss_pdf = pdf_factory(mrange=MRANGE, key="signal")
    dcbgauss_cdf = cdf_factory(mrange=MRANGE, key="signal")

    # minimise the cost [https://iminuit.readthedocs.io/en/stable/notebooks/cost_functions.html#Binned-Fit]
    nh, xe = np.histogram(MC["B_M"].to_numpy(), bins=FIT_BINS, range=MRANGE)
    cost_obj = cost.BinnedNLL(nh, xe, dcbgauss_cdf)
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
    mi.limits["mug"] = (5200, 5400)
    mi.limits["sgg"] = (0, 50)
    mi.limits["sgl"] = (0, 50)
    mi.limits["sgr"] = (0, 50)
    mi.limits["al"] = (0, 10)
    mi.limits["ar"] = (0, 10)
    mi.limits["nl"] = (0, 20)
    mi.limits["nr"] = (0, 20)

    # minimise and error estimation
    mi.migrad()
    mi.hesse()

    # sanity checks
    SanityChecks(mi)()

    # viz results
    # -----------
    fig, ax = simple_ax()

    # fit model; handle the norm of the pdf
    nhp, _ = np.histogram(MC.B_M, range=MRANGE, bins=PLOT_BINS)
    _norm = np.sum(nhp) * (MRANGE[1] - MRANGE[0]) / PLOT_BINS
    _x = np.linspace(*MRANGE, PLOT_BINS)
    viz_signal(x=_x, y=_norm * dcbgauss_pdf(_x, *mi.values), ax=ax)(label="Fit")

    # observation
    plot_data(
        data=MC.B_M,
        range=MRANGE,
        bins=PLOT_BINS,
        ax=ax,
        label="Simulation 2016",
    )
    # cosmetics & save
    ax.legend()
    ax.set_yscale("log")
    ax.set_xlabel(r"$m(J/\psi \, K^+)$ [MeV$/c^2$]")
    save_to(outdir="test_plots", name="test")

    # save the mc fit pars
    mc_postfit_pars = {}
    for p in mi.parameters:
        mc_postfit_pars[p] = ufloat(mi.params[p].value, mi.params[p].error)

    # RFE: saving here -> fill later

    # # ghost-prob and track-match chi2 study
    # # -------------------------------------

    # # relax truth conditions
    # EVT_SELECTION = "((muplus_L0MuonDecision_TOS==1) | (muminus_L0MuonDecision_TOS==1)) \
    #         & (K_PIDK>3) \
    #         & (B_M>5200) \
    #         & (B_M<5400)"

    # MC = load_ntuple(
    #     file_path=f"root://eoslhcb.cern.ch/{args.sim}",
    #     key="B2DsMuNuTuple",
    #     tree_name="DecayTree",
    #     max_entries=-1,
    #     branches=BRANCHES,
    #     cut=EVT_SELECTION,
    # )

    # # having sourced the MC, can now define our study samples
    # ghost_matched_b_mu = MC[(MC.K_TRUEID == 0) & (MC.K_MC_MOTHER_ID == 0)]
    # ghost_matched_mu = MC[MC.K_TRUEID == 0]
    # hi_match_chi2 = MC[MC.K_TRACK_MatchCHI2 > 30]
    # hi_match_ghostprob = MC[MC.K_TRACK_GhostProb > 0.3]

    # # plots
    # # -----
    # # is the track-match chi2 different amongs the various definition of ghosts?
    # for var in (
    #     "K_TRACK_MatchCHI2",
    #     "K_TRACK_GhostProb",
    # ):
    #     fig, ax = simple_ax()
    #     match var:
    #         case "K_TRACK_MatchCHI2":
    #             _range = (0, 50)
    #         case "K_TRACK_GhostProb":
    #             _range = (0, 1)

    #     plt_config = {"range": _range, "ax": ax, "bins": 25, "norm": True}
    #     hist_step_fill(
    #         data=ghost_matched_b_mu[var],
    #         label=r"$\{B, K\}$ \texttt{TRUEID==0}",
    #         **plt_config,
    #     )
    #     plot_data(
    #         data=ghost_matched_mu[var],
    #         label=r"$K$ \texttt{TRUEID==0}",
    #         **plt_config,
    #     )
    #     hist_step_fill(
    #         data=MC[var],
    #         label=r"$B_c^+ \to J/\psi \, K^+$ MC",
    #         **plt_config,
    #     )

    #     make_legend(ax)
    #     ax.set_xlabel(r"\texttt{%s}" % var)
    #     save_to(outdir="test_plots", name=f"{var}")

    # # compare distributions
    # for var in [
    #     "K_PT",
    #     "K_P",
    #     "K_PIDmu",
    #     "K_PIDK",
    #     "K_ProbNNghost",
    # ]:

    #     # plot-range config
    #     fig, ax = simple_ax()
    #     match var:
    #         case "K_PT":
    #             _range = (0, 1e4)
    #         case "K_P":
    #             _range = (0, 1e5)
    #         case "K_PIDmu" | "K_PIDK":
    #             _range = (-15, 15)
    #         case "K_ProbNNghost":
    #             _range = (0, 1)
    #     plt_config = {"range": _range, "ax": ax, "bins": 20, "norm": True}

    #     plot_data(
    #         data=ghost_matched_mu[var],
    #         label=r"$K^+$ \texttt{TRUEID==0}",
    #         **plt_config,
    #     )
    #     hist_step_fill(
    #         data=hi_match_chi2[var],
    #         label=r"Track-match $\chi^2 > 30$",
    #         **plt_config,
    #     )
    #     hist_step_fill(
    #         data=hi_match_ghostprob[var],
    #         label=r"Ghost prob. $> 0.3$",
    #         **plt_config,
    #     )

    #     make_legend(ax)
    #     ax.set_xlabel(r"\texttt{%s}" % var)
    #     save_to(outdir="test_plots", name=f"{var}")

    # =================================================================================================
    #                                        DATA SECTION
    # =================================================================================================
    EVT_SELECTION = "((muplus_L0MuonDecision_TOS==1) | (muminus_L0MuonDecision_TOS==1)) \
            & (K_PIDK>3) \
            & (B_M>5200) \
            & (B_M<5400)"
    DATA = load_root(
        file_path=f"root://eoslhcb.cern.ch/{args.data}",
        key="B2DsMuNuTuple",
        tree_name="DecayTree",
        library="pd",
        max_entries=-1,
        branches=BRANCHES,
        cut=EVT_SELECTION,
    )
    # generate a composite model for the fit:
    # - signal: 2CB + gauss
    # - combinatorial: exp

    # fix paramters extracted from MC fit
    data_model = composite_pdf_factory(
        mrange=MRANGE,
        key="twoclass",
    )

    nh, xe = np.histogram(DATA["B_M"], bins=FIT_BINS, range=MRANGE)
    cost_obj = cost.ExtendedBinnedNLL(nh, xe, data_model)
    mi = Minuit(
        cost_obj,
        f1=mc_postfit_pars["f1"].n,
        f2=mc_postfit_pars["f2"].n,
        mug=6275,  # shared mean
        sgg=20,
        sgl=20,
        sgr=20,
        al=mc_postfit_pars["al"].n,
        ar=mc_postfit_pars["ar"].n,
        nl=mc_postfit_pars["nl"].n,
        nr=mc_postfit_pars["nr"].n,
        lb=200,
        sig_yield=0.98 * len(DATA["B_M"]),
        comb_yield=0.2 * len(DATA["B_M"]),
    )
    breakpoint()

    # define the parameter ranges
    mi.limits["sig_yield"] = (0, len(DATA["B_M"]))
    mi.limits["comb_yield"] = (0, len(DATA["B_M"]))
    mi.limits["mug"] = (5200, 5400)
    mi.limits["sgg"] = (0, 50)
    mi.limits["sgl"] = (0, 50)
    mi.limits["sgr"] = (0, 50)
    mi.limits["lb"] = (0, 500)
    mi.fixed["f1"] = True
    mi.fixed["f1"] = True
    mi.fixed["al"] = True
    mi.fixed["ar"] = True
    mi.fixed["nl"] = True
    mi.fixed["nr"] = True

    # minimise and error estimation
    mi.migrad()
    mi.hesse()

    # sanity checks
    SanityChecks(mi)()
