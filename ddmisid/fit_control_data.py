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
    composite_cdf_factory,
    twoclass_cdf,
    twoclass_pdf,
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
    eff_plot,
)  # plotting
import numpy as np
from typing import Callable, Any
from uncertainties import ufloat
from scipy.stats import expon
from numpy.typing import ArrayLike
from functools import partial

# efficiency error calculation
def eff_err(n_pass: int, n_tot: int) -> float:
    """Efficiency error calculation.
    Source: https://lss.fnal.gov/archive/test-tm/2000/fermilab-tm-2286-cd.pdf (eq. 2)

    Parameters
    ----------
    n_pass: int
        Number of events passing the selection

    n_tot: int
        Total number of events

    Returns
    -------
    float
        Efficiency error
    """
    return np.sqrt(n_pass * (1 - n_pass / n_tot)) / n_tot


# template for the event selection: simpe
data_selection = (
    lambda m_min, m_max: f"((muplus_L0MuonDecision_TOS==1) | (muminus_L0MuonDecision_TOS==1)) \
            & (B_M>{m_min}) \
            & (B_M<{m_max})"
)

# pep8 for immutables
# number of mins considered in the NLL minimisation in the fit
FIT_BINS = 40
PLOT_BINS = 100
MASS_MIN = 5200
MASS_MAX = 5400

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

    # load the mc, with appropriate truthmatching
    mc_sel = f"{data_selection(MASS_MIN, MASS_MAX)} \
        & (abs(B_TRUEID)==521) & (Jpsi_TRUEID==443) & (abs(Jpsi_MC_MOTHER_ID)==521) & (abs(muplus_TRUEID)==13) & (muplus_MC_MOTHER_ID==443) & (abs(muminus_TRUEID)==13) & (muminus_MC_MOTHER_ID==443) & (abs(K_TRUEID)==321) & (abs(K_MC_MOTHER_ID)==521)"
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
        cut=mc_sel,
    )
    # mass rage considered - remove low-mass bkg for simplicity, taking advantage of the high-mass Bc peak
    MRANGE = (MASS_MIN, MASS_MAX)

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
    for scale in ("linear", "log"):
        ax.set_yscale(scale)
        ax.set_xlabel(r"$m(J/\psi \, K^+)$ [MeV$/c^2$]")
        save_to(outdir="test_plots", name=f"test_y{scale}")

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

    # Prefit the upper-mass sideband to get an estimate for lb in the exp pdf & then fix it in data fit
    # -------------------------------------------------------------------------------------------------
    _upper_mass_range = (5.5e3, 5.75e3)
    upper_mass = load_root(
        file_path=f"root://eoslhcb.cern.ch/{args.data}",
        key="B2DsMuNuTuple",
        tree_name="DecayTree",
        library="pd",
        max_entries=-1,
        branches=["B_M"],
        cut=data_selection(*_upper_mass_range),  # make sure we don't include Bc decays
    )

    # fit with decayin exponential
    expon_pdf = pdf_factory(mrange=_upper_mass_range, key="combinatorial")
    expon_cdf = cdf_factory(mrange=_upper_mass_range, key="combinatorial")

    # minimise the cost [https://iminuit.readthedocs.io/en/stable/notebooks/cost_functions.html#Binned-Fit]
    nh, xe = np.histogram(upper_mass.B_M, bins=FIT_BINS, range=_upper_mass_range)
    cost_obj = cost.BinnedNLL(nh, xe, expon_cdf)
    mi = Minuit(cost_obj, lb=100)

    # define the parameter ranges
    mi.limits["lb"] = (0, 1_000)

    # minimise and error estimation
    mi.migrad()
    mi.hesse()

    # sanity checks
    SanityChecks(mi)()

    # viz results
    # -----------
    fig, ax = simple_ax()

    # fit model; handle the norm of the pdf
    nhp, _ = np.histogram(upper_mass.B_M, range=_upper_mass_range, bins=PLOT_BINS)
    _norm = np.sum(nhp) * (_upper_mass_range[1] - _upper_mass_range[0]) / PLOT_BINS
    _x = np.linspace(*_upper_mass_range, PLOT_BINS)
    viz_signal(x=_x, y=_norm * expon_pdf(_x, *mi.values), ax=ax)(label="Fit")

    # observation
    plot_data(
        data=upper_mass.B_M,
        range=_upper_mass_range,
        bins=PLOT_BINS,
        ax=ax,
        label="Data 2016",
    )
    # cosmetics & save
    ax.legend()
    ax.set_xlabel(r"$m(J/\psi \, K^+)$ [MeV$/c^2$]")
    save_to(outdir="test_plots", name="test_upper_mass")

    # save the mc fit pars
    himass_postfit_pars = {}
    for p in mi.parameters:
        himass_postfit_pars[p] = ufloat(mi.params[p].value, mi.params[p].error)

    # fit to B+ -> J/psi K+ mass spectrum
    # -----------------------------------
    DATA = load_root(
        file_path=f"root://eoslhcb.cern.ch/{args.data}",
        key="B2DsMuNuTuple",
        tree_name="DecayTree",
        library="pd",
        max_entries=-1,
        branches=BRANCHES,
        cut=data_selection(*MRANGE),  # make sure we don't include Bc decays
    )

    # generate a composite model for the fit:
    # - signal: 2CB + gauss
    # - combinatorial: exp
    # data_model = composite_pdf_factory(
    #     mrange=MRANGE,
    #     key="twoclass",
    # )
    data_model_pdf = pdf_factory(mrange=MRANGE, key="signal")
    data_model_cdf = cdf_factory(mrange=MRANGE, key="signal")

    nh, xe = np.histogram(DATA.B_M, bins=FIT_BINS, range=MRANGE)
    cost_obj = cost.BinnedNLL(nh, xe, data_model_cdf)
    mi = Minuit(
        cost_obj,
        f1=mc_postfit_pars["f1"].n,
        f2=mc_postfit_pars["f2"].n,
        mug=5270,
        sgg=mc_postfit_pars["sgg"].n,
        sgl=mc_postfit_pars["sgl"].n,
        sgr=mc_postfit_pars["sgl"].n,
        al=mc_postfit_pars["al"].n,
        ar=mc_postfit_pars["ar"].n,
        nl=mc_postfit_pars["nl"].n,
        nr=mc_postfit_pars["nr"].n,
        # lb=himass_postfit_pars["lb"].n,
        # sig_yield=0.98 * len(DATA["B_M"]),
        # comb_yield=0.2 * len(DATA["B_M"]),
    )

    # define the parameter ranges
    # mi.limits["sig_yield"] = (0, len(DATA.B_M))
    # mi.limits["comb_yield"] = (0, len(DATA.B_M))
    mi.limits["mug"] = (5200, 5400)
    mi.limits["sgg"] = (0, 100)
    mi.limits["sgl"] = (0, 100)
    mi.limits["sgr"] = (0, 100)
    mi.limits["nl"] = (0, 10)
    mi.limits["nr"] = (0, 10)
    mi.limits["al"] = (0, 10)
    mi.limits["ar"] = (0, 10)
    # mi.limits["lb"] = (0, 500)

    # fix the parameters taken from MC/upper-mass fits
    _fixed = (
        "f1",
        "f2",
        "al",
        "ar",
        "nl",
        "nr",
        # "sgl",
        # "sgr",
        # "sgg",
        # "lb",
    )
    for fixpar in _fixed:
        mi.fixed[fixpar] = True

    # minimise and error estimation
    mi.migrad()
    mi.hesse()
    print(mi.params)

    # sanity checks
    SanityChecks(mi)()

    # viz results
    # -----------
    fig, ax = simple_ax()

    # observation
    plot_data(
        data=DATA.B_M,
        range=MRANGE,
        bins=PLOT_BINS,
        ax=ax,
        label="Data 2016",
    )

    # # plot the total pdf
    # _norm = (MRANGE[1] - MRANGE[0]) / PLOT_BINS
    # _x = np.linspace(*MRANGE, PLOT_BINS)
    # viz_signal(x=_x, y=_norm * data_model(_x, *mi.values), ax=ax)(label="Fit") # # plot the total pdf
    # _norm = (MRANGE[1] - MRANGE[0]) / PLOT_BINS
    # _x = np.linspace(*MRANGE, PLOT_BINS)
    # viz_signal(x=_x, y=_norm * data_model(_x, *mi.values), ax=ax)(label="Fit")

    # fit model; handle the norm of the pdf
    nhp, _ = np.histogram(DATA.B_M, range=MRANGE, bins=PLOT_BINS)
    _norm = np.sum(nhp) * (MRANGE[1] - MRANGE[0]) / PLOT_BINS
    _x = np.linspace(*MRANGE, PLOT_BINS)
    viz_signal(x=_x, y=_norm * data_model_pdf(_x, *mi.values), ax=ax)(label="Fit")

    # cosmetics & save
    ax.legend()
    for scale in ("linear", "log"):
        ax.set_yscale(scale)
        ax.set_xlabel(r"$m(J/\psi \, K^+)$ [MeV$/c^2$]")
        save_to(outdir="test_plots", name=f"test_data_y{scale}")

    # ==============================================================================
    #                               sFit / COWs section
    # ==============================================================================


    # ==============================================================================
    #                               Efficiency section
    # ==============================================================================

    # ------------------------------------------------------------------------------
    #                               Preliminaries
    # a) verify that the data and MC are consistent in Track-match CHI2 on "reak" K
    # ------------------------------------------------------------------------------
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
        cut=mc_sel,
    ) # true B2JpsiK

    # assume under the peak is all signal -> restrict mass range
    DATA = load_root(
        file_path=f"root://eoslhcb.cern.ch/{args.data}",
        key="B2DsMuNuTuple",
        tree_name="DecayTree",
        library="pd",
        max_entries=-1,
        branches=BRANCHES,
        cut=data_selection(
            5279 - 50, 5279 + 50
        ),  # tight mass window about nominal B+ mass
    )

    fig, ax = simple_ax()
    hist_step_fill(
        data = MC.K_TRACK_MatchCHI2,
        bins = 50,
        range= (0, 50),
        ax = ax,   
        label=r"$B^+ \to J/\psi K^+$ MC",
    )    
    plot_data(
        data = DATA.K_TRACK_MatchCHI2,
        bins = 50,
        range= (0, 50),
        ax = ax,   
        label=r"$B^+ \to J/\psi\,K^+$ Data (tight mass window)",
    )    
    make_legend(ax=ax, on_plot=False, ycoord=-0.4)
    ax.set_xlabel(r"Track-match $\chi^2$")
    for scale in ("linear", "log"):
        ax.set_yscale(f"{scale}")
        save_to(outdir="test_plots", name=f"track_match_chi2_consistency_y{scale}")



    fig, ax = simple_ax()
    hist_step_fill(
        data = MC.K_ProbNNghost,
        bins = 10,
        range= (0, 1),
        ax = ax,   
        label=r"$B^+ \to J/\psi K^+$ MC",
    )    
    plot_data(
        data = DATA.K_ProbNNghost,
        bins = 10,
        range= (0, 1),
        ax = ax,   
        label=r"$B^+ \to J/\psi\,K^+$ Data (tight mass window)",
    )    
    make_legend(ax=ax, on_plot=False, ycoord=-0.4)
    ax.set_xlabel(r"Track-match $\chi^2$")
    for scale in ("linear", "log"):
        ax.set_yscale(f"{scale}")
        save_to(outdir="test_plots", name=f"probnnghost_consistency_y{scale}")
    # ~ END OF PRELIMINARIES ~

    
    # ghosts in MC
    # ------------
    # load the mc, with appropriate truthmatching
    mc_ghost_sel = f"{data_selection(MASS_MIN, MASS_MAX)} \
        & (K_TRUEID == 0)"
    mc_ghosts = load_ntuple(
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
        cut=mc_ghost_sel,
    )
    sden_nh, sden_xe, sden_xc, _ = fill_hist_w(
        data=mc_ghosts.K_P,
        bins=10,
        range=(0, 1e5),
        weights=None,
    )
    snum_nh, _, _, _ = fill_hist_w(
        data=mc_ghosts[mc_ghosts.K_PIDmu > 3].K_P,
        bins=10,
        range=(0, 1e5),
        weights=None,
    )

    mc_eff = snum_nh / sden_nh
    mc_eff_err = eff_err(n_pass=snum_nh, n_tot=sden_nh)

    # data
    # ----
    
    dden_nh, dden_xe, dden_xc, dden_nh_err = fill_hist_w(
        data=DATA[DATA.K_TRACK_MatchCHI2 > 30].K_P,
        bins=10,
        range=(0, 1e5),
        weights=None,
    )
    dnum_nh, dnum_xe, dnum_xc, dnum_nh_err = fill_hist_w(
        data=DATA[(DATA.K_TRACK_MatchCHI2 > 30) & (DATA.K_PIDmu > 3)].K_P,
        bins=10,
        range=(0, 1e5),
        weights=None,
    )

    data_eff = dnum_nh / dden_nh
    data_eff_err = eff_err(n_pass=dnum_nh, n_tot=dden_nh)

    # eff plot
    # --------
    fig, ax = simple_ax()
    eff_plot(
        x=sden_xc,
        eff=mc_eff,
        eff_err=mc_eff_err,
        xerr=(sden_xe[1] - sden_xe[0]) / 2,
        label="MC TRUEID$(K)=0$",
        ax=ax,
    )
    eff_plot(
        x=dden_xc,
        eff=data_eff,
        eff_err=data_eff_err,
        xerr=(dden_xe[1] - dden_xe[0]) / 2,
        label="Data track-match $\chi^2>30$",
        ax=ax,
    )
    make_legend(ax=ax, on_plot=False, ycoord=-0.4)
    ax.set_ylabel(
        r"$\varepsilon($DLL$_{\mu}>3)~[\%]$"
    )  # supersedes the label in simple_ax
    ax.set_xlabel(r"$p$ [MeV/$c$]")
    save_to(outdir="test_plots", name=f"data_eff_test")


    