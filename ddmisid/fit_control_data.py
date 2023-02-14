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
from ddmisid import config_simple_model, SanityChecks  # fitting
from ddmisid import load_ntuple  # I/O
from ddmisid import simple_ax, plot_data, save_to  # plotting


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fit control data")
    parser.add_argument(
        "-d", "--data", help="Path to data considered in the fit", required=True
    )
    args = parser.parse_args()

    # load the data
    # pep8 convention to signify immutable variables
    DATA = load_ntuple(
        file_path=f"root://eoslhcb.cern.ch/{args.data}",
        key="B2JpsiMuNuXTuple",
        tree_name="DecayTree",
        max_entries=1e4,
        branches=["B_plus_M"],
    )

    # mass rage considered - remove low-mass bkg for simplicity, taking advantage of the high-mass Bc peak
    MRANGE = (6000, 6600)

    # execute the extended unbinned ML fit
    # ------------------------------------
    # define appropriate cost function
    twoclass_model = config_simple_model(mrange=MRANGE)

    cost_obj = cost.ExtendedUnbinnedNLL(
        DATA["B_plus_M"].to_numpy(),  # fit the beauty invariant-mass spectrum
        twoclass_model,
    )
    mi = Minuit(cost_obj, ns=80e3, nb=20e3, mu=6280, sg=20, lb=400)

    # define the parameter ranges
    mi.limits["ns"] = (0, 100e3)
    mi.limits["nb"] = (0, 100e3)
    mi.limits["mu"] = (6200, 6340)
    mi.limits["sg"] = (0, 50)
    mi.limits["lb"] = (0, 5000)

    # minimise and error estimation
    mi.migrad()
    mi.hesse()

    # sanity checks
    SanityChecks(mi)

    # viz results
    # -----------
    fig, ax = simple_ax()
    plot_data(
        data=DATA["B_plus_M"].to_numpy(),
        range=MRANGE,
        ax=ax,
    )
    save_to(outdir="test", name="test")
