"""Execute binned maximum-likelihood (BML) fit to reco partiions of a nominal bin of kinematics and occupancy"""

__authors__ = ["Blaise Delaney", "Kevin Kurashima"]
__email__ = "blaise.delaney at cern.ch"

import argparse
import pathlib
from ddmisid.fitter import BMLFitter

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(
        description="BML fit machinery to reco partitions."
    )

    # Add the arguments
    parser.add_argument(
        "-d",
        "--data",
        type=str,
        required=True,
        help="Path to observations considered in the fit.",
    )

    parser.add_argument(
        "-t",
        "--template_path",
        type=str,
        required=True,
        help="The path to the binned templates",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="The path to the output file containing the fit results",
    )

    # Parse the arguments
    args = parser.parse_args()

    # Instantiate the fitter
    fitter = BMLFitter(
        args.template_path, args.data
    )  # FIXME: should interface path and config
    fitter.fit()
    fitter.visualize()

    # Save the results
    fitter.save(filename=args.output)
