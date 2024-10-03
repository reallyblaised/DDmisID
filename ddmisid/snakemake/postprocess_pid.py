"""
Handle the cases where effic
"""

from ddmisid import load_hist, _pid_eff_tolerance
import numpy as np
import warnings
import pickle
import argparse


def process_eff(hist):
    """Squash to zero efficiencies that are compatible with zero within 1 sigma"""
    # Access the storage (weights and variances) of the histogram
    view = hist.view()

    # Iterate over multi-dimensional indices in the histogram
    for index in np.ndindex(view.shape):
        # Extract the value and variance from the weighted storage
        cval = view[index].value
        variance = view[index].variance

        # case by case handling of the PID efficiency values
        if cval < 0.0:
            if abs(cval) < variance**0.5:  # <0 and within 1 sigma from 0
                warnings.warn(
                    f"PID efficiency value below 0.0 detected at index {index} within 1 sigma of 0.0: {view[index]}. Setting to 0.0"
                )
                view[index] = (
                    0.0,
                    abs(cval),  # FIXME: ascertain this is a sensible choice
                )
            # within tolerance
            if (abs(cval) > variance**0.5) and (
                abs(cval) < _pid_eff_tolerance
            ):  # <0 and more than 1 sigma from 0 but within tolerance
                warnings.warn(
                    f"PID efficiency value below 0.0 detected at index {index} with >1 sigma of 0.0 but below tolerance of {_pid_eff_tolerance}: {view[index]}. Setting to 0.0"
                )
                view[index] = (
                    0.0,
                    abs(cval),  # FIXME: ascertain this is a sensible choice
                )
            # kill the process if the PID efficiency is below 0.0 and outside tolerance
            if (abs(cval) > variance**0.5) and (abs(cval) > _pid_eff_tolerance):
                raise ValueError(
                    f"PID efficiency value below 0.0 detected at index {index} outside 1 sigma of 0.0 and above tolerance: {view[index]}. Aborting."
                )
        else:
            pass

    return hist


# write me argparse that has input and output file
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Load the PID efficiency histograms and set entries to zero if within tolerance"
    )
    parser.add_argument("-i", "--input", help="Input PID efficiency file")
    parser.add_argument("-o", "--output", help="Output PID efficiency file")
    opts = parser.parse_args()

    # load PID histograms
    pid_eff = process_eff(load_hist(opts.input))

    # set entries to zero if within tolerance
    with open(opts.output, "wb") as f_out:
        pickle.dump(pid_eff, f_out)
