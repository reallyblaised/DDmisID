"""
Load the PID efficiency histograms and set entries to zero if within tolerance
"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

from ddmisid import load_hist, _pid_eff_tolerance
import argparse

# write me argparse that has input and output file
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Load the PID efficiency histograms and set entries to zero if within tolerance"
    )
    parser.add_argument("-i", "--input", help="Input PID efficiency file")
    parser.add_argument("-o", "--output", help="Output PID efficiency file")
    opts = parser.parse_args()

    # load PID histograms
    pid_eff = load_hist(opts.input)

    breakpoint()
    # set entries to zero if within tolerance