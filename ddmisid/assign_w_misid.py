"""
Assign the w_misid weight following the main equation
   
    w_misid = sum_i Ni/Nred eff_i(h->mu)/eff_i(!mu)

see LHCb-ANA-2021-052
"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

import polars as pl
from ddmisid.utils import read_config, simple_load
import argparse



if __name__ == "__main__":
    # get path to data file
    parser = argparse.ArgumentParser(
        description="tally all effs, yields and observations for the assignment of w_misid"
    )
    opts = parser.parse_args()

    # binning
    binning = read_config("config/main.yml", key="pid")["binning"]

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
