"""Generate the PIDCalib2 binning spec based on the config YAML file."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import json


class PIDCalibAliasFactory:
    """Factory to process variable names based on year and species."""

    @staticmethod
    def process_variable(var: str, year: str, spc: str) -> str:
        """
        Apply year-based processing to variable names.
        Source: https://twiki.cern.ch/twiki/bin/view/LHCb/PIDCalibPackage
        """
        # Run 1 vs Run 2 aliases
        if year in ["2016", "2017", "2018"]:
            if var in ["P", "PT"]:
                var = f"Brunel_{var}"
            elif var == "nTracks":
                var = f"{var}_Brunel"
        elif year not in ["2011", "2012", "2015"]:
            raise ValueError(f"Year {year} not recognized.")

        # bespoke handling of electrons
        if spc == "e_B_Jpsi" and year == "2016":
            if var in ["Brunel_P", "Brunel_PT"]:
                var = var.replace("Brunel_", "")
            elif var == "nTracks_Brunel":
                var = var.replace("_Brunel", "")

        return var


class BinningGeneratorBase(ABC):
    """Abstract base class for binning generators."""

    @abstractmethod
    def fetch_info(self) -> Tuple[Dict[str, str], Dict[str, Dict[str, List[float]]]]:
        """Fetch the user-defined binning information."""
        pass

    @abstractmethod
    def build(self, year: str, outdir: str = ".data", verbose: bool = True) -> None:
        """Build the binning spec JSON file."""
        pass


@dataclass(frozen=True)
class DefaultBinningGenerator(BinningGeneratorBase):
    """Default binning generator for PIDCalib2."""

    species: Dict[str, str]
    binning: Dict[str, Dict[str, List[float]]]
    binning_alias: Optional[str] = None

    def fetch_info(self) -> Tuple[Dict[str, str], Dict[str, Dict[str, List[float]]]]:
        """Fetch the user-defined binning information for the current configuration."""
        return self.species, self.binning

    def build(self, year: str, outdir: str = ".data", verbose: bool = True) -> None:
        """Build the JSON spec file for the supplied years and species."""
        bins, species = self.fetch_info()

        for spc in species.values():
            binning2json = {spc: {}}

            for pidc_alias, edges in bins.items():
                # ensure correct PIDCalib2 variable-alias compatibility
                pidc_alias = PIDCalibAliasFactory.process_variable(
                    pidc_alias, year, spc
                )
                binning2json[spc][pidc_alias] = edges

            # persist to JSON file
            if self.binning_alias:
                outfile_path = Path(
                    f"{outdir}/binning_{year}/{self.binning_alias}.json"
                )
            else:
                outfile_path = Path(f"{outdir}/binning_{year}/{spc}.json")
            outfile_path.parent.mkdir(parents=True, exist_ok=True)
            with open(outfile_path, "w") as f:
                json.dump(binning2json, f, indent=4)

            if verbose:
                print(f"Binning file for {spc} written to {outfile_path}")
