from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import json


class PIDCalibAliasFactory:
    """Factory to process variable names based on year and species."""

    @staticmethod
    def process_variable(var: str, year: str, spc: str | None = None) -> str:
        """
        Apply year-based processing to variable names.
        Source: https://twiki.cern.ch/twiki/bin/view/LHCb/PIDCalibPackage

        Parameters
        ----------
        var : str
            The binning variable to process (e.g., 'P', 'PT', 'ETA', 'nTracks').
        year : str
            The year of data taking (e.g., '2018', '2016').

        Returns
        -------
        str
            The processed variable name based on the year.
        """
        if var in ["P", "PT", "ETA"]:
            if year in ["2016", "2017", "2018"]:
                return f"Brunel_{var}"
            elif year in ["2011", "2012", "2015"]:
                return var
        elif var == "nTracks":
            if year in ["2016", "2017", "2018"]:
                return f"{var}_Brunel"
            elif year in ["2011", "2012", "2015"]:
                return var
        else:
            raise ValueError(f"Variable {var} not supported for year {year}.")
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
    """
    Default binning generator for PIDCalib2.

    Attributes:
    ----------
    species : dict
        Mapping of species (e.g., 'kaon', 'pion', etc.) to their aliases.
    binning : dict
        Binning information with variable names as keys and bin edges as values.
    binning_alias : str, optional
        An optional alias for the binning file.
    """

    species: Dict[str, str]
    binning: Dict[str, Dict[str, List[float]]]
    binning_alias: Optional[str] = None

    def fetch_info(self) -> Tuple[Dict[str, str], Dict[str, Dict[str, List[float]]]]:
        """Fetch the user-defined binning information for the current configuration."""
        return self.species, self.binning

    def get_binning_variables(self, year: str) -> List[str]:
        """
        Return the list of binning variables with the correct PIDCalib2 variable aliases.

        Parameters:
        ----------
        year : str
            The year for which to generate the binning variable aliases.

        Returns:
        -------
        List[str]
            A list of binning variables processed with the correct aliases.
        """
        return [
            PIDCalibAliasFactory.process_variable(var, year)
            for var in self.binning.keys()
        ]

    def build(self, year: str, outdir: str = ".data", verbose: bool = True) -> None:
        """
        Build the JSON spec file for the supplied years and species.

        Parameters:
        ----------
        year : str
            The year of data taking (e.g., '2018').
        outdir : str
            The output directory for storing the binning JSON files.
        verbose : bool
            If True, prints out the path to the written binning files.

        Returns:
        -------
        None
        """
        species, bins = self.fetch_info()

        for spc in species.values():
            binning_json = {spc: {}}

            for var, edges in bins.items():
                # Ensure correct PIDCalib2 variable-alias compatibility
                processed_var = PIDCalibAliasFactory.process_variable(var, year)
                binning_json[spc][processed_var] = edges

            # Create the output file path
            output_file = self.generate_output_filepath(spc, year, outdir)

            # Write the binning data to the JSON file
            output_file.parent.mkdir(parents=True, exist_ok=True)
            with open(output_file, "w") as f:
                json.dump(binning_json, f, indent=4)

            if verbose:
                print(f"Binning file for {spc} written to {output_file}")

    def generate_output_filepath(self, species: str, year: str, outdir: str) -> Path:
        """
        Generate the file path for the binning JSON file.

        Parameters:
        ----------
        species : str
            The species (e.g., 'kaon', 'pion') to generate the file path for.
        year : str
            The year of data taking (e.g., '2018').
        outdir : str
            The output directory for storing the binning JSON files.

        Returns:
        -------
        Path
            The path for the binning JSON file.
        """
        if self.binning_alias:
            return Path(f"{outdir}/binning_{year}/{self.binning_alias}.json")
        else:
            return Path(f"{outdir}/binning_{year}/{species}.json")
