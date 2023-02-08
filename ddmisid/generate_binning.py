"""generate and write the json binning file for pidcalib2

Use config/main.yml to generate binning for each variable for pidcalib2

__author__ = Blaise Delaney
__email__  = blaise.delaney at cern.ch
"""

from dataclasses import dataclass
from ddmisid.utils import read_config
from typing import Dict, Tuple
import json
from pathlib import Path
from termcolor2 import c as tc


@dataclass(frozen=True)
class BinningGenerator:
    path: str

    def gather_info(
        self,
        _key: str = "pid",
    ) -> Tuple[Dict, Dict]:
        """Generate the json binning file for pidcalib2

        Parameters
        ----------
        key : str
            The key in the config file to read from (default: "pid")

        Returns
        -------
        Tuple[Dict, Dict]
            The edges and species for binning entry
        """
        config = read_config(self.path, key=_key)
        edges = config["binning"]
        species = config["species"]

        return edges, species

    def __post_init__(self) -> None:
        """Attribute detailing the sourced binning parameters"""
        _edges, _species = self.gather_info()
        print(
            tc(f"\nUser-defined binning scheme:").underline.blue,
            f"\n* species: {_species}\n* variables: {_edges}\n",
        )

    def build(
        self, year: str, outdir: str = "data", print_outpath: bool = True, **kwargs
    ) -> None | str:
        """Generate the json binning file for pidcalib2

        Parameters
        ----------
        outdir : str
            The output directory to write the binning file to

        year: str
            Year of data taking identifier

        print_outpath : bool
            If true, print the output path to the binning file (default: True)

        Returns
        -------
        None
        """
        bins, species = self.gather_info(**kwargs)

        binning2json = {}
        for spc in species.values():
            binning2json[spc] = {}
            for var, edges in bins.items():
                match year:
                    case "2016" | "2017" | "2018":
                        match var:
                            case "P" | "PT" as alias:
                                var = f"Brunel_{alias}"
                            case "nTracks" as alias:
                                var = f"{alias}_Brunel"
                    case "2011" | "2012" | "2015":
                        pass
                    case other:
                        raise ValueError("Year identifier not recognised")

                # HACK: 2016 electrons dont follow the convention
                if spc == "e_B_Jpsi" and year == "2016":
                    match var:
                        case "Brunel_P" | "Brunel_PT" as alias:
                            var = alias.replace("Brunel_", "")
                        case "nTracks_Brunel" as alias:
                            var = alias.replace("_Brunel", "")
                binning2json[spc][var] = edges

        Path(outdir).mkdir(parents=True, exist_ok=True)
        outfile_path = f"{outdir}/binning_{year}.json"
        with open(f"{outfile_path}", "w") as f:
            json.dump(binning2json, f)

        if print_outpath:
            return outfile_path
