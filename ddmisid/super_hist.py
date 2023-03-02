"""Generate multidimensional efficiency histograms accounting for:

- calibration species (true hadrons, ghosts, electrons) [`calib_species`]
- reco partictions of the hadron-enriched data [`reco_partitions`]
- efficiencies, obtained via pidcalib2 via calibration samples, plus ghosts 
"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

import pickle
import hist
from ddmisid.utils import read_config


class PIDEffMixin:
    """Mixin for sanity checks for PID(calib) eff hists.

    Methods
    -------
    check_h_type
        Verify the hist.Histogram type
    open_hist_file
        Simple open of pkl file
    chech_h_axes
        Verity the consistency of the axes and the user-defined binning config
    """

    def open_hist_file(self) -> Any:
        """Open the pkl file to obtain the Histogram efficiency object"""
        with open(self.hist_path, "rb") as f_in:
            return pickle.load(f_in)

    @staticmethod
    def check_h_type(hpid: Any) -> None:
        """Verify that a hist-like object is of type Histogram"""
        if not isinstance(hpid, hist.Histogram):
            raise TypeError("The hpid object must be of type Histogram")

    @staticmethod
    def check_h_axes(
        hpid: hist.Histogram,
        config_path: str = "config/main.yml",
        config_key: str = "pid",
    ) -> None:
        """Check the axes of the hpid histogram are compatible with the user-defined config

        Parameters
        ----------
        hpid: hist.Histogram
            PID eff histogram
        config_path: str
            The path to the user-defined config file
        config_key: str
            The key in the config file to the binning config

        Returns
        -------
        None
            Checks compatibility of the axes of the hpid histogram with the user-defined config
        """
        # open custom binning config
        binning = read_config(config_path, key=config_key)["binning"]

        # sanity check: same number of axes
        assert len(binning.keys()) == len(
            hpid.axes
        ), "The binning dimensions in the config does not match the number of axes in the histogram"

        _axes = hpid.axes
        for ax in _axes:
            if ax.name in binning.keys():
                if ax.edges != binning[ax.name]:
                    raise ValueError(
                        f"The edges of the {ax.name} axis do not match the config"
                    )
            else:
                raise ValueError(f"The {ax.name} axis is not in the config")


class SuperHistFactory:
    """Book the super eff container, without filling ir"""

    def __init__(
        self,
        config_path: str = "config/main.yml",
    ) -> None:
        self._config_path = config_path

    @property
    def config_path(self) -> str:
        return self._config_path

    @staticmethod
    def book_binning_axes(binning_dict: dict[str, list[float]]) -> list | tuple:
        """Match binning vars to hist axes

        Parameters
        ----------
        binning: dict[str, list[float]]
            The binning config

        Returns
        -------
        list | tuple
            The axes of the histogram
        """
        axes = []
        match binvar:
            case [*_, "P"]:
                axes.append(
                    hist.axis.Variable(
                        binning_dict[binvar], name="P", label=r"$p$ [MeV$/c$]"
                    )
                )
            case [*_, "PT"]:
                axes.append(
                    hist.axis.Variable(
                        binning_dict[binvar], name="PT", label=r"$p_T$ [MeV$/c$]"
                    )
                )
            case [*_, "ETA"]:
                axes.append(
                    hist.axis.Variable(
                        binning_dict[binvar], name="ETA", label=r"$\eta$"
                    )
                )
            case [*_, "nTracks"]:
                axes.append(
                    hist.axis.Variable(
                        binning_dict[binvar], name="nTracks", label=r"\texttt{nTracks}"
                    )
                )
            case [*_, "nSPDhits"]:
                axes.append(
                    hist.axis.Variable(
                        binning_dict[binvar], name="nSPDhits", label=r"\texttt{nSPD}"
                    )
                )

        # sanity check
        assert len(axes) == len(
            binning_dict.keys()
        ), "The number of axes does not match the number of binning variables"

        return axes

    def make(self, config_key: str = "pid") -> hist.Hist:
        """Generate the super efficiency histogram

        Parameters
        ----------
        config_path: str
            The path to the user-defined config file [default: "config/main.yml"]
        config_key: str
            The key in the config file to the binning config [default: "pid"]

        Returns
        -------
        hist.Hist
            The super efficiency histogram
        """

        # based on the config, get calib species & reco partitions axes
        config = read_config(self.config_path, key=config_key)
        binning = config["binning"]
        calib_species = config["species"].keys()
        reco_partitions = config["reco_cuts"].keys()

        # initialise the super efficiency histogram, unpacking the binning axes
        super_hpid = hist.Hist(
            hist.axis.StrCategory(
                [
                    "test_calib",
                ],
                name="calib",
                label="True Category",
            ),
            hist.axis.StrCategory(
                [
                    "test_reco",
                ],
                name="reco",
                label="Reco Category",
            ),
            *self.book_binning_axes(binning),
            storage=hist.storage.Weight(),
        )
