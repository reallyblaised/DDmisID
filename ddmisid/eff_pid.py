"""Generate multidimensional efficiency histograms accounting for:

- calibration species (true hadrons, ghosts, electrons) [`calib_species`]
- reco partictions of the hadron-enriched data [`reco_partitions`]
- efficiencies, obtained via pidcalib2 via calibration samples, plus ghosts 
"""

__author__ = "Blaise Delaney"
__email__  = "blaise.delaney at cern.ch"

import pickle
import hist
from ddmisid.utils import read_config

class PIDEffMixin: 
    """Mixin for sanity checks for PID eff hists.
    
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
        config_path: str = "config/main.yml"   
    ) -> None:
        """Check the axes of the hpid histogram are compatible with the user-defined config
        
        Parameters
        ----------
        hpid: hist.Histogram
            PID eff histogram
        config_path: str
            The path to the user-defined config file
        
        Returns
        -------
        None
            Checks compatibility of the axes of the hpid histogram with the user-defined config
        """
        # open custom binning config
        binning = read_config(config_path, key="pid")["binning"]
        
        # sanity check: same number of axes
        assert len(binning.keys()) == len(hpid.axes), "The binning dimensions in the config does not match the number of axes in the histogram"
        
        _axes = hpid.axes
        for ax in _axes:
            if ax.name in binning.keys():
                if ax.edges != binning[ax.name]:
                    raise ValueError(f"The edges of the {ax.name} axis do not match the config")
            else:
                raise ValueError(f"The {ax.name} axis is not in the config")        


class PIDEffExtractor(PIDEffMixin): # inherit sanity checks from mixin
    """Open the PID eff pkl file and generate the super eff container"""
    def __init__(
        self, 
        hist_path: str,
    ) -> None:
        self._hist_path = hist_path
        
    @property
    def hist_path(self) -> str:
        return self._hist_path
    
    def super_hist(self) -> hist.Hist:
        """Generate the super efficiency histogram"""
        pass
        
        # implement opening, checking, reading in calib and reco specie, and assembing this into a super hist
        