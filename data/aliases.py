"""Alias dictionaries required by pidcalib2.

Note: binning vars, calib samples, selection cuts and tunings
are defined as dataclass attributes post-init, ie mutable.

__authors__ = Blaise Delaney
__emails__  = blaise.delaney at cern.ch
"""
from dataclasses import dataclass
from dataclasses import field
from ddmisid.utils import read_config

_config_path = "config/main.yml"

def update_selstr_run2(sel: str, year: str) -> str:
    """Update a nominal selection string with appropriate aliases"""
    
    match year:
        case "_2016" | "_2017" | "_2018":
            # would love to use match case here
            if "PT" in sel:
                sel = sel.replace("PT", "Brunel_PT")
            if "P" in sel:
                sel = sel.replace("P", "Brunel_P")
            if "nTracks" in sel:   
                sel = sel.replace("nTracks", "nTracks_Brunel")
        
        case "_2011" | "_2012" | "_2015":
            pass
        
        case other:
            raise ValueError("Year identifier not recognised")  

    return sel


@dataclass(frozen=True, slots=True)
class CalibSamples:
    """Container dataclass for the calibration samples of interest

    Attributes
    ----------
    {hadron/e}_{year}: str
        The calibration sample for the given year, diverging between electrons and hadrons
    """

    # hadrons, from calibration samples
    hadron_2018: str = "Turbo18"
    hadron_2017: str = "Turbo17"
    hadron_2016: str = "Turbo16"
    hadron_2015: str = "Turbo15"
    hadron_2012: str = "21"
    hadron_2011: str = "21r1"

    # electrons, from calibration samples
    e_2018: str = "Electron18"
    e_2017: str = "Electron17"
    e_2016: str = "Electron16"
    e_2015: str = "Electron15"
    e_2012: str = "21"
    e_2011: str = "21r1"


@dataclass(frozen=True, slots=True)
class MCTunings:
    """Container dataclass for the MC tunings

    Attributes
    ----------
    {year_conditions}: str
    """

    _2018: str = "MC15TuneV1"
    _2017: str = "MC15TuneV1"
    _2016: str = "MC15TuneV1"
    _2015: str = "MC15TuneV1"
    _2012: str = "MC12TuneV2"
    _2011: str = "MC12TuneV2"


@dataclass(frozen=True, slots=True)
class BinningVars:
    """Container for binning vars
   
    Attributes
    ----------
        {year_conditions}: str
    
    Methods
    -------
    build
        Builds the pidcalib2 binning dictionary
    """
    _2011: str = field(init=False)
    _2012: str = field(init=False)
    _2015: str = field(init=False)
    _2016: str = field(init=False)
    _2017: str = field(init=False)
    _2018: str = field(init=False)


    @staticmethod
    def alias_upater_run2(varlist: list[str,...], year: str) -> list[str,...]:
        """Update the list of variables with appropriate aliases"""
        
        match year:
            case "_2016" | "_2017" | "_2018":
                # would love to use match case here
                if "PT" in varlist:
                        varlist = list(map(lambda x: x.replace("PT", "Brunel_PT"), varlist))    
                if "P" in varlist:
                        varlist = list(map(lambda x: x.replace("P", "Brunel_P"), varlist))    
                if "nTracks" in varlist:   
                        varlist = list(map(lambda x: x.replace("nTracks", "nTracks_Brunel"), varlist))    
            
            case "_2011" | "_2012" | "_2015":
                pass
            
            case other:
                raise ValueError("Year identified in BinningVars not recognised")       
        
        return varlist 


    def build(self) -> None:
        """Populate the dictionary pidcalib binning vars
        appropriately labelled
        """
        years = ("_2011", "_2012", "_2015", "_2016", "_2017", "_2018")

        # read the user config options
        bin_vars = list(read_config(_config_path, key="pid")["binning"].keys())
       
        # if require, add pre-/suf-fix to variable aliases
        for y in years:
            self.alias_upater_run2(varlist=bin_vars, year=y)
        
            # post alias-update, set as attributed of the class
            object.__setattr__(self, y, bin_vars)

    def __post_init__(self) -> None:
        self.build()
        
        
@dataclass(frozen=True, slots=True)
class CommonCuts:
    """Container dataclass for evt sel cuts

    Attributes
    ----------
    year_conditions: dict[str, str]

    Methods
    -------
    build
        Builds the dictionary with PID, kinematic and occupancy criteria
        imposed in the evt selection
    """
    _2011: str = field(init=False)
    _2012: str = field(init=False)
    _2015: str = field(init=False)
    _2016: str = field(init=False)
    _2017: str = field(init=False)
    _2018: str = field(init=False)


    def __post_init__(self) -> None:
        """Populate the per-year evt sel cuts
        labelled appropriately depending on the year
        """ 
        # read the user config options for the acceptance and kinematics
        # cuts shared between HE and signal samples
        common_sel = read_config(_config_path, key="pid")["common_sel"] 
        
        years = ("_2011", "_2012", "_2015", "_2016", "_2017", "_2018")
        for y in years:
            # establish the common cuts for the given year & assign as class attribute
            object.__setattr__(self, y, update_selstr_run2(sel=common_sel, year=f"{y}"))