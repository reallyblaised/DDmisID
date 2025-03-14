"""Config model validation of the user-defined DDmisID configuration."""

from pydantic import BaseModel, Field, validator
from typing import List, Dict, Optional
from pydantic import BaseModel, Field, validator
from typing import List, Dict, Optional, Union


class PIDConfig(BaseModel):
    """PIDCalib2 configuration model."""

    sweight_binning: Dict[str, List[float]] = Field(
        ...,
        description="Coarser binning for the BML fits to extract particle-species sWeights",
    )
    pid_extrap_binning: Dict[str, List[float]] = Field(
        ...,
        description="Binning in unfolding eff(!mu) and folding in eff(mu), from PIDCalib2",
    )
    species: Dict[str, str] = Field(
        ...,
        description="User-defined species alias for use within DDmisID : PIDCalib2 species name (e.g., electron: e_B_Jpsi)",
    )
    year: str = Field(..., description="Data-taking years")
    magpol: str = Field(..., description="Magnet polarities")
    control: str = Field(..., description="PID selection for control region")
    target: str = Field(..., description="PID selection for target region")
    common_selection: Optional[str] = Field(
        None,
        description="Common cuts between the control and target region (e.g., InMuonAcc==1.0)",
    )
    reco_partitions: Dict[str, str] = Field(
        ...,
        description="DLL and ProbNN cuts defining the reco partitions within the control region",
    )
    ghost_config: Optional[Dict[str, str]] = Field(
        ..., description="Configuration for the ghost candidates"
    )

    # year either of Run 1 or Run 2
    @validator("year")
    def validate_year(cls, value):
        allowed_years = ("2011", "2012", "2015", "2016", "2017", "2018")
        if value not in allowed_years:
            raise ValueError(
                f"Invalid year of data taking: {value}. Allowed values are '2011', '2012', '2015', '2016', '2017', '2018'."
            )
        return value

    # magpol in ["up", "down"]
    @validator("magpol")
    def validate_magpol(cls, value):
        allowed_magpols = ("up", "down")
        if value not in allowed_magpols:
            raise ValueError(
                f"Invalid magnet polarity: {value}. Allowed values are 'up' or 'down'."
            )
        return value

    # full ghost MC spec validation, if provided [Optional]
    @validator("ghost_config", always=True)
    def validate_ghost_config(cls, value):
        if value is None:
            return value  # No validation needed if ghost_config is None
        required_keys = {"path", "key", "tree", "hadron_enriched_sel", "branch_prefix"}
        missing_keys = required_keys - value.keys()
        if missing_keys:
            raise ValueError(f"ghost_config is missing required keys: {missing_keys}")
        return value


class DataConfig(BaseModel):
    """Data configuration for DDmisID."""

    input_path: str = Field(..., description="Path to the input data file")
    data_key: Optional[str] = Field(
        None, description="Key for the data tree in the input file"
    )
    data_tree: str = Field(..., description="Name of the data tree in the input file")
    data_prefixes: Dict[str, str] = Field(
        ..., description="Variable: alias prefix in the data tree"
    )
    data_reco_partitions: Dict[str, str] = Field(
        ...,
        description="Reconstruction partitions, including cuts defining the control region",
    )
    output_path: str = Field(..., description="Path to the output data file")


class DDmisIDConfig(BaseModel):
    """Main configuration model for DDmisID."""

    user_id: str = Field(..., description="CERN user ID")
    max_calib_files: int = Field(
        ...,
        description="Must be a non-negative integer specifying the number of calibration files, or -1 for full set of calibration files",
    )
    verbose: bool = Field(..., description="Whether to produce verbose output")
    pid: PIDConfig = Field(..., description="PID configuration")
    data: DataConfig = Field(..., description="Data configuration")

    class Config:
        frozen = False  # ensure immutability of the configuration post-validation

    # non-negative calib files read in from eos
    @validator("max_calib_files")
    def validate_max_calib_files(cls, v):
        if v < 0 and v != -1:  # allow -1 to signal full suite of calibration files
            raise ValueError("max_calib_files must be non-negative")
        return v

    # user_id cannot be empty - must be a valid CERN user ID for Kerberos authentication to access EOS
    @validator("user_id")
    def validate_user_id(cls, v):
        if not v:
            raise ValueError("user_id cannot be empty")
        return v
