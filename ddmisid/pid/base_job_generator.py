"""
Template pattern for generating PID efficiency extraction jobs jobs.
This differentiates between paticle (PIDCalib2) and ghost (MC) job PID efficiency extraction jobs.
"""

from abc import ABC, abstractmethod
from ddmisid.engine import config


class BaseJobGenerator(ABC):
    """Abstract base class for PID-efficiency job creation, accounting ofr particle-species strategy, magpol, year, and PID region."""

    def __init__(self, config):
        self._sweight_binning = config.pid.sweight_binning
        self._pid_extrap_binning = config.pid.pid_extrap_binning
        self._control_selection = config.pid.control
        self._target_selection = config.pid.target
        self._common_selection = config.pid.common_selection
        self._reco_partitions = config.pid.reco_partitions

    @property
    def sweight_binning(self):
        return self._sweight_binning

    @property
    def pid_extrap_binning(self):
        return self._pid_extrap_binning

    @property
    def control_region(self):
        return self._control_selection

    @property
    def target_region(self):
        return self._target_selection

    @property
    def common_selection(self):
        return self._common_selection

    @property
    def reco_partitions(self):
        return self._reco_partitions

    @abstractmethod
    def generate_jobs(
        self,
        year: str,
        magpol: str,
        region_id: list,
        outdir_dir: str,
        verbose: bool,
        test: bool,
    ) -> None:
        pass
