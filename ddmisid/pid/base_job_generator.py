"""
Template pattern for generating PID efficiency extraction jobs.
This differentiates between particle (PIDCalib2) and ghost (MC) job PID efficiency extraction jobs.
"""

from abc import ABC, abstractmethod


class BaseJobGenerator(ABC):
    """Abstract base class for PID-efficiency job creation, accounting for particle-species strategy, magpol, year, and PID region."""

    def __init__(self, config):
        self._sweight_binning = config.pid.sweight_binning
        self._pid_extrap_binning = config.pid.pid_extrap_binning
        self._reco_partitions = config.pid.reco_partitions
        self._control_pid_selection = config.pid.control
        self._target_pid_selection = config.pid.target
        self._common_selection = config.pid.common_selection
        self._verbose = config.verbose

    @property
    def sweight_binning(self):
        return self._sweight_binning

    @property
    def pid_extrap_binning(self):
        return self._pid_extrap_binning

    @property
    def control_pid_selection(self):
        return self._control_pid_selection

    @property
    def target_pid_selection(self):
        return self._target_pid_selection

    @property
    def common_selection(self):
        return self._common_selection

    @property
    def reco_partitions(self):
        return self._reco_partitions

    @property
    def verbose(self):
        return self._verbose

    @abstractmethod
    def generate_jobs(
        self,
        year: str,
        magpol: str,
        region_id: list,
        outdir_dir: str,
    ) -> None:
        pass
