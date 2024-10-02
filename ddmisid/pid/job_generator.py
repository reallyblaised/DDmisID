"""
Specialised job creator classes setting up query protocols to PIDCalib2 for particle species,
and a bespoke MC-based protocol for ghost efficiencies.
"""

from pathlib import Path
from ddmisid.pid.species_strategy import ParticleStrategy, GhostStrategy
from ddmisid.pid.base_job_generator import BaseJobGenerator
from loguru import logger
import os
from config.calib_tuning import CalibSamples, MCTunings
from ddmisid.utils.binning import DefaultBinningGenerator


class ParticleJobGenerator(BaseJobGenerator):
    """
    Concrete implementation for generating PIDCalib2 jobs for particles.
    Inherits from the BaseJobCreator to access the relevant config parameters.
    """

    def __init__(self, config, strategy: ParticleStrategy):
        super().__init__(config)
        self.strategy = strategy
        self.max_calib_files = config.max_calib_files  # Picked up from the config

    def generate_jobs(
        self,
        year: str,
        magpol: str,
        region_id: list,  # PIDCalibJobFactory sets this by default to ["control", "target"]
        verbose: bool,
        output_dir: str = "bin",  # Indicating where bash files will be written to
    ) -> None:
        """
        Generate the PIDCalib2 jobs for the given year, magnetic polarity, and regions.

        Parameters
        ----------
        year : str
            The year of data taking (e.g., "2018").
        magpol : str
            The magnetic polarity (allowed values are "up" and "down").
        region_id : list
            List of regions to apply (allowed values are "control" and "target").
        verbose : bool
            If True, include verbose output.
        output_dir : str
            Directory where the output will be saved [by default set to `bin`].
        """
        species = self.strategy.get_species_name()  # e.g., "kaon", "pion", etc.
        alias = self.strategy.get_species_alias()  # e.g., "K", "Pi", etc.

        # Fetch external inputs: calibration sample and MCTuning version
        calib_sample = CalibSamples().fetch(year=year, species=species)
        mc_tuning = MCTunings().fetch(year=year)  # Prefix for the ProbNNghost selection

        # Iterate through control and target regions to dynamically create PIDCalib2 jobs
        for region in region_id:

            # Establish PIDCalib-specific alias for the species
            species_pidcalib_alias = self.strategy.get_species_alias()

            # Establish the path for the executable (run.sh)
            namespace = f"{species}_to_{region}_like"
            script_path = Path(
                f"{output_dir}/{year}/{magpol}/{region}/{species}/{namespace}/run.sh"
            )
            script_path.parent.mkdir(parents=True, exist_ok=True)

            # Set the appropriate PID cuts for each region
            pid_cut = (
                f"{mc_tuning}{self.control_pid_selection}"
                if region == "control"
                else f"{mc_tuning}{self.target_pid_selection}"
            )

            # Establish the binning coordinates
            Binning = DefaultBinningGenerator(species, alias, self.pid_extrap_binning)
            binning_vars = Binning.get_binning_variables(
                year
            )  # List of binning variables
            binning_path = Binning.generate_output_filepath(
                species, year, output_dir
            )  # Path to binning JSON
            Binning.build(year, output_dir, verbose)  # Write the binning JSON file

            # Write the job executable
            self._write_job_executable(
                script_path,
                calib_sample,
                magpol,
                species_pidcalib_alias,
                pid_cut,
                binning_vars,
                binning_path,
            )

    def _write_job_executable(
        self,
        script_path: Path,
        calib_sample: str,
        magpol: str,
        species_pidcalib_alias: str,
        pid_cut: str,
        binning_vars: list,
        binning_path: Path,
    ) -> None:
        """
        Write the job executable (run.sh) for PIDCalib2.

        Parameters
        ----------
        script_path : Path
            The path to the script file (run.sh).
        calib_sample : str
            The calibration sample (e.g., 'Turbo18').
        magpol : str
            The magnetic polarity (e.g., 'up' or 'down').
        species_pidcalib_alias : str
            The species alias for PIDCalib (e.g., 'Pi', 'K').
        pid_cut : str
            The PID cut applied to the numerator in the efficiency ratio.
        binning_vars : list
            List of binning variables to pass to PIDCalib2.
        binning_path : Path
            The path to the binning JSON file.
        """
        # Construct the job command
        job_conf = (
            f"source /cvmfs/lhcb.cern.ch/lib/LbEnv &&\n"
            f"lb-conda pidcalib pidcalib2.make_eff_hists --sample {calib_sample} "
            f"--magnet {magpol} --particle {species_pidcalib_alias} --pid-cut '{pid_cut}' "
            f"--cut '{self.common_selection}' --binning-file {binning_path} --output-dir {script_path.parent}"
        )

        # Add binning variables
        for binning_var in binning_vars:
            job_conf += f" --bin-var {binning_var}"

        # Apply the max file limit based on the config
        if self.max_calib_files > 0:
            job_conf += f" --max-files {self.max_calib_files}"

        # Final step: Touch a done file after completion
        job_conf += (
            f" &&\ntouch {script_path.parent}/pidcalib2.make_eff_hists.done\n"
            f"for f in {script_path.parent}/*.pkl; do\n"
            f'mv "$f" {script_path.parent}/perf.pkl\ndone\n'
        )

        # Write the job script to the file
        with open(script_path, "w") as f:
            f.write(job_conf)

        logger.info(f"Job script written to {script_path}")


class GhostJobGenerator(BaseJobGenerator):
    """
    Concrete implementation for generating PIDCalib2 jobs for ghosts.
    Inherits from the BaseJobCreator to access the relevant config parameters.
    """

    def __init__(self, config, strategy: GhostStrategy):
        super().__init__(config)
        self.strategy = strategy
        self.max_calib_files = config.max_calib_files  # Picked up from the config

    def generate_control_target_jobs(self):
        pass
