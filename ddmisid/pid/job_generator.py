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


class JobWriterMixin:
    """Mixin to handle the writing and permissions of bash scripts."""

    def _construct_job_command(
        self,
        calib_sample: str,
        magpol: str,
        species_pidcalib_alias: str,
        pid_cut: str,
        binning_vars: list,
        binning_path: Path,
        common_selection: str,
        max_calib_files: int,
        output_dir: Path,
    ) -> str:
        """Construct the bash command for PID-efficiency extraction job execution."""
        # Build the basic job command
        job_conf = (
            f"source /cvmfs/lhcb.cern.ch/lib/LbEnv &&\n"
            f"lb-conda pidcalib pidcalib2.make_eff_hists --sample {calib_sample} "
            f"--magnet {magpol} --particle {species_pidcalib_alias} --pid-cut '{pid_cut}' "
            f"--cut '{common_selection}' --binning-file {binning_path} --output-dir {output_dir}"
        )

        # Add binning variables
        for binning_var in binning_vars:
            job_conf += f" --bin-var {binning_var}"

        # Apply max file limit if specified [default is -1 for all files]
        if max_calib_files > 0:
            job_conf += f" --max-files {max_calib_files}"

        # Add verbose flag if specified
        if self.verbose:
            job_conf += " --verbose"

        # Final step: touch a done file after completion
        job_conf += (
            f" &&\ntouch {output_dir}/pidcalib2.make_eff_hists.done\n"
            f"for f in {output_dir}/*.pkl; do\n"
            f'mv "$f" {output_dir}/perf.pkl\ndone\n'
        )

        return job_conf

    def _write_job_script(self, script_path: Path, job_conf: str) -> None:
        """Write the bash script to the given path and make it executable."""
        with open(script_path, "w") as f:
            f.write(job_conf)

        # Make the script executable
        os.chmod(script_path, 0o755)
        logger.info(f"Job script written and made executable: {script_path}")


class JobSetter(BaseJobGenerator, JobWriterMixin):
    """Abstract base class for both Particle and Ghost job generation."""

    def __init__(self, config, strategy):
        super().__init__(config)
        self.strategy = strategy
        self._max_calib_files = config.max_calib_files

    @property
    def max_calib_files(self):
        return self._max_calib_files

    def _setup_binning(
        self,
        species: str,
        species_alias: str,
        year: str,
        output_dir: str,
        binning_dict: dict[str, list[float]],
        binning_alias: str,
    ):
        """Helper method to setup binning for the jobs."""
        binning = DefaultBinningGenerator(
            species, species_alias, binning_dict, binning_alias
        )
        binning_vars = binning.get_binning_variables(year)
        binning_path = binning.generate_output_filepath(species, year, output_dir)
        binning.build(year, output_dir)
        return binning_vars, binning_path

    def _generate_job_script_paths(
        self,
        species: str,
        pid_category: str,
        year: str,
        magpol: str,
        output_dir: str,
        partition_label: str | None = None,
    ) -> Path:
        """Helper to generate the script path for the given parameters."""
        # establish a directory tree refleting the PID maps to be derived
        if partition_label:
            namespace = f"{species}_to_{partition_label}"
        else:
            namespace = f"{species}_to_{pid_category}_like"

        # compose the path to the executable script
        script_path = Path(
            f"{output_dir}/{year}/{magpol}/{pid_category}/{species}/{namespace}/run.sh"
        )
        script_path.parent.mkdir(parents=True, exist_ok=True)
        return script_path

    def _fetch_calib_and_mctuning(self, species: str, year: str) -> tuple:
        """Helper to fetch calibration sample and MC tuning."""
        calib_sample = CalibSamples().fetch(year=year, species=species)
        mc_tuning = MCTunings().fetch(year=year)
        return calib_sample, mc_tuning


class ParticleJobGenerator(JobSetter):
    """
    concrete implementation of job script generator to extract efficiencies related to the pid selections defining control and target.
    specific to particles, hence the compliance with pidcalib2 directives.
    """

    def generate_jobs(
        self,
        year: str,
        magpol: str,
        region_id: list,
        output_dir: str = "bin",
    ) -> None:
        """Generate PIDCalib2 jobs for particles."""
        species = self.strategy.get_species_name()
        species_alias = self.strategy.get_species_alias()

        # Fetch calibration sample and MC tuning
        calib_sample, mc_tuning = self._fetch_calib_and_mctuning(species, year)

        # Iterate over regions (control/target) to generate jobs
        for region in region_id:

            # Set PID cuts based on the region
            pid_cut = (
                f"{mc_tuning}{self.control_pid_selection}"
                if region == "control"
                else f"{mc_tuning}{self.target_pid_selection}"
            )

            # Setup binning variables
            binning_vars, binning_path = self._setup_binning(
                species,
                species_alias,
                year,
                output_dir,
                self.pid_extrap_binning,
                "extrapolation_binning",
            )

            # Construct and write the executable jobsbash script
            script_path = self._generate_job_script_paths(
                species, region, year, magpol, output_dir
            )
            job_conf = self._construct_job_command(
                calib_sample=calib_sample,
                magpol=magpol,
                species_pidcalib_alias=species_alias,
                pid_cut=pid_cut,
                binning_vars=binning_vars,
                binning_path=binning_path,
                common_selection=self.common_selection,
                max_calib_files=self.max_calib_files,
                output_dir=script_path.parent,
            )
            self._write_job_script(script_path, job_conf)


class ParticleRecoPartitionJobGenerator(JobSetter):
    """
    Cconcrete implementation of job script generator to extract efficiencies related to the reco partitions within the control sample, eff(reco partition | control PID, kinematics, topology, occupancy).
    specific to particles, hence the compliance with pidcalib2 directives.
    """

    def generate_jobs(
        self,
        year: str,
        magpol: str,
        output_dir: str = "bin",
    ) -> None:
        """Generate PIDCalib2 jobs with reconstruction partitions."""
        species = self.strategy.get_species_name()
        species_alias = self.strategy.get_species_alias()

        # Fetch calibration sample and MC tuning
        calib_sample, mc_tuning = self._fetch_calib_and_mctuning(species, year)

        # Iterate over reco partitions
        for partition_label, partition_pid_criteria in self.reco_partitions.items():

            # Set PID cut for each partition
            pid_cut = partition_pid_criteria
            hadron_enriched_selection = (
                f"{mc_tuning}{self.control_pid_selection} & {self.common_selection}"
            )

            # Setup binning variables
            binning_vars, binning_path = self._setup_binning(
                species,
                species_alias,
                year,
                output_dir,
                self.sweight_binning,
                "reco_partition_binning",
            )

            # Construct and write the executable jobsbash script
            script_path = self._generate_job_script_paths(
                species, "control", year, magpol, output_dir, partition_label
            )
            job_conf = self._construct_job_command(
                calib_sample=calib_sample,
                magpol=magpol,
                species_pidcalib_alias=species_alias,
                pid_cut=pid_cut,
                binning_vars=binning_vars,
                binning_path=binning_path,
                common_selection=hadron_enriched_selection,
                max_calib_files=self.max_calib_files,
                output_dir=script_path.parent,
            )
            self._write_job_script(script_path, job_conf)


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


class GhostRecoPartitionJobGenerator(BaseJobGenerator, JobWriterMixin):
    """
    Ghost-specific generator for jobs to extract eff(partittion|hadron-enriched sample).
    This too relies on suitably truth-matched user-specified simulation samples.
    """

    pass
