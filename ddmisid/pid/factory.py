"""Factory interface to engage PIDCalib2 job generation depending on the species specified in the YAML configuration file."""

from ddmisid.pid.job_generator import ParticleJobGenerator, GhostJobGenerator, ParticleRecoPartitionJobGenerator, GhostRecoPartitionJobGenerator
from ddmisid.pid.species_strategy import (
    KaonStrategy,
    PionStrategy,
    ProtonStrategy,
    ElectronStrategy,
    GhostStrategy,
)
from ddmisid.engine import config
from loguru import logger

class PIDEffXJobFactory:
    """
    Factory interface to engage PIDCalib2 job generation depending on the species specified in the YAML configuration file.
    This defines the job creation strategy for DDmisID whilst allowing the user to specify a (sub)set of species ({π,p,K,e,g}) to run over.
    """

    def __init__(self):
        self.years = config.pid.years
        self.magpols = config.pid.magpols
        self.species = config.pid.species
        self._validate_species()

    def generate_jobs(
        self, output_dir: str, region_ids: list = ["control", "target"], verbose: bool = False
    ) -> None:
        """Factory interface to spawn all the requisite PIDCalib2 jobs."""
        for species_id, species_alias in self.species.items():  # e.g., {"electron": "e_B_Jpsi"}
            logger.info(f"Generating PIDCalib2 jobs for species: {species_id}")

            # Fetch the appropriate strategy and job generators for each species (control, target, reco partition | control)
            strategy, control_target_job_generator, reco_partition_job_generator = self._get_strategy_and_generator(species_id) 
            
            # Dynamically generate the jobs, looping through the years and magnetic polarities
            for year in self.years:
                for magpol in self.magpols:
                    # first, control and target pid-efficiency-extra jobs
                    logger.info(f"  -> Generating PIDCalib2 jobs for year {year} with polarity {magpol} targeting {region_ids}.")
                    control_target_job_generator(config, strategy).generate_jobs(
                        year=year,
                        magpol=magpol,
                        region_id=region_ids,
                        output_dir=output_dir,
                    )
                    # second, reco partition pid-efficiency-extra jobs
                    logger.info(f"  -> Generating PIDCalib2 jobs for year {year} with polarity {magpol} targeting reco partition.")
                    reco_partition_job_generator(config, strategy).generate_jobs(
                        year=year,
                        magpol=magpol,
                        output_dir=output_dir,
                    )

    def _get_strategy_and_generator(self, species: str) -> tuple:
        """Fetch the appropriate strategy and generator class for the species.
        species identifier; pid eff for control and target; pid eff for reco partitions | control selection.
        """
        match species:
            case "pion": return PionStrategy(), ParticleJobGenerator, ParticleRecoPartitionJobGenerator
            case "proton": return ProtonStrategy(), ParticleJobGenerator, ParticleRecoPartitionJobGenerator
            case "kaon": return KaonStrategy(), ParticleJobGenerator, ParticleRecoPartitionJobGenerator
            case "electron": return ElectronStrategy(), ParticleJobGenerator, ParticleRecoPartitionJobGenerator
            case "ghost": return GhostStrategy(), GhostJobGenerator, GhostRecoPartitionJobGenerator
            case _: raise ValueError(f"Species {species} PIDCalib2 strategy missing or not recognised. Allowed values: ['pion', 'proton', 'kaon', 'electron', 'ghost']")

    def _validate_species(self):
        """Ensure species in the configuration file are valid."""
        allowed_species = {"kaon", "pion", "proton", "electron", "ghost"}
        for species_id in self.species.keys():
            if species_id not in allowed_species:
                raise ValueError(f"Invalid species: {species_id}. Must be one of {allowed_species}.")