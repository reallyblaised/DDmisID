"""Factory interface to engage PIDCalib2 job generation depending on the species specified in the YAML configuration file."""

from ddmisid.pidcalib.job_generator import ParticleJobGenerator, GhostJobGenerator
from ddmisid.pidcalib.species_strategy import (
    KaonStrategy,
    PionStrategy,
    ProtonStrategy,
    ElectronStrategy,
    GhostStrategy,
)
from ddmisid.engine import config
from loguru import logger

class PIDCalibJobFactory:
    """
    Factory interface to engage PIDCalib2 job generation depending on the species specified in the YAML configuration file.
    This defines the job creation strategy for DDmisID whilst allowing the user to specify a (sub)set of species ({π,p,K,e,g}) to run over.
    """

    def __init__(self):
        self.years = config.pid.years
        self.magpols = config.pid.magpols
        self.species = config.pid.species
        self.max_calib_files = config.max_calib_files  # Max files to process (-1 for all, else user-defined)
        self._validate_species()

    def generate_jobs(
        self, output_dir: str, region_ids: list = ["control", "target"], verbose: bool = False
    ) -> None:
        """Factory interface to spawn all the requisite PIDCalib2 jobs."""
        for species_id, species_alias in self.species.items():  # e.g., {"electron": "e_B_Jpsi"}
            logger.info(f"Generating PIDCalib2 jobs for species: {species_id}")

            # Fetch the appropriate strategy and job generator class for the species
            strategy, generator_class = self._get_strategy_and_generator(species_id)
            job_generator = generator_class(config, strategy)

            # Dynamically generate the jobs, looping through the years and magnetic polarities
            for year in self.years:
                for magpol in self.magpols:
                    logger.info(f"  -> Generating PIDCalib2 jobs for year {year} with polarity {magpol}")
                    job_generator.generate_jobs(
                        year=year,
                        magpol=magpol,
                        region_id=region_ids,
                        output_dir=output_dir,
                        verbose=verbose,
                        max_calib_files=self.max_calib_files  # Pass the correct max file limit
                    )

    def _get_strategy_and_generator(self, species_id: str) -> tuple:
        """Fetch the appropriate strategy and generator class for the species."""
        match species_id:
            case "pion": return PionStrategy(), ParticleJobGenerator
            case "proton": return ProtonStrategy(), ParticleJobGenerator
            case "kaon": return KaonStrategy(), ParticleJobGenerator
            case "electron": return ElectronStrategy(), ParticleJobGenerator
            case "ghost": return GhostStrategy(), GhostJobGenerator
            case _: raise ValueError(f"Species {species_id} PIDCalib2 strategy missing or not recognised. Allowed values: ['pion', 'proton', 'kaon', 'electron', 'ghost']")

    def _validate_species(self):
        """Ensure species in the configuration file are valid."""
        allowed_species = {"kaon", "pion", "proton", "electron", "ghost"}
        for species_id in self.species.keys():
            if species_id not in allowed_species:
                raise ValueError(f"Invalid species: {species_id}. Must be one of {allowed_species}.")