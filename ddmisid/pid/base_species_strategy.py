"""
Strategy pattern base class to specifically retrieve the correct species name and alias for the PID efficiency jobs, and other relevant configuration parameters.
Specialising the strategy to each species, inheriting from a common base class, allows for a bespoke treatment (e.g., electrons and ghosts).
This is crucial for general application within DDmisID and its potential application within LHCb.
"""

from abc import ABC, abstractmethod
from ddmisid.engine import config


class SpeciesValidatorMixin:
    """Mixin to validate species provided in the YAML configuration file."""

    def validate_species(self):
        """
        Validate that the species specified in the YAML config file are compliant with the PIDCalib2 aliases.
        Raises an error if the species is not found in the config file.
        """
        if self.species not in self.species_alias_map.keys():
            raise ValueError(
                f"Species '{self.species}' not found in the configuration file or alias map."
            )


class SpeciesStrategyBase(ABC):
    """Abstract base class for species strategies."""

    def __init__(self, config):
        self._species_alias_map = config.pid.species

    @property
    def species_alias_map(self):
        """Return the map of species aliases from the configuration."""
        return self._species_alias_map

    @abstractmethod
    def get_species_name(self) -> str:
        """
        Return the internal species identifier within DDmisID, as specified by the user in the config YAML file (i.e., the keys of the `species` field in the YAML config file).
        This identifier is used within DDmisID for species-specific operations except direct PIDCalib2 queries.
        """
        pass

    @abstractmethod
    def get_species_alias(self) -> str:
        """
        Fetch the PIDCalib2-compliant alias for each species from the user-specified YAML files (i.e., the values of the `species` field in the YAML config file).
        This alias is exclusively used for PIDCalib2 queries.
        """
        pass


class ParticleStrategy(SpeciesStrategyBase, SpeciesValidatorMixin):
    """
    Strategy for particle species (e.g., kaon, pion, proton, electron) that rely on PIDCalib2 to extract PID efficiencies.
    This strategy fetches aliases from the config YAML file and validates them against PIDCalib2 requirements.
    """

    def __init__(self, config, species):
        super().__init__(config)
        self._species = species
        self.validate_species()

    @property
    def species(self):
        """Return the species specified."""
        return self._species

    def get_species_name(self) -> str:
        """Return the DDmisID internal species identifier."""
        return self.species

    def get_species_alias(self) -> str:
        """Return the PIDCalib2-compliant species alias from the config."""
        return self.species_alias_map[self.species]
