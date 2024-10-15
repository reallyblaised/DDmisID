"""Species-specific strategies for PID-efficency extraction jobs."""

from ddmisid.pid.base_species_strategy import BaseSpeciesStrategy, ParticleStrategy
from ddmisid.engine import config


class PionStrategy(ParticleStrategy):
    def __init__(self, config=config):
        super().__init__(config=config, species="pion")


class KaonStrategy(ParticleStrategy):
    def __init__(self, config=config):
        super().__init__(config=config, species="kaon")


class ProtonStrategy(ParticleStrategy):
    def __init__(self, config=config):
        super().__init__(config=config, species="proton")


class ElectronStrategy(ParticleStrategy):
    def __init__(self, config=config):
        super().__init__(config=config, species="electron")


class MuonStrategy(ParticleStrategy):
    def __init__(self, config=config):
        super().__init__(config=config, species="muon")


class GhostStrategy(BaseSpeciesStrategy):
    """Bespoke treatment for ghosts, using suitably-truthmatched user-specified simulation samples."""

    def __init__(self, config=config):
        self._species = "ghost"
        self._hadron_enriched_sel = config.pid.ghost_config.hadron_enriched_selection

    @property
    def species(self):
        return self._species

    @property
    def hadron_enriched_sel(self):
        return self._hadron_enriched_sel

    def get_species_name(self) -> str:
        return "ghost"

    def get_species_alias(self) -> str:
        return "ghost"
