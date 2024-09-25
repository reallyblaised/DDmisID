"""Engine build and runs for the DDmisID pipeline."""

from ddmisid.utils import read_config
from loguru import logger


def build_schema(config_path: str = "config/main.yml") -> dict:
    """Builds commonly referenced objects from the configuration."""
    # Read configuration
    config = read_config(config_path)
    logger.info(f"Configuration loaded successfully from {config_path}")

    # Example: Build objects based on config
    # You can add logic here to create objects, initialize components, etc.
    logger.info(f"Building objects with config: {config}")

    # Returning some object for example
    return config


def run_workflow() -> None:
    """Runs the Snakemake backend pipeline."""
    pass
