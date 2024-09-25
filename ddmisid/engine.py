from ddmisid.utils import read_config
from loguru import logger


def build_objects(config_path):
    """Builds commonly referenced objects from the configuration."""
    # Read configuration
    config = read_config(config_path)
    logger.info(f"Configuration loaded successfully from {config_path}")

    # Example: Build objects based on config
    # You can add logic here to create objects, initialize components, etc.
    logger.info(f"Building objects with config: {config}")

    # Returning some object for example
    return config
