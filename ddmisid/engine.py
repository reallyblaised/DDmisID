"""Build the engine config and run the Snakemake backend pipeline."""

from pydantic import ValidationError  # Correct the typo
from .pydantic_config_model import DDmisIDConfig
from ddmisid import read_config  # Assuming this is the correct path
import subprocess
from loguru import logger

# Global variable to store the validated config
config = None


def _load_config(config_path: str):
    """Load and validate the configuration file."""
    global config
    config_data = read_config(config_path)
    config = DDmisIDConfig(**config_data)  # Validate the config with Pydantic


def _run_snakemake(**kwargs):
    """Wrapper for running the Snakemake pipeline with dynamic flags."""
    try:
        # Build the Snakemake command as a list, starting with "snakemake"
        cmd = ["snakemake"]

        # Convert the passed keyword arguments (**kwargs) into command-line flags
        for key, value in kwargs.items():
            flag = f"--{key.replace('_', '-')}"  # Convert underscores to hyphens (standard CLI)
            if isinstance(value, bool) and value:
                cmd.append(
                    flag
                )  # For boolean flags, just add the flag (e.g., --dry-run)
            else:
                cmd.append(
                    f"{flag}={value}"
                )  # For others, add the flag with its value (e.g., --cores=4)

        # Execute the Snakemake command
        logger.info(f"Running Snakemake with command: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        logger.info("Snakemake backend pipeline ran successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Snakemake backend fail: {e}")
