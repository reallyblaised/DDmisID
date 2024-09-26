"""Build the engine config and run the Snakemake backend pipeline."""

from pydantic import ValidationError  # Correct the typo
from .pydantic_config_model import DDmisIDConfig
from ddmisid import read_config  # Assuming this is the correct path
import subprocess
from loguru import logger
from tabulate import tabulate
from pathlib import Path
import json

# Global variable to store the validated config
config = None
config_path_json = Path(
    ".ddmisid/validated_config.json"
)  # Hidden directory for the validated config


def _load_config(config_path: str):
    """Load and validate the configuration file."""
    global config
    config_data = read_config(config_path)
    config = DDmisIDConfig(**config_data)

    # Persist the validated config to a hidden directory
    config_path_json.parent.mkdir(exist_ok=True)
    with config_path_json.open("w") as f:
        f.write(config.json(indent=4))


def get_config():
    """Return the current value of the config."""
    return config


def _load_validated_config():
    """Load the pre-validated configuration from the JSON file."""
    global config
    if not config_path_json.exists():
        logger.error(f"Run 'ddmisid-engine build' first.")
        raise FileNotFoundError(
            f"Config build missing. Run 'ddmisid-engine build' first."
        )

    with config_path_json.open("r") as f:
        config_data = json.load(f)
        config = DDmisIDConfig(**config_data)


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
