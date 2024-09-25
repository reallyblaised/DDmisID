"""Wrapper for running the Snakemake pipeline within the DDmisID engine with Snakemake's dynamic flags."""

import subprocess
from loguru import logger


def run_snakemake(**kwargs):
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
