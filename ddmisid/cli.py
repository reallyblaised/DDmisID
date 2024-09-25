"""Command-line interface to build the DDmisID engine"""

import click
from ddmisid.engine import build_objects
from ddmisid.snakemake_wrapper import run_snakemake
from loguru import logger
from pathlib import Path

# setup logging
logdir = Path("logs/engine")
logdir.mkdir(parents=True, exist_ok=True)
logger.add(
    logdir / "ddmisid_log_{time}.log",
    rotation="1 day",
    retention="7 days",
    level="INFO",
)


@click.group()
def cli():
    """Command-line interface to build the DDmisID engine"""
    pass


@cli.command()
@click.option(
    "-c",
    "--config",
    help="Path to the YAML user-defined configuration file",
    required=True,
)
def build(config):
    """Build the DDmisID configuration and objects.

    Example:
        ddmisid-engine build --config config/main.yml
    """
    logger.info(f"Building the DDmisID engine from: {config}")
    build_objects(config)


@cli.command()
@click.argument("snakemake_args", nargs=-1)
def run(snakemake_args):
    """Run the Snakemake backend pipeline with dynamic flags.

    This command passes any Snakemake flag to the workflow.

    Example:
        ddmisid-engine run --cores=4 --dry-run
    """
    logger.info(f"Running the Snakemake pipeline with arguments: {snakemake_args}")

    # Convert the list of arguments into a dictionary of keyword arguments for Snakemake
    kwargs = {}
    for arg in snakemake_args:
        if "=" in arg:
            key, value = arg.split("=", 1)
            kwargs[key] = value
        else:
            kwargs[arg] = True  # Treat as a boolean flag (e.g., --dry-run)

    # Pass the keyword arguments to run_snakemake
    run_snakemake(**kwargs)


if __name__ == "__main__":
    cli()
