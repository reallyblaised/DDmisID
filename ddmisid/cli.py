"""Command-line interface to build the DDmisID engine"""

import click
from .engine import build_schema, run_workflow
from .snakemake_wrapper import run_snakemake
from ddmisid.auth import kinit
from loguru import logger
from pathlib import Path


# Setup logging
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
    default="config/main.yml",
    required=False,
)
def build(config):
    """Build the DDmisID configuration and objects.

    Example:
        ddmisid-engine build --config config/main.yml
    """
    logger.info(f"Building the DDmisID engine from: {config}")
    try:
        build_schema(config)
        logger.info(f"DDmisID engine built successfully from: {config}")
    except Exception as e:
        logger.error(f"Failed to build DDmisID engine: {e}")


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
    i = 0
    while i < len(snakemake_args):
        arg = snakemake_args[i]
        if "=" in arg:
            key, value = arg.split("=", 1)
            kwargs[key] = value
        elif (
            arg.startswith("--")
            and i + 1 < len(snakemake_args)
            and not snakemake_args[i + 1].startswith("--")
        ):
            # Handle flags with space-separated values (e.g., --cores 4)
            kwargs[arg] = snakemake_args[i + 1]
            i += 1
        else:
            kwargs[arg] = True
        i += 1

    # Pass the keyword arguments to run_snakemake
    try:
        # First, validate Kerberos ticket

        run_snakemake(**kwargs)
        logger.info(
            f"Snakemake pipeline ran successfully with arguments: {snakemake_args}"
        )
    except Exception as e:
        logger.error(f"Snakemake execution failed: {e}")


if __name__ == "__main__":
    cli()
