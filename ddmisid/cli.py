"""Command-line interface to build the DDmisID engine"""

import click
from .engine import _run_snakemake, _load_config, config
from ddmisid.auth import kinit
from pydantic import ValidationError
from loguru import logger
from pathlib import Path
from .pydantic_config_model import DDmisIDConfig  # Import the model for validation

# Setup logging
logdir = Path("logs/engine")
logdir.mkdir(parents=True, exist_ok=True)
logger.add(
    logdir / "ddmisid_log_{time}.log",
    rotation="1 day",
    retention="7 days",
    level="INFO",
)


@click.group()  # Add this to ensure cli group is registered
def cli():
    """Command-line interface to build the DDmisID engine"""
    pass


@cli.command()
@click.option(
    "-c",
    "--config-path",  # This is the CLI flag (external option for the user)
    help="Path to the YAML user-defined configuration file",
    default="config/main.yml",
    required=False,
)
def build(config_path):  # This is the corresponding Python variable for the option
    """Build the DDmisID configuration and objects.

    Example:
        ddmisid-engine build --config-path config/main.yml
    """
    logger.info(f"Building the DDmisID engine from: {config_path}")

    # Load and validate the configuration
    try:
        _load_config(config_path)  # Pass config_path to _load_config
        logger.info(f"DDmisID engine built successfully from: {config_path}")

        # Tabulate the configuration for a cleaner display
        config_dict = config.dict()  # Convert to a dictionary
        table_data = []

        # Flatten and format the dictionary into a table-friendly structure
        for section, subsection in config_dict.items():
            if isinstance(subsection, dict):  # Handle nested dicts
                for key, value in subsection.items():
                    table_data.append([f"{section}.{key}", value])
            else:
                table_data.append([section, subsection])

        # Print the table
        table_report = tabulate(
            table_data, headers=["Config Key", "Value"], tablefmt="grid"
        )
        click.echo(table_report)

    except ValidationError as e:
        logger.error(f"Validation failed for configuration: {e}")
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

    kwargs = parse_snakemake_args(snakemake_args)

    # Pass the keyword arguments to _run_snakemake
    try:
        _run_snakemake(**kwargs)
        logger.info(
            f"Snakemake pipeline ran successfully with arguments: {snakemake_args}"
        )
    except Exception as e:
        logger.error(f"Snakemake execution failed: {e}")


def parse_snakemake_args(snakemake_args):
    """Parse Snakemake CLI arguments into keyword arguments."""
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
            kwargs[arg] = snakemake_args[i + 1]
            i += 1
        else:
            kwargs[arg] = True
        i += 1
    return kwargs


if __name__ == "__main__":
    cli()
