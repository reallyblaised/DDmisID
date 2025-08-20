"""Hack my way through the generation of configuration files for the misID study"""

from pathlib import Path
import yaml
from jinja2 import Template
import os


def load_yaml_template(file_path):
    with open(file_path, "r") as f:
        return yaml.safe_load(f)


def save_yaml(data, output_path):
    with open(output_path, "w") as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)


def get_charm_modes(stream):
    # Define charm decay modes for each stream
    modes = {
        "D0MuNu": ["D02KPi", "D02K3Pi"],
        "DpMuNu": ["DpKPiPi"],  # Add more Dp modes if needed
    }
    return modes[stream]


def generate_configs():
    # Configuration parameters
    streams = ["D0MuNu", "DpMuNu"]
    years = ["2011", "2012", "2015", "2016", "2017", "2018"]
    magpols = ["MU", "MD"]

    # Load template
    template_config = load_yaml_template("config/main.yml")

    for stream in streams:
        charm_modes = get_charm_modes(stream)
        for charm_mode in charm_modes:
            for year in years:
                for magpol in magpols:
                    # Create a deep copy of the template
                    new_config = yaml.safe_load(yaml.dump(template_config))

                    # Update the year and magpol in PID section
                    new_config["pid"]["year"] = year
                    new_config["pid"]["magpol"] = "up" if magpol == "MU" else "down"

                    # Update the data_key with charm mode
                    new_config["data"][
                        "data_key"
                    ] = f"B2DMuNuX_{charm_mode}_FakeMuonTuple"

                    # Update the file paths
                    base_path = (
                        f"/ceph/submit/data/user/b/blaised/bc2dmunu/fakemuon/DATA"
                    )

                    # Input path
                    new_config["data"][
                        "input_path"
                    ] = f"{base_path}/merged_sj_data/{stream}/{year}/{magpol}/Bc2{stream}X_{charm_mode}.root"

                    # Output path
                    new_config["data"][
                        "output_path"
                    ] = f"{base_path}/merged_sj_data_misid_w/{stream}/{year}/{magpol}/Bc2{stream}X_{charm_mode}.root"

                    # Generate output filename
                    output_filename = f"{stream}_{charm_mode}_{year}_{magpol}.yml"

                    # Create output directory if it doesn't exist
                    os.makedirs("generated_configs", exist_ok=True)
                    output_path = Path("generated_configs") / output_filename

                    # Save the configuration
                    save_yaml(new_config, output_path)
                    print(f"Generated: {output_filename}")


if __name__ == "__main__":
    generate_configs()
