"""Objects for fitting the reco categories in each partition of the hadron-enriched control sample."""

__authors__ = ["Blaise Delaney", "Kevin Kurashima"]
__email__ = "blaise.delaney at cern.ch"

import pyhf
import yaml
import numpy as np
import pickle
from typing import List
import json
from typing import Any
import matplotlib.pyplot as plt
import scienceplots
from numpy.typing import ArrayLike
import hist
import cabinetry

plt.style.use(["science", "no-latex"])


class BMLFitter:
    """Class for fitting the reco categories in each partition of the hadron-enriched control sample.

    Args:
        data_file (str): Path to the data file.
        config_file (str): Path to the config file.

    Attributes:
        config (dict): The configuration dictionary.
        data (np.ndarray): The data.
        spec (dict): The workspace specification.
        model (pyhf.Model): The pyhf model.
        result (np.ndarray): The result of the fit.
        result_dict (dict): The result of the fit as a dictionary.
    """

    def __init__(
        self, data_file: str, template_dir: str, config_file: str = "config/main.yml"
    ) -> None:
        # Load the config to identify species under scrutiny
        with open(config_file, "r") as file:
            self.config = yaml.safe_load(file)

        # Load the data
        with open(data_file, "rb") as file:
            self.data = pickle.load(file)

        self.template_dir = template_dir

        # Set the template attributes from the config file
        self._set_templates()

    def __str__(self):
        attributes = []
        for attr in self.__dict__:
            value = getattr(self, attr)
            if isinstance(value, np.ndarray):
                value = value.tolist()
            attributes.append(f"{attr}: {value}")
        return "\n".join(attributes)

    def __repr__(self):
        return f"BMLFitter(config_file={self.config_file}, data_file={self.data_file})"

    def _normalise_templates(self):
        """Normalise the templates to unit area; the normfactor modifies is now the yield."""
        for template_name in self.config.keys():
            attr_name = f"{template_name}_template"
            if hasattr(self, attr_name):
                template = getattr(self, attr_name)
                template_sum = template.sum()
                if template_sum != 0:
                    normalized_template = template / template_sum
                    setattr(self, attr_name, normalized_template)
                else:  # Handle empty template
                    # Create a null-populated Hist object with the same axes
                    null_template = hist.Hist(template.axes)
                    setattr(self, template_name, null_template)

    def _set_templates(self) -> None:
        """Set the templates, normalised to unity, as attributes of the class."""
        for template_name in self.config.keys():
            # load the pkl template, as identified by the key
            with open(f"{self.template_dir}/{template_name}.pkl", "rb") as f:
                template = pickle.load(
                    f
                ).view()  # load only the bin content of a regular - unweighted - histogram

            if isinstance(
                template, (list, np.ndarray)
            ):  # Check if it is a list or ndarray
                setattr(self, f"{template_name}_template", np.array(template))
            else:
                raise ValueError(
                    f"Template {template_name} is not a list or ndarray."
                )  # Provide a clear error message

        # Normalise the templates to attain unit area
        self._normalise_templates()

    def build_workspace(self) -> dict[str, Any]:
        """Build the workspace spec."""
        channel = {"name": "singlechannel", "samples": []}

        for true_species in self.config.keys():
            if hasattr(self, f"{true_species}_template"):
                template = getattr(self, f"{true_species}_template")
                sample = {
                    "name": f"{true_species}_template",
                    "data": template.tolist(),
                    "modifiers": [
                        {
                            "name": f"{true_species}_mu",
                            "type": "normfactor",
                            "data": None,
                        }
                    ],
                }
                channel["samples"].append(sample)

        spec = {
            "channels": [channel],
            "observations": [{"name": "singlechannel", "data": self.data.tolist()}],
            "measurements": [
                {
                    "name": "Measurement",
                    "config": {
                        "poi": f"{list(self.config.keys())[0]}_mu",  # first species is the POI -> pyhf does not yet support multiple POIs
                        "parameters": [],
                    },
                }
            ],
            "version": "1.0.0",
        }

        # set spec as attribute
        setattr(self, "spec", spec)

        return spec

    def export_workspace(self, filename: str = "workspace.json") -> None:
        """Export the workspace to a JSON file."""
        if not hasattr(self, "spec"):
            raise ValueError("Workspace must be built before it can be exported.")

        with open(filename, "w") as file:
            json.dump(self.spec, file)

        print(f"Workspace exported to {filename}")

    def fit(self) -> dict[str, float]:
        """Execute the NLL minimisation, and persiste the result to dict"""

        # NOTE: interface with cabinetry; not sure how to avoid writing and building the workspace yet
        self.build_workspace()
        workspace_path = "workspace.json"
        self.export_workspace(filename=workspace_path)
        ws = cabinetry.workspace.load(workspace_path)
        model, data = cabinetry.model_utils.model_and_data(ws)

        # accordingly, create the pyhf model
        self.model = model

        # minimize the negative log likelihood
        self.result = cabinetry.fit.fit(model, data)

        return self.result

    # def visualize(self):
    #     fig, ax = plt.subplots(2, 1, figsize=(10, 8))
    #     ax[0].plot(self.data, label="Data")
    #     ax[0].plot(self.result, label="Fit Result")
    #     ax[0].set_xlabel("Bin number")
    #     ax[0].set_ylabel("Counts")
    #     ax[0].legend()
    #     ax[1].plot(
    #         (self.data - self.result) / np.sqrt(self.data)
    #     )  # compute the pulls, assume poisson error
    #     ax[1].set_xlabel("Bin number")
    #     ax[1].set_ylabel("Pull")
    #     plt.show()

    # def save_results(self, save: bool = False, filename: str = "results.json"):
    #     if save:
    #         with open(filename, "w") as file:
    #             json.dump(self.result_dict, file)
    #     return self.result_dict
