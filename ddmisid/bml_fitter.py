"""Objects for fitting the reco categories in each partition of the hadron-enriched control sample."""

__authors__ = ["Blaise Delaney", "Kevin Kurashima"]
__email__ = "blaise.delaney at cern.ch"

import pyhf
import yaml
import numpy as np
import pickle
import matplotlib.pyplot as plt
from typing import List
import json
from typing import Any
import matplotlib.pyplot as plt
import scienceplots

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

    def __init__(self, data_file: str, config_file: str = "config/main.yml"):
        # Load the configuration
        with open(config_file, "r") as file:
            self.config = yaml.safe_load(file)

        # Load the data
        with open(data_file, "rb") as file:
            self.data = pickle.load(file)

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
            if hasattr(self, template_name):
                template = getattr(self, template_name)
                template_sum = np.sum(template.view().value)
                if template_sum != 0:
                    normalized_template = template / template_sum
                    setattr(self, template_name, normalized_template)
                else:  # Handle empty template
                    # Setting a null-populated template
                    null_template = np.zeros_like(template)
                    setattr(self, template_name, null_template)

    def _set_templates(self) -> None:
        """Set the templates, normalised to unity, as attributes of the class."""
        for template_name, template in self.config.items():
            if isinstance(template, List):
                setattr(self, template_name, np.array(template))
            else:
                raise ValueError(f"Template {template_name} is not a list.")

        # Normalise the templates to attain unit area
        self._normalise_templates()

    def build_workspace(self) -> dict[str, Any]:
        """Build the workspace spec."""
        channel = {"name": "singlechannel", "samples": [], "data": self.data.tolist()}

        for template_name in self.config.keys():
            if hasattr(self, template_name):
                template = getattr(self, template_name)
                sample = {
                    "name": template_name,
                    "data": template.tolist(),
                    "modifiers": [{"name": "mu", "type": "normfactor", "data": None}],
                }
                channel["samples"].append(sample)

        spec = {"channels": [channel]}

        # return the workspace
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
        # build the workspace spec
        self.spec = self.build_workspace()

        # accordingly, create the pyhf model
        self.model = pyhf.Model(self.spec)

        # minimize the negative log likelihood
        self.result = pyhf.infer.mle.fit(self.data, self.model)

        # store result to dict
        self.result_dict = {"fit_result": self.result.tolist()}

        return self.result

    def visualize(self):
        fig, ax = plt.subplots(2, 1, figsize=(10, 8))
        ax[0].plot(self.data, label="Data")
        ax[0].plot(self.result, label="Fit Result")
        ax[0].set_xlabel("Bin number")
        ax[0].set_ylabel("Counts")
        ax[0].legend()
        ax[1].plot(
            (self.data - self.result) / np.sqrt(self.data)
        )  # compute the pulls, assume poisson error
        ax[1].set_xlabel("Bin number")
        ax[1].set_ylabel("Pull")
        plt.show()

    def save_results(self, save: bool = False, filename: str = "results.json"):
        if save:
            with open(filename, "w") as file:
                json.dump(self.result_dict, file)
        return self.result_dict
