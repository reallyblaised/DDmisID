"""Objects for fitting the misID data."""

__authors__ = ["Blaise Delaney", "Kevin Kurashima"]
__email__ = "blaise.delaney at cern.ch"

import pyhf
import yaml
import numpy as np
import pickle
import matplotlib.pyplot as plt
from typing import List


# RFE:
# - [ ] Add documentation
# - [ ] normalise templates to unity (else mu != the yield)
# - [ ] verify path to templates
# - [ ] verify path to data
class BMLFitter:
    def __init__(self, config_file: str, data_file: str):
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
        return f"MaxLikelihoodFitter(config_file={self.config_file}, data_file={self.data_file})"

    def _set_templates(self):
        for template_name, template in self.config.items():
            if isinstance(template, List):
                setattr(self, template_name, np.array(template))
            else:
                raise ValueError(f"Template {template_name} is not a list.")

    def fit(self):
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
        self.model = pyhf.Model(spec)
        self.result = pyhf.infer.mle.fit(self.data, self.model)
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
            (self.data - self.result) / np.sqrt(self.result)
        )  # compute the pulls
        ax[1].set_xlabel("Bin number")
        ax[1].set_ylabel("Pull")
        plt.show()

    def save_results(self, save: bool = False, filename: str = "results.json"):
        if save:
            with open(filename, "w") as file:
                json.dump(self.result_dict, file)
        return self.result_dict
