"""Objects for fitting the reco categories in each partition of the hadron-enriched control sample."""

__authors__ = ["Blaise Delaney", "Kevin Kurashima"]
__email__ = "blaise.delaney at cern.ch"

import pyhf
import yaml
import numpy as np
import pickle
from typing import List
from typing import Any, TypeAlias
import matplotlib.pyplot as plt
import scienceplots
import hist
import json, requests, jsonschema
from tabulate import tabulate
from typing import Callable
from functools import wraps
from pathlib import Path
from numpy.typing import ArrayLike
from functools import namedtuple
import logging

pyhf.set_backend(
    # https://scikit-hep.org/pyhf/_generated/pyhf.optimize.opt_minuit.minuit_optimizer.html
    # 0.5 errordef for NLL
    "numpy",
    pyhf.optimize.minuit_optimizer(verbose=1, errordef=0.5),
)  # get uncertainties and covariance

plt.style.use(["science", "no-latex"])

# Define a type alias for the complex return type
FitResultType: TypeAlias = tuple[ArrayLike, ...] | tuple[namedtuple, ...]


# decorator to verify the correctness of the likelihood minimisation
def validate_nllmin(func: Callable, verbose: bool = False) -> Callable:
    """Exploit the fact that we are using iminuit in the backend to verify the correcteness
    of the NLL minimisation

    if verbose: print all the parameters and the relevant limit
    """

    @wraps(func)
    def wrapping(*args, **kwargs):
        fit_result, result_obj = func(*args, **kwargs)

        # firstly, sanity check: the iminuit errordef if 0.5
        # ref: https://scikit-hep.org/pyhf/_generated/pyhf.optimize.opt_minuit.minuit_optimizer.html
        assert (
            pyhf.optimizer.errordef == 0.5
        ), "MinimiserError: iminuit errordef is not 0.5"

        # validity checks
        assert result_obj.minuit.valid == True, "MinuitError: NLL minimum not valid"
        assert (
            result_obj.minuit.accurate == True
        ), "MinuitError: NLL minimum not accurate"
        assert (
            result_obj.minuit.fmin.has_made_posdef_covar == False
        ), "MinuitError: Hessian matrix forced pos def"
        if result_obj.minuit.fmin.has_parameters_at_limit == True:
            print("Fitted parameters report:")
            if verbose:
                for p in result_obj.minuit.params:
                    print(
                        f"{p.name}, value = {p.value} +/- {p.error}; allowed range = [{p.lower_limit}, {p.upper_limit}]",
                    )
        return fit_result, result_obj

    return wrapping


# RFE: depending on integration with snakemake, we may want to pass the data counts instead of the path
# - [ ] make the data a property
# - [ ] make the templates a property
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
        self,
        data_file: str,
        template_dir: str,
        config_file: str = "config/main.yml",
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
        return f"BMLFitter(data_file={repr(self.data_file)}, template_dir={repr(self.template_dir)}, config_file={repr(self.config)})"

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
                    logging.warning(
                        f"Template {template_name} is empty and cannot be normalized. Normalization is skipped for this template."
                    )
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

        # convert and type-hint to Path
        filename: Path = Path(filename)
        filename.parent.mkdir(parents=True, exist_ok=True)

        # JSON formatting for workspace criteria
        formatted_spec = (
            json.dumps(self.spec, indent=4)
            .replace("None", "null")
            .replace("True", "true")
        )
        with open(filename, "w") as file:
            file.write(formatted_spec)

        print(f"Workspace exported to {filename}")

    def build_validate_workspace(self, filename: str) -> pyhf.Workspace:
        """Build and validate the model and observation"""
        # validate workspace
        workspace = json.load(open(filename))
        schema = requests.get(
            "https://scikit-hep.org/pyhf/schemas/1.0.0/workspace.json"
        ).json()
        # If no exception is raised by validate(), the instance is valid.
        jsonschema.validate(instance=workspace, schema=schema)
        workspace = pyhf.Workspace(workspace)

        return workspace

    @validate_nllmin
    @staticmethod
    def run_minimisation(
        data: object, pdf: object, verbose: bool = True
    ) -> FitResultType:
        """wrapper of the pyhf minisation to make `return_uncertainties` attribute private and True by default.
        Additionally, this enables the converge-checks decorator.
        If verbose is True, print the iminuit fmin object
        """
        fit_result, result_obj = pyhf.infer.mle.fit(
            data, pdf, return_result_obj=True, return_uncertainties=True
        )

        if verbose:
            print(result_obj.minuit.fmin)

        return fit_result, result_obj

    def fit(self, **kwargs) -> dict[str, float]:
        """Execute the NLL minimisation, and persiste the result to dict"""

        # NOTE: interface with cabinetry; not sure how to avoid writing and building the workspace yet
        self.build_workspace()
        workspace_path = "workspace.json"
        breakpoint()
        self.export_workspace(filename=workspace_path)
        self.workspace = self.build_validate_workspace(filename=workspace_path)

        # accordingly, create the pyhf model
        model = self.workspace.model(poi_name=list(self.config.keys())[0] + "_mu")
        data = self.workspace.data(model)

        # minimize the negative log likelihood
        self.result, _ = self.run_minimisation(data, model, **kwargs)

        # print results
        self._to_table()

        return self.result

    def _to_table(self, tex: bool = False) -> None:
        """Print yield results. If tex is True, print in LaTeX format."""
        table = []
        for k in self.result.keys():
            table.append([k, self.result[k].n, self.result[k].s])

        print(
            tabulate(table, headers=["Parameter", "cval", "err"], tablefmt="fancy_grid")
        )
        if tex:
            print(
                tabulate(table, headers=["Parameter", "cval", "err"], tablefmt="latex")
            )

    def prefit_plot(self, **kwargs) -> None:
        """Plot the templates and data, pre-fit."""
        pass

    def postfit_plot(self, **kwargs) -> None:
        """Plot the templates and data, post-fit."""
        pass

    # def save_results(self, save: bool = False, filename: str = "results.json"):
    #     if save:
    #         with open(filename, "w") as file:
    #             json.dump(self.result_dict, file)
    #     return self.result_dict
