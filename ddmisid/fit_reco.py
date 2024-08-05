import pyhf
import pickle
import argparse
from typing import Optional, Callable, Any
import pyhf.workspace
from ddmisid import load_hist, data_pull_plot, make_binning
from pathlib import Path
import numpy as np
import cabinetry
from pathlib import Path
import matplotlib.pyplot as plt
from uncertainties import ufloat
import json

plt.style.use(["science", "no-latex"])
cabinetry.set_logging()


def assign_reco_labels(categories: list) -> list:
    """Assign reco labels to the categories"""
    labels = []
    for i in categories:
        match i:
            case "kaon": labels.append(r"kaon")
            case "pion": labels.append(r"pion")
            case "proton": labels.append(r"proton")
            case "electron": labels.append(r"electron")
            case "ghost": labels.append(r"ghost")
    
    return labels


def load_hist(path: str) -> list:
    """Load the data from the file"""
    with open(f"{path}", "rb") as f:
        return pickle.load(f)


def build_template_spec(
    template_name: str,
    template_data: list
) -> dict:
    """Build the schema for each template - assuming same likelihood modifiers for all templates
    """
    return {
        "name": template_name,
        "data": template_data,
        "modifiers": [
            {
                "name": f"{template_name}",
                "type": "normfactor",
                "data": None,
            }
        ]
    }
        

def build_sample_spec(
    proton_template: Path | str | None = None,
    pion_template:  Path | str | None = None,
    kaon_template:  Path | str | None = None,
    electron_template:  Path | str | None = None,
    ghost_template:  Path | str | None = None,
)->list:
    """Build the schema for the BML fit to hadron-enriched data"""
    samples = []
    for template in [proton_template, pion_template, kaon_template, electron_template, ghost_template]:
        if template is not None:
            # assign label
            if template==proton_template:
                template_name = r"proton"
            elif template==pion_template:
                template_name = f"pion"
            elif template==kaon_template:
                template_name = r"kaon"
            elif template==electron_template:
                template_name = r"electron"
            elif template==ghost_template:
                template_name = r"ghost"
            # if there is something to load, do it
            _template_h = load_hist(template)
            # build the samples
            samples.append(
               build_template_spec(template_name, 
                    list(_template_h.view().value / np.sum(_template_h.view().value)) # normalise to unity                   
               )
            )
    return samples


def build_channel_spec(
    obs: list,
    proton_template: Optional[list] = None,
    pion_template: Optional[list] = None,
    kaon_template: Optional[list] = None,
    electron_template: Optional[list] = None,
    ghost_template: Optional[list] = None,
)->object:
    """Build the schema for the channel"""
    schema = {
        "channels": [
            {
                "name": "Hadron-enriched data",
                "samples": build_sample_spec(
                    proton_template,
                    pion_template,
                    kaon_template,
                    electron_template,
                    ghost_template,
                ),
            }
        ],
        "observations":[
            {"name": "Hadron-enriched data", "data": list(load_hist(obs).view()[:-1])} # FIXME: neglect ghost temporarily
        ],
        "measurements": [
            {
            "name": "True !mu abundance extraction",
            "config": {"poi": 'pion', "parameters": []}
            }
        ],
        "version": "1.0.0"
    }
    # build and validate according to pyhf specification
    ws = pyhf.workspace.Workspace(schema, validate=True)
    return ws 



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Template and observations for binned ML fit")
    parser.add_argument(
        "--obs", help="path to observations", required=True,
    )
    parser.add_argument(
        "--proton", help="path to proton template", default=None
    )
    parser.add_argument(
        "--pion", help="path to pion template", default=None
    )
    parser.add_argument(
        "--kaon", help="path to kaon template", default=None
    )
    parser.add_argument(
        "--electron", help="path to muon electron", default=None
    )
    parser.add_argument(
        "--ghost", help="path to muon ghost", default=None
    )
    opts = parser.parse_args()

    # binned maximum likelihood fit
    # -----------------------------
    workspace = build_channel_spec(
        obs = opts.obs,
        proton_template = opts.proton,
        pion_template = opts.pion,
        kaon_template = opts.kaon,
        electron_template = opts.electron,
        ghost_template = opts.ghost,
    ) 
    
    # having built the workspace in pyhf, use cabinetry to fit and viz
    model, data = cabinetry.model_utils.model_and_data(workspace)
    fit_results = cabinetry.fit.fit(model, data, custom_fit=True) # use iminuit directly
    model_pred_postfit = cabinetry.model_utils.prediction(model, fit_results=fit_results)

    # persist the results
    # -------------------
    results_info = {}    
    for i in range(len(fit_results.labels)):
        results_info[fit_results.labels[i]] = ufloat(
            fit_results.bestfit[i], fit_results.uncertainty[i]
        ) / np.array([ufloat(d, np.sqrt(d)).n for d in data]).sum() # relative abundance of species wrt the observation    

    # persist the results to pickle file
    _outpath = Path(str(Path(opts.obs).parent).replace('obs', 'postfit'))
    _outpath.mkdir(parents=True, exist_ok=True)
    with open(f"{_outpath}/yields.pkl", "wb") as f:
        pickle.dump(results_info, f)
    
    # viz results 
    # -----------
    colors = {
        "proton": "#2166ac",
        "pion": "#92c5de",
        "kaon": "#b2182b",
        "electron": "#f4a582",
        "ghost": "#878787",
    }
    figures = cabinetry.visualize.data_mc(model_pred_postfit, data, colors=colors, save_figure=False)
    fig = figures[0]["figure"]
    ratio_panel = fig.get_axes()[1]
    ratio_panel.set_xlabel(r"Hadron-enriched $reco$ categories")
    ratio_panel.set_ylabel(r"Data/Model")
    hist_plot = fig.get_axes()[0]
    hist_plot.set_ylabel(r"Candidates [Arbitrary Units]") 


    # NOTE: custom viz of results
    # obs = load_hist(opts.obs)
    #_ = cabinetry.tabulate.yields(model_pred_postfit, data)
    # fig, ax, axp = data_pull_plot(
    #     axp_xlabel=r"Hadron-enriched $reco$ categories",
    #     annotation=None,
    # )
    # breakpoint()
    # # plot the data
    # ax.errorbar(
    #     assign_reco_labels([obs.axes['reco'][i] for i in range(obs.axes['reco'].size)]),
    #     obs.view(),
    #     yerr = np.sqrt(obs.view()), 
    #     color="black",
    #     fmt=".",
    #     markersize=3,
    #     label="Data",
    # )

    # breakpoint()
    # # plot the model
    # ax.plot(
    #     assign_reco_labels([obs.axes['reco'][i] for i in range(obs.axes['reco'].size)]),
    #     model_pred_postfit[0],
    #     color="red",
    #     label="Model",
    # )

    # go onto saving all the info 
    # figure

    [fig.savefig(f"{_outpath}/projection.{ext}") for ext in ("pdf", "png", "eps")]
