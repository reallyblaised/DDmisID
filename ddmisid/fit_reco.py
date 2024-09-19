import pyhf
import pickle
import argparse
from typing import Optional, Callable, Any
import pyhf.workspace
from ddmisid import load_hist, data_pull_plot, make_binning, make_legend
from pathlib import Path
import numpy as np
import cabinetry
from pathlib import Path
import matplotlib.pyplot as plt
from uncertainties import ufloat, ufloat_fromstr
import json
import matplotlib.patches as patches

plt.style.use(["science", "no-latex"])
cabinetry.set_logging()


def assign_reco_labels(categories: list, include_ghosts: bool = True) -> list:
    """Assign reco labels to the categories"""
    labels = []
    for i in categories:
        match i:
            case "kaon": labels.append(r"$K$")
            case "pion": labels.append(r"$\pi$")
            case "proton": labels.append(r"$p$")
            case "electron": labels.append(r"$e$")
            case "ghost": 
                if include_ghosts:
                    labels.append(r"g")
    
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
                "name": f"{template_name}_yield",
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
            if np.any(_template_h.view().value): # at least one non-zero bin entry
                samples.append(
                    build_template_spec(template_name, 
                        list(_template_h.view().value / np.sum(_template_h.view().value)) # normalise to unity                   
                )
            )
            else: # empty - fill with placeholder small template bin contents
                samples.append(
                    build_template_spec(template_name, 
                        [1e-6 for val in _template_h.view().value] # dummy value to avoid 0.0, aiding the NLL minimisation                    
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
    observations = list(load_hist(obs).view())
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
            {"name": "Hadron-enriched data", "data": observations} 
        ],
        "measurements": [
            {
            "name": "True !mu abundance extraction",
            "config": {
                "poi": 'pion_yield', "parameters": [
                    # bounds on floating normalisations
                    {
                        "name": "proton_yield",
                        "bounds": [[0.0, np.sum(observations)]],
                        "inits": [np.sum(observations)/10.0] # initialise at 10% of the data   
                    },
                    {
                        "name": "kaon_yield",
                        "bounds": [[0.0, np.sum(observations)]],
                        "inits": [np.sum(observations)/10.0] # initialise at 10% of the data   
                    },
                    {
                        "name": "pion_yield",
                        "bounds": [[0.0, np.sum(observations)]],
                        "inits": [np.sum(observations)/10.0] # initialise at 10% of the data   
                    },
                    {
                        "name": "electron_yield",
                        "bounds": [[0.0, np.sum(observations)]],
                        "inits": [np.sum(observations)/10.0] # initialise at 10% of the data   
                    },
                    {
                        "name": "ghost_yield",
                        "bounds": [[0.0, np.sum(observations)]],
                        "inits": [np.sum(observations)/10.0] # initialise at 10% of the data   
                    },
                ]}
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
        "proton": "#225ea8",
        "pion": "#1d91c0",
        "kaon": "#7fcdbb",
        "electron": "#edf8b1",
        "ghost": "#081d58",
    }
    figures = cabinetry.visualize.data_mc(model_pred_postfit, data, colors=colors, save_figure=False)
    fig = figures[0]["figure"]
    ratio_panel = fig.get_axes()[1]
    ratio_panel.set_xlabel(r"Hadron-enriched $reco$ categories")
    ratio_panel.set_ylabel(r"Data/Model")
    hist_plot = fig.get_axes()[0]
    hist_plot.set_ylabel(r"Candidates [Arbitrary Units]") 



    # NOTE: custom viz of results
    postfit_info = cabinetry.tabulate.yields(model_pred_postfit, data) 
    obs = np.array(
        list(postfit_info['yields_per_bin'][-1].values())[1:], 
        dtype=float
    ) # unweighted

    # sanity checks - inspect the compatibility (assume ghosts are the last entry)
    _obs = load_hist(opts.obs)
    assert _obs.view().all() == obs.all()

    # book canvas
    fig_pull, ax, axp = data_pull_plot(
        axp_xlabel=r"Hadron-enriched $reco$ categories",
        annotation=None,
        ylabel="Candidates [A.U.]",
        is_pull=False, # data/model
        axp_ylabel = r"$\frac{\mathrm{Data}}{\mathrm{Model}}$",
    )
    # plot the data
    ax.errorbar(
        assign_reco_labels([_obs.axes['reco'][i] for i in range(_obs.axes['reco'].size)]),
        obs,
        yerr = np.sqrt(obs.view()), # unweighted observations
        color="black",
        fmt=".",
        markersize=3,
        label="Data",
    )

    include_ghosts = False
    bottom = np.zeros(len(obs))
    for s_idx in range(len(postfit_info['yields_per_bin'])):
        category = postfit_info['yields_per_bin'][s_idx]['sample']
        match category: 
            case 'proton': cosmetics = {'color' : colors['proton'], 'edgecolor' : colors['proton'], 'label' : r"$p$"}
            case 'pion': cosmetics = {'color' : colors['pion'], 'edgecolor' : colors['pion'], 'label' : r"$\pi$"}
            case 'kaon': cosmetics = {'color' : colors['kaon'], 'edgecolor' : colors['kaon'], 'label' : r"$K$"}
            case 'electron': cosmetics = {'color' : colors['electron'], 'edgecolor' : colors['electron'], 'label' : r"$e$"}
            case 'ghost': cosmetics = {'color' : colors['ghost'], 'edgecolor' : colors['ghost'], 'label' : r"g"}
            case 'total': cosmetics = {'edgecolor': '#b10026', 'label' : 'Model', 'facecolor' : None}
            case _:
                break
                #raise KeyError(f"Unexpected sample: {postfit_info['yields_per_bin'][s_idx]['sample']}")
        


        # stacked per-species yields
        if category not in ('data', 'total'):
            y = np.array(
                [ufloat_fromstr(i).n for i in list(postfit_info['yields_per_bin'][s_idx].values())[1:]],
                dtype=int
            )   
            # plot stacked hist
            ax.bar(
                height = y,
                x = assign_reco_labels([_obs.axes['reco'][i] for i in range(_obs.axes['reco'].size)]),
                width=1.0,
                bottom=bottom,
                **cosmetics,
            )
            bottom += y       
        
        # total fit model
        if category=='total':
            y = np.array(
                [ufloat_fromstr(i) for i in list(postfit_info['yields_per_bin'][s_idx].values())[1:]],
            )   
            
            # fit model uncertainty
            ax.bar(
                x = assign_reco_labels([_obs.axes['reco'][i] for i in range(_obs.axes['reco'].size)]),
                height = [b.s*2 for b in y],
                bottom = [b.n - b.s for b in y],
                facecolor='none',
                hatch='///////',
                linewidth=0.0000,
                width=1.0,
                label="Uncertainty",
                edgecolor='tab:grey',
                alpha=0.5
            )                

            # data/model subplot
            # HACK: replace absolute zero yield with small values
            y[y==0.0] = ufloat(1e-6, 1e-6)
            dom = obs/y
            uobs = np.array([
                ufloat(o, o**.5) for o in obs
            ])
            udom = uobs/y # data over model, accounting for poison errors
            axp.axhline(1.0, color="black", lw=0.33, ls=":")
            axp.errorbar(
                x = assign_reco_labels([_obs.axes['reco'][i] for i in range(_obs.axes['reco'].size)]),
                y = [r.n for r in dom],
                yerr = [r.s for r in dom],
                color="black",
                fmt=".",
                markersize=3,
            )
            # hatched error total error budget 
            axp.bar(
                x = assign_reco_labels([_obs.axes['reco'][i] for i in range(_obs.axes['reco'].size)]),
                height = [2*r.s for r in udom], # all uncertainties folded in
                bottom = [1.0-r.s for r in udom],
                width=1.0,
                edgecolor='tab:grey',
                facecolor='none',
                alpha=0.5,
                hatch='///////',
                lw=0.0,
            ) 
            axp.set_ylim(0.5, 1.5)
            axp.set_yticks([0.75, 1.0, 1.25]) 
            axp.tick_params(
                axis="y", labelsize="small"
            )  # Set y-axis tick label font size to small

    # legend
    handles, labels = ax.get_legend_handles_labels()
    order =[0, 1, 3, 4, 5, 2, 6] # order of legend
    ax.legend(
        [handles[idx] for idx in order], [labels[idx] for idx in order],    
        loc='center left', bbox_to_anchor=(1, 0.5)
    )
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))  # Force scientific notation

    # save figs
    # [fig.savefig(f"{_outpath}/postfit_cabinetry.{ext}") for ext in ("pdf", "png", "eps")] # backup
    [fig_pull.savefig(f"{_outpath}/postfit.{ext}") for ext in ("pdf", "png", "eps")]
