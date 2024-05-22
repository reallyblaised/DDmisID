"""
Snakemake routine to orchestrate the requisite DDmisID stages:
    - PIDCalib2 setup
    - Spawning of PIDCalib2 jobs
    - Running of PIDCalib2 jobs -> efficiency histograms (anti-muon, signal)
    - Binned template creation for the species in the hadron-enriched makeup
    - Discretisation of the data into hadron-enriched, species-specific bins
    - Fit to the reco data bins with !mu PID templates; account for cross-contamination between high-purity data partitions
    - Extraction of the true abundance of each !mu species
    - Weight assignment yielding the effective per-species count in the single-muon misID composition
"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

import shutil
import os
from snakemake.io import glob_wildcards     


# cleanup
onsuccess: 
    shutil.rmtree(".snakemake/metadata")


rule all:
    input:
        "antimuon_abundances.done"


# ==================================================================================
#                          PID Efficiency Map Section
# ----------------------------------------------------------------------------------
#   - spawn the executables to run PIDCalib2 to generate the PID maps (!mu, mu)
#   - run the executables -> histograms.pkl with binning directory sub-structure
#   - signal completion of all PIDCalib2 jobs 
# ==================================================================================
checkpoint gen_sh_files:
    """
    Spawn the PIDCalib2 jobs as instructed in config/main.yml; executables store in /bin/
    """
    input:
        config_executable = "ddmisid/setup_pidcalib_jobs.py",
    params:
        is_test = f"--test",
        reco = f"--pid_regime reco",
        he_all = f"--pid_regime he_all",
    output:
        dir = directory("bin/"),
    shell:
        "python {input.config_executable} {params.is_test} {params.reco} {params.he_all}"


rule run_pidcalib:
    """
    Run the PIDCalib2 jobs as spawned in the previous step; store the efficiency histograms in /bin/
    """
    input:
        bash_file = "bin/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/run.sh",
    output:
        pidcalib_hists = "bin/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/perf.pkl",
    run:
        shell("bash {input.bash_file}")


def aggregate_pidcalib(wildcards):
    '''
    Collect the variable number of sub-directories produced by setup.py that live in /exec/,
    and expand to use perfHists.root files produced by PIDCalib2
    '''
    checkpoint_output = checkpoints.gen_sh_files.get(**wildcards).output[0]
    year, magpol, dllmu_bin, species, eff_dirs = glob_wildcards(os.path.join(checkpoint_output, '{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/run.sh'))

    return expand(checkpoint_output+"/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/perf.pkl", zip, year=year, magpol=magpol, dllmu_bin=dllmu_bin, species=species, eff_dirs=eff_dirs)


rule collect_pidcalib:
    """
    Signal completion of all PIDCalib2 jobs
    """
    input:
        aggregate_pidcalib
    output:
        combined = temp("pid_efficiency_maps.done"),
    shell:
        "echo {input} > {output.combined}"


# # =========================================================================================
# #                     Hadron-Enriched Data Partitioning Section
# # -----------------------------------------------------------------------------------------
# #   - partition HE data into species-specific bins, within kinematic & occupancy partitions 
# # =========================================================================================
checkpoint discretize_data:
    """
    Spawn the data bins in the form bespoke histograms (/obs/*/*.pkl)
    These have all kinematic & occupancy cuts applied, as well has the PID cuts for the hadron-enriched bins specified in the config file
    """
    input:
        data_path = "scratch/99_2021_01_443970_443970746_Bc2D0MuNuXSlim.root", # TODO: retrieve from config file
    output:
        directory("obs/"),
    shell:
        "python ddmisid/data_discretizer.py {input.data_path}"


def aggregate_discretizer(wildcards):
    """
    Collect the variable number of sub-directories produced by partictioning the data
    """
    checkpoint_output = checkpoints.discretize_data.get(**wildcards).output[0]
    p, eta, ntracks = glob_wildcards(os.path.join(checkpoint_output, '{p}/{eta}/{ntracks}/obs.pkl'))
    return expand(checkpoint_output+"/{p}/{eta}/{ntracks}/obs.pkl", zip, p=p, eta=eta, ntracks=ntracks)


rule collect_discretizer:
    """
    Signal completion of all data discretization jobs
    """
    input:
        aggregate_discretizer
    output:
        combined = temp("data_partitions.done"),
    shell:
        "echo {input} > {output.combined}"


# ==================================================================================
#                          Binned Maximum Likelihood Section
# ----------------------------------------------------------------------------------
#   - build hadron-enriched templates for each species across all reco partitions
#   - binned ML fit to the hadron-enriched data in bins of kinematics and occupancy
#   - extract the true abundance of each species, accounting for cross-contamination 
# ==================================================================================
checkpoint make_templates:
    input:
        pidcalib_done = "pid_efficiency_maps.done" # required to build the templates
    params:
        script = "ddmisid/make_templates.py",
        path_prefix = r"bin/2018/down/antimu_id", # TODO: retrieve from config, and deal with magpol
    output: 
        directory("templates/"),
    shell:
        "python {params.script} {params.path_prefix}" 


rule bml_fit: 
    input:
        # dummy -> correct workflow
        data_partitions = "data_partitions.done", # signal spawning of hadron-enriched data partitions
        # templates
        pion_template = lambda wildcards: [
            "templates/{p}/{eta}/{ntracks}/{species}.pkl".\
            format(p=wildcards.p, eta=wildcards.eta, ntracks=wildcards.ntracks, species="pion")
        ],
        kaon_template = lambda wildcards: [
            "templates/{p}/{eta}/{ntracks}/{species}.pkl".\
            format(p=wildcards.p, eta=wildcards.eta, ntracks=wildcards.ntracks, species="kaon")
        ],
        proton_template = lambda wildcards: [
            "templates/{p}/{eta}/{ntracks}/{species}.pkl".\
            format(p=wildcards.p, eta=wildcards.eta, ntracks=wildcards.ntracks, species="proton")
        ],
        electron_template = lambda wildcards: [
            "templates/{p}/{eta}/{ntracks}/{species}.pkl".\
            format(p=wildcards.p, eta=wildcards.eta, ntracks=wildcards.ntracks, species="electron")
        ],
        # observations
        obs = lambda wildcards: [
            "obs/{p}/{eta}/{ntracks}/obs.pkl".\
            format(p=wildcards.p, eta=wildcards.eta, ntracks=wildcards.ntracks)
        ],
    output:
        # true abundance of each species, accounting for cross-contamination between reco bins
        "hadron_enriched_yields/{p}/{eta}/{ntracks}/proton.json",
        "hadron_enriched_yields/{p}/{eta}/{ntracks}/kaon.json",
        "hadron_enriched_yields/{p}/{eta}/{ntracks}/pion.json",
        "hadron_enriched_yields/{p}/{eta}/{ntracks}/electron.json",
    run:
        shell("touch {output}")


def hadron_enriched_abundances(wildcards):
    """Trick: book the desired process outcome out of the checkpoint, so to be able to use wildcards to enact the fit"""
    checkpoint_output = checkpoints.make_templates.get(**wildcards).output[0]
    p = list(set(glob_wildcards(os.path.join(checkpoint_output, '{p}/{eta}/{ntracks}/{species}.pkl')).p))
    eta = list(set(glob_wildcards(os.path.join(checkpoint_output, '{p}/{eta}/{ntracks}/{species}.pkl')).eta))
    ntracks = list(set(glob_wildcards(os.path.join(checkpoint_output, '{p}/{eta}/{ntracks}/{species}.pkl')).ntracks))
    species = list(set(glob_wildcards(os.path.join(checkpoint_output, '{p}/{eta}/{ntracks}/{species}.pkl')).species))
    
    return expand("hadron_enriched_yields/{p}/{eta}/{ntracks}/{species}.json", p=p, eta=eta, ntracks=ntracks, species=species)


rule aggregate_antimuon_abundances:
    input:
        hadron_enriched_abundances
    output:
        temp("antimuon_abundances.done")
    run:
        shell("touch {output}")
        





# checkpoint consistency_check:
#     input:
#         pidcalib_done = "pidcalib_full_run.done",
#         templates_done = "templates.done",
#     params:
#         script = "ddmisid/consistency_checker.py",
#         eff_hist_path_prefix = r"bin/2018/down/antimu_id",
#         template_path_prefix = r"templates",
#     output:
#         directory("checks/"),
#     shell:
#         "python {params.script} {params.eff_hist_path_prefix} {params.template_path_prefix}"


# def aggregate_checks(wildcards):
#     checkpoint_output = checkpoints.consistency_check.get(**wildcards).output[0]
#     species = glob_wildcards(os.path.join(checkpoint_output, '{species}.{ext}'))
#     return expand(checkpoint_output+"/{species}.{ext}", zip, species=species, ext=ext)


# rule collect_checks:
#     input:
#         aggregate_checks,
#     output:
#         combined = "consistency_checks.done",
#     shell:
#         "echo {input} > {output.combined}"