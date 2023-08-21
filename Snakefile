"""snakefile to run PIDCalib jobs locally or on Condor"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

import shutil
import os
from snakemake.io import glob_wildcards


# full report of all the efficiencies
rule all:
    input:
        "pidcalib_full_run.done",
        "discretizer.done",
        "templates.done",
        "consistency_checks.done",

onsuccess:
    shutil.rmtree(".snakemake/metadata")


# RFE: test runs 

# checkpoint for unknown number of reco-dir/sh files to run pidcalib2; stored in directory bin by default 
checkpoint gen_sh_files:
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
    input:
        bash_file = "bin/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/run.sh",
    output:
        pkl_file = "bin/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/perf.pkl",
    params:
        output_dir = "logs/pidcalib/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}", # FIXEME: need to adjust this path + maybe make the same name for all hists.pkl
    log:
        log_file = "logs/pidcalib/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/pidcalib2.make_eff_hists.log",
    shell:
        "bash {input.bash_file}"


def aggregate_pidcalib(wildcards):
    '''
    Collect the variable number of sub-directories produced by setup.py that live in /exec/,
    and expand to use perfHists.root files produced by PIDCalib
    '''
    checkpoint_output = checkpoints.gen_sh_files.get(**wildcards).output[0]
    year, magpol, dllmu_bin, species, eff_dirs = glob_wildcards(os.path.join(checkpoint_output, '{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/run.sh'))

    return expand(checkpoint_output+"/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/perf.pkl", zip, year=year, magpol=magpol, dllmu_bin=dllmu_bin, species=species, eff_dirs=eff_dirs)


rule collect_pidcalib:
    input:
        aggregate_pidcalib,
    output:
        combined = "pidcalib_full_run.done",
    shell:
        "echo {input} > {output.combined}"


checkpoint discretize_data:
    input:
        data_path = "99_2021_01_443970_443970746_Bc2D0MuNuXSlim.root", # TODO: retrieve from config file
    params:
        script = "ddmisid/data_discretizer.py",
    output:
        directory("obs/"),
    shell:
        "python {params.script} {input.data_path}"


def aggregate_discretizer(wildcards):
    checkpoint_output = checkpoints.discretize_data.get(**wildcards).output[0]
    p, eta, ntracks = glob_wildcards(os.path.join(checkpoint_output, '{p}/{eta}/{ntracks}/obs.pkl'))
    return expand(checkpoint_output+"/{p}/{eta}/{ntracks}/obs.pkl", zip, p=p, eta=eta, ntracks=ntracks)


rule collect_discretizer:
    input:
        aggregate_discretizer,
    output:
        combined = "discretizer.done",
    shell:
        "echo {input} > {output.combined}"


checkpoint make_templates:
    input:
        pidcalib_done = "pidcalib_full_run.done",
    params:
        script = "ddmisid/make_templates.py",
        path_prefix = r"bin/2018/down/antimu_id", # TODO: retrieve from config, and deal with magpol
    output: 
        directory("templates/"),
    shell:
        "python {params.script} {params.path_prefix}" 


def aggregate_templates(wildcards):
    checkpoint_output = checkpoints.make_templates.get(**wildcards).output[0]
    p, eta, ntracks, species = glob_wildcards(os.path.join(checkpoint_output, '{p}/{eta}/{ntracks}/{species}.pkl'))
    return expand(checkpoint_output+"/{p}/{eta}/{ntracks}/{species}.pkl", zip, p=p, eta=eta, ntracks=ntracks, species=species)


rule collect_templates:
    input:
        aggregate_templates,
    output:
        combined = "templates.done",
    shell:
        "echo {input} > {output.combined}"


checkpoint consistency_check:
    input:
        pidcalib_done = "pidcalib_full_run.done",
        templates_done = "templates.done",
    params:
        script = "ddmisid/consistency_checker.py",
        eff_hist_path_prefix = r"bin/2018/down/antimu_id",
        template_path_prefix = r"templates",
    output:
        directory("checks/"),
    shell:
        "python {params.script} {params.eff_hist_path_prefix} {params.template_path_prefix}"


def aggregate_checks(wildcards):
    checkpoint_output = checkpoints.consistency_check.get(**wildcards).output[0]
    species = glob_wildcards(os.path.join(checkpoint_output, '{species}.{ext}'))
    return expand(checkpoint_output+"/{species}.{ext}", zip, species=species, ext=ext)


rule collect_checks:
    input:
        aggregate_checks,
    output:
        combined = "consistency_checks.done",
    shell:
        "echo {input} > {output.combined}"


# rule run_bml_fit: # ADDITION OF THIS CHECKPOINT WILL REQUIRE MODIFICATION OF AGGREGATE_TEMPLATE_WILDCARDS FUNCTION, SIMILAR TO PIDCALIB
#     # additionally, collect_templates will also be altered by this so that bml fit and then results are aggregated
#     input:
#         script = "ddmisid/bml_fitter.py",
#         templates = "templates/{p}/{eta}/{ntracks}/{species}.results.json"
#     output:
#         pass
#     shell:
#         pass


# rule aggregate_weight_components:
#     input:

# data = lambda wildcards: [
#             data_storage+"/{stream}/{filetype}/merged_magpols/{channel}/{year}/{mode}.root".\
#             format(stream=wildcards.stream, filetype=wildcards.filetype, channel=wildcards.channel, year=wildcards.year, mode=<>e)
#         ] if wildcards.filetype!="all" else [],