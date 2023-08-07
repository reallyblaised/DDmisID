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
        # "templates.done",

onsuccess:
    shutil.rmtree(".snakemake/metadata")


# RFE: test runs 

# checkpoint for unknown number of reco-dir/sh files to run pidcalib2; stored in directory bin by default 
checkpoint gen_sh_files:
    input:
        config_executable = "ddmisid/setup_pidcalib_jobs.py",
    params:
        is_test = f"--test"
    output:
        directory("bin/")
    shell:
        "python {input.config_executable} {params.is_test}"


rule run_pidcalib:
    input:
        bash_file = "bin/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/run.sh"
    output:
        pkl_file = "bin/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/perf.pkl",
    shell:
        "bash {input.bash_file}"


def aggregate(wildcards):
    '''
    Collect the variable number of sub-directories produced by setup.py that live in /exec/,
    and expand to use perfHists.root files produced by PIDCalib
    '''
    checkpoint_output = checkpoints.gen_sh_files.get(**wildcards).output[0]
    year, magpol, dllmu_bin, species, eff_dirs = glob_wildcards(os.path.join(checkpoint_output, '{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/run.sh'))
    return expand(checkpoint_output+"/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/perf.pkl", zip, year=year, magpol=magpol, dllmu_bin=dllmu_bin, species=species, eff_dirs=eff_dirs)


rule collect:
    input:
        aggregate,
    output:
        combined = "pidcalib_full_run.done",
    shell:
        '''
        echo {input} > {output.combined}
        '''


checkpoint discretize_data:
    input:
        script = "ddmisid/data_discretizer.py",
        data_path = "99_2021_01_443970_443970746_Bc2D0MuNuXSlim.root" # TODO: retrieve from config file
    output:
        directory("obs/")
    shell:
        "python {input.script} {input.data_path}"


def aggregate_discretizer_wildcards(wildcards):
    checkpoint_output = checkpoints.discretize_data.get(**wildcards).output[0]
    p, eta, ntracks = glob_wildcards(os.path.join(checkpoint_output, '{p}/{eta}/{ntracks}/obs.pkl'))
    return expand(checkpoint_output+"/{p}/{eta}/{ntracks}/obs.pkl", zip, p=p, eta=eta, ntracks=ntracks)


rule collect_discretizer:
    input:
        aggregate_discretizer_wildcards,
    output:
        combined = "discretizer.done"
    shell:
        "echo {input} > {output.combined}"


checkpoint make_templates:
    input:
        pidcalib_done = "pidcalib_full_run.done",
        script = "ddmisid/make_templates.py",
    params:
        path_prefix = r"bin/2018/down/antimu_id"
    output: 
        directory("templates/")
    shell:
        "python {input.script} {params.path_prefix}" 


def aggregate_template_wildcards(wildcards):
    checkpoint_output = checkpoints.make_templates.get(**wildcards).output[0]
    p, eta, ntracks, species = glob_wildcards(os.path.join(checkpoint_output, '{p}/{eta}/{ntracks}/{species}.pkl'))
    return expand(checkpoint_output+"/{p}/{eta}/{ntracks}/{species}.pkl", zip, p=p, eta=eta, ntracks=ntracks, species=species)


rule collect_templates:
    input:
        aggregate_template_wildcards,
    output:
        combined = "templates.done"
    shell:
        "echo {input} > {output.combined}"
