"""snakefile to run PIDCalib jobs locally or on Condor"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

import shutil
import os
from snakemake.io import glob_wildcards

# full report of all the efficiencies
rule all:
    input:
        "pidcalib_full_run.done"

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
        # bash_file = "bin/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/run.sh"
    shell:
        "python {input.config_executable} {params.is_test}"


# run pidcalib2 on each of the sub-directories produced by setup.py
rule run_pidcalib:
    input:
        bash_file = "bin/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/run.sh"
    output:
        pkl_file = "bin/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/perf.pkl",
    params:
        output_dir = "logs/pidcalib/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}", # FIXEME: need to adjust this path + maybe make the same name for all hists.pkl
    log:
        log_file = "logs/pidcalib/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/pidcalib2.make_eff_hists.log",
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