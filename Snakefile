"""snakefile to run PIDCalib jobs locally or on Condor"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

import shutil

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
        config_executable = "setup_pidcalib_jobs.py",
    params:
        is_test = f"--test"
    output:
        directory("bin")
    shell:
        "python {input.config_executable} {params.is_test}"

# rule run_pidcalib:
#     input:        
#         "exec/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/run.sh"
#     output:
#         protected("exec/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/perfHist.root")
#     log:
#         "logs/{year}_{magpol}_{dllmu_bin}_{species}_{eff_dirs}.log"
#     shell:
#         '''
#         sh {input} &> {log}
#         '''

# rule tally_effs:
#     input:        
#         "exec/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/perfHist.root"
#     output:
#         "exec/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/eff.tex"
#     log:
#         "logs/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}_eff.log"
#     params:
#         yml_file = "exec/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/preserv.yaml"    
#     shell:
#         "python /usera/delaney/private/Bc2D0MuNuX/MisID/python/read_calib_file.py --input {input} --preserv {params.yml_file} &> {log}"

def aggregate(wildcards):
    '''
    Collect the variable number of sub-directories produced by setup.py that live in /exec/,
    and expand to use perfHists.root files produced by PIDCalib
    '''
    checkpoint_output = checkpoints.gen_sh_files.get(**wildcards).output[0]
    year, magpol, dllmu_bin, species, eff_dirs = glob_wildcards(os.path.join(checkpoint_output, '{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/run.sh'))
    return expand(checkpoint_output+"/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/run.sh", zip, year=year, magpol=magpol, dllmu_bin=dllmu_bin, species=species, eff_dirs=eff_dirs)


rule collect:
    input:
        aggregate
    output:
        combined = "pidcalib_full_run.done",
    shell:
        '''
        echo {input} > {output.combined}
        '''