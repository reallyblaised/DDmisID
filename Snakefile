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
import pyfiglet

# user-defined config, accessible in-scope via `config[key]`
configfile: "config/main.yml"

# cleanup
onsuccess: 
    shutil.rmtree(".snakemake/metadata")

rule all:
    """Target of the entire DDmisID workflow"""
    input:
        "ddmisid.done"
    run: 
        print(pyfiglet.figlet_format("misID", font = "isometric1" ))
        print(pyfiglet.figlet_format("done", font = "isometric1" ))



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
        test = config['test'],
        verbose = config['verbose']
    output:
        dir = directory("bin/"),
    params: 
        user_id = config['user_id']
    priority: 100
    run:
        command = "python {input.config_executable} --pid_regime reco --pid_regime he_all" 
        
        # run [nominal, test, verbose]
        if params.verbose:
            print(">>> WARNING: running PIDCalib2 in `verbose` mode")
            command+=" --verbose"
        if params.test: 
            print(">>> WARNING: running PIDCalib2 in `test` mode")
            command+=" --test"

        shell(command)


rule run_pidcalib:
    """
    Run the PIDCalib2 jobs as spawned in the previous step; store the efficiency histograms in /bin/
    """
    input:
        bash_file = "bin/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/run.sh",
    output:
        pidcalib_hists = "bin/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/perf.pkl",
    log: 
        "log/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/pidcalib.log"
    run:
        shell("bash {input.bash_file}") #" &> {log}")


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
        data_path = config["data"]["input_path"] 
    output:
        directory("obs/"),
    log: 
       "log/data_discretizer.log" 
    shell:
        "python ddmisid/data_discretizer.py {input.data_path} &> {log}"


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


# ================================================================================================
#                               Binned Maximum Likelihood Section
# ------------------------------------------------------------------------------------------------
#   - build hadron-enriched templates for each species across all reco partitions
#   - binned ML fit to the hadron-enriched data in bins of kinematics and occupancy
#   - extract the true abundance of each species, accounting for cross-contamination 
#   - build a histogram of the N_i/N_ref, per species, across all bins of kinematics and occupancy
# ================================================================================================
checkpoint make_templates:
    input:
        pidcalib_done = "pid_efficiency_maps.done", # NOTE: required to build the templates
        data_partitions_done = "data_partitions.done" # avoid checkpoint conflict with data discretisation
    params:
        executable = "ddmisid/make_templates.py",
        path_prefix = r"bin/2018/up/antimu_id", # TODO: retrieve from config, and deal with magpol
    output: 
        directory("templates/"),
    shell:
        "python {params.executable} {params.path_prefix}" 


rule bml_fit: 
    input:
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
    params:
        executable = "ddmisid/fit_reco.py"
    log: 
        "log/{p}/{eta}/{ntracks}/bml_fit.log"
    output:
        # true abundance of each species, accounting for cross-contamination between reco bins
        "postfit/{p}/{eta}/{ntracks}/yields.pkl",
    run:
        shell("python {params.executable}\
            --pion {input.pion_template}\
            --kaon {input.kaon_template}\
            --proton {input.proton_template}\
            --electron {input.electron_template}\
            --obs {input.obs}" #&> {log}" # FIXME: ghosts missing 
        )


def hadron_enriched_abundances(wildcards):
    """Book the desired process outcome out of the checkpoint, so to be able to use wildcards to enact the fit"""
    checkpoint_output = checkpoints.make_templates.get(**wildcards).output[0]
    p = list(set(glob_wildcards(os.path.join(checkpoint_output, '{p}/{eta}/{ntracks}/{species}.pkl')).p))
    eta = list(set(glob_wildcards(os.path.join(checkpoint_output, '{p}/{eta}/{ntracks}/{species}.pkl')).eta))
    ntracks = list(set(glob_wildcards(os.path.join(checkpoint_output, '{p}/{eta}/{ntracks}/{species}.pkl')).ntracks))
    
    return expand("postfit/{p}/{eta}/{ntracks}/yields.pkl", p=p, eta=eta, ntracks=ntracks)


rule build_relative_abundance_lookup_table:
    """Build a histogram of the N_i/N_ref, per species, across all bins of kinematics and occupancy"""
    input:
        hadron_enriched_abundances
    params: 
        executable = "ddmisid/collect_yields.py"
    output:
        "postfit/rel_abundances_hist.pkl", # global histogram; a lookup table in {species, kinematics, occupancy} -> N_i/N_ref
    run:
        shell("python {params.executable} --input {input} --output {output}")


# =====================================================================================================================
#                                       Extrapolation to Signal Region Section
# ---------------------------------------------------------------------------------------------------------------------
#   - Source the relative true abundance of each species in the control sample  
#   - Extract misID weight through the linear combination of the above and eff(h->mu), unfolding eff(!isMuon & PIDmu<0)
# =====================================================================================================================
def fetch_pid_eff_sp_id(_dllmu_bin, _species, partition=None, _magpol="up"): # FIXME: expand to both magpol
    def fetch_eff(wildcards):
        '''
        Collect the variable number of sub-directories produced by setup.py that live in /exec/,
        and expand to use perfHists.root files produced by PIDCalib2
        '''
        checkpoint_output = checkpoints.gen_sh_files.get(**wildcards).output[0]
        
        year, eff_dirs = glob_wildcards(os.path.join(checkpoint_output, '{year}/'+_magpol+'/'+_dllmu_bin+'/'+_species+'/{eff_dirs}/run.sh'))
        
        # care in handling the inclusive !muon species-specific efficiency histogram
        match partition:
            case None:
                return expand(checkpoint_output+"/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/perf.pkl", zip, year=year, magpol=_magpol, dllmu_bin=_dllmu_bin, species=_species, eff_dirs=eff_dirs)
            case _: 
                return expand(checkpoint_output+"/{year}/{magpol}/{dllmu_bin}/{species}/{eff_dirs}/perf.pkl", zip, year=year, magpol=_magpol, dllmu_bin=_dllmu_bin, species=_species, eff_dirs=partition) # user-defined `eff_dirs`
    
    return fetch_eff


# for magpol in ["up", "down"]: FIXME: need to expand to include magpol variety
rule extract_misid_weights: 
    """
    In each bin of kinematics and occupancy, perform the weighted sum of the c_i * eff(h_i->mu)/eff(!isMuon & PIDmu<0), 
    where c_i is the true abundance of each species in that bin of hadron-enriched control data.
    That is, unfold the efficiency related to isolating the hadron-enriched control sample as a whole, and fold back in the 
    species-specific misID efficiency driving the h->mu pollution in the signal channel.
    """
    input:
        # $N_i/N_{\mathrm{ref}}$ for each species 
        # ---------------------------------------
        relative_true_abundances = "postfit/rel_abundances_hist.pkl",
        # eff(h->mu) for each species
        # -----------------------------
        pion_to_mu_eff = fetch_pid_eff_sp_id("mu_id", "pion"),
        electron_to_mu_eff = fetch_pid_eff_sp_id("mu_id", "electron"),
        kaon_to_mu_eff = fetch_pid_eff_sp_id("mu_id", "kaon"),
        proton_to_mu_eff = fetch_pid_eff_sp_id("mu_id", "proton"),
        # FIXME: ghosts
        # species-specific PID efficiencies for h \in hadron-enriched-dataset [no partitions required]
        # --------------------------------------------------------------------------------------------
        proton_to_antimu_eff = fetch_pid_eff_sp_id("antimu_id", "proton", "all"), # all signifies in {!isMuon & DLLmu<0}
        kaon_to_antimu_eff = fetch_pid_eff_sp_id("antimu_id", "kaon", "all"),
        pion_to_antimu_eff = fetch_pid_eff_sp_id("antimu_id", "pion", "all"),
        electron_to_antimu_eff = fetch_pid_eff_sp_id("antimu_id", "electron", "all"),
        # FIXME: global_antimuon_eff_ghosts = fetch_pid_eff_sp_id("antimu_id", "proton", "all"),
        # hadron-enriched data, to which we want to assign weights
        # --------------------------------------------------------
        hadron_enriched_data = config["data"]["input_path"]
    params:
        executable = "ddmisid/assign_w_misid.py",
        root_key = config['data']['root_config']['root_key'],
        root_tree_name = config['data']['root_config']['root_tree_name'],
        outfile_path = config['data']['output_path']
    output:
        temp("ddmisid.done") # signal completion
    run:
        shell("python {params.executable} --rel_abundances {input.relative_true_abundances} --obs {input.hadron_enriched_data}\
        --proton_to_mu {input.proton_to_mu_eff} --pion_to_mu {input.pion_to_mu_eff} --kaon_to_mu {input.kaon_to_mu_eff} --electron_to_mu {input.electron_to_mu_eff}\
        --proton_to_antimu {input.proton_to_antimu_eff} --pion_to_antimu {input.pion_to_antimu_eff} --kaon_to_antimu {input.kaon_to_antimu_eff} --electron_to_antimu {input.electron_to_antimu_eff}\
        --output {params.outfile_path} --key {params.root_key} --tree {params.root_tree_name} \
        && touch {output}") # temporary signal file only if main executable has run successfully