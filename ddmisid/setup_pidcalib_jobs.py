"""Generate PIDCalib jobs for the DDmisID project.

This script generates a bash file which is mostly a wrapper for pidcalib2, 
generated from the user-defined config file.
"""

from ddmisid import BinningGenerator
from ddmisid.utils import read_config, debug, timing
from data.aliases import CalibSamples, MCTunings, BinningVars, CommonCuts
from termcolor2 import c as tc
from pathlib import Path
import argparse
import os
import functools
from typing import Callable, TypeVar
from typing_extensions import ParamSpec

T = TypeVar("T")
P = ParamSpec("P")


def check_mu_region(
    func: Callable[P, T], valid_ids: "list[str]" = ["antimu_id", "mu_id"]
) -> Callable[P, T]:
    """Decorator to check whether muon id is valid"""

    @functools.wraps(func)
    def id_checker(*args: P.args, **kwargs: P.kwargs) -> T:
        id = kwargs["region_id"]
        assert id in valid_ids, f"ValueError: region_id must be in {valid_ids}"
        return func(*args, **kwargs)

    return id_checker


def match_muid_criteria(
    region_id: str,
    pid_config: dict,
) -> "list[str]":
    """Establish the pid categories from which to extract efficiencies via pidcalib2

    Parameters
    ----------
    region_id: str
        The region identifier (hadron-enriched or signal). Allowed values: "antimu_id", "mu_id"

    pid_config: dict
        The pid configuration dictionary, user-defined

    Returns
    -------
    list[str]:
        The pid category [{k, pi, p, e}_like, muon_like]
    """
    if (
        region_id == "antimu_id"
    ):  # partition the hadron-enriched sample into reco categories
        return list(pid_config["reco_cuts"].keys())
    elif region_id == "mu_id":  # signal selection
        return ["muon_like"]
    else:
        raise ValueError("Invalid region_id; must be either 'antimu_id' or 'mu_id'")


@timing
@check_mu_region  # check whether the region_id is valid
def generate_jobs(
    pid_config: dict,
    region_id: str,  # allowed values: "antimu_id", "mu_id"
    parent_outdir: str = "bin",
    verbose: bool = False,
    test: bool = True,
    summarise: bool = True,
) -> None:
    """Generate the executable for pidcalib2 jobs.

    Parameters
    ----------
    pid_config : dict
        Read-in user-defined config

    region_id : str
        The region identifier (hadron-enriched or signal). Allowed values: "antimu_id", "mu_id"

    verbose: bool

    Returns
    -------
    None
        Generates the appropriate bash file for pidcalib2 jobs
    """
    # establish the reco criteria: whether reco partitions of HE, or signal
    reco_set = match_muid_criteria(region_id, pid_config)

    # having read the config, write a suitable json file
    for y in pid_config["years"]:

        # generate the all binning files, per year
        BinningGenerator(path="config/main.yml").build(year=y)

        for true_sp_id, true_sp_alias in pid_config[
            "species"
        ].items():  # true hadrons, ghosts, electons whose abundance must be extracted
            binning_path = Path(f"data/binning_{y}/{true_sp_alias}.json")

            for reco_sp in reco_set:  # reco partition identifiers
                for magpol in pid_config["magpols"]:
                    # differentiate between electron and hadrons for calibration samples
                    if true_sp_id != "electron":
                        CALIBRATION_SAMPLE = getattr(CalibSamples(), f"hadron_{y}")
                    if true_sp_id == "electron":
                        CALIBRATION_SAMPLE = getattr(CalibSamples(), f"e_{y}")

                    # book the binning variables, with per-year appropriate
                    BINNING_VARS = getattr(BinningVars(), f"_{y}")

                    # based on the desidered category, the pid, occupancy & kinematic selection criteria are different
                    if region_id == "antimu_id":
                        PID_SEL = f"{getattr(MCTunings(), f'_{y}')}{pid_config['reco_cuts'][reco_sp]} & {pid_config[region_id]['pid_cut']}"
                        RECO_SEL = f"{pid_config['common_sel']}"  # HACK: absorb the `--cut` directive here
                    if region_id == "mu_id":
                        PID_SEL = f"{getattr(MCTunings(), f'_{y}')}{pid_config['mu_id']['pid_cut']}"
                        RECO_SEL = f"{pid_config['common_sel']}"  # HACK: absorb the `--cut` directive here

                    # HACK: there is an exception for e 2016
                    if true_sp_id == "e_B_Jpsi" and y == "2016":
                        BINNING_VARS = getattr(
                            BinningVars(), f"e_2012"
                        )  # apparently Brunel prefix does not apply in this case

                    # produce summary
                    if summarise:
                        print(tc("Pidcalib jobs config summary:").underline.yellow)
                        print(
                            f"""\n* year: {y}\n* polarity: {magpol}\n* calib species: {true_sp_alias}\n* calibration sample: {CALIBRATION_SAMPLE }\
                            \n* binning variables: {BINNING_VARS}\n* reco selection: {RECO_SEL}\n* pid selection: {PID_SEL}\n
                        
                        """
                        )

                    # proceed with writing job executable
                    # -----------------------------------
                    # book a suitable directory
                    sp_outdir = f"{parent_outdir}/{y}/{magpol}/{region_id}/{true_sp_id}"
                    namespace = f"{true_sp_id}_to_{reco_sp}"

                    # clean up the directory to remove relic jobs
                    try:
                        os.system(f"rm -rf {parent_outdir}/{namespace}")
                        Path(f"{sp_outdir}/{namespace}").mkdir(
                            parents=True, exist_ok=True
                        )
                    except:
                        Path(f"{sp_outdir}/{namespace}").mkdir(
                            parents=True, exist_ok=True
                        )

                    # establish the pidcalib2 command and relative args
                    job_conf = f'source /cvmfs/lhcb.cern.ch/lib/LbEnv &&\nlb-conda pidcalib pidcalib2.make_eff_hists --sample {CALIBRATION_SAMPLE} --magnet {magpol} --particle {true_sp_alias} --pid-cut "{PID_SEL}" --cut "{RECO_SEL}" --binning-file {binning_path} --output-dir {sp_outdir}/{namespace}'

                    # job_conf = f'lb-conda pidcalib pidcalib2.make_eff_hists --sample {CALIBRATION_SAMPLE} --magnet {magpol} --particle {true_sp_alias} --pid-cut "{RECO_SEL}" --binning-file {binning_path} --output-dir {sp_outdir}/{namespace}'

                    for bv in BINNING_VARS:
                        job_conf += f" --bin-var {bv}"

                    # if booked test, run over one calibration file only
                    if test:
                        job_conf += " --max-files 25"  # ensure sufficient statistics

                    # verbose pidcalib output
                    if verbose:
                        job_conf += " --verbose"

                    # need to preserve the selection info to read in the eff histograms from root fike
                    preserv_file = open(
                        f"{sp_outdir}/{namespace}/preserv.yaml",
                        "w",
                    )
                    preserv_file.write(f'real     : "{true_sp_alias}"\n')
                    preserv_file.write(f'reco     : "{reco_sp}"\n')
                    preserv_file.write(f'muid     : "{region_id}"\n')
                    preserv_file.write(f'eff_hist : "{true_sp_alias}_{reco_sp}_all"')
                    preserv_file.close()

                    # bash file to run PIDCalib
                    sh_file = open(f"{sp_outdir}/{namespace}/run.sh", "w")

                    # pipe command to bashfile
                    sh_file.write(f"{job_conf} &&\n")

                    sh_file.write(
                        f"touch {sp_outdir}/{namespace}/pidcalib2.make_eff_hists.done"
                    )  # placeholder filler to signal complete execution of pidcalib2

                    # # HACK: rename so that all root tuples are called the same - helps with the snake pipeline
                    epilogue = f"""\nfor f in {sp_outdir}/{namespace}/*.pkl; do\nmv \"$f\" {sp_outdir}/{namespace}/perf.pkl\ndone
                    """

                    sh_file.write(epilogue + "\n")

                    sh_file.close()

                    # make the bash file executable
                    os.chmod(f"{sp_outdir}/{namespace}/run.sh", 0o0777)

                    # report successul run
                    if verbose:
                        print(
                            tc(
                                f"Success: generated {sp_outdir}/{namespace}/run.sh"
                            ).green
                        )


@timing
@check_mu_region
def generate_he_jobs(
    pid_config: dict,
    region_id: str,  # allowed values: "antimu_id", "mu_id"
    parent_outdir: str = "bin",
    verbose: bool = True,
    test: bool = True,
    summarise: bool = True,
) -> None:
    """Generate the executable for pidcalib2 jobs.

    Parameters
    ----------
    pid_config : dict
        Read-in user-defined config

    region_id : str
        The region identifier (hadron-enriched or signal). Allowed values: "antimu_id", "mu_id"

    verbose: bool

    Returns
    -------
    None
        Generates the appropriate bash file for pidcalib2 jobs
    """
    # establish the reco criteria: whether reco partitions of HE, or signal
    reco_set = match_muid_criteria(region_id, pid_config)

    # having read the config, write a suitable json file
    for y in pid_config["years"]:

        # generate the all binning files, per year
        BinningGenerator(path="config/main.yml").build(year=y)

        for true_sp_id, true_sp_alias in pid_config[
            "species"
        ].items():  # true hadrons, ghosts, electons whose abundance must be extracted
            binning_path = Path(f"data/binning_{y}/{true_sp_alias}.json")

            for magpol in pid_config["magpols"]:
                # differentiate between electron and hadrons for calibration samples
                if true_sp_id != "electron":
                    CALIBRATION_SAMPLE = getattr(CalibSamples(), f"hadron_{y}")
                if true_sp_id == "electron":
                    CALIBRATION_SAMPLE = getattr(CalibSamples(), f"e_{y}")

                # book the binning variables, with per-year appropriate
                BINNING_VARS = getattr(BinningVars(), f"_{y}")

                # based on the desidered category, the pid, occupancy & kinematic selection criteria are different
                if region_id == "antimu_id":
                    PID_SEL = f"{getattr(MCTunings(), f'_{y}')}{pid_config['he_all']} & {pid_config[region_id]['pid_cut']}"
                    RECO_SEL = f"{pid_config['common_sel']}"

                # HACK: there is an exception for e 2016
                if true_sp_id == "e_B_Jpsi" and y == "2016":
                    BINNING_VARS = getattr(
                        BinningVars(), f"e_2012"
                    )  # apparently Brunel prefix does not apply in this case

                # produce summary
                if summarise:
                    print(tc("Pidcalib jobs config summary:").underline.yellow)
                    print(
                        f"""\n* year: {y}\n* polarity: {magpol}\n* calib species: {true_sp_alias}\n* calibration sample: {CALIBRATION_SAMPLE }\
                        \n* binning variables: {BINNING_VARS}\n* reco selection: {RECO_SEL}\n* pid selection: {PID_SEL}\n
                    
                    """
                    )

                # proceed with writing job executable
                # -----------------------------------
                # book a suitable directory
                sp_outdir = f"{parent_outdir}/{y}/{magpol}/{region_id}/{true_sp_id}"
                namespace = f"all"

                # clean up the directory to remove relic jobs
                try:
                    os.system(f"rm -rf {parent_outdir}/{namespace}")
                    Path(f"{sp_outdir}/{namespace}").mkdir(parents=True, exist_ok=True)
                except:
                    Path(f"{sp_outdir}/{namespace}").mkdir(parents=True, exist_ok=True)

                # establish the pidcalib2 command and relative args
                job_conf = f'source /cvmfs/lhcb.cern.ch/lib/LbEnv &&\nlb-conda pidcalib pidcalib2.make_eff_hists --sample {CALIBRATION_SAMPLE} --magnet {magpol} --particle {true_sp_alias} --pid-cut "{PID_SEL}" --cut "{RECO_SEL}" --binning-file {binning_path} --output-dir {sp_outdir}/{namespace}'

                # job_conf = f'lb-conda pidcalib pidcalib2.make_eff_hists --sample {CALIBRATION_SAMPLE} --magnet {magpol} --particle {true_sp_alias} --pid-cut "{RECO_SEL}" --binning-file {binning_path} --output-dir {sp_outdir}/{namespace}'

                for bv in BINNING_VARS:
                    job_conf += f" --bin-var {bv}"

                # if booked test, run over one calibration file only
                if test:
                    job_conf += " --max-files 25"  # ensure sufficient statistics

                # verbose pidcalib output
                if verbose:
                    job_conf += " --verbose"

                # need to preserve the selection info to read in the eff histograms from root fike
                preserv_file = open(
                    f"{sp_outdir}/{namespace}/preserv.yaml",
                    "w",
                )
                preserv_file.write(f'real     : "{true_sp_alias}"\n')
                preserv_file.write(f'reco     : "all"\n')
                preserv_file.write(f'muid     : "{region_id}"\n')
                preserv_file.write(f'eff_hist : "{true_sp_alias}_all_all"')
                preserv_file.close()

                # bash file to run PIDCalib
                sh_file = open(f"{sp_outdir}/{namespace}/run.sh", "w")

                # pipe command to bashfile
                sh_file.write(f"{job_conf} &&\n")

                sh_file.write(
                    f"touch {sp_outdir}/{namespace}/pidcalib2.make_eff_hists.done"
                )  # placeholder filler to signal complete execution of pidcalib2

                # # HACK: rename so that all root tuples are called the same - helps with the snake pipeline
                epilogue = f"""\nfor f in {sp_outdir}/{namespace}/*.pkl; do\nmv \"$f\" {sp_outdir}/{namespace}/perf.pkl\ndone
                """

                sh_file.write(epilogue + "\n")

                sh_file.close()

                # make the bash file executable
                os.chmod(f"{sp_outdir}/{namespace}/run.sh", 0o0777)

                # report successul run
                if verbose:
                    print(
                        tc(f"Success: generated {sp_outdir}/{namespace}/run.sh").green
                    )


def pidcalib_jobs_selector(
    job_type: str, **kwargs
) -> None:  # execute pertinent function
    if job_type == "reco":
        for id in ("antimu_id", "mu_id"):
            generate_jobs(
                pid_config=read_config("config/main.yml", key="pid"),
                region_id=id,
                **kwargs,
            )
    elif job_type == "he_all":
        generate_he_jobs(
            pid_config=read_config("config/main.yml", key="pid"),
            region_id="antimu_id",
            **kwargs,
        )
    else:
        raise KeyError("PID mode identifier not allowed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PIDCalib jobs generator")
    parser.add_argument(
        "-t", "--test", action="store_true", help="Run pidcalib2 with --max-files==25"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Run pidcalib2 with --verbose"
    )
    parser.add_argument(
        "-p",
        "--pid_regime",
        type=str,
        choices=["reco", "he_all"],
        action="append",
        help="Run pidcalib2 in region of interect. Valid arguments: `reco`, `he_all`, with the former applying reco partitions to the sample of interest.",
    )
    args = parser.parse_args()

    # establish the credentials to access eos
    os.system(f"kinit {read_config('config/main.yml', key='user_id')}@CERN.CH")

    # generate the pidcalib2 jobs in the hadron-enriched region and extrapolation to signal
    for pr in args.pid_regime:
        pidcalib_jobs_selector(pr, test=args.test, verbose=args.verbose)
