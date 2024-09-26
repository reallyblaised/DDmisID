<div style="text-align:left;">
  <img src="assets/DDmisID_logo.png" alt="Example postfit plot" style="width:40%; max-width:40%;">
</div>

# DDMisID

Python library for **D**ata**D**riven**MisID** modelling. This derives a single-track misID template and event yield from a control region in data. Such a task is executed via the assignment of per-event weights obtained from particle identification efficiency maps. In turn, these are (mostly) produced from bespoke high-purity, high-statistics calibration data samples via [PIDCalib2](https://pypi.org/project/pidcalib2/).

## Methodology

$$w_{\mathrm{misID}} = \sum_{i\in\{p,K,\pi,e,g\}} \frac{N_i}{N_{\mathrm{ref}}} \frac{1}{\varepsilon_{\mathrm{PID}}{(i \to !\mu)}}\varepsilon_{\mathrm{PID}}(i\to \mu)$$

where the true abundance of each species in the reference control sample, $N_i$, is given by unfolding the observed abundance of each species in high-purity partitions of the control sample, accounting for cross-contamination due to imperfect particle identification (PID):

$$N_i^{\mathrm{obs}} = \sum_i \sum_j N_i^{\mathrm{}}\varepsilon(i\to j)$$

Here, $j$ indexes the high-purity partitions in the reference sample, and $i \in \{p, K, \pi, e, g\}$ (with $g$ denoting _ghosts_ in the LHCb reconstruction jargon).

The unfolding is executed by means of binned maximum-likelihood fits within. In turn, each fit is executed in bins of kinematics and occupaancy, to account for the variation of PID responses with momentum, pseudo-rapidity, and detector occupancy. 

As an example, `DDmisID` extracts the true abundance of each species, in each bin of kinematics and occupancy, as yields extracted in fits such as this one:

<div style="text-align:center;">
  <img src="assets/ddmisid_postfit.png" alt="Example postfit plot" style="width:100%; max-width:100%;">
</div>

*Binned maximum likelihood fit to orthogonal, high-purity partitions of the hadron-enriched data. The filled coloured histograms illustrate the post-fit extracted abundance of each species, accounting of cross-contamination between the partitions due to imperfect PID. Generated with in-house pseudo-data mimicking the LHCb reconstruction*.

## Installation

### 1. Clone the Repository 
```bash
git clone git@github.com:reallyblaised/DDmisID.git
cd DDmisID
```
### 2. Install Dependencies

For regular use: 

```bash
pip install -r requirements.txt
```

For developement (editable mode, and recommended until further notice):

```bash
pip install -e . 
```

## Running DDmisID

DDmisID runs through an engine that 
  1. Builds and validates a user-specified configuration YAML file [by default `config/main.yml`].
  2. Orchestrates the DDmisID pipeline, powered by a [Snakemake](https://snakemake.readthedocs.io/en/stable/#) backend.

### Instructions
1. #### Edit YAML main configuration file
Open and modify `config/main.yml` (or a custom path to a YAML file with the same key structure), following the in-inline field descriptios.

2. #### Build the DDmisID engine
Parse the YAML configuration file and ensure it complies with the DDmisID engine specification:
```bash
$ ddmisid-engine build # sources config/main.yml by default
```
or, for custom config-file locations: 
```bash 
$ ddmisid-engine build --config-path=<custom_YAML_config_path>
```
3. #### Verify the engine spec has been built correctly
Run 
```bash
$ python -c "from ddmisid.engine import config; print(config)"
```
to verify that the configuration report is compatible with the user's directives.

4. #### Run the DDmisID engine
Execute the pipeline, passing Snakemake’s dynamic flags as needed, for example:
```bash
ddmisid-engine run -- <snakemake options>
```
for example,
```
ddmisid-engine run -- --dryrun --cores 1
```
where, as a tip, `--cores all` can be used to exploit all the resources available on your machine or cluster.

## License

This project is licensed under the [MIT License](https://opensource.org/license/MIT). See [LICENSE](LICENSE) for more information.

![License](https://img.shields.io/badge/license-MIT-blue.svg)
