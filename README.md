# nextstrain.org/lassa dmsa_pred build

This repository is forked from [nextstrain.org/lassa](https://nextstrain.org/lassa) and modified to include the analysis of deep mutational scanning data. The analysis is performed using the [dmsa_pred](https://github.com/matsengrp/dmsa_pred) package. This repo is currently set up to run the data as presented in the manuscript [Deep mutational scanning reveals functional constraints and antigenic variability of Lassa virus glycoprotein complex](https://www.biorxiv.org/content/10.1101/2024.02.05.579020v1).

To run the pipeline, 
install snakemake following their 
[documentation](https://snakemake.readthedocs.io/en/v8.4.11/getting_started/installation.html), 
This pipeline has been tested with versions >= 8.4.11 

then clone this repository and the submodules
```
git clone https://github.com/matsengrp/lassa-dmsa.git --recurse-submodules
```
Run the pipeline with the following command
```
snakemake --use-conda --cores 2 --configfile my_profiles/configs/manuscript.yaml
```
Each rule currently uses the environment specified in 
[my_profiles/dmsa_pred/dmsa_env.yaml](my_profiles/dmsa_pred/dmsa_env.yaml).
The configuration file specified above provides paths to the 
escape data, reference sequence, annotation file, sequences, and metadata within the 
[LASV_Josiah_GP_DMS](https://github.com/dms-vep/LASV_Josiah_GP_DMS.git) submodule. 
The pipeline will output a JSON file(s) under the "auspice/" directory
that can be visualized with [auspice](https://auspice.us/) software.
CSV's of phenotype predictions can be found under "results/dmsa-phenotype/"


## Parent repository README


[![Build Status](https://github.com/nextstrain/lassa/actions/workflows/ci.yaml/badge.svg?branch=master)](https://github.com/nextstrain/lassa/actions/workflows/ci.yaml)

This is the [Nextstrain](https://nextstrain.org) build for Lassa, visible at
[nextstrain.org/lassa](https://nextstrain.org/lassa).

The build encompasses fetching data, preparing it for analysis, doing quality
control, performing analyses, and saving the results in a format suitable for
visualization (with [auspice][]).  This involves running components of
Nextstrain such as [fauna][] and [augur][].

All Lassa-specific steps and functionality for the Nextstrain pipeline should be
housed in this repository.


## Usage

If you're unfamiliar with Nextstrain builds, you may want to follow our
[quickstart guide][] first and then come back here.

The easiest way to run this pathogen build is using the [Nextstrain
command-line tool][nextstrain-cli]:

    nextstrain build .

See the [nextstrain-cli README][] for how to install the `nextstrain` command.

Alternatively, you should be able to run the build using `snakemake` within an
suitably-configured local environment.  Details of setting that up are not yet
well-documented, but will be in the future.

Build output goes into the directories `data/`, `results/` and `auspice/`.

Once you've run the build, you can view the results in auspice:

    nextstrain view auspice/


## Configuration

Configuration takes place entirely with the `Snakefile`. This can be read top-to-bottom, each rule
specifies its file inputs and output and also its parameters. There is little redirection and each
rule should be able to be reasoned with on its own.


### fauna / RethinkDB credentials

This build starts by pulling sequences from our live [fauna][] database (a RethinkDB instance). This
requires environment variables `RETHINK_HOST` and `RETHINK_AUTH_KEY` to be set.

If you don't have access to our database, you can run the build using the
example data provided in this repository.  Before running the build, copy the
example sequences into the `data/` directory like so:

    mkdir -p data/
    cp example_data/lassa_*.fasta data/


[Nextstrain]: https://nextstrain.org
[fauna]: https://github.com/nextstrain/fauna
[augur]: https://github.com/nextstrain/augur
[auspice]: https://github.com/nextstrain/auspice
[snakemake cli]: https://snakemake.readthedocs.io/en/stable/executable.html#all-options
[nextstrain-cli]: https://github.com/nextstrain/cli
[nextstrain-cli README]: https://github.com/nextstrain/cli/blob/master/README.md
[quickstart guide]: https://nextstrain.org/docs/getting-started/quickstart
