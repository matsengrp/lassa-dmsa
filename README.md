# nextstrain.org/lassa dmsa_pred build

For those unfamiliar with Nextstrain, we recommend checking out the [Nextstrain documentation](https://docs.nextstrain.org/en/latest/) before reading further.

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
snakemake --use-conda --cores 2 --configfile config/config.yaml
```
Each rule currently uses the environment specified in 
[my_profiles/dmsa_pred/dmsa_env.yaml](my_profiles/dmsa_pred/dmsa_env.yaml).
The configuration file specified above provides paths to the 
escape data within the 
[LASV_Josiah_GP_DMS](https://github.com/dms-vep/LASV_Josiah_GP_DMS.git) submodule. 
The pipeline will output a JSON file(s) under the "auspice/" directory
that can be visualized with [auspice](https://auspice.us/) software.
CSV's of phenotype predictions can be found under "results/dmsa-phenotype/"

## Genbank accessions

To download the accessions, go to [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/) and click *Search by virus*. In the *Search by virus name or taxonomy* box, enter *Mammarenavirus lassaense, taxid:3052310* and hit enter. Then click the  *Download* option, select *Accession List* and *Nucleotide* options and hit *Next*. On the next page, select *Download All Records* and hit *Next*. On the next page, select *Accession with version* and click *Download*. Sequences are downloaded from the list of accessions because more information is extracted from the genbank file during the download process. The current accession list was downloaded on August 10, 2023. 


## Configuration

Configuration takes place within the `Snakefile` and the `config/config.yaml` files. The `Snakefile` can be read top-to-bottom, each rule
specifies its file inputs and output and also its parameters. There is little redirection and each
rule should be able to be reasoned with on its own. The `config/config.yaml` is important for configuring the DMSA phenotype prediction component of this build (e.g., paths to the deep mutational scanning directories for escape data). 


[Nextstrain]: https://nextstrain.org
[augur]: https://github.com/nextstrain/augur
[auspice]: https://github.com/nextstrain/auspice
[snakemake cli]: https://snakemake.readthedocs.io/en/stable/executable.html#all-options
[nextstrain-cli]: https://github.com/nextstrain/cli
[nextstrain-cli README]: https://github.com/nextstrain/cli/blob/master/README.md
[quickstart guide]: https://nextstrain.org/docs/getting-started/quickstart
