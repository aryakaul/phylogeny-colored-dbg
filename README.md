# Phylogeny-Colored de Bruijn Graphs

<p>
Workflow for building phylogeny-colored de Bruijn graphs (pcDBGs)
from microbial genomes.

pcDBG first generates the classical compacted colored-de Bruijn graph with
`Cuttlefish` and then calculates the parsimonious explanation of each unitig
given a phylogeny of the user-provided samples.

<!--For more information, see the <a href="LOLSOON">associated paper</a>.-->
</p><br/>

<!--[![Paper DOI](ZENODOLINK)](DOI4PAPER)-->
<!--[![MOF-Compress test](https://github.com/karel-brinda/mof-compress/actions/workflows/main.yaml/badge.svg)](https://github.com/karel-brinda/mof-compress/actions/)-->
[![GitHub release](https://img.shields.io/github/release/aryakaul/phylogeny-colored-dbg.svg)](https://github.com/aryakaul/phylogeny-colored-dbg/releases/)
[![Tests](https://github.com/aryakaul/phylogeny-colored-dbg/actions/workflows/main.yaml/badge.svg)](https://github.com/aryakaul/phylogeny-colored-dbg/actions/workflows/main.yaml)
<h2>Contents</h2>

<!-- vim-markdown-toc GFM -->

* [1. Introduction](#1-introduction)
* [2. Dependencies](#2-dependencies)
* [3. Installation](#3-installation)
* [4. Usage](#4-usage)
    * [4a. Basic example](#4a-basic-example)
    * [4b. Adjusting configuration](#4b-adjusting-configuration)
* [5. Citation](#5-citation)
* [6. Issues](#6-issues)
* [7. Changelog](#7-changelog)
* [8. License](#8-license)
* [9. Contacts](#9-contacts)
* [10. Acknowledgements](#10-acknowledgements)

<!-- vim-markdown-toc -->


## 1. Introduction

The user provides files of files for individual batches in the `input/` directory and specifies the requested k-mer length in the
[configuration file](config.yaml). By running `make`, the pipeline builds a phylogeny-colored de Bruijn graph for each list of files
in `input/`. The GFA file is outputted in `output/` and will have the same basename as the input list of files file.

In addition, the minimal cuts corresponding to each unitig in the compacted de Bruijn graph will be outputted in `output/`.


## 2. Dependencies

* [Conda](https://docs.conda.io/en/latest/miniconda.html) (unless the use of Conda is switched off in the configuration) and ideally also [Mamba](https://mamba.readthedocs.io/) (>= 0.20.0)
* [GNU Make](https://www.gnu.org/software/make/)
* [Python](https://www.python.org/) (>=3.7)
* [Snakemake](https://snakemake.github.io) (>=6.2.0)

These can be installed by Conda by
```bash
bash conda install -c conda-forge -c bioconda -c defaults \
  make "python>=3.7" "snakemake>=6.2.0" "mamba>=0.20.0"
```

Other dependencies are installed automatically by
Snakemake when they are requested. The specifications of individual environments can be found in [`workflow/envs/`](workflow/envs/),
and they contain:
- [ETE 3](http://etetoolkit.org/),
- [Pandas](https://pandas.pydata.org/),
- [Cuttlefish](https://github.com/COMBINE-lab/cuttlefish),
- [BubbleGun](https://github.com/fawaz-dabbaghieh/bubble_gun),


All dependencies across all protocols can also be
installed at once by `make conda`.


## 3. Installation

Clone and enter the repository by

```bash
git clone https://github.com/aryakaul/phylogeny-colored-dbg
cd phylogeny-colored-dbg
```

Alternatively, the repository can also be installed using cURL by
```bash
mkdir phylogeny-colored-dbg
cd phylogeny-colored-dbg
curl -L https://github.com/aryakaul/phylogeny-colored-dbg/tarball/main \
    | tar xvf - --strip-components=1
```


## 4. Usage

### 4a. Basic example

* ***Step 1: Provide lists of input files.*** \
  For every batch, create a txt list of input files in the `input/`
  directory (i.e., as `input/{batch_name}.txt`. Use either absolute paths (recommended),
  or paths relative to the root of the Github repository (not relative to the txt files).

  Such a list can be generated, for instance, by `find` by
  ```bash
  find ~/dir_with_my_genomes -name '*.fa' > input/my_first_batch.txt
  ```
  The supported input file formats include FASTA and FASTQ (possibly compressed by GZip).

* ***Step 2: Provide corresponding phylogenies.*** \
  The tree files should be named `input/{batch_name}.nwk`,
  and the leave names inside should correspond
  to FASTA filenames (without FASTA suffixes).

* ***Step 3 (optional): Adjust configuration.*** \
  By editing [`config.yaml`](config.yaml) it is possible to specify
  value of `k` and other parameters.

* ***Step 4: Run the pipeline.*** \
  Run the pipeline by `make`; this is run by
  Snakemake with the corresponding parameters.

* ***Step 5: Retrieve the output files.*** \
  All output files will be located in `output/`.
  

### 4b. Adjusting configuration

The workflow can be configured via the [`config.yaml`](./config.yaml) file, and
all options are documented directly there. The configurable functionality includes:
* switching off Conda,
* *k* for de Bruijn graph 


### 4c. List of workflow commands

phylogeny-colored-dbg is executed via [GNU Make](https://www.gnu.org/software/make/), which handles all parameters and passes them to Snakemake.
Here's a list of all implemented commands (to be executed as `make {command}`):


```yaml
######################
## General commands ##
######################
    all                  Run everything
    help                 Print help messages
    conda                Create the conda environments
    clean                Clean all output archives and files with statistics
    cleanall             Clean everything but Conda, Snakemake, and input files
    cleanallall          Clean completely everything
###############
## Reporting ##
###############
    viewconf             View configuration without comments
    reports              Create html report
####################
## For developers ##
####################
    test                 Run the workflow on test data
    format               Reformat all source code
    checkformat          Check source code format
```


### 4d. Troubleshooting

Tests can be run by `make test`.


## 5. Citation

TODO


## 6. Issues

Please use [Github issues](https://github.com/aryakaul/phylogeny-colored-dbg/issues).



## 7. Changelog

See [Releases](https://github.com/aryakaul/phylogeny-colored-dbg/releases).



## 8. License

[GPL3](https://github.com/aryakaul/phylogeny-colored-dbg/blob/main/LICENSE)



## 9. Contacts

* [Arya Kaul](https://arya.casa) \<arya_kaul@g.harvard.edu\>
* [Karel Brinda](http://karel-brinda.github.io) \<karel.brinda@inria.fr\>


## 10. Acknowledgements

Structure and format for this pipeline, documentation, and some code was heavily inspired 
and modeled after [Miniphy](https://github.com/karel-brinda/Miniphy)! Check it out!
