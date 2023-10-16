# RGS (Replicator-Gene Symbiosis)

R implementation of the replication-prone genomic region detection method described in the BioRxiv preprint:

**Replicators, genes, and the C-value enigma: High-quality genome assembly of barley provides direct evidence that self-replicating DNA forms 'cooperative' associations with genes in arms races**

M. Timothy Rabanus-Wallace, Thomas Wicker, and Nils Stein

doi: https://doi.org/10.1101/2023.10.01.560391

## System Requirements

These scripts have been tested running on a Red Hat Enterprise Linux v9.2 OS, R version 4.2.2, via Rstudio Server 2022.07.2+576. R running on almost any linux OS should be fine.

## Installation

These scripts are run interactively in an R session. Installation is as simple as copying the github repo using (e.g.) `git clone https://github.com/mtrw/RGS`, then open the scripts in an Rstudio session.

The scripts call lastz (https://www.bx.psu.edu/~rsharris/lastz/; tested on lastz v1.04.03), you will need a version running on your system.

The scripts rely on the functions and packages of https://github.com/mtrw/tim_r_functions, which are loaded in the second line of ldprIdent_RUN.R.

To install the necessary package dependencies easily, after `source`ing tim_r_functions, just run `install_load_Tims_packages()`.

## Running

1) Set up your desired parameters in 0_ldprIdent_setup.R. Importantly, your reference genome file, and the file to which Long-duplication-prone regions flagged will be saved (in .rds format). You will have to provide a path to the lastz binary in 0_ldprIdent_setup.R in the variable `lastzBin`.

2) Run the steps in ldprIdent_RUN.R. This will call first 0_ldprIdent_setup.R, then 1_llprIdent_selfAlignments.R, then 2_0_ldprIdent_setup.R.

To get tricky with the details, dig into the scripts that are called.

## Demo

A demo dataset to run the scripts on is available in data/test_data. This is the first ~6 Mb of the Barley cv 'Morex' reference genome used in the release publication. The scripts left as-is will run on the demo dataset, the user only has to set the working directory to wherever the repo was cloned in the first line of ldprIdent_RUN.R, and the lastz binary.

## Performance

The run time of these scripts is almost entirely dependent on the complexity of the lastz alignment. On an allocation of 16 CPUs and 48 Gb RAM, the demo ran in around 30 minutes.

To do a full genome, it is highly recommended to use parallel processing with `mclappy`. The lines of 1_llprIdent_selfAlignments.R to [un]comment to achieve this are indicated in comments, around line 6.