# Discussion of PC selection methods for scRNA-seq data

This repository contains some scripts to assess different methods of choosing the number of PCs to retain.
The `text` directory contains LaTeX files for the report, a compiled PDF of which can be found [here](https://jmlab-gitlab.cruk.cam.ac.uk/aaron/technical-reports/raw/master/pc-selection.pdf).
The `simulations` directory contains R scripts for performing the basic simulations:

- `functions.R`, a central R script containing definitions of useful functions for the simulations.
- `sim_gaussclust.R`, a template for simulations of clusters with Gaussian noise.
- `sim_trajectory.R`, a tempalte for simulations of trajectories between multiple nodes.
- `submitter.sh`, a Bash script for SLURM job submission of the simulations.
- `plot_results.R`, an R script to generate the plots.
- `simulate_noise.R`, an R script examining the effect of removing biological noise.

The `real` directory contains R scripts for performing the real data-based simulations:

- `proc_kolod.R`, an R script for pre-processing the mESC data set.
- `proc_pbmc4k.R`, an R script for pre-processing the PBMC data set.
- `run_kolod.R`, a template for performing simulations based on the mESC data set.
- `run_pbmc4k.R`, a template for performing simulations based on the PBMC data set.
- `submitter.sh`, a Bash script for SLURM job submission of the simulations.
- `plot_results.R`, an R script to generate the plots.

In addition, `batching/batching.Rmd` contains an example of how batch removal in the presence of zeroes can distort the PCA results.
