#################################################################
# This script submits all simulation jobs to the SLURM scheduler.
# The use of separate files makes it trivial to set up multiple jobs on different cores.
# This allows timely execution of all of the scripts, yet with easily altered parameters.
#################################################################

mkdir -p logs
R=/Users/lun01/software/R/devel/bin/R

# Submitting the PBMC 4K simulations.

sbatch << EOT
#!/bin/bash
#SBATCH -o logs/out-pbmc4k
#SBATCH -e logs/err-pbmc4k
#SBATCH -n 1
#SBATCH --mem 32000

echo 'source("run_pbmc4k.R")' | ${R} --slave --vanilla
EOT

# Submitting the Kolodziejczyk simulations.

sbatch << EOT
#!/bin/bash
#SBATCH -o logs/out-kolod
#SBATCH -e logs/err-kolod
#SBATCH -n 1
#SBATCH --mem 16000

echo 'source("run_kolod.R")' | ${R} --slave --vanilla
EOT
