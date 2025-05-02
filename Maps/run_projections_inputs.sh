#!/bin/bash

#SBATCH --output log/mtp-oncho-projections-prep.out
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00

#####################
# Your task is here #
#####################
# Add extra commands here to load a recent version of R
module purge
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.2.1

stdbuf -i0 -o0 -e0 command
####################
# End of your task #
####################

# Now that you have loaded R above, we can run our R script
cd ../EPIONCHO-IBM
Rscript ../Maps/multipletimepoints_projections_inputs.R

