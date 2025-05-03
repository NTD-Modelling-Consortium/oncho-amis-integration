#!/bin/bash

#SBATCH --output outputs/log/mtp-oncho.out-%A_%a
#SBATCH --array 1
#SBATCH --nodes=6
#SBATCH --cpus-per-task 12
#SBATCH --time=30:00:00

# add required batch numbers to SBATCH --array line

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
unset RETICULATE_PYTHON
cd EPIONCHO-IBM
Rscript ../oncho-endgame-multipletimepts-sigma0.025.R ${SLURM_ARRAY_TASK_ID} 

