#!/bin/bash

#SBATCH --output outputs/log/re-endgame-mtp-proj_until_2025.out-%A_%a
#SBATCH --array=1001-1178
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=48:00:00

####################
# Your task is here #
#####################
# Add extra commands here to load a recent version of python
module purge
module load GCCcore/12.2.0 Python/3.10.8 

stdbuf -i0 -o0 -e0 command
####################
# End of your task #
####################
# Now that you have loaded python  above, we can run our script
cd oncho-mtp-skylake/EPIONCHO-IBM
poetry install
poetry run python wrappersimulationsmultipletimepoints.py ${SLURM_ARRAY_TASK_ID}

