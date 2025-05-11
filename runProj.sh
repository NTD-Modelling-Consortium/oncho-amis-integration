#!/bin/bash

#SBATCH --output outputs/log/mtp-proj_until_2025.out-%A_%a
#SBATCH --array=1-707
#SBATCH --time=24:00:00

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
cd EPIONCHO-IBM
source .venv/bin/activate
poetry install
poetry run python wrappersimulationsmultipletimepoints.py ${SLURM_ARRAY_TASK_ID}

