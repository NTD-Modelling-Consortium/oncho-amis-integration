#!/bin/bash

export SLURM_ARRAY_TASK_ID=$1

# Activate the oncho model virtual environment
cd ./model/EPIONCHO-IBM
source $(poetry env info --path)/bin/activate

# Get back to the oncho-amis-integration directory
cd $ONCHO_AMIS_DIR

Rscript oncho-endgame-multipletimepts.R ${SLURM_ARRAY_TASK_ID}