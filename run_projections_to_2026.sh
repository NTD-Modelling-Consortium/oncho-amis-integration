#!/bin/bash

export SLURM_ARRAY_TASK_ID=$1

# Activate the oncho model virtual environment if not already done
# Check if the virtual environment is already activated
if [ -z "$VIRTUAL_ENV" ]; then
    # If not, activate it
    cd ./model/EPIONCHO-IBM
    source $(poetry env info --path)/bin/activate
else
    echo "Virtual environment is already activated."
fi

# Get back to the oncho-amis-integration directory
# Check if current directory is not ONCHO_AMIS_DIR
if [ "$(basename "$PWD")" != "oncho-amis-integration" ]; then
    # If not, change to the ONCHO_AMIS_DIR
    cd $ONCHO_AMIS_DIR
else
    echo "Already in the oncho-amis-integration directory."
fi

python wrapper-simulations-multipletimepoints.py