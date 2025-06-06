#!/bin/bash

# Activate the oncho model virtual environment if not already done
# Check if the virtual environment is already activated
if [ -z "$VIRTUAL_ENV" ]; then
    # If not, activate it
    cd $ONCHO_MODEL_DIR
    source $(poetry env info --path)/bin/activate
else
    echo "Virtual environment is already activated."
fi

Rscript $PATH_TO_FITTING_SCRIPTS/oncho-endgame-multipletimepts.R $@
