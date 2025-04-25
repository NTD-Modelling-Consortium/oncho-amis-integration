#!/bin/bash

export SLURM_ARRAY_TASK_ID=$1
export ONCHO_PYVENV=/root/.cache/pypoetry/virtualenvs/epioncho-ibm-U5jSXRn5-py3.10

source $ONCHO_PYVENV/bin/activate
export RETICULATE_PYTHON=$(which python)

Rscript oncho-endgame-multipletimepts.R ${SLURM_ARRAY_TASK_ID}