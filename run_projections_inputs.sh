#!/bin/bash

# Parse named arguments
while [ $# -gt 0 ]; do
    case "$1" in
        --slurm-task-id=*)
            SLURM_ARRAY_TASK_ID="${1#*=}"
            ;;
        *)
            echo "Error: Invalid argument $1"
            echo "Usage: $0 --slurm-task-id=<id> "
            exit 1
            ;;
    esac
    shift
done

# Check if slurm-task-id is provided
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
echo "Error: slurm-task-id is required"
echo "Usage: $0 --slurm-task-id=<id>"
        exit 1
fi

export SLURM_ARRAY_TASK_ID
Rscript "$PATH_TO_MTP_PREPROCESS_PROJ_SCRIPTS/multipletimepoints_projections_inputs.R"