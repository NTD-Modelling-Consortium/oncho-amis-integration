Rscript "$PATH_TO_MTP_PREPROCESS_PROJ_SCRIPTS/multipletimepoints_projections_inputs.R"

# Copying the artifacts to a temporary location makes them available to the host, outside the container
cp -r "$PATH_TO_OUTPUTS" /tmp/outputs