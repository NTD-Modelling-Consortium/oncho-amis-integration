Rscript "$PATH_TO_MTP_PREPROCESS_PROJ_SCRIPTS/multipletimepoints_preprocess_map_and_histories.R"

# Copying the artifacts to a temporary location makes them available to the host, outside the container
cp -r "$PATH_TO_MAPS" /tmp/Maps
cp -r "$PATH_TO_MTP_PREPROCESS_PROJ_SCRIPTS/model_output" /tmp/model_output