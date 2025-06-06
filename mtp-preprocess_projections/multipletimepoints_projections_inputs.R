#### This file samples 200 parameters from the fitting, generates the histories for each IU up to 2025 and reassigns batches for the projections (IU's are treated individually rather than in batches)
#### Relies on outputs from multipletimepoints_preprocess_maps.R and the fitting

library(dplyr)
library(tidyr)
library(magrittr)
library(optparse)
library(AMISforInfectiousDiseases)

set.seed(1234)

kPathToFittingPrep <- Sys.getenv("PATH_TO_FITTING_PREP")
kPathToFittingPrepInputs <- Sys.getenv("PATH_TO_FITTING_PREP_INPUTS")
kPathToFittingPrepArtefacts <- Sys.getenv("PATH_TO_FITTING_PREP_ARTEFACTS")
kPathToFittingArtefacts <- Sys.getenv("PATH_TO_FITTING_ARTEFACTS")
kPathToOutputs <- Sys.getenv("PATH_TO_PROJECTIONS_PREP_ARTEFACTS")
kPathToModelOutput <- file.path(kPathToOutputs, "model_output")

if (!dir.exists(kPathToModelOutput)) {
  dir.create(kPathToModelOutput, recursive = TRUE)
}

option_list <- list(
  make_option(c("-i", "--id"),
    type = "integer",
    help = "Single ID to process. If not provided, all IDs will be processed."
  ),
  make_option(c("-f", "--failed-ids"),
    type = "character",
    default = "",
    help = paste(
      "Comma-separated list of failed IDs to skip.",
      "Needs to be filled in once we know which batch-ids failed.",
      "Will be ignored when a single ID is specified with the -i/--id option."
    )
  ),
  make_option(c("--amis-n-samples"),
    type = "integer",
    default = 200,
    help = "Number of AMIS samples (default: 200)"
  ),
  make_option(c("--amis-sigma"),
    type = "numeric",
    default = 0.0025,
    help = "AMIS 'sigma' parameter (default: 0.0025)"
  ),
  make_option(c("--ess-threshold"),
    type = "numeric",
    default = 200,
    help = "ESS threshold parameter (default: 200)",
  )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

# ius that failed AMIS so we have no results
failed_ids <- if (!is.null(opts$"failed-ids")) {
  as.numeric(strsplit(opts$"failed-ids", ",")[[1]])
} else {
  c() # empty vector if no failed IDs provided
}

# load data and histories
mda_file <- readRDS(file.path(kPathToFittingPrepInputs, "Full_histories_df_popinfo_ALL_minimal_070425_listlabels.rds"))
mda_file[mda_file == "NA"] <- NA

load(file.path(kPathToFittingPrepArtefacts, "Maps", "ALL_prevalence_map_multipletimespoints.rds"))
load(file.path(kPathToFittingPrepArtefacts, "Maps", "iu_task_lookup.rds"))

# Save MDA files

# Load table that linked IU to batch id in the fitting process
table_iu_idx <- map_all_mtp[[1]][, c("IUID", "TaskID")] %>%
  mutate(country = substr(IUID, 1, 3))
table_iu_idx$TaskID <- as.integer(table_iu_idx$TaskID)
write.csv(table_iu_idx, file = file.path(kPathToOutputs, "table_iu_idx.csv"), row.names = F)


process_batch <- function(id) {
  # Load the output for the given ID
  # this is just to get ESS initial run with sigma=0.0025
  output_file <- file.path(kPathToFittingArtefacts, paste0("output_", id, ".Rdata"))
  load(output_file) # loads amis_output

  # colnames of weight_matrix should reflect IUID
  iu_list <- colnames(output$weight_matrix)

  # find ESS with ESS<200
  ess <- output$ess
  iu_names_lt_ess_threshold <- iu_list[ess < opts$"ess-threshold"]
  iu_names_ge_ess_threshold <- iu_list[!ess < opts$"ess-threshold"]

  sampled_params_id <- c()
  for (iu_idx in 1:length(iu_list)) {
    iu <- iu_list[iu_idx]

    # use sigma=0.025 results if ESS<ess_threshold when sigma=0.0025
    if (iu %in% iu_names_lt_ess_threshold) {
      output_file <- file.path(
        kPathToFittingArtefacts,
        paste0("output_", id, "_sigma", opts$"amis-sigma", ".Rdata")
      )
    } else {
      output_file <- file.path(
        kPathToFittingArtefacts,
        paste0("output_", id, ".Rdata")
      )
    }

    load(output_file) # loads amis_output

    # sample posterior parameters
    sampled_params <- sample_parameters(output,
      n_samples = opts$"amis-n-samples",
      locations = which(iu_list == iu)
    )
    sampled_params <- cbind(IUID = iu, sampled_params)
    # colnames for prevalence so manually change
    colnames(sampled_params) <- c(colnames(sampled_params[1:5]), "prev_t1", "prev_t2", "prev_t3")

    cat(sprintf("Producing files InputPars_MTP_proj[%d] and InputVC_MTP_proj[%d] for the fitted batch '%d' \n", id, id, id))
    input_file <- file.path(kPathToModelOutput, paste0("InputPars_MTP_proj_", iu, ".csv"))
    write.csv(cbind(sampled_params, input_file), file = input_file, row.names = F) # write input parameter file
    # save to all params object
    sampled_params_id <- rbind(sampled_params_id, sampled_params)

    # save MDA
    mda_file_iu <- mda_file %>%
      filter(IUID == iu) %>%
      filter(number_rnds > 0 & Cov.in2 > 0) %>%
      dplyr::select(Year, number_rnds, Cov.in2, adherence_par) %>%
      mutate(treatment_interval = 1 / number_rnds)
    colnames(mda_file_iu) <- c("Year", "number_rnds", "ModelledCoverage", "adherence_par", "treatment_interval")
    mda_path <- file.path(kPathToModelOutput, paste0("InputMDA_MTP_proj_", iu, ".csv"))
    write.csv(mda_file_iu, file = mda_path, row.names = F) # write input MDA file

    # save VC
    vc_file_iu <- mda_file %>%
      filter(IUID == iu) %>%
      filter(vector_control > 0) %>%
      dplyr::select(Year, vector_control)
    vc_all_years <- data.frame(Year = 1975:2025) %>%
      left_join(vc_file_iu, by = c("Year"))
    vc_all_years$vector_control[is.na(vc_all_years$vector_control)] <- 0
    vc_all_years %<>%
      mutate(abr_multiplier = case_when(
        vector_control == 1 ~ 0.2,
        vector_control == 2 ~ 0,
        TRUE ~ 1
      ))
    vc_path <- file.path(kPathToModelOutput, paste0("InputVC_MTP_proj_", iu, ".csv"))
    write.csv(vc_all_years, file = vc_path, row.names = F) # write input vector control file
  }

  sampled_params_id
}

process_all_batches <- function(id_vec) {
  # # Loop over all the batches -------------------
  sampled_params_all <- c()

  # Produce
  for (id in id_vec) {
    sampled_params <- process_batch(id)
    sampled_params_all <- rbind(sampled_params_all, sampled_params)
  }
  save(sampled_params_all, file = file.path(kPathToOutputs, "InputPars_MTP_allIUS.rds"))

  wh <- lapply(id_vec, function(id) which(table_iu_idx$TaskID == id))
  IUs <- lapply(wh, function(x) table_iu_idx[x, "IUID"])

  for (i in seq_along(IUs)) {
    iu_file <- file.path(kPathToModelOutput, paste0("IUs_MTP_proj_", id_vec[i], ".csv"))
    cat(sprintf("Writing %s\n", iu_file))
    write.table(IUs[[i]], file = iu_file, row.names = F, col.names = F, quote = F, sep = ",") # write input parameter file
  }

  cat(paste0("Produced samples for all IUs \n"))
}


if (!is.null(opts$"id")) {
  # If a single ID is specified, use it
  ids <- c(as.numeric(opts$"id"))
} else {
  # Otherwise, use all IDs except the failed ones
  ids <- setdiff(sort(unique(table_iu_idx$TaskID)), failed_ids)
}

cat(sprintf("Processing IDs: %s\n", paste(ids, collapse = ", ")))
process_all_batches(ids)
