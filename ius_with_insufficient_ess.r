# assuming working directory is post_AMIS_analysis

library("AMISforInfectiousDiseases")
library("dplyr")
library("readxl")
library("tidyr")

kPathToMaps <- Sys.getenv("PATH_TO_FITTING_PREP_ARTEFACTS")
kPathToOutputs <- Sys.getenv("PATH_TO_FITTING_ARTEFACTS")
kPathToModelOutput <- file.path(
  Sys.getenv("PATH_TO_PROJECTIONS_PREP_ARTEFACTS"),
  "model_output"
)

load(file.path(kPathToMaps, "iu_task_lookup.rds"))
num_batches <- max(iu_task_lookup$TaskID)

option_list <- list(
  make_option(c("-f", "--failed-ids"),
    type = "character",
    default = "",
    help = paste(
      "Comma-separated list of failed IDs to skip.",
      "Needs to be filled in once we know which batch-ids failed."
    )
  ),
  make_option(c("--amis-sigma"),
    type = "numeric",
    default = 0.0025,
    help = "AMIS 'sigma' parameter (default: 0.0025)"
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

cat(paste0("disease: oncho \n"))
cat(paste0("num_batches: ", num_batches, " \n"))

ids_sample_pars <- setdiff(1:num_batches, failed_ids)

# sample parameters and save draws
insufficient_ess <- c()

for (id in ids_sample_pars) {
  #### Load AMIS output
  output_file <- file.path(kPathToModelOutput, paste0("output_", id, "_sigma", opts$"amis-sigma", ".Rdata"))
  if (!file.exists(output_file)) {
    load(output_file)

    ess <- output$ess
    iu_names <- rownames(output$prevalence_map[[1]]$data)
    iu_names_lt200 <- iu_names[ess < 200]
    iu_names_ge200 <- iu_names[!ess < 200]

    insufficient_ess <- c(insufficient_ess, iu_names_lt200)

    if (id %% 100 == 0) {
      cat(paste0("id=", id, "; "))
    }
  }
}

insufficient_ess_mat <- data.frame(IU = insufficient_ess)
write.csv(insufficient_ess_mat, file = file.path(kPathToModelOutput, "IUsWithInsufficientESS.csv"), row.names = F)

print("Produced IUsWithInsufficientESS.csv")
