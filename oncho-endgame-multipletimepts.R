library(reticulate)
library(truncnorm)
library(dplyr)
library(magrittr)
library(optparse)
library(invgamma)
library(tidyr)
library(AMISforInfectiousDiseases)

kPathToWorkingDir <- Sys.getenv("ONCHO_AMIS_DIR")
kPathToModel <- file.path(kPathToWorkingDir, "model", "EPIONCHO-IBM")
kPathToMaps <- file.path(Sys.getenv("PATH_TO_FITTING_PREP_ARTEFACTS"), "Maps")

# Load python version of oncho model
# cpu <- "skylake"
# use_virtualenv(paste0("./oncho-mtp-",cpu),required=T)
# pathToOnchoPyvenv <- Sys.getenv("ONCHO_PYVENV")
# use_virtualenv(pathToOnchoPyvenv,required=T)
# setwd("/ntdmc/EPIONCHO-IBM")
reticulate::py_config()
# Define python modules
module <- reticulate::import_from_path(module = "epioncho_ibm", path = kPathToModel)
wrapper_fitting <- reticulate::import_from_path(module = "r_wrapper_endgame_fitting_multipletimepts", path = Sys.getenv("PATH_TO_FITTING_SCRIPTS"))
print("Loaded python")

# outputs directory for the batch
kPathToOutputs <- file.path(Sys.getenv("PATH_TO_FITTING_ARTEFACTS"))
if (!dir.exists(kPathToOutputs)) {
  dir.create(kPathToOutputs, recursive = TRUE)
}

# Load data and extract IUs by taskID
load(file.path(kPathToMaps, "ALL_prevalence_map_multipletimespoints.rds"))


option_list <- list(
  make_option(c("-i", "--id"),
    type = "integer",
    help = "Single ID to process."
  ),
  make_option(c("--amis-sigma"),
    type = "numeric",
    default = 0.0025,
    help = "AMIS 'sigma' parameter (default: 0.0025)"
  ),
  make_option(
    c("--amis-max-iters"),
    type = "integer",
    default = 50,
    help = "Maximum number of AMIS iterations (default: 50)"
  ),
  make_option(c("--amis-n-samples"),
    type = "integer",
    default = 500,
    help = "Number of AMIS samples (default: 500)"
  ),
  make_option(c("--amis-target-ess"),
    type = "integer",
    default = 500,
    help = "Target ESS parameter for AMIS (default: 500)"
  ),
  make_option(c("--amis-delete-induced-prior"),
    type = "logical",
    default = FALSE,
    help = paste(
      "Flag controlling whether the induced prior density is to be deleted when updating weights.",
      "Is taken to be FALSE if the flag is not specified."
    )
  )
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)


prevalence_map <- lapply(1:length(map_all_mtp), function(t) {
  map_t <- map_all_mtp[[t]]
  wh <- which(map_t$TaskID == opts$"id")
  map_t <- map_t[wh, ]
  rownames(map_t) <- map_t$IUID
  map_t <- as.matrix(map_t[, (ncol(map_t) - 1000 + 1):ncol(map_t), drop = FALSE])
  return(list(data = map_t))
})

# Priors
k_alpha <- 10
k_beta <- 3
k_min <- 0.01
k_max <- 4
abr_mean <- log(1000)
abr_sd <- log(5000) - log(1000)
abr_min <- log(0.01)
abr_max <- log(200000)

dinvgammatrunc <- function(x, shape, rate, lower, upper, log = F) {
  res <- ifelse(x < lower | x > upper, 0, dinvgamma(x, shape, rate) / diff(pinvgamma(c(lower, upper), shape, rate)))
  if (log == TRUE) {
    res <- log(res)
  }
  return(res)
}
rinvgammatrunc <- function(n, shape, rate, lower, upper) {
  samples <- rep(NA, n)
  for (i in 1:n) {
    samples[i] <- rinvgamma(1, shape, rate)
    while (samples[i] < lower | samples[i] > upper) {
      samples[i] <- rinvgamma(1, shape, rate)
    }
  }
  return(samples)
}

rprior <- function(n) {
  params <- matrix(NA, n, 2)
  colnames(params) <- c("individual_exposure", "bite_rate_per_person_per_year")
  params[, 1] <- rinvgammatrunc(n, shape = k_alpha, rate = k_beta, lower = k_min, upper = k_max)
  params[, 2] <- rtruncnorm(n, a = abr_min, b = abr_max, mean = abr_mean, sd = abr_sd)
  return(params)
}

dprior <- function(x, log = FALSE) {
  if (!is.matrix(x)) {
    if (log) {
      return(sum(
        dinvgammatrunc(x[1], shape = k_alpha, rate = k_beta, lower = k_min, upper = k_max, log = TRUE),
        log(dtruncnorm(x[2], a = abr_min, b = abr_max, mean = abr_mean, sd = abr_max))
      ))
    } else {
      return(prod(
        dinvgammatrunc(x[1], shape = k_alpha, rate = k_beta, lower = k_min, upper = k_max),
        dtruncnorm(x[2], a = abr_min, b = abr_max, mean = abr_mean, sd = abr_max)
      ))
    }
  } else {
    if (log) {
      return(sum(
        dinvgammatrunc(x[, 1], shape = k_alpha, rate = k_beta, lower = k_min, upper = k_max, log = TRUE),
        log(dtruncnorm(x[, 2], a = abr_min, b = abr_max, mean = abr_mean, sd = abr_max))
      ))
    } else {
      return(prod(
        dinvgammatrunc(x[, 1], shape = k_alpha, rate = k_beta, lower = k_min, upper = k_max),
        dtruncnorm(x[, 2], a = abr_min, b = abr_max, mean = abr_mean, sd = abr_max)
      ))
    }
  }
}
prior <- list(rprior = rprior, dprior = dprior)

# Define transmission model
sigma_suffix <- ifelse(opts$"amis-sigma" == 0.0025, "", paste0("_sigma", opts$"amis-sigma"))
trajectories <- c() # save simulated trajectories as code is running
path_to_trajectories_id <- file.path(kPathToOutputs, paste0("trajectories_", opts$"id", sigma_suffix, ".Rdata"))
save(trajectories, file = path_to_trajectories_id)

transmission_model <- function(seeds, parameters, n_tims = length(map_all_mtp)) {
  parameters[, 2] <- exp(parameters[, 2])
  # Wrapper function for python model
  output <- wrapper_fitting$wrapped_parameters(parameters = cbind(seeds, parameters))

  # save trajectories
  load(file = path_to_trajectories_id)
  trajectories <- rbind(trajectories, t(sapply(1:length(output), function(i) output[[i]][[2]])))
  save(trajectories, file = path_to_trajectories_id)

  output_prev <- t(sapply(1:length(output), function(i) output[[i]][[1]]))
  return(output_prev)
}

# Algorithm parameters
amis_params <- default_amis_params()
amis_params$max_iters <- opts$"amis-max-iters" # 50
amis_params$n_samples <- opts$"amis-n-samples" # 500
amis_params$target_ess <- opts$"amis-target-ess" # 500
amis_params$sigma <- opts$"amis-sigma" # 0.0025
amis_params$delete_induced_prior <- opts$"amis-delete-induced-prior" # FALSE

# Run AMIS
set.seed(NULL)
st <- Sys.time()
output <- amis(prevalence_map, transmission_model, prior, amis_params, seed = opts$"id")
en <- Sys.time()
dur_amis <- as.numeric(difftime(en, st, units = "mins"))
print(dur_amis)

## Save AMIS output
save(output, file = file.path(kPathToOutputs, paste0("output_", opts$"id", sigma_suffix, ".Rdata")))

# Output to summary file
ess <- output$ess
n_success <- length(which(ess >= amis_params[["target_ess"]]))
failures <- which(ess < amis_params[["target_ess"]])
n_failure <- length(failures)
if (n_failure > 0) {
  cat(paste(rownames(prevalence_map[[1]][[1]])[failures], opts$"id", ess[failures]),
    file = file.path(kPathToOutputs, paste0("ESS_NOT_REACHED", sigma_suffix, ".txt")),
    sep = "\n", append = TRUE
  )
}
n_sim <- nrow(output$weight_matrix)

summary_file <- file.path(kPathToOutputs, paste0("summary", sigma_suffix, ".csv"))
if (!file.exists(summary_file)) {
  cat("ID,n_failure,n_success,n_sim,min_ess,duration_amis\n", file = summary_file)
}
cat(opts$"id", n_failure, n_success, n_sim, min(ess), dur_amis, "\n",
  sep = ",", file = summary_file, append = TRUE
)
