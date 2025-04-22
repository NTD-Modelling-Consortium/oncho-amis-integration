id<-as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(id)
# id = 25

#args <- commandArgs(trailingOnly = TRUE)
#options(echo = TRUE)
#print(args)
#id <- as.numeric(args[1])
#print(id)

Sys.setenv(RETICULATE_PYTHON = "/gpfs/home/maths/stsnmt/NTDs/oncho-mtp-EK/oncho-mtp-skylake/bin/python")
library(reticulate)
library(truncnorm)
library(dplyr)
library(magrittr)
library(invgamma)
library(tidyr)
library(AMISforInfectiousDiseases, lib.loc = '/gpfs/home/maths/stsnmt/R/x86_64-pc-linux-gnu-library/4.2')

# Load python version of oncho model
cpu <- "skylake"
use_virtualenv(paste0("./oncho-mtp-",cpu),required=T)
setwd(paste0("./oncho-mtp-",cpu,"/EPIONCHO-IBM"))
reticulate::py_config()

# Define python modules
module = reticulate::import_from_path('epioncho_ibm','.')
wrapper_fitting = reticulate::import('r_wrapper_endgame_fitting_multipletimepts')
print('Loaded python')

# Load data and extract IUs by taskID
# mda_file =  read.csv('Full_histories_df_popinfo_210324.csv',header=T)
load('../../Maps/ALL_prevalence_map_multipletimespoints.rds')

prevalence_map = lapply(1:length(map_all_mtp), function(t) {
  map_t = map_all_mtp[[t]]
  wh<-which(map_t$TaskID==id)
  map_t<-map_t[wh,]
  rownames(map_t)<-map_t$IU_CODE
  map_t = as.matrix(map_t[,(ncol(map_t)-1000+1):ncol(map_t),drop=FALSE])
  return(list(data=map_t))
})

# Priors
k_alpha = 10
k_beta = 3
k_min = 0.01
k_max = 4
abr_mean = log(1000)
abr_sd = log(5000)-log(1000)
abr_min = log(0.01)
abr_max = log(200000)

dinvgammatrunc <- function(x, shape, rate, lower, upper, log=F) {
  res=ifelse(x < lower | x > upper, 0, dinvgamma(x, shape, rate) / diff(pinvgamma(c(lower,upper), shape, rate)))
  if (log==TRUE){
    res=log(res)
  }
  return(res)
}
rinvgammatrunc <- function(n,shape, rate, lower, upper){
  samples=rep(NA,n)
  for (i in 1:n){
    samples[i]=rinvgamma(1, shape, rate)
    while (samples[i] < lower | samples[i] > upper) {samples[i]=rinvgamma(1, shape, rate)}
  }
  return(samples)
}

rprior <- function(n) {
  params<-matrix(NA,n,2)
  colnames(params)<-c("individual_exposure","bite_rate_per_person_per_year")
  params[,1]<-rinvgammatrunc(n, shape=k_alpha, rate=k_beta, lower=k_min, upper=k_max)
  params[,2]<-rtruncnorm(n,a=abr_min,b=abr_max,mean=abr_mean,sd=abr_sd)
  return(params)
}

dprior <- function(x,log=FALSE) {
  if (!is.matrix(x)){
    if (log) {
      return(sum(dinvgammatrunc(x[1],shape=k_alpha, rate=k_beta, lower=k_min, upper=k_max,log=TRUE),
                 log(dtruncnorm(x[2],a=abr_min,b=abr_max,mean=abr_mean,sd=abr_max))))
    } else {
      return(prod(dinvgammatrunc(x[1],shape=k_alpha, rate=k_beta, lower=k_min, upper=k_max),
                  dtruncnorm(x[2],a=abr_min,b=abr_max,mean=abr_mean,sd=abr_max)))
    }

  } else {
    if (log) {
      return(sum(dinvgammatrunc(x[,1],shape=k_alpha, rate=k_beta, lower=k_min, upper=k_max,log=TRUE),
                 log(dtruncnorm(x[,2],a=abr_min,b=abr_max,mean=abr_mean,sd=abr_max))))
    } else {
      return(prod(dinvgammatrunc(x[,1],shape=k_alpha, rate=k_beta, lower=k_min, upper=k_max),
                  dtruncnorm(x[,2],a=abr_min,b=abr_max,mean=abr_mean,sd=abr_max)))
    }
  }
}
prior<-list(rprior=rprior,dprior=dprior)


# outputs directory for the batch
outputs_path_id <- paste0("../../outputs/outputs_batch_",id)
if (!dir.exists(outputs_path_id)) {dir.create(outputs_path_id)}

# Define transmission model
trajectories = c() # save simulated trajectories as code is running
save(trajectories,file=paste0(outputs_path_id,"/trajectories_",id,".Rdata"))

transmission_model=function(seeds,parameters,n_tims=length(map_all_mtp)) {
  parameters[,2] = exp(parameters[,2])
  # Wrapper function for python model
  #output = as.matrix(wrapper_fitting$wrapped_parameters(parameters=cbind(seeds,parameters)))
  output = wrapper_fitting$wrapped_parameters(parameters=cbind(seeds,parameters))
  
  # save trajectories
  load(paste0(outputs_path_id,"/trajectories_",id,".Rdata"))
  trajectories =  rbind(trajectories,t(sapply(1:length(output), function(i) output[[i]][[2]])))
  save(trajectories, file=paste0(outputs_path_id,"/trajectories_",id,".Rdata"))
  
  output_prev = t(sapply(1:length(output), function(i) output[[i]][[1]]))
  return(output_prev)
}


# Algorithm parameters
amis_params<-default_amis_params()
amis_params$max_iters <- 50    #
amis_params$n_samples <- 500   #
amis_params$target_ess <- 500  #
amis_params$sigma <- 0.0025

# Run AMIS
set.seed(NULL)
st = Sys.time()
output = amis(prevalence_map, transmission_model, prior, amis_params, seed=id)
en<-Sys.time()
dur_amis<-as.numeric(difftime(en,st,units="mins"))
print(dur_amis)

## Save AMIS output
save(output,file=paste0(outputs_path_id,"/output_",id,".Rdata"))

# Output to summary file
ess <- output$ess
n_success <- length(which(ess>=amis_params[["target_ess"]]))
failures <- which(ess<amis_params[["target_ess"]])
n_failure <- length(failures)
if (n_failure>0) {cat(paste(rownames(prevalence_map[[1]][[1]])[failures],id,ess[failures]),
                      file = "ESS_NOT_REACHED.txt",sep = "\n", append = TRUE)}
n_sim = nrow(output$weight_matrix)

if (!file.exists("../../outputs/summary.csv")) {cat("ID,n_failure,n_success,n_sim,min_ess,duration_amis\n",file="../../outputs/summary.csv")}
cat(id,n_failure,n_success,n_sim,min(ess),dur_amis,"\n",sep=",",file="../../outputs/summary.csv",append=TRUE)

