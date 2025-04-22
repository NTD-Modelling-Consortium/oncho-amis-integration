
#### This file creates the MDA and VC histories for each unique treatment history
#### Relies on outputs from multipletimepoints_preprocess_maps.R

# cpu <- Sys.getenv("MODULE_CPU_TYPE")
# setwd(paste0("../oncho-mtp-",cpu,"/EPIONCHO-IBM"))
# setwd("~/Desktop/Oncho/oncho-mtp-EK/oncho-mtp-skylake/EPIONCHO-IBM/")
cpu <- "skylake"
setwd(paste0("../oncho-mtp-",cpu,"/EPIONCHO-IBM"))
library(dplyr)
library(tidyr)
library(magrittr)
library(AMISforInfectiousDiseases)
amis_params <- default_amis_params()
resamples <- 200

#load data and histories
mda_file_raw <- read.csv('../../Maps/Full_histories_df_popinfo_210324.csv',header=T) 
colnames(mda_file_raw)[1] <- "X"
mda_file_raw$adherence_par[which(mda_file_raw$Cov.in2 == 0.52)] = 0.3 # correct error in histories csv
load("../../Maps/ALL_prevalence_map_multipletimespoints.rds")
load("../../Maps/iu_task_lookup.rds")
iu_task_lookup$TaskID = as.integer(iu_task_lookup$TaskID)

# append histories for new IUs
mda_file_updates = read.csv("../../Maps/Full_histories_df_popinfo_4countries_230524.csv",header=T) %>%
  select(-new_MDA_IUs_2022)
mda_file_updates$adherence_par[which(mda_file_updates$adherence_par %in% c(0.5,0.6))] = 0.3 # use 0.3 for endgame
mda_file = rbind(mda_file_raw %>% filter(!IU_CODE_MAPPING %in% mda_file_updates$IU_CODE_MAPPING),
                 mda_file_updates)

# Subsampling function
sample_init_values <- function(params, weights, seeds, nsamples, replace=F) {
  sampled_idx <- sample.int(
    nrow(params),
    nsamples,
    replace = replace,
    prob = weights
  )
  sampled <- cbind(seeds[sampled_idx], params[sampled_idx,,drop=FALSE])
  colnames(sampled) <- c("seeds", "individual_exposure","bite_rate_per_person_per_year")
  return(sampled)
}

# Save MDA files 
if (!dir.exists("model_output")) {dir.create("model_output")}

# Load table that linked IU to batch id in the fitting process
table_iu_idx <- read.csv("../../outputs/table_iu_idx.csv")

# Define new batches for projections
countries <- sort(unique(table_iu_idx$country))
for (country_id in seq_along(countries)){
  country <- countries[country_id]
  wh <- which(table_iu_idx$country%in%country)
  IUs <- unique(table_iu_idx[wh,"IU_CODE"])  
  IUs <- matrix(IUs,ncol=1)
  iu_file <- file.path("model_output",paste0("IUs_MTP_proj_",country_id,".csv"))
  # write.table(IUs, file=iu_file, row.names=F, col.names = F, quote=F, sep=",")# write input parameter file
}


# Splitting COD
country_id <- which(countries=="COD")
country <- countries[country_id]
wh <- which(table_iu_idx$country%in%country)
IUs <- unique(table_iu_idx[wh,"IU_CODE"])
idxs <- seq(1, 359, by=72)
parts <- vector("list", length(idxs))
for(i in 1:length(idxs)){
  parts[[i]] <- idxs[i]:(idxs[i]+71)
}
parts[[5]] <- parts[[5]][-length(parts[[5]])]
for(part in 1:length(parts)){
  IUs_part <- matrix(IUs[ parts[[part]] ],ncol=1)
  iu_file <- file.path("model_output",paste0("IUs_MTP_proj_",country_id,part,".csv"))
  # write.table(IUs_part, file=iu_file, row.names=F, col.names = F, quote=F, sep=",")# write input parameter file
}


# Splitting ETH
country_id <- which(countries=="ETH")
country <- countries[country_id]
wh <- which(table_iu_idx$country%in%country)
IUs <- unique(table_iu_idx[wh,"IU_CODE"])  
parts <- list(c(1:74), c(75:148), c(149:222))
for(part in 1:length(parts)){
  IUs_part <- matrix(IUs[ parts[[part]] ],ncol=1)
  iu_file <- file.path("model_output",paste0("IUs_MTP_proj_",country_id,part,".csv"))
  # write.table(IUs_part, file=iu_file, row.names=F, col.names = F, quote=F, sep=",")# write input parameter file
}


# Splitting NGA
country_id <- which(countries=="NGA")
country <- countries[country_id]
wh <- which(table_iu_idx$country%in%country)
IUs <- unique(table_iu_idx[wh,"IU_CODE"])
idxs <- seq(1, length(IUs), by=72)
parts <- vector("list", length(idxs))
for(i in 1:length(idxs)){
  parts[[i]] <- idxs[i]:(idxs[i]+71)
  if(i==length(idxs)){
    parts[[i]] <- idxs[i]:(length(IUs))  
  }
}
for(part in 1:length(parts)){
  IUs_part <- matrix(IUs[ parts[[part]] ],ncol=1)
  iu_file <- file.path("model_output",paste0("IUs_MTP_proj_",country_id,part,".csv"))
  write.table(IUs_part, file=iu_file, row.names=F, col.names = F, quote=F, sep=",")# write input parameter file
}

# Splitting 404 batches that were refitted due to an error
batchesIDtoRefit <- read.csv("../../Maps/taskIDs_with_VC.csv")

TaskID <- sort(unique(batchesIDtoRefit$TaskID))
wh <- (iu_task_lookup$TaskID)%in%TaskID
IUs <- sort(unique(iu_task_lookup[wh,"IU_CODE"]))

idxs <- seq(1, length(IUs), by=72)
parts <- vector("list", length(idxs))
for(i in 1:length(idxs)){
  if(i<length(idxs)){
    parts[[i]] <- idxs[i]:(idxs[i]+71)
  }else{
    parts[[i]] <- idxs[i]:length(IUs)
  }
}
for(part in 1:length(parts)){
  IUs_part <- matrix(IUs[ parts[[part]] ],ncol=1)
  iu_file <- file.path("model_output",paste0("IUs_MTP_proj_",90,part,".csv"))
  print(iu_file)
  print(IUs_part)
  write.table(IUs_part, file=iu_file, row.names=F, col.names = F, quote=F, sep=",")
}


# # Loop over all the batches -------------------
# ids <- 1:max(as.integer(iu_task_lookup$TaskID))
# # Loop over specific countries ----------------
# this will run over all IUs in the batches where these countries are found

# wh <- which(table_iu_idx$country%in%c("SSD"))
# # wh <- which(table_iu_idx$country%in%c("CMR","ETH","GHA","SSD"))
# ids <- unique(table_iu_idx[wh,"TaskID"])

if(TRUE){
  
ids <- sort(unique(table_iu_idx$TaskID))
  
# Produce 
cat(paste0("Producing files InputPars_MTP_proj[id] and InputVC_MTP_proj[id] for the fitted batch 'id' \n"))
for (id in ids){
  
  cat(paste0("id = ", id, "; "))
  
  load(paste0("../../outputs/outputs_batch_",id,"/output_",id,".Rdata"))
  
  for (iu in colnames(output$weight_matrix)){
    
    if(amis_params$log==TRUE){weights_iu = exp(as.data.frame(output$weight_matrix)[[iu]])}else{weights_iu=as.data.frame(output$weight_matrix)[[iu]]}
    sampled_params <- sample_init_values(params=cbind(as.data.frame(output$param)[["individual_exposure"]],
                                                      exp(as.data.frame(output$param)[["bite_rate_per_person_per_year"]])),
                                         weights=weights_iu, seeds=output$seeds,
                                         nsamples = resamples,
                                         replace=T)
    sampled_params = cbind(IU_CODE=iu,sampled_params,output$simulated_prevalences[sampled_params[,1],])
    input_file<-file.path("model_output",paste0("InputPars_MTP_proj_",iu,".csv"))
    write.csv(cbind(sampled_params,input_file),file=input_file, row.names = F)# write input parameter file
    
    
    mda_file_iu = mda_file %>%
      filter(IU_CODE_MAPPING == iu ) %>%
      filter(number_rnds>0 & Cov.in2 > 0) %>%
      dplyr::select(Year,number_rnds,Cov.in2,adherence_par) %>%
      mutate(treatment_interval = 1/number_rnds)
    colnames(mda_file_iu) = c("Year","number_rnds","ModelledCoverage","adherence_par","treatment_interval")
    mda_path<-file.path("model_output",paste0("InputMDA_MTP_proj_",iu,".csv"))
    write.csv(mda_file_iu,file=mda_path, row.names = F) # write input MDA file
    
    vc_file_iu = mda_file %>%
      filter(IU_CODE_MAPPING == iu) %>%
      filter(vector_control>0) %>%
      dplyr::select(Year,vector_control)
    vc_all_years = data.frame(Year=1975:2025) %>%
      left_join(vc_file_iu, by=c("Year"))
    vc_all_years$vector_control[is.na(vc_all_years$vector_control)] = 0
    vc_all_years %<>%
      mutate(abr_multiplier = case_when(vector_control==1 ~ 0.2,
                                        vector_control==2 ~ 0,
                                        TRUE ~ 1))
    vc_path<-file.path("model_output",paste0("InputVC_MTP_proj_",iu,".csv"))
    write.csv(vc_all_years,file=vc_path, row.names = F) # write input vector control file
    
  }
}

}
