
#### This file samples 200 parameters from the fitting, generates the histories for each IU up to 2025 and reassigns batches for the projections (IU's are treated individually rather than in batches)
#### Relies on outputs from multipletimepoints_preprocess_maps.R and the fitting

# assume working directory is EPIONCHO-IBM

library(readxl)
library(dplyr)
library(tidyr)
library(magrittr)
library(AMISforInfectiousDiseases)
resamples <- 200

failed_ids=c() #ius that failed AMIS so we have no results
# # Loop over all the batches -------------------
#ids <- sort(unique(table_iu_idx$TaskID))
#just loop over batches that were rerun
ids=setdiff(675:721,failed_ids)

set.seed(1234)

#load data and histories
mda_file = read_excel('../Maps/Full_histories_df_popinfo_ALL_minimal_070425_listlabels.xlsx')
mda_file[mda_file=="NA"] = NA

load("../Maps/ALL_prevalence_map_multipletimespoints.rds")
load("../Maps/iu_task_lookup.rds")

# Save MDA files 
if (!dir.exists("model_output")) {dir.create("model_output")}

# Load table that linked IU and country to batch id in the fitting process
table_iu_idx <- map_all_mtp[[1]][,c("IUID","TaskID")] %>%
  mutate(country = substr(IUID,1,3))
table_iu_idx$TaskID <- as.integer(table_iu_idx$TaskID)
write.csv(table_iu_idx,file=paste0("../outputs/table_iu_idx.csv"),row.names=F)


if(TRUE){
  sampled_params_all = c()
 
  # Produce 
  cat(paste0("Producing files InputPars_MTP_proj[id] and InputVC_MTP_proj[id] for the fitted batch 'id' \n"))
  for (id in ids){
    
    cat(paste0("id = ", id, "; "))
    load(paste0("../outputs/output_",id,".Rdata")) # this is just to get ESS initial run with sigma=0.0025
    
    # IUIDs weren't added to prevalance_map so manually add them here
    wh=which(map_all_mtp[[1]]$TaskID==id)
    iu_list=map_all_mtp[[1]]$IUID[wh]
    # # use this line if rerunning (code updated so colnames of weight_matrix should reflect IUID in future runs)
    # iu_list = colnames(output$weight_matrix)
    
    # find ESS with ESS<200
    ess = output$ess
    iu_names_lt200 = iu_list[ess<200]
    iu_names_ge200 = iu_list[!ess<200]
      
    for (iu_idx in 1:length(iu_list)){
      
      iu = iu_list[iu_idx]
      
      # use sigma=0.025 results if ESS < 200 when sigma=0.0025
      if(iu %in% iu_names_ge200){
        load(paste0("../outputs/output_",id,".Rdata")) # loads amis_output
      } else {
        load(paste0("../outputs/output_",id,"_sigma0.025.Rdata")) # loads amis_output
      }
      
      # sample posterior parameters
      sampled_params <- sample_parameters(output,
                                          n_samples = resamples,
                                          locations=which(iu_list==iu))
      sampled_params = cbind(IUID=iu,sampled_params)
      # colnames for prevalence weird so manually change
      colnames(sampled_params) = c(colnames(sampled_params[1:5]),"prev_t1","prev_t2","prev_t3")
      
      input_file<-file.path("model_output",paste0("InputPars_MTP_proj_",iu,".csv"))
      write.csv(cbind(sampled_params,input_file),file=input_file, row.names = F)# write input parameter file
      #save to all params object
      sampled_params_all = rbind(sampled_params_all,sampled_params)
      
      #save MDA 
      mda_file_iu = mda_file %>%
        filter(IUID == iu ) %>%
        filter(number_rnds>0 & Cov.in2 > 0) %>%
        dplyr::select(Year,number_rnds,Cov.in2,adherence_par) %>%
        mutate(treatment_interval = 1/number_rnds)
      colnames(mda_file_iu) = c("Year","number_rnds","ModelledCoverage","adherence_par","treatment_interval")
      mda_path<-file.path("model_output",paste0("InputMDA_MTP_proj_",iu,".csv"))
      write.csv(mda_file_iu,file=mda_path, row.names = F) # write input MDA file
      
      #save VC
      vc_file_iu = mda_file %>%
        filter(IUID == iu) %>%
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
  save(sampled_params_all,file= paste0("../outputs/InputPars_MTP_allIUs.rds"))
}


# Define new batches for projections
# remove the 2 instances of '%>% filter(TaskID %in% ids)' when doing all batches
num_IUs_per_batch <- 12
num_batches <- ceiling(nrow(table_iu_idx %>% filter(TaskID %in% ids))/num_IUs_per_batch)
for(id in 1:num_batches){
  wh <- 12*(id-1) + 1:12
  IUs <- (table_iu_idx %>% filter(TaskID %in% ids))[wh,"IUID"]
  IUs <- matrix(IUs[which(!is.na(IUs))],ncol=1)
  iu_file <- file.path("model_output",paste0("IUs_MTP_proj_",id,".csv"))
  write.table(IUs, file=iu_file, row.names=F, col.names = F, quote=F, sep=",")# write input parameter file
}

cat(paste0("Produced samples for all IUs \n"))
