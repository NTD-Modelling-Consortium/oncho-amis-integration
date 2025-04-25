
#### This file creates the MDA and VC histories for each unique treatment history
#### Relies on outputs from multipletimepoints_preprocess_maps.R

# cpu <- Sys.getenv("MODULE_CPU_TYPE")
# setwd(paste0("../oncho-mtp-",cpu,"/EPIONCHO-IBM"))
# setwd("~/Desktop/Oncho/oncho-mtp-EK/oncho-mtp-skylake/EPIONCHO-IBM/")
# cpu <- "skylake"
# setwd(paste0("../oncho-mtp-",cpu,"/EPIONCHO-IBM"))
# setwd("/ntdmc/EPIONCHO-IBM")

library(dplyr)
library(tidyr)
library(magrittr)

# pathToMaps <- "../oncho-amis-integration/Maps"
pathToMaps <- "Maps"

#load data and histories
mda_file = read.csv(paste0(pathToMaps, '/Full_histories_df_popinfo_210324.csv'),header=T)
colnames(mda_file)[1] = "X"
mda_file$adherence_par[which(mda_file$Cov.in2 == 0.52)] = 0.3 # correct error in histories csv
# append histories for new IUs
newIUs = c( "ETH0194318576", "ETH0194318658", "ETH0194318537","ETH0194318591", "ETH0194318592", "ETH0194318651", "ETH0194818833","ETH0194818871", "ETH0194818905", "ETH0195019197","ETH0195019253", "ETH0195019566", "ETH0195019567")
mda_file_newIUs = read.csv(paste0(pathToMaps, "/Full_histories_df_popinfo_4countries_230524.csv"),header=T) %>%
  filter(IU_CODE_MAPPING %in% newIUs) %>%
  select(-new_MDA_IUs_2022)
mda_file_newIUs$adherence_par[which(mda_file_newIUs$adherence_par %in% c(0.5,0.6))] = 0.3 # use 0.3 for endgame
mda_file = rbind(mda_file %>% filter(!IU_CODE_MAPPING %in% newIUs),
                 mda_file_newIUs)

load(paste0(pathToMaps, "/ALL_prevalence_map_multipletimespoints.rds"))
load(paste0(pathToMaps, "/iu_task_lookup.rds"))

pathToModelOutput <- "mtp-preprocess_projections/model_output"

# Save MDA files 
if (!dir.exists(pathToModelOutput)) {dir.create(pathToModelOutput, recursive = T)}
for (id in 1:max(as.integer(iu_task_lookup$TaskID))){
  
  iu_file<-file.path(pathToModelOutput,paste0("IUs_MTP_",id,".csv"))
  map_all_mtp_id = map_all_mtp[[1]] %>%
    filter(TaskID == id)
  ius = matrix(map_all_mtp_id$IU_CODE,ncol=1)
  write.table(ius,file=iu_file, row.names=F, col.names = F, quote=F, sep=",")# write input parameter file
  
  for (iu in ius[1]) { # histories are all the same within a batch so only need to generate once
    mda_file_iu = mda_file %>%
      filter(IU_CODE_MAPPING == iu & Year <= 2018) %>%
      filter(number_rnds>0 & Cov.in2 > 0) %>%
      dplyr::select(Year,number_rnds,Cov.in2,adherence_par) %>%
      mutate(treatment_interval = 1/number_rnds)
    colnames(mda_file_iu) = c("Year","number_rnds","ModelledCoverage","adherence_par","treatment_interval")
    mda_path<-file.path(pathToModelOutput,paste0("InputMDA_MTP_",id,".csv"))
    write.csv(mda_file_iu,file=mda_path, row.names = F) # write input MDA file
    
    vc_file_iu = mda_file %>%
      filter(IU_CODE_MAPPING == iu & Year <= 2018) %>%
      filter(vector_control>0) %>%
      dplyr::select(Year,vector_control)
    vc_all_years = data.frame(Year=1975:2025) %>%
      left_join(vc_file_iu, by=c("Year"))
    vc_all_years$vector_control[is.na(vc_all_years$vector_control)] = 0
    vc_all_years %<>%
      mutate(abr_multiplier = case_when(vector_control==1 ~ 0.2,
                                        vector_control==2 ~ 0,
                                        TRUE ~ 1))
    vc_path<-file.path(pathToModelOutput,paste0("InputVC_MTP_",id,".csv"))
    write.csv(vc_all_years,file=vc_path, row.names = F) # write input vector control file
  }
}



