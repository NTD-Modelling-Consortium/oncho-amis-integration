
#### This file assigns batch IDs by finding IUs with the same treatment histories and generates the map with multiple time points

# cpu <- Sys.getenv("MODULE_CPU_TYPE")
# setwd(paste0("../oncho-mtp-",cpu,"/EPIONCHO-IBM"))
# setwd("~/Desktop/Oncho/oncho-mtp-EK/oncho-mtp-skylake/EPIONCHO-IBM/")
cpu <- "skylake"
setwd(paste0("../oncho-mtp-",cpu,"/EPIONCHO-IBM"))
library(dplyr)
library(tidyr)
library(magrittr)

#load data and histories
mda_file = read.csv('../../Maps/Full_histories_df_popinfo_210324.csv',header=T) 
mda_file$adherence_par[which(mda_file$Cov.in2 == 0.52)] = 0.3 # correct error in histories csv


load('../../Maps/ALL_prevalence_map.rds')
wh<-which(map_all$IU_CODE %in% mda_file$IU_CODE)
map_all<-map_all[wh,]

prevalence_map = as.matrix(map_all[,(ncol(map_all)-1000+1):ncol(map_all),drop=FALSE])
rownames(prevalence_map)<-map_all$IU_CODE

#find unique mda histories
mda_history = mda_file %>%
  filter(Year <= 2018) %>%
  filter(IU_CODE_MAPPING %in% rownames(prevalence_map)) %>%
  select(IU_CODE_MAPPING,Year,number_rnds) %>%
  pivot_wider(names_from="Year",values_from="number_rnds")

unique_histories = mda_history %>% 
  group_by(across(c(-IU_CODE_MAPPING))) %>% 
  summarise(ius= list(IU_CODE_MAPPING)) %>%
  ungroup() %>%
  mutate(index = row_number())

histories_indexes = lapply(1:nrow(unique_histories), function(i) data.frame(unique_histories[["ius"]][[i]],unique_histories[["index"]][i]))
histories_indexes = do.call(rbind,histories_indexes)
colnames(histories_indexes) = c("IU","histories_index")

# find unique coverage histories
coverage_history = mda_file %>%
  filter(Year <= 2018) %>%
  filter(IU_CODE_MAPPING %in% rownames(prevalence_map)) %>%
  mutate(Cov.in2 = ifelse(number_rnds>0,Cov.in2,NA)) %>% 
  select(IU_CODE_MAPPING,Year,Cov.in2) %>%
  pivot_wider(names_from="Year",values_from="Cov.in2")

unique_coverage = coverage_history %>% 
  group_by(across(c(-IU_CODE_MAPPING))) %>% 
  summarise(ius= list(IU_CODE_MAPPING)) %>%
  ungroup() %>%
  mutate(index = row_number())

coverage_indexes = lapply(1:nrow(unique_coverage), function(i) data.frame(unique_coverage[["ius"]][[i]],unique_coverage[["index"]][i]))
coverage_indexes = do.call(rbind,coverage_indexes)
colnames(coverage_indexes) = c("IU","coverage_index")

# find unique adherence histories
adherence_history = mda_file %>%
  filter(Year <= 2018) %>%
  filter(IU_CODE_MAPPING %in% rownames(prevalence_map)) %>%
  select(IU_CODE_MAPPING,Year,adherence_par) %>%
  pivot_wider(names_from="Year",values_from="adherence_par")

unique_adherence = adherence_history %>% 
  group_by(across(c(-IU_CODE_MAPPING))) %>% 
  summarise(ius= list(IU_CODE_MAPPING)) %>%
  ungroup() %>%
  mutate(index = row_number())

adherence_indexes = lapply(1:nrow(unique_adherence), function(i) data.frame(unique_adherence[["ius"]][[i]],unique_adherence[["index"]][i]))
adherence_indexes = do.call(rbind,adherence_indexes)
colnames(adherence_indexes) = c("IU","adherence_index")

# find unique vector control histories
vc_history = mda_file %>%
  filter(Year <= 2018) %>%
  filter(IU_CODE_MAPPING %in% rownames(prevalence_map)) %>%
  select(IU_CODE_MAPPING,Year,vector_control) %>%
  pivot_wider(names_from="Year",values_from="vector_control") 

unique_vc = vc_history %>% 
  group_by(across(c(-IU_CODE_MAPPING))) %>% 
  summarise(ius= list(IU_CODE_MAPPING)) %>%
  ungroup() %>%
  mutate(index = row_number())

vc_indexes = lapply(1:nrow(unique_vc), function(i) data.frame(unique_vc[["ius"]][[i]],unique_vc[["index"]][i]))
vc_indexes = do.call(rbind,vc_indexes)
colnames(vc_indexes) = c("IU","vc_index")

# join all and find unique combinations
all_treatments = data.frame(IU=rownames(prevalence_map)) %>%
  left_join(histories_indexes) %>%
  left_join(coverage_indexes) %>%
  left_join(adherence_indexes) %>%
  left_join(vc_indexes)

unique_all = all_treatments %>%
  group_by(across(c(-IU))) %>% 
  summarise(ius= list(IU)) %>%
  ungroup() %>%
  mutate(index = row_number())

iu_task_lookup = lapply(1:nrow(unique_all), function(i) data.frame(unique_all[["ius"]][[i]],unique_all[["index"]][i]))
iu_task_lookup = do.call(rbind,iu_task_lookup)
colnames(iu_task_lookup) = c("IU_CODE","TaskID")
iu_task_lookup=iu_task_lookup[which(!iu_task_lookup$IU_CODE %in% c("ETH0194318576", "ETH0194318658")),]
# add new IUs
newIUs = c( "ETH0194318576", "ETH0194318658", "ETH0194318537","ETH0194318591", "ETH0194318592", "ETH0194318651", "ETH0194818833","ETH0194818871", "ETH0194818905", "ETH0195019197","ETH0195019253", "ETH0195019566", "ETH0195019567")
iu_task_lookup = rbind(iu_task_lookup,cbind(IU_CODE=newIUs,TaskID=(max(as.integer(iu_task_lookup$TaskID))+1)))
save(iu_task_lookup, file="../../Maps/iu_task_lookup.rds")

# reload map to add in new IUs
load('../../Maps/ALL_prevalence_map.rds')
wh<-which(map_all$IU_CODE %in% iu_task_lookup$IU_CODE)
map_all<-map_all[wh,]

map_all_mtp = iu_task_lookup %>% 
  left_join(map_all %>% select(-TaskID))
map_all_mtp = list(map_all_mtp)

# add 2000 map
map_2000 = read.csv('../../Maps/maps_joint/samples_2000.csv',header=T) 
map_2000 = map_2000[,c(2,6:1005)]
map_all_mtp[[2]] = map_all_mtp[[1]] %>%
  select(IU_CODE,TaskID) %>%
  left_join(map_2000,by=c("IU_CODE"="IU_code"))

# add 2018 map
map_2018 = read.csv('../../Maps/maps_joint/samples_2018.csv',header=T) 
map_2018 = map_2018[,c(2,6:1005)]
map_all_mtp[[3]] = map_all_mtp[[1]] %>%
  select(IU_CODE,TaskID) %>%
  left_join(map_2018,by=c("IU_CODE"="IU_code"))

save(map_all_mtp, file="../../Maps/ALL_prevalence_map_multipletimespoints.rds")

## CAREFUL with lines below. maps_100pc_rankcorr are being used

# # make map with all years
# map_mtp_all_years = list(map_all_mtp[[1]])
# start_year = 2000
# for (map_t_ind in 2:20){
#   year = start_year + (map_t_ind-2)
#   print(year)
#   map = read.csv(paste0("../maps_100pc_rankcorr/samples_ordered_",year,".csv"),header=T) 
#   
#   map = map[,c(2,6:1005)]
#   map_mtp_all_years[[map_t_ind]] = map_mtp_all_years[[1]] %>%
#     select(IU_CODE,TaskID) %>%
#     left_join(map,by=c("IU_CODE"="IU_code"))
# }
# save(map_all_mtp, file="../../Maps/ALL_prevalence_map_multipletimespoints_allyears.rds")


