# assuming working directory is post_AMIS_analysis

library("AMISforInfectiousDiseases")
library("dplyr")
library("readxl")
library("tidyr")


load("Maps/iu_task_lookup.rds")
num_batches <- max(iu_task_lookup$TaskID)
failed_ids = c()

cat(paste0("disease: oncho \n"))
cat(paste0("num_batches: ",num_batches, " \n"))

ids_sample_pars = setdiff(1:num_batches,failed_ids)

# sample parameters and save draws
insufficient_ess = c()

for(id in ids_sample_pars){

  #### Load AMIS output
  if(file.exists(paste0("outputs/output_",id,"_sigma0.025.Rdata"))){
  
  load(paste0("outputs/output_",id,"_sigma0.025.Rdata"))
	  
  ess = output$ess
  iu_names <- rownames(output$prevalence_map[[1]]$data)
  iu_names_lt200 = iu_names[ess<200]
  iu_names_ge200 = iu_names[!ess<200]

  insufficient_ess = c(insufficient_ess,iu_names_lt200)

  if(id%%100==0){cat(paste0("id=",id, "; "))}
  }
}

insufficient_ess_mat = data.frame(IU=insufficient_ess)
write.csv(insufficient_ess_mat, file=paste0("EPIONCHO-IBM/model_output/IUsWithInsufficientESS_oncho.csv"),row.names=F)

print("Produced IUsWithInsufficientESS.csv")


