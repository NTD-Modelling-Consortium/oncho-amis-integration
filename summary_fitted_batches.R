

# setwd("~/Desktop/Oncho/oncho-mtp-EK/outputs")
setwd(paste0("./outputs"))

countries_batches <- read.csv("countries_batches.csv")

all_batches <- sort(unique(countries_batches$batch))

num_batches <- length(all_batches)
summary_fitted <- NULL
summary_not_fitted <- NULL
for(id in all_batches){
  
  outputs_batch_id <- paste0("outputs_batch_",id,"/output_",id,".Rdata")
  
  if (file.exists(outputs_batch_id)) {
    load(outputs_batch_id) # load AMIS output
    amis_params <- output$amis_params
    ess <- output$ess
    n_success <- length(which(ess>=amis_params[["target_ess"]]))
    n_failure <- sum(ess<amis_params[["target_ess"]])
    n_sim <- nrow(output$weight_matrix)
    n_iters <- ncol(output$ess_per_iteration)
    n_locs <- ncol(output$weight_matrix)
    summary_id <- c(id, n_locs, n_failure, n_success, n_sim, n_iters)
    summary_fitted <- rbind(summary_fitted, summary_id)
  }else{
    summary_not_fitted <- rbind(summary_not_fitted, id)
  }
  
  
}
summary_fitted <- as.data.frame(summary_fitted)
colnames(summary_fitted) <- c("id","n_locs","n_failure","n_success","n_sim","n_iters")
rownames(summary_fitted) <- NULL
summary_fitted
write.csv(summary_fitted, file = "summary_fitted.csv", row.names = F)

rownames(summary_not_fitted) <- NULL
write.csv(summary_not_fitted, file = "summary_not_fitted.csv", row.names = F)
summary_not_fitted

