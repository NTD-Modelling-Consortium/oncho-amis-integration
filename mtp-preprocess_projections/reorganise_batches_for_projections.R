
#### This file creates the MDA and VC histories for each unique treatment history
#### Relies on outputs from multipletimepoints_preprocess_maps.R

# cpu <- Sys.getenv("MODULE_CPU_TYPE")
# setwd(paste0("../oncho-mtp-",cpu,"/EPIONCHO-IBM"))
# setwd("~/Desktop/Oncho/oncho-mtp-EK/oncho-mtp-skylake/EPIONCHO-IBM/")
cpu <- "skylake"
setwd(paste0("../oncho-mtp-",cpu,"/EPIONCHO-IBM"))
library(dplyr)

# Save MDA files 
if (!dir.exists("model_output")) {dir.create("model_output")}

# Load table that linked IU to batch id in the fitting process
table_iu_idx <- read.csv("../../outputs/table_iu_idx.csv")
ord <- order(table_iu_idx$IU_CODE)
table_iu_idx <- table_iu_idx[ord, ]
# Define new batches for projections

num_IUs_per_batch <- 12
num_batches <- ceiling(nrow(table_iu_idx)/num_IUs_per_batch)
for(id in 1:num_batches){
  wh <- 12*(id-1) + 1:12
  IUs <- unique(table_iu_idx[wh,"IU_CODE"])
  IUs <- matrix(IUs,ncol=1)
  iu_file <- file.path("model_output",paste0("IUs_MTP_proj_",1000+id,".csv"))
  write.table(IUs, file=iu_file, row.names=F, col.names = F, quote=F, sep=",")# write input parameter file
}
