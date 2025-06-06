#### This file assigns batch IDs by finding IUs with the same treatment histories and generates the map with multiple time points and saves the history for each batch

library(dplyr)
library(tidyr)
library(magrittr)

kPathToInputs <- Sys.getenv("PATH_TO_FITTING_PREP_INPUTS")
kPathToMaps <- file.path(Sys.getenv("PATH_TO_FITTING_PREP_ARTEFACTS"), "Maps")

# load data and histories
mda_file <- readRDS(file.path(kPathToInputs, "Full_histories_df_popinfo_ALL_minimal_070425_listlabels.rds"))
mda_file[mda_file == "NA"] <- NA

# around 2000 IUs weren't included in the endgame (multiple time point fitting), and now there are around 600 additional IUs being included (still ~1500 unfitted to multiple time points)
load(file.path(kPathToInputs, "ALL_prevalence_map.rds")) # load baseline map that has all IUs (endgame and non-endgame)
# convert IU_CODE to IUID
map_all <- map_all %>%
        select(-TaskID) %>% # remove old taskID
        rename(IUID = IU_CODE) %>%
        mutate(IUID = paste0(substr(IUID, 1, 3), substr(IUID, (nchar(IUID) - 4), nchar(IUID)))) %>%
        filter(IUID %in% mda_file$IUID)

prevalence_map <- as.matrix(map_all[, (ncol(map_all) - 1000 + 1):ncol(map_all), drop = FALSE])
rownames(prevalence_map) <- map_all$IUID

# find unique mda histories
mda_history <- mda_file %>%
        filter(Year <= 2018) %>%
        filter(IUID %in% rownames(prevalence_map)) %>%
        select(IUID, Year, number_rnds) %>%
        pivot_wider(names_from = "Year", values_from = "number_rnds")

unique_histories <- mda_history %>%
        group_by(across(c(-IUID))) %>%
        summarise(ius = list(IUID)) %>%
        ungroup() %>%
        mutate(index = row_number())

histories_indexes <- lapply(1:nrow(unique_histories), function(i) data.frame(unique_histories[["ius"]][[i]], unique_histories[["index"]][i]))
histories_indexes <- do.call(rbind, histories_indexes)
colnames(histories_indexes) <- c("IU", "histories_index")

# find unique coverage histories
coverage_history <- mda_file %>%
        filter(Year <= 2018) %>%
        filter(IUID %in% rownames(prevalence_map)) %>%
        mutate(Cov.in2 = ifelse(number_rnds > 0, Cov.in2, NA)) %>%
        select(IUID, Year, Cov.in2) %>%
        pivot_wider(names_from = "Year", values_from = "Cov.in2")

unique_coverage <- coverage_history %>%
        group_by(across(c(-IUID))) %>%
        summarise(ius = list(IUID)) %>%
        ungroup() %>%
        mutate(index = row_number())

coverage_indexes <- lapply(1:nrow(unique_coverage), function(i) data.frame(unique_coverage[["ius"]][[i]], unique_coverage[["index"]][i]))
coverage_indexes <- do.call(rbind, coverage_indexes)
colnames(coverage_indexes) <- c("IU", "coverage_index")

# find unique adherence histories
adherence_history <- mda_file %>%
        filter(Year <= 2018) %>%
        filter(IUID %in% rownames(prevalence_map)) %>%
        select(IUID, Year, adherence_par) %>%
        pivot_wider(names_from = "Year", values_from = "adherence_par")

unique_adherence <- adherence_history %>%
        group_by(across(c(-IUID))) %>%
        summarise(ius = list(IUID)) %>%
        ungroup() %>%
        mutate(index = row_number())

adherence_indexes <- lapply(1:nrow(unique_adherence), function(i) data.frame(unique_adherence[["ius"]][[i]], unique_adherence[["index"]][i]))
adherence_indexes <- do.call(rbind, adherence_indexes)
colnames(adherence_indexes) <- c("IU", "adherence_index")

# find unique vector control histories
vc_history <- mda_file %>%
        filter(Year <= 2018) %>%
        filter(IUID %in% rownames(prevalence_map)) %>%
        select(IUID, Year, vector_control) %>%
        pivot_wider(names_from = "Year", values_from = "vector_control")

unique_vc <- vc_history %>%
        group_by(across(c(-IUID))) %>%
        summarise(ius = list(IUID)) %>%
        ungroup() %>%
        mutate(index = row_number())

vc_indexes <- lapply(1:nrow(unique_vc), function(i) data.frame(unique_vc[["ius"]][[i]], unique_vc[["index"]][i]))
vc_indexes <- do.call(rbind, vc_indexes)
colnames(vc_indexes) <- c("IU", "vc_index")

# join all and find unique combinations
all_treatments <- data.frame(IU = rownames(prevalence_map)) %>%
        left_join(histories_indexes) %>%
        left_join(coverage_indexes) %>%
        left_join(adherence_indexes) %>%
        left_join(vc_indexes)

unique_all <- all_treatments %>%
        group_by(across(c(-IU))) %>%
        summarise(ius = list(IU)) %>%
        ungroup() %>%
        mutate(index = row_number())

# make lookup table between IUs and batches
iu_task_lookup <- lapply(1:nrow(unique_all), function(i) data.frame(unique_all[["ius"]][[i]], unique_all[["index"]][i]))
iu_task_lookup <- do.call(rbind, iu_task_lookup)
colnames(iu_task_lookup) <- c("IUID", "TaskID")

# split up batch 1 (IUs with no treatment) into 3 groups based on baseline prevalence
batch_with_no_treatment <- 1
batch_to_split_mapped_prev <- prevalence_map[which((rownames(prevalence_map) %in% iu_task_lookup$IUID[iu_task_lookup$TaskID == batch_with_no_treatment])), ]
batch_to_split_mapped_prev_avg <- apply(batch_to_split_mapped_prev, 1, mean)
batch_to_split_mapped_prev_avg_quantiles <- quantile(batch_to_split_mapped_prev_avg, probs = c(0.25, 0.75))
lower_prev_batch <- names(batch_to_split_mapped_prev_avg)[which(batch_to_split_mapped_prev_avg < batch_to_split_mapped_prev_avg_quantiles[1])]
upper_prev_batch <- names(batch_to_split_mapped_prev_avg)[which(batch_to_split_mapped_prev_avg > batch_to_split_mapped_prev_avg_quantiles[2])]

# assign new ID to IUs being refitted
new_taskID_middlegroup <- max(iu_task_lookup$TaskID)
iu_task_lookup$TaskID[which(iu_task_lookup$IUID %in% lower_prev_batch)] <- new_taskID_middlegroup + 1
iu_task_lookup$TaskID[which(iu_task_lookup$IUID %in% upper_prev_batch)] <- new_taskID_middlegroup + 2

# note that these batches don't necessarily match the batch IDs from previous multiple time point fitting
save(iu_task_lookup, file = file.path(kPathToMaps, "iu_task_lookup.rds"))

# join batch IDs to map
map_all_mtp <- iu_task_lookup %>%
        left_join(map_all)
map_all_mtp <- list(map_all_mtp)


# add 2000 map
map_2000 <- read.csv(file.path(kPathToMaps, "maps_joint/samples_2000.csv"), header = T) %>%
        mutate(IU_ID = paste0(substr(IU_code, 1, 3), substr(IU_code, (nchar(IU_code) - 4), nchar(IU_code)))) %>%
        select(-X) %>%
        select(IU_ID, starts_with("X"))
map_all_mtp[[2]] <- map_all_mtp[[1]] %>%
        select(IUID, TaskID) %>%
        left_join(map_2000, by = c("IUID" = "IU_ID"))

# add 2018 map
map_2018 <- read.csv(file.path(kPathToMaps, "maps_joint/samples_2018.csv"), header = T) %>%
        mutate(IU_ID = paste0(substr(IU_code, 1, 3), substr(IU_code, (nchar(IU_code) - 4), nchar(IU_code)))) %>%
        select(-X) %>%
        select(IU_ID, starts_with("X"))
map_all_mtp[[3]] <- map_all_mtp[[1]] %>%
        select(IUID, TaskID) %>%
        left_join(map_2018, by = c("IUID" = "IU_ID"))

save(map_all_mtp, file = file.path(kPathToMaps, "ALL_prevalence_map_multipletimespoints.rds"))

kPathToModelOutput <- file.path(Sys.getenv("PATH_TO_FITTING_PREP_ARTEFACTS"), "model_output")

# Save MDA files
if (!dir.exists(kPathToModelOutput)) {
        dir.create(kPathToModelOutput, recursive = TRUE)
}
for (id in 1:max(as.integer(iu_task_lookup$TaskID))) {
        iu_file <- file.path(kPathToModelOutput, paste0("IUs_MTP_", id, ".csv"))
        map_all_mtp_id <- map_all_mtp[[1]] %>%
                filter(TaskID == id)
        ius <- matrix(map_all_mtp_id$IUID, ncol = 1)
        write.table(ius, file = iu_file, row.names = F, col.names = F, quote = F, sep = ",") # write input parameter file

        for (iu in ius[1]) { # histories are all the same within a batch so only need to generate once
                mda_file_iu <- mda_file %>%
                        filter(IUID == iu & Year <= 2018) %>%
                        filter(number_rnds > 0 & Cov.in2 > 0) %>%
                        dplyr::select(Year, number_rnds, Cov.in2, adherence_par) %>%
                        mutate(treatment_interval = 1 / number_rnds)
                colnames(mda_file_iu) <- c("Year", "number_rnds", "ModelledCoverage", "adherence_par", "treatment_interval")
                mda_path <- file.path(kPathToModelOutput, paste0("InputMDA_MTP_", id, ".csv"))
                write.csv(mda_file_iu, file = mda_path, row.names = F) # write input MDA file

                vc_file_iu <- mda_file %>%
                        filter(IUID == iu & Year <= 2018) %>%
                        filter(vector_control > 0) %>%
                        dplyr::select(Year, vector_control)
                vc_all_years <- data.frame(Year = 1975:2025) %>%
                        left_join(vc_file_iu, by = c("Year"))
                vc_all_years$vector_control[is.na(vc_all_years$vector_control)] <- 0
                vc_all_years %<>%
                        mutate(abr_multiplier = case_when(
                                vector_control == 1 ~ 0.2,
                                vector_control == 2 ~ 0,
                                TRUE ~ 1
                        ))
                vc_path <- file.path(kPathToModelOutput, paste0("InputVC_MTP_", id, ".csv"))
                write.csv(vc_all_years, file = vc_path, row.names = F) # write input vector control file
        }
}
