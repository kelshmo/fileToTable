# CMC merged metadata
# Author: Kelsey Montgomery
# Date: 2019.07.24
library(synapser)
library(dplyr)
library(tidyverse)
library(purrr)
library(tibble)
library(knitr)
require(githubr)

source("process_functions.R")

SYN_LIST <- tibble(id = c("syn2279441","syn2279442","syn17114462"))

data <- get_files(SYN_LIST)

clinical <- data$filecontents[grep(c("clinical"),data$filename)] %>%
  purrr::map(~ dplyr::rename_all(.,funs(gsub(" ","_",.)))) %>%
  reduce(bind_rows)

brainRegion <- rename_bind(data,
                           datatype = "brainRegion",
                           label = "ChIPSeq_dissection",
                           join_columns = c("Individual ID", "Institution Dissection ID"))

ChIPSeq <- data$filecontents[[3]]

full <- ordered_join(clinical = clinical,
                     brainRegion = brainRegion,
                     assay = ChIPSeq, 
                     datatype = c("clinical","brainRegion","ChIPSeq"), 
                     seq_1 = c("clinical", "brainRegion"))

full <- select_if(full, ~!all(is.na(.)))

nameFile <- c("../files/CMC_Human_ChIPSeq_mergedMetadata.csv")

write.csv(full, nameFile , row.names = FALSE)

# get github history for Provenance 
thisFilename <- c("mergedMetadata_script.R")
thisRepo <- githubr::getRepo(repository = "kelshmo/fileToTable")
thisFile <- githubr::getPermlink(repository = thisRepo, repositoryPath = thisFilename)

#describe activity for Provenance
activityName = "Join metadata"
activityDescription = "changes to metadata - new merge"

# push file to Synapse
file <- File(path = nameFile , parent = "syn16809995", versionComment = paste0("Merged ", date()))

ids <- paste0(SYN_LIST$id, collapse = "','")

synStore(file,
         activityName = activityName,
         activityDescription = activityDescription,
         used = c("syn2279441","syn2279442","syn17114462"),
         executed = thisFile)
         



