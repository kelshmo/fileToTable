# CMC merged metadata
# Author: Kelsey Montgomery
# Date: 2019.07.24
library(synapser)
library(dplyr)
library(tidyverse)
library(purrr)
library(tibble)
library(knitr)
library(githubr)

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

write.csv(full, "../files/CMC_Human_ChIPSeq_mergedMetadata.csv", row.names = FALSE)

