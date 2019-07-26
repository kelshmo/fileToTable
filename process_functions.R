# Merge Synapse files into a Synapse table
# Author: Kelsey Montgomery
# Date: 2019-04-10

#' get_files
#' 
#' This function uses synIds to download files from Synapse
#' 
#' @param synFile a tibble of synIds in column id tibble(id = c("syn####","syn####"))
#' @return a tibble of synIds in column id, file property metadata in column thefile, .csv filename in column filename and nested tibbles of .csv data in column filecontents. 

#NOTES: could modify this function to pull filename from the synGetChildren() R6 
#consider getting version number   mutate(version = purrr::map(Sync, function(x) x$get('versionNumber'))) to check if there should be a change.

get_files <- function(synFile = tibble()){
  files <- synFile %>% 
    mutate(thefile = purrr::map(id,synGet)) %>% 
    mutate(filename = purrr::map_chr(thefile, function(x) x$get('name'))) %>% #map_chr to force character values
    mutate(filecontents = purrr::map(thefile, function(x) readr::read_csv(x$path,
                                                                        col_types = cols(.default = "c"),
                                                                        na = c("","NA")))) 
  files
}

#' rename_bind
#' 
#' This function uses dataype = c() to parse the file of interest from get_files and rename headers appropriately for a join.
#' 
#' @param data a tibble of nested files.
#' @param datatype a string to grep() the filename.
#' @param label a string will be appended to column headers, except the "join by" column. 
#' @param join_columns a string matching the "join by" column name.
#' @return an unnested tibble of .csv with edits to column headers 
#' @examples rename_bind(data, datatype = "brainRegion", label = "RNASeq_report", join_columns = "Assay_Sample_ID")

rename_bind <- function(data, datatype = c(), label = c(), join_columns = c()) {
  
  files_of_interest <- data$filecontents[grep(c(datatype),data$filename)]
  
  temp_drop <- purrr::map(files_of_interest, ~ dplyr::select(.x,one_of(join_columns))) %>% 
    purrr::map(~ dplyr::rename_all(.,funs(gsub(" ","_",.))))
  
  new <- files_of_interest %>% 
    purrr::map(~ dplyr::select(.x,-one_of(join_columns))) %>% 
    purrr::map(~ dplyr::rename_all(.,funs(gsub(" ","_",.)))) %>% 
    purrr::map(~ dplyr::rename_all(.,funs(paste0(label,":",.)))) %>% 
    purrr::map2(temp_drop, ~ cbind(.y,.x)) %>%
    reduce(bind_rows)
  
  as_tibble(new)
}
#'ordered_join
#'
#'Data frames ordered for a join where the assay data defines the final dataset dimensions
#'
#'@param clinical 
#'@param brainRegion
#'@param isolation 
#'@param assay 
#'@param datatype
#'@param seq_1
#'@return a single dataframe
ordered_join <- function(clinical, brainRegion = NULL, isolation = NULL , assay = NULL, datatype = c(), seq_1 = c()){
  join_list <- tibble(datatype)
  tojoin <- join_list %>% 
    mutate(dtable = purrr::map(datatype, function(x) get(x)))
  firstJoin <- purrr::map(seq_1, ~ tojoin$dtable[grep(.,tojoin$datatype)]) %>%
    setNames(seq_1) %>% 
    flatten() %>% 
    reduce(full_join, by = "Individual_ID") %>% 
    select(Individual_Notes, Individual_ID, Institution, Brain_ID, SCZ_Pair, BP_Pair, everything())
  new_join_list <- c("firstJoin", datatype[!(datatype %in% seq_1)])
  join_list <- tibble(new_join_list)
  tojoin <- join_list %>% 
    mutate(dtable = purrr::map(new_join_list, function(x) get(x)))
  next_iteration <- reduce(tojoin$dtable, right_join)
  next_iteration
}
