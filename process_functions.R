# Merge Synapse files into a Synapse table
# Author: Kelsey Montgomery
# Date: 2019-04-10

#' get_files
#' 
#' This function uses synIds to download files from Synapse
#' 
#' @param synFile a tibble of synIds in column id 
#' 
#' @return a tibble of synIds in column id, file property metadata in column thefile, .csv filename in column filename and nested tibbles of .csv data in column filecontents. 



#could modify this function to skip making S3: File its own column in thefile and instead attach that to filecontents and pull filename from the synGetChildren bit 
#consider getting version number   mutate(version = purrr::map(Sync, function(x) x$get('versionNumber'))) to check if there should be a change.

getFiles <- function(synFile = tibble()){
  files <- synFile %>% 
    mutate(thefile=purrr::map(id,synGet)) %>% 
    mutate(filename=purrr::map_chr(thefile, function(x) x$get('name'))) %>% #map_chr to force character values
    mutate(filecontents=purrr::map(thefile, function(x) readr::read_csv(x$path,
                                                                        col_types = cols (.default = "c"),
                                                                        na=c("","NA")))) 
  files
}
