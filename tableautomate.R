#if there are no changes don't push a new version
#remove edit permissions from table
#write to file and save new table version at the same time
#use primary key Sample ID or Data ID

#syncFromSynapse()
#check for version updates (possible to check against stored versions? How to account for new files added?)
#if there is a new version, readr::read_csv
#

library(tibble)
library(synapser)
library(dplyr)
library(tidyverse)
library(synapserutils)

#check for new files 
check_newFiles <- function(x){}

testSync <- syncFromSynapse("syn2946054")

#testSync[[1]]$path

testSyncTibble<- as.tibble(seq(testSync))

files <- testSyncTibble %>% 
  mutate(filename = purrr::map_chr(testSync, function(x) x$get('name'))) %>% 
  mutate(filecontents = purrr::map(testSync, function(x) readr::read_csv(x$path,
                                                                       col_types = cols (.default = "c"),
                                                                       na=c("","NA")))) %>% 
  mutate(version = purrr::map(testSync, function(x) x$get('versionNumber')))


###TEMP
get_files <- function(SYN_LIST){
  files <- SYN_LIST %>% 
    mutate(thefile=purrr::map(id,synGet)) %>% 
    mutate(filename=purrr::map_chr(thefile, function(x) x$get('name'))) %>% #map_chr to force character values
    mutate(filecontents=purrr::map(thefile, function(x) readr::read_csv(x$path,
                                                                        col_types = cols (.default = "c"),
                                                                        na=c("","NA")))) %>% 
    mutate(version=purrr::map(thefile, function(x) x$get('versionNumber')))
  files
}
###


###assertr
brainbank <- c("HBCC","MSSM","PENN","PITT")
dx <- c("AFF", "BP", "Control","MDD","SCZ")
funding <- c("BX002395","Generated at HBCC","Grant Supplement-Sklar") #not a complete list 
###

SYN_LIST<- tibble(id=c("syn2279441","syn5908858","syn5666154","syn8048017","syn2929048","syn2279442","syn3153298","syn2279444","syn5666155","syn5666156","syn2279445","syn2929052","syn2279443","syn2929053","syn5650452","syn5666157","syn5666158","syn5666159","syn2279446","syn5698227","syn5698228"))

masterTable <- get_files(SYN_LIST)
masterTable$version #how to store previous versions? 
masterCols <- c("Individual_ID","Brain_Bank","HBCC_Brain_ID", "NDA_GUID","SCZ_Pair","BP_Pair","Gender","Ethnicity","Age_of_Death","Dx","Funding","ASSAY","ASSAY_Target","Tissue","Cell_Type","Dissection_ID","Sample_ID","Data_ID","Exclude","Exclude_Reason", "QC_Metric","QC_Value","Biomaterials_Available")

selectCols <-c("Individual ID","Institution","Brain ID","SCZ Pair", "BP Pair","Ethnicity","Age of Death","Dx","Institution Dissection ID", "Institution Source ID","Sample DNA ID","Brain Region", "Cell Type","Exclude?","Exclude Reason","Sample RNA ID")


test <- masterTable
#need to specify join columns, need to join sequentially


datatype = c("Clinical|Brain_Region")
#parse into set1 and set2 
parse <- function(data, set1 = c("Clinical|Brain_Region")) {
        first <- data[grep(set1,data$filename),]
        other <- data$filecontents[!grep(set1,data$filename),]

        #add list bit 
}





selectCols <-c("Individual ID","Institution","Brain ID","SCZ Pair", "BP Pair","Ethnicity","Age of Death","Dx","Institution Dissection ID", "Institution Source ID","Sample DNA ID","Brain Region", "Cell Type","Exclude?","Exclude Reason","Sample RNA ID")

bind_like <- function(data, selectCols = c()) {
        bound <- data$filecontents %>% 
                reduce(bind_rows) %>% 
                select(.,one_of(selectCols))
        bound
}

brain_region <- test[grepl("Brain_Region", test$filename),]
isolation <- test[grepl("Isolation",test$filename),]
rnaseq <- test[grepl("RNAseq", test$filename),]
microarray <- test[grepl("MicroArray", test$filename),]
atacseq <- test[grepl("ATACseq", test$filename),]

#all like types need to be row bound first
bind_br <- bind_like(brain_region, selectCols)
bind_isolation <- bind_like(isolation, selectCols)
bind_rnaseq <- bind_like(rnaseq, selectCols)
bind_microarray <- bind_like(microarray, selectCols)
bind_atacseq <- bind_like(atacseq, selectCols)




firstjoin <-testnest[[1]]$filecontents %>% 
        purrr::map(., ~ dplyr::select(.,one_of(selectCols))) %>% 
        purrr::map(., ~ dplyr::rename_all(.,funs(gsub(" ","_",.)))) %>% 
        reduce(full_join)


#column naming convention is new_name = current_name
#cannot select and rename in the same function. In fact, select() on multiple columns only works in the following format:
subset_cols <- test %>% 
        purrr::map(test$filecontents, ~ dplyr::select(.,one_of(selectCols))) %>% 
        purrr::map(test$filecontents, ~ dplyr::rename_all(.,funs(gsub(" ","_",.))))

firstJoin <- seq_join(subset_cols$filecontents, datatype = c("Clinical","Brain_Region"), join_colums = c("Individual_ID, SCZ_Pair, BP_Pair"))



#joins must be sequential based on type 

datatype = c("Clinical|Brain_Region")
join_colums = c("Individual_ID, SCZ_Pair, BP_Pair")
files_of_interest<-subset_cols[grep(c(datatype),data$filename)]

new <- files_of_interest %>%
        reduce(full_join, by = join_columns)

seq_join <- function(data, datatype = c(), join_columns = c()) {
  
  files_of_interest<-data$filecontents[grep(c(datatype),data$filename)]
  
  new <- files_of_interest %>%
          purrr::map(test$filecontents, ~ dplyr::select(.,one_of(selectCols))) %>% 
          purrr::map(test$filecontents, ~ dplyr::rename_all(.,funs(gsub(" ","_",.))))
    reduce(full_join, by = join_columns)
  
  as.tibble(new)
}



to_table <- function(data, masterCols = c() ){
  drop <- purrr::map(data, ~ dplyr::select(.,one_of(selectCols))) %>% 
    purrr::map(~ dplyr::rename_all(.,funs(gsub(" ","_",.)))) %>% 
    reduce(full_join) %>% 
    group_by(Individual_ID,Institution,Brain_ID, SCZ_Pair, BP_Pair,Ethnicity,Age_of_Death,Dx,Institution_Dissection_ID, Institution_Source_ID,Sample_DNA_ID,Brain_Region, Cell_Type,`Exclude?`,Exclude_Reason,Sample_RNA_ID) %>% 
    filter(row_number() == 1)

  
  new <- files_of_interest %>% 
    purrr::map(~ dplyr::select(.x,-one_of(join_columns))) %>% 
    purrr::map(~ dplyr::rename_all(.,funs(gsub(" ","_",.)))) %>% 
    purrr::map(~ dplyr::rename_all(.,funs(paste0(label,":",.)))) %>% 
    purrr::map2(temp_drop, ~ cbind(.y,.x)) %>%
    reduce(bind_rows)
}
 

for (i in 1:dim(masterTable)[1]) {
  
  print(i)
  print (colnames(masterTable$filecontents[[i]]))
}

aggregate_excludes <- function (x){}





##################test how row versioning works when files to table could change underlying data structure 
####use case 1 rows added
clinical <- masterTable$filecontents[[1]]
samp_clinical <- clinical[1:100,]
tissue <- masterTable$filecontents[[5]]
samp_tissue <- tissue[1:100,]

collapse <- samp_clinical %>% 
  full_join(samp_tissue) %>% 
  select(`Individual ID`,`Brain ID`,`Institution Dissection ID`)

newrow <- samp_tissue[5,]
newrow[,2]<-c("CMC_MSSM_dummy")
samp_tissue <- rbind(samp_tissue,newrow)

currentTable <- synTableQuery("SELECT * FROM syn15665630")
ct<-as.data.frame(currentTable)
ct<-as.tibble(ct)

collapse_w_change <- samp_clinical %>% 
  full_join(samp_tissue) %>% 
  select(`Individual ID`,`Brain ID`,`Institution Dissection ID`)

synStore(Table("syn15665630",collapse_w_change))

changeTable <- synTableQuery("SELECT * FROM syn15665630")
changet <- as.data.frame(changeTable)
changet <- as.tibble(changet)


bttest<-masterTable$filecontents[[5]]
#standardized naming happens once values are joined.....how to account for sample RNA ID and sample DNA ID
stand_naming <-function(data){
        data <- data %>%
                purrr::map(~ dplyr::select(., 
                                           HBCC_Brain_ID = "Brain ID")) %>% 
                purrr::map(~ dplyr::select(.,
                                           scz = "SCZ_Pair"))
        data
}

query<-synTableQuery("SELECT * FROM syn11801978")
sm<-as.data.frame(query)
sm<-as_tibble(sm)

