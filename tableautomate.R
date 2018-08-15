#if there are no changes don't push a new version
#remove edit permissions from table
#write to file and save new table version at the same time
#use primary key Sample ID or Data ID 

library(tibble)
library(synapser)
library(dplyr)
library(tidyverse)
library(synapserutils)



testSync <- syncFromSynapse("syn1978914")

#testSync[[1]]$path

testSyncTibble<- as.tibble(seq(testSync))

files <- testSyncTibble %>% 
  mutate(filename = purrr::map_chr(testSync, function(x) x$get('name'))) %>% 
  mutate(filecontents = purrr::map(testSync, function(x) readr::read_csv(x$path,
                                                                       col_types = cols (.default = "c"),
                                                                       na=c("","NA")))) %>% 
  mutate(version = purrr::map(testSync, function(x) x$get('versionNumber')))

SYN_LIST<- tibble(id=c("syn2279441","syn5908858"))


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


###GENERAL REGEX
#var mapObj = {Institution:"Brain_Bank", Brain_ID:"HBCC_Brain_ID",Sample_RNA_ID:"Sample_ID"}
#once all metadata files are in list of tibbles, select variable for table 
stand_naming <-function(data){
  data <- data %>%
    purrr::map(~ dplyr::select(., 
                               HBCC_Brain_ID = "Brain ID")) %>% 
    purrr::map(~ dplyr::select(.,
                               scz = "SCZ_Pair"))
  data
}


purrr::map(~ dplyr::rename_all(.,funs(gsub(" ","_",.))))
  


#column naming convention is new_name = current_name
#perhaps do not rename until everything is bound/joined. 




#purrr::map(~ dplyr::rename_all(.,funs(sub('Brain_ID', 'HBCC_Brain_ID',.)))) 
cols_current = quos(Brain_ID, SCZ_Pair)
cols_new = quos(HBCC_Brain_ID, scz)
cols_list = list(cols_current,cols_new)
df1<-all$filecontents[[1]]

df <- df1 %>% rename_all(., funs(gsub(" ","_",.)))
df <- df %>% dplyr::select(HBCC_Brain_ID = "Brain_ID", scz = "SCZ_Pair", everything()) 



###assertr
brainbank <- c("HBCC","MSSM","PENN","PITT")
dx <- c("AFF", "BP", "Control","MDD","SCZ")
funding <- c("BX002395","Generated at HBCC","Grant Supplement-Sklar") #not a complete list 


###merged files to test automate(don't want to use merged files due to column names )
file1<-File(path='/Users/kelsey/Documents/CMC/table_automate/Metadata_merge_RNAseq.csv',parentId= 'syn15660098')
file2<-File(path='/Users/kelsey/Documents/CMC/table_automate/Metadata_merge_SNP.csv',parentId= 'syn15660098')
mergeRNA <- readr::read_csv("Metadata_merge_RNAseq.csv")

SYN_LIST<- tibble(id=c("syn2279441","syn5908858","syn5666154","syn8048017","syn2929048","syn2279442","syn3153298","syn2279444","syn5666155","syn5666156","syn2279445","syn2929052","syn2279443","syn2929053","syn5650452","syn5666157","syn5666158","syn5666159","syn2279446","syn5698227","syn5698228"))

masterTable <- get_files(SYN_LIST)
masterTable$version #how to store previous versions? 
masterCols <- c("Individual_ID","Brain_Bank","HBCC_Brain_ID", "NDA_GUID","SCZ_Pair","BP_Pair","Gender","Ethnicity","Age_of_Death","Dx","Funding","ASSAY","ASSAY_Target","Tissue","Cell_Type","Dissection_ID","Sample_ID","Data_ID","Exclude","Exclude_Reason", "QC_Metric","QC_Value","Biomaterials_Available")

selectCols <-c("Individual ID","Institution","Brain ID","SCZ Pair", "BP Pair","Ethnicity","Age of Death","Dx","Institution Dissection ID", "Institution Source ID","Sample DNA ID","Brain Region", "Cell Type","Exclude?","Exclude Reason","Sample RNA ID")
data<-masterTable$filecontents

#need to specify join columns, need to join sequentially

data<-masterTable[1:7,]
data<-data$filecontents


subset_cols <- purrr::map(data, ~ dplyr::select(.,one_of(selectCols))) %>% 
  purrr::map(~ dplyr::rename_all(.,funs(gsub(" ","_",.))))

seq_join(subset_cols, datatype = c("Clinical","Brain_Region"))



#joins must be sequential based on type 

seq_join <- function(data, datatype = c(), join_columns = c()) {
  
  files_of_interest<-data$filecontents[grep(c(datatype),data$filename)]
  
  new <- files_of_interest %>%
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
