query<-synTableQuery("SELECT * FROM syn11801978")
sm<-as.data.frame(query)
sm<-as_tibble(sm)
#if there are no changes don't push a new version
#remove edit permissions from table
#write to file and save new table version at the same time
#use primary key Sample ID or Data ID

#syncFromSynapse()
#check for version updates (possible to check against stored versions? How to store previous versions? How to account for new files added?)
#if there is a new version, readr::read_csv
library(tibble)
library(synapser)
library(dplyr)
library(tidyverse)
library(synapserutils)
library(purrr)

##check for new files 
check_newFiles <- function(x){}

##synFromSynapse to get files
#Sync <- syncFromSynapse("syn2946054")
#Specify parentsIDs not interested in using, problem is then have to go to bottom of folder structure 
test <- tempSync[lapply(tempSync, function(x) x$get('parentId')) != c("syn4565094")]

files <- tibble(
  filename = purrr::map_chr(Sync, function(x) x$get('name'))
) %>% 
  mutate(filecontents = purrr::map(Sync, function(x) readr::read_csv(x$path,
                                                                       col_types = cols (.default = "c"),
                                                                       na=c("","NA")))) %>% 
  mutate(version = purrr::map(Sync, function(x) x$get('versionNumber')))



##TEMP get_files
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

##select relevant columns
colSubset <- function(data, selectCols = c()){
  data$filecontents <- data$filecontents %>% 
    purrr::map(., ~ dplyr::select(.,one_of(selectCols))) %>% 
    purrr::map(., ~ dplyr::rename_all(.,funs(gsub(" ","_",.)))) %>% 
    purrr::map(., ~ dplyr::rename_all(.,funs(gsub("Sample_RNA_ID", "Sample_ID",.)))) %>% 
    purrr::map(., ~ dplyr::rename_all(.,funs(gsub("Sample_DNA_ID", "Sample_ID",.)))) %>% 
    purrr::map(., ~ dplyr::rename_all(.,funs(gsub("Assay_Sample_ID", "Sample_ID",.))))
  data
}

##assertr
brainbank <- c("HBCC","MSSM","PENN","PITT")
dx <- c("AFF", "BP", "Control","MDD","SCZ")
funding <- c("BX002395","Generated at HBCC","Grant Supplement-Sklar") #not a complete list 



SYN_LIST<- tibble(id=c("syn2279441","syn5908858","syn5666154","syn8048017","syn2929048","syn2279442","syn3153298","syn2279444","syn5666155","syn5666156","syn2279445","syn2929052","syn2279443","syn2929053","syn5650452","syn5666157","syn5666158","syn5666159","syn2279446","syn5698227","syn5698228"))
masterTable <- get_files(SYN_LIST)
masterTable$version
masterCols <- c("Individual_ID","Brain_Bank","HBCC_Brain_ID", "NDA_GUID","SCZ_Pair","BP_Pair","Gender","Ethnicity","Age_of_Death","Dx","Funding","ASSAY","ASSAY_Target","Tissue","Cell_Type","Dissection_ID","Sample_ID","Data_ID","Exclude","Exclude_Reason", "QC_Metric","QC_Value","Biomaterials_Available")
selectCols <-c("Individual ID","Institution","Brain ID","SCZ Pair", "BP Pair","Ethnicity","Age of Death","Dx","Institution Dissection ID", "Institution Source ID","Sample DNA ID","Brain Region", "Cell Type","Exclude?","Exclude Reason","Sample RNA ID", "Assay_Sample_ID")
datatype = c("Clinical", "Brain", "Isolation", "Genotyping", "WGS", "MicroArray", "RNAseq","ATACseq")
seq_1 = c("Clinical","Brain")
seq_2 = c("Isolation")





####with consistent naming can rely on datatype in position 3, therefore code doesnt need to be modified as types of data are increased**************
collapseRows <- function(data, datatype = c(), seq_1 = c()) {
  next_iteration <- purrr::map(datatype, ~ data$filecontents[grep(.,data$filename)]) %>%
    setNames(datatype) %>% 
    purrr::map(., ~ dplyr::bind_rows(.))
  first <- next_iteration[seq_1] %>% 
    reduce(full_join)
  second <- next_iteration[-which(names(next_iteration) %in% seq_1)] %>% 
    reduce(full_join)
  finalJoin <- full_join(first,second, by = "Institution_Dissection_ID", "Sample_ID")
}


temp <- first %>% 
  full_join(next_iteration$Isolation) %>% 
  full_join(next_iteration$Genotyping) %>% 
  full_join(next_iteration$WGS) %>% 
  full_join(next_iteration$MicroArray) %>% 
  full_join(next_iteration$RNAseq)


###concatenate exclude reason horizontally across isolation and assays; If I proceed with this need to do a check that exclude values of 1 populate every reason 
#select only Sample_ID,
#try to take out institution dissection ID

#try coallese 
concatenateCols <- c("Sample_ID", "Exclude_Reason")
concatenate_excludes <- next_iteration[-which(names(next_iteration) %in% seq_1)] %>% 
  purrr::map(., ~ dplyr::select(.,one_of(concatenateCols))) %>% 
  reduce(full_join, by= c("Sample_ID")) %>% 
  unite(Exclude_Reason_Summary, -one_of(c("Sample_ID")), sep = "_") #can use case_when in this situation?? 

aggregate_excludes <- function (x){}






##Call functions
working <- colSubset(masterTable, selectCols)
finalJoin <- collapseRows(working, datatype, seq_1)
data <- working


#check for lost IDs: seems like the left join is functioning properly, though it drops sample IDs that are not mapped to a isolation ID. Could enforce a check to make sure 
#those IDs are associated with an exclusion value and if so, still push the table. For IDs that do not have an exclusion value, an error will be thrown. 
checkEx <- next_iteration[-which(names(next_iteration) %in% seq_1)] %>% 
  reduce(bind_rows)
View(checkEx[(checkEx$Sample_ID %in% missingDID),])
missingDID <- test$Sample_ID[is.na(test$Individual_ID) & is.na(test$Institution_Dissection_ID)]
tempsamplelist<-c(sample1,sample2,sample3,sample4)
check <- as.tibble(tempsamplelist)
names(check) <- c("Sample_ID")
test <- check %>% 
  left_join(finalJoin)

#also check whether isolation IDs with no Sample ID are present in finalJoin: they are, so do we want these to stay included? 
checkSID <- finalJoin %>% 
  filter(is.na(Sample_ID))

#check duplicate IDS in RNAseq ACC, all duplicate ACC values have 1 of 2 entries marked as exclude
rnaseqacc<-masterTable$filecontents[[14]]
rnaseqpfc<-masterTable$filecontents[[19]]
dupACC<- rnaseqacc %>% 
  count(`Sample RNA ID`) %>% 
  filter(n>1)

duplicated_rnaseqacc<-rnaseqacc[(rnaseqacc$`Sample RNA ID` %in% dupACC$`Sample RNA ID`),]
check <- duplicated_rnaseqacc %>% 
  filter(is.na(`Exclude?`)) %>% 
  count(`Sample RNA ID`)

#check redundant sampleIDs 
test <- finalJoin %>% 
  filter(is.na(`Exclude?`)) %>% 
  group_by(Sample_ID) %>% 
  summarize(num_SampleID = n()) %>% 
  filter(num_SampleID != 1)


x<-output
for (i in 1:length(x)) {
  
  print(i)
  print (dim(x[[i]]))
}

#http://zevross.com/blog/2014/04/30/mini-post-for-large-tables-in-r-dplyrs-function-inner_join-is-much-faster-than-merge/
#https://blog.dominodatalab.com/a-quick-benchmark-of-hashtable-implementations-in-r/
#https://codereview.stackexchange.com/questions/94253/identify-changes-between-two-data-frames-explain-deltas-for-x-columns
#test method for tracking new and old changes
df.new <-head(mtcars, n=10)
newrow <- mtcars[30,]
df.old <-head(mtcars, n=10) %>% 
  mutate(mpg = replace(mpg, which(qsec == 18.61), 100)) %>% 
  bind_rows(newrow)

df.changes <- function(df.old, df.new, 
                       KEYS = c("id"), 
                       VAL = NULL, 
                       retain.columns = NULL) {
N <- transform(df.new, is = TRUE)
O <- transform(df.old, is = TRUE)

M <- merge(N,O, by = KEYS, all = TRUE, suffixes = c(".new", ".old"))
M$is.new <- !is.na(M$is.new)
M$is.old <- !is.na(M$is.old)

O <- M[KEYS]

O$row.changed <- with(M, ifelse(is.old & is.new, "10.Retained",
                                ifelse(is.old, "05.Lost", 
                                       "00.New")))
original.vars <- setdiff(names(df.new), KEYS)



}
