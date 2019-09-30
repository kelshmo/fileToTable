# query<-synTableQuery("SELECT * FROM syn11801978")
# sm<-as.data.frame(query)
# sm<-as_tibble(sm)
#if there are no changes don't push a new version
#remove edit permissions from table
#write to file and save new table version at the same time
#use primary key Sample ID or Data ID

#check for version updates (possible to check against stored versions? How to store previous versions? How to account for new files added?)
#if there is a new version, readr::read_csv
library(tibble)
library(synapser)
library(dplyr)
library(tidyverse)
library(synapserutils)
library(purrr)

source("process_functions.R")

#' Get files from a folder in Synapse 
#' 
#' `synGetChildren` returns the metadata (properties) for each file in a folder when `includeTypes = list("file)`
#' 
#' This is a synapser function: https://r-docs.synapse.org/reference/synGetChildren.html

folder <- synGetChildren(parent = c("syn16809995"), includeTypes = list("file"))$asList()

#' Get file synIds
synFiles <- tibble(id = unlist(lapply(folder, function(x) x$id)))

selectCols <- c("Individual_ID", "Brain_ID", "SCZ_Pair", "BP_Pair", "Sex", "Ethnicity", "Age_of_Death", "Dx", "Institution_Dissection_ID", "Sample_ID", "assayType", "organism", "Exclude")
testCols <- paste0(selectCols, collapse = "|")
merged <- get_files(synFiles)

withType <- merged %>% 
  mutate(assay = purrr::map_chr(merged$filename, ~ stringr::str_split_fixed(., "[_\\.]", 5)[3])) %>% 
  mutate(species = purrr::map_chr(merged$filename, ~ stringr::str_split_fixed(., "[_\\.]", 5)[2]))

withType$filecontents <- purrr::map2(withType$filecontents, withType$assay, ~ dplyr::mutate(.x, assayType = .y))
withType$filecontents <- purrr::map2(withType$filecontents, withType$species, ~ dplyr::mutate(.x, organism = .y))

test_SM <- withType$filecontents %>%
  purrr::map(., ~ dplyr::rename_all(., funs(gsub("Sample_RNA_ID", "Sample_ID",.)))) %>% 
  purrr::map(., ~ dplyr::rename_all(., funs(gsub("Sample_DNA_ID", "Sample_ID",.)))) %>% 
  purrr::map(., ~ dplyr::rename_all(., funs(gsub("Assay_Sample_ID", "Sample_ID",.)))) %>%
  purrr::map(., ~ dplyr::select(., one_of(selectCols))) %>% 
  bind_rows %>% 
  select("organism", "Individual_ID","SCZ_Pair", "BP_Pair","Brain_ID", "Institution_Dissection_ID", everything())

# tb <- synBuildTable("All Characterized Assays", "syn1867134", test_SM)
# synStore(tb)

##assertr
#check every file individually with assertr before putting into sample master table join
not.empty.p <- function(x) if (x=="") return(FALSE)
indID.regex <- function(x) grepl('^CMC_\\w\\w\\w\\w_\\d\\d\\d',x)
four.digit.ID <- function(x) nchar(x)==4
two.digit <- function(x) nchar(x)==2

testassertr <- masterTable$filecontents[[1]] %>%
        assert(indID.regex, `Individual ID`) %>% 
        assert(four.digit.ID, Brain_ID) %>% 
        assert(two.digit, Age_of_Death)
        assert(in_set(gender), Gender)
        assert(is_uniq(Data_ID))
        assert(not.empty.p, Individual_ID)
        assert(in_set(brainbank), Institution)

brainbank <- c("NIMH-HBCC","HBCC","MSSM","Penn","Pitt")
gender <- c("Male","Female")
funding <- c("BX002395","Generated at HBCC","Grant Supplement-Sklar") #not a complete list
assay_schema <- c("ATACSeq", "HI-C", "ChIPSeq", "rnaSeq","wholeGenomeSeq", "rnaArray", "mmPCRSeq", "snpArray")
assayTarget_schema <- c("H3K27ac","H3K27me3","H3K4me3","input","")
celltype_schema <- c("GLUtamatergic neurons", "NeuN+", "NeuN-","GABAergic neurons","oligodendrocyte")
exclude_schema <- c(1, "")


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



masterCols <- c("Individual_ID","Brain_Bank","HBCC_Brain_ID", "NDA_GUID","SCZ_Pair","BP_Pair","Gender","Ethnicity","Age_of_Death","Dx","Funding","ASSAY","ASSAY_Target","Tissue","Cell_Type","Dissection_ID","Sample_ID","Data_ID","Exclude","Exclude_Reason", "QC_Metric","QC_Value","Biomaterials_Available")
selectCols <-c("Individual ID","Institution","Brain ID","SCZ Pair", "BP Pair","Ethnicity","Age of Death","Dx","Institution Dissection ID", "Institution Source ID","Sample DNA ID","Brain Region", "Cell Type","Exclude?","Exclude Reason","Sample RNA ID", "Assay_Sample_ID")
datatype = c("Clinical", "brainRegion", "Isolation", "Genotyping", "WGS", "MicroArray", "RNAseq","ATACseq")
seq_1 = c("Clinical","brainRegion")
seq_2 = c("Isolation")





####with consistent naming can rely on datatype in position 3, therefore code doesnt need to be modified as types of data are increased**************
collapseRows <- function(data, datatype = c(), seq_1 = c()) {
  next_iteration <- purrr::map(datatype, ~ data$filecontents[grep(.,data$filename)]) %>%
    setNames(datatype) %>% 
    purrr::map(., ~ dplyr::bind_rows(.))
  first <- next_iteration[seq_1] %>% 
    reduce(full_join, by = "Individual_ID")
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
