#' This script generates a clean compound / MOA table. Compounds with identical
#' structure are mapped to a single unique compound ID. The file 
#' `data/Unique_Screened_Compound_MOA.csv` lists names of compounds that are
#' equivalent with respect to structure. Unique compound IDs are then merged 
#' with the MOA table `data/Screened_Compound_ID_MOA_added_20200420.csv`, 
#' which lists MOA by (non-unique) compound IDs. Non-unique compound IDs
#' are taken from the raw data (e.g. by vendor). Compounds with no MOA
#' annotation are dropped from analysis.
#' 
#' Last updated 10/4/2022, Karl Kumbier
library(data.table)
library(tidyverse)

ccl.dir <- Sys.getenv('CCL_DIR')
output.file <- str_c(ccl.dir, 'data/Table1_compound_library.csv')

# Load compound ID table - compounds grouped by structure
xcpd <- fread(str_c(ccl.dir, 'data/Unique_Screened_Compound_MOA.csv')) %>%
  dplyr::rename(Compound_ID=Compound_ID_1) %>%
  dplyr::select(-MOA)

# Load compound moa table
xmoa <- str_c(ccl.dir, 'data/Screened_Compound_ID_MOA_added_20200420.csv') %>%
  fread() %>%
  dplyr::select(-InchIKeys_All) %>%
  dplyr::rename(Category=MOA) %>%
  dplyr::select(-Compound_Category, -Pathway, -Target) %>%
  dplyr::select(-Accuracy, -nBioactive_Cell) 

# Map MOA compounds to unique ID based on structure
xmoa <- mutate(xmoa, Unique_Compound_ID=Compound_ID)

id2 <- match(xmoa$Compound_ID, xcpd$Compound_ID_2)
xmoa$Unique_Compound_ID[!is.na(id2)] <- xcpd$Compound_ID[na.omit(id2)]

id3 <- match(xmoa$Compound_ID, xcpd$Compound_ID_3)
xmoa$Unique_Compound_ID[!is.na(id3)] <- xcpd$Compound_ID[na.omit(id3)]

# Clean compound MOA for unique compound IDs,drop compounds missing MOA
xmoa.unique <- dplyr::select(xmoa, Unique_Compound_ID, Category) %>%
  group_by(Unique_Compound_ID) %>%
  summarize(Category=list(unique(Category))) %>%
  mutate(Category=lapply(Category, function(z) setdiff(z, ''))) %>%
  mutate(NCat=sapply(Category, length))

#########################**************************#############################
######################### TEMPORARY TAKE FIRST MOA #############################
xmoa.unique <- mutate(xmoa.unique, Category=sapply(Category, function(z) z[1]))
################################################################################
#########################*************************##############################

xmoa <- dplyr::rename(xmoa, CategoryOrig=Category) %>%
  left_join(xmoa.unique, by='Unique_Compound_ID') %>%
  dplyr::select(-NCat, -CategoryOrig) %>%
  filter(!is.na(Category))

write.csv(file=output.file, xmoa, row.names=FALSE)