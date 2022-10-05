#' This script merges KS profile data with annotated compound library. Compound
#' IDs from the raw data (e.g. from vendor) are replaced with compound IDs that
#' have been grouped by structure. 
#' 
#' - KS profiles for compounds with no MOA annotation are dropped from analysis.
#' - KS profiles are filtered to samples from highest dose of compound.
#' - Median KS profiles are taken in case of compound replicates.
#' 
#' Last updated 10/4/2022, Karl Kumbier
library(data.table)
library(tidyverse)
library(testthat)

################################################################################
# Load phenotypic profiling dataset
################################################################################
ccl.dir <- Sys.getenv('CCL_DIR')
compound.file <- str_c(ccl.dir, 'data/Table1_compound_library.csv')
ks.file <- str_c(ccl.dir, 'data/screens/LH_CDC_1/profiles_normalized.Rdata')

source(str_c(ccl.dir, 'scripts/utilities.R'))

# Load compound library and KS datasets
xmoa <- fread(compound.file)
load(ks.file)

# Clean fibroblast line name
x <- dplyr::select(x, -matches('^PC')) %>%
  mutate(Cell_Line=ifelse(Cell_Line == 'ALS-WT', 'FB', Cell_Line))

# Merge MOA labels with phenotypic profiles, replace ID with unique ID
x <- left_join(x, xmoa, by='Compound_ID') %>%
  mutate(Unique_Compound_ID=ifelse(Compound_ID == 'DMSO', 'DMSO', Unique_Compound_ID)) %>%
  mutate(Category=ifelse(Compound_ID == 'DMSO', 'DMSO', Category)) %>%
  mutate(Compound_ID=Unique_Compound_ID)

# Filter KS profiles from unannotated compounds
x <- filter(x, !is.na(Category))

# Filter to highest dose
x <- group_by(x, Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup()

################################################################################
# Merge phenotypic profiles for replicates - average over wells
################################################################################
xmeta <- dplyr::select(x, Cell_Line, Compound_ID, Category, Dose)
x <- dplyr::select(x, matches('^nonborder'))
x <- cbind(xmeta, x)

# Initialize DMSO set
id.dmso <- xmeta$Category == 'DMSO'
xdmso <- x[id.dmso,]

# Merge replicates in non-dmso set
xtreat <- group_by(x[!id.dmso,], Cell_Line, Compound_ID, Category, Dose) %>%
  summarize_if(is.numeric, median) %>%
  ungroup()

# Set positive/negative control labels
x <- rbind(xdmso, xtreat) %>%
  mutate(Usage=ifelse(Compound_ID == 'DMSO', 'negative_ctrl_cpd', 'query_cpd')) %>%
  mutate(Usage=ifelse(Compound_ID == 'Gemcitabine', 'positive_ctrl_cpd', Usage)) %>%
  mutate(Usage=ifelse(Compound_ID == 'Bortezomib', 'positive_ctrl_cpd', Usage))

################################################################################
# Tests for final output dataset
################################################################################
test_that("All compounds appear once for each cell line", {
  cpd.table <- table(x$Compound_ID)
  cpd.table <- cpd.table[names(cpd.table) != 'DMSO']
  ncells <- length(unique(x$Cell_Line))
  expect_true(all(cpd.table == ncells))
})

test_that("Compounds map to unique Category", {
  xcat.tab <- group_by(x, Compound_ID) %>%
    summarize(Category=list(unique(Category))) %>%
    mutate(NCat=sapply(Category, length))
  expect_true(all(xcat.tab$NCat == 1))
})

test_that("Single dose per compound", {
  xdose.tab <- group_by(x, Compound_ID) %>%
    summarize(Dose=list(unique(Dose))) %>%
    mutate(NDose=sapply(Dose, length))
  expect_true(all(xdose.tab$NDose == 1))
})



