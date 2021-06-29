library(data.table)
library(tidyverse)
library(rprofiler)

setwd('~/github/cancer_translator/')
data.dir <- 'data/screens/20210616/'

nfeat <- 93
plates <- list.dirs(data.dir, recursive=FALSE)
n.core <- 16

read_plate <- function(plate.dir, nfeat) {
  
  # Set evaluation directory to be read
  eval.dir <- str_c(plate.dir, '/Evaluation1')
  if (length(list.dirs(plate.dir, recursive=FALSE)) > 1) {
    eval.dir <- check_eval_dirs(plate.dir, nfeat)
  }
  
  x <- fread(str_c(eval.dir, '/Objects_Population - nonborder.txt'), skip=9)
  x <- mutate(x, PlateID=sapply(str_split(plate.dir, '/'), tail, 1))
              
  # Check that output will contain specified number of features
  if (ncol(x) != (nfeat + 1)) {
    warning(str_c('Feature size mismatch: ', plate.dir))
    return(NULL)
  }
  
  return(x)
}

check_eval_dirs <- function(plate.dir, nfeat) {
  # Check whether feature files in eval dirs contain correct number of features
  eval.dirs <- list.dirs(plate.dir) %>% str_subset('Evaluation')
  
  nfeat.eval <- sapply(eval.dirs, function(d) {
    x <- fread(str_c(d, '/Objects_Population - nonborder.txt'), skip=9, nrows=1)
    return(ncol(x))
  })
  
  return(eval.dirs[nfeat.eval == nfeat][1])
}

select_features <- function(x, fselect) {
  # Select features from data frame matching regex
  if (is.null(x)) return(NULL)
  x <- dplyr::select(x, matches(fselect))
}

generate_ks_profiles <- function(x, prop=1, n.core=1) {
  # Generate KS profiles for select plate
  # NOTE: controls are defined by colums (2, 23)
  xcontrol <- filter(x, Column %in% c(2, 23)) %>% sample_frac(prop)
  
  # Note: features are defined by nonborder key
  id.feat <- which(str_detect(colnames(x), '^nonborder'))
  features <- colnames(x)[id.feat]
  plate <- unique(x$PlateID)
  
  # Set well ids
  x <- mutate(x, Well=str_c(Row, '_', Column))
  wells <- unique(x$Well)

  # Generate KS profiles for each well
  x <- lapply(wells, function(w) filter(x, Well == w))
  xks <- mclapply(x, wellKS, xcontrol=xcontrol, id.feat=id.feat, mc.cores=n.core)
  xks <- do.call(rbind, xks)
  
  colnames(xks)[1:(ncol(xks) - 1)] <- features
  xks <- data.frame(xks) %>% mutate(PlateID=plate, WellID=wells)
  return(xks)
}

################################################################################
# Load in data for each plate
################################################################################
x <- mclapply(plates, read_plate, nfeat=nfeat, mc.cores=n.core)

# Subset to features and well metadata
fselect <- '(PlateID|Row|Column|^nonborder.*)'
x <- lapply(x, select_features, fselect=fselect)

# Generate KS profiles for each plate
xks <- lapply(x, generate_ks_profiles, n.core=16, prop=0.1)
xks <- rbindlist(xks)
save(file='data/screens/20210616/ks_profiles.Rdata', xks)
