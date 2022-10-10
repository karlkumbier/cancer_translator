#' This script maps cell-level features extracted raw imaging data to 
#' well-level KS profiles. To run this script, you must first install the R 
#' package `rprofiler`, which contains utility functions for processing imaging 
#' data. A version of the `rprofiler` package is contained with this repository.
#' 
#' Last updated 10/5/2022, Karl Kumbier
library(data.table)
library(tidyverse)
library(rprofiler)
library(parallel)

ccl.dir <- Sys.getenv('CCL_DIR')
plate.dir <- str_c(ccl.dir, '/data/screens/LH_CDC_1/raw_data/')
platemap.dir <- str_c(ccl.dir, '/data/screens/LH_CDC_1/platemaps/')

plate.dirs <- list.dirs(plate.dir, recursive=TRUE)
platemaps <- str_c(platemap.dir, list.files(platemap.dir))

nfeat <- 93
n.core <- 24

read_plate <- function(plate.dir, nfeat) {
  
  # Set evaluation directory to be read
  eval.dir <- list.dirs(plate.dir) %>% str_subset('Evaluation')
  if (length(list.dirs(plate.dir, recursive=FALSE)) > 1) {
    eval.dir <- check_eval_dirs(plate.dir, nfeat)
  }
  
  x <- fread(str_c(eval.dir, '/Objects_Population - nonborder.txt'), skip=9)
  x <- mutate(x, PlateID=sapply(str_split(plate.dir, '/'), tail, 2)[1])
              
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

generate_ks_profiles <- function(x, id.control, prop=1, n.core=1) {
  # Generate KS profiles for select plate
  xcontrol <- filter(x, id.control) %>%  sample_frac(prop)
  
  # Note: features are defined by nonborder key
  id.feat <- which(str_detect(colnames(x), '^nonborder'))
  features <- colnames(x)[id.feat]
  plate <- unique(x$PlateID)
  
  # Set well ids
  wells <- unique(x$WellID)

  # Generate KS profiles for each well
  x <- lapply(wells, function(w) filter(x, WellID == w))
  xks <- mclapply(x, wellKS, xcontrol=xcontrol, id.feat=id.feat, mc.cores=n.core)
  xks <- do.call(rbind, xks)
  
  colnames(xks)[1:(ncol(xks) - 1)] <- features
  xks <- data.frame(xks) %>% mutate(PlateID=plate, WellID=wells)
  return(xks)
}


################################################################################
# Load in metadata for each plate
################################################################################
xmeta <- lapply(platemaps, function(p) {
  plate.id <- str_remove_all(p, "^.*platemaps/") %>% str_remove_all('\\.xlsx')
  out <- loadMeta(p) %>% mutate(PlateID=plate.id)
})

# Subset to metadata features in all files
features <- lapply(xmeta, colnames)
features <- Reduce(intersect, features)
xmeta <- lapply(xmeta, function(x) dplyr::select(x, one_of(features)))

# Filter out empty wells
xmeta <- rbindlist(xmeta) %>% filter(!is.na(Compound_ID))

################################################################################
# Load in data for each plate
################################################################################
plates <- unique(xmeta$PlateID)

x <- mclapply(plates, function(p) {
  set.seed(47)
  tryCatch({
    xmeta.p <- filter(xmeta, PlateID == p)
    
    # Match metadata/data
    plate.dir.p <- str_subset(plate.dirs, p) %>% str_subset('Evaluation')
    plate.dir.p <- tail(plate.dir.p, 1)

    out <- read_plate(plate.dir.p, nfeat=nfeat) %>% mutate(WellID=str_c(Row, '-', Column))
    xmeta.p <- xmeta.p[match(out$WellID, xmeta.p$WellID),]
    
    out <- dplyr::select(out, matches('(PlateID|WellID|^nonborder.*)'))
    id.control <- xmeta.p$Compound_ID == 'DMSO'
    xks <- generate_ks_profiles(out, id.control=id.control, n.core=n.core)
    
    # Match metadata with KS profiles
    xmeta.p <- distinct(xmeta.p)
    xmeta.p <- xmeta.p[match(xks$WellID, xmeta.p$WellID),]
    return(list(xks=xks, xmeta=xmeta.p))
  }, error=function(e) {
    print(p)
    return(NULL)
  })
}, mc.cores=1)


xks <- rbindlist(lapply(x, function(z) return(z$xks)))
xmeta <- rbindlist(lapply(x, function(z) return(z$xmeta)))
save(file='./ks_profiles.Rdata', xks, xmeta)
