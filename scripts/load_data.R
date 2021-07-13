library(data.table)
library(tidyverse)
library(rprofiler)
library(parallel)

setwd('/awlab/Lab_temp/Louise/opfeatures/')
plate.dir <- 'LH_CDC_1/plates/'
platemap.dir <- 'LH_CDC_1/platemaps/'

plates <- list.dirs(plate.dir, recursive=FALSE) %>% 
  str_subset('reimaged', negate=TRUE)

platemaps <- str_c(platemap.dir, list.files(platemap.dir))

nfeat <- 93
n.core <- 16
prop.dmso <- 0.25


read_plate <- function(plate.dir, nfeat) {
  
  # Set evaluation directory to be read
  eval.dir <- list.dirs(plate.dir) %>% str_subset('Evaluation')
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

x <- mclapply(plates, function(p) {
  set.seed(47)
  tryCatch({
    plate.id <- str_remove_all(p, "^.*plates//")
    xmeta.p <- filter(xmeta, PlateID == plate.id)
    
    # Match metadata/data
    out <- read_plate(p, nfeat=nfeat) %>% mutate(WellID=str_c(Row, '-', Column))
    xmeta.p <- xmeta.p[match(out$WellID, xmeta.p$WellID),]
    
    out <- dplyr::select(out, matches('(PlateID|WellID|^nonborder.*)'))
    id.control <- xmeta.p$Compound_ID == 'DMSO'
    xks <- generate_ks_profiles(out, id.control=id.control, n.core=16, prop=prop.dmso)
    
    # Match metadata with KS profiles
    xmeta.p <- distinct(xmeta.p)
    xmeta.p <- xmeta.p[match(xks$WellID, xmeta.p$WellID),]
    return(list(xks=xks, xmeta=xmeta))
  }, error=function(e) {
    print(p)
    return(NULL)
  })
}, mc.cores=1)


xks <- rbindlist(lapply(x, function(z) return(z$xks)))
xmeta <- rbindlist(lapply(x, function(z) return(z$xmeta)))
save(file='./ks_profiles.Rdata', xks, xmeta)
