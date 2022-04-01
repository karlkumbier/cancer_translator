################################################################################
# distance functions in phenotypic space
################################################################################
category_distance <- function(x, categories) {
  
  dist.c <- sapply(unique(categories), function(ctg) {
    median(dist(x[categories == ctg,]))
  })
  
  
}

################################################################################
# Feature transformation functions
################################################################################
correlation_weight <- function(x) {
  
  # Apply correlation weighting within each cell line
  x <- lapply(unique(x$Cell_Line), function(cl) {
    
    center <- function(z) z - mean(z)
    
    # Compute feature correlation
    xfeat <- dplyr::select(x, matches('(^non|Cell_Line)'))
    xfeat.c <- filter(xfeat, Cell_Line == cl) %>% dplyr::select(-Cell_Line)
    xnfeat.c <- filter(x, Cell_Line == cl) %>% dplyr::select(-matches('(^non)'))
    
    # Center features
    xfeat.c <- mutate_if(xfeat.c, is.numeric, center)
    
    # Compute feature weight
    xcor <- cor(xfeat.c)
    xcor.wt <- rowSums(1 - abs(xcor))
    xcor.wt <- xcor.wt / max(xcor.wt)
    
    # Weight features
    xfeat.c <- data.frame(t(t(xfeat.c) * xcor.wt))
    
    configs <- umap.defaults
    configs$n_neighbors <- 5
    xfeat.c <- umap(xfeat.c)$layout
    
    return(cbind(xnfeat.c, xfeat.c))
  })

  return(rbindlist(x))
}

################################################################################
# General utility functions
################################################################################
median_impute <- function(x) {
  x[is.na(x)] <- median(x, na.rm=TRUE)
  return(x)
}

select_category <- function(x) {
  # Select majority value from vector of categories
  x <- na.omit(x)
  if (length(x) == 0) return(NA)
  xtab <- table(x)
  return(names(xtab)[which.max(xtab)])
}

count_categories <- function(x) {
  # Count number of unique values in vector of categories
  x <- na.omit(x)
  return(length(unique(x)))
}


################################################################################
# Normalization functions
################################################################################
intensity_normalize <- function(x) {
  # Wrapper function to normalize all ntensity features
  col <- as.numeric(x$Col)
  features <- colnames(x) %>% str_subset('^nonborder')

  for (f in features) {
    x[[f]] <- intensity_normalize_(col, x[[f]])
  }
  
  return(x)
}

intensity_normalize_ <- function(col, y) {
  # Normalize intensity as residuals from lm using column as predictor
  out <- rep(NA, length(y))
  names(out) <- 1:length(y)
  
  fit <- lm(y ~ as.numeric(col))
  out[names(fit$residuals)] <- fit$residuals
  return(out)
}


################################################################################
# Functions for evaluating bioactivity of compounds
################################################################################
bioactivity <- function(x, null_summary=max) {
  # Wrapper function to compute distance for each well to DMSO

  # Feature-center data
  center <- function(z) z - mean(z)
  features <- '(^PC|^non)'
  xwell <- select(x, matches(features))
  
  # Initialize DMSO point cloud center
  xctl <- apply(filter(xwell, x$Compound_ID == 'DMSO'), MAR=2, median)

  # Compute well to DMSO centroid distance
  well.dist <- sqrt(colMeans((t(xwell) - xctl) ^ 2))
  
  # Generate null distribution of maximum distance DMSO to DMSO center
  null <- well.dist[x$Compound_ID == 'DMSO']
  null.thresh <- null_summary(null)
  
  xdist <- data.frame(Dist=well.dist) %>% mutate(DistNorm=Dist / null.thresh)
  xmeta <- dplyr::select(x, !matches('(^PC|^nonborder)'))
  
  return(cbind(xdist, xmeta))
}

null_dist <- function(x) {
  # Wrapper function to compute null distribution of DMSO to centroid distances
  
  # Feature-center data
  center <- function(z) z - mean(z)
  features <- '(^PC|^non)'
  xwell <- select(x, matches(features))
  
  # Initialize DMSO point cloud center
  xctl <- apply(filter(xwell, x$Compound_ID == 'DMSO'), MAR=2, median)

  # Compute mahalanobis distance
  well.dist <- sqrt(colMeans((t(xwell) - xctl) ^ 2))
  return(well.dist[x$Compound_ID == 'DMSO'])
}

opt_bioactive <- function(x, categories, weights, k=1) {
  # Evaluate cell line sets with optimal bioactivity
  
  # Initialize cell line combinations of size k
  cell.combs <- combn(rownames(x), k, simplify=FALSE)
  
  # Compute proportion of bioactivity within each compound category
  unq.category <- unique(categories)
  unq.category <- unq.category[unq.category != '']
  
  bioactive.cat <- lapply(cell.combs, function(cl) {
    sapply(unq.category, function(ct) {
      mean(colSums(x[cl, categories == ct, drop=FALSE]) > 0)
    })
  })
  
  # Filter DMSO from weight vector
  weights <- weights[names(weights) != '']
  
  # Compute weighted bioactivity
  bioactive.score <- sapply(bioactive.cat, function(z) {
    sum(z[names(weights)] * c(weights)) / sum(weights)
  })
  
  names(bioactive.score) <- sapply(cell.combs, str_c, collapse=', ')
  opt <- bioactive.score[which.max(bioactive.score)]
  return(list(opt=opt, score=bioactive.score))
}
