################################################################################
# distance functions in phenotypic space
################################################################################
category_distance <- function(x, categories, summary_fun) {
  
  # Compute pairwise distances
  dist.mat <- dist(x) %>% as.matrix
  
  # Initialize DMSO median distance
  id.dmso <- categories == 'DMSO'
  dmso.dist <- dist.mat[id.dmso, id.dmso]
  dmso.dist <- dmso.dist[upper.tri(dmso.dist)]
  
  # Compute within category distances
  dist.summary <- sapply(unique(categories), function(ctg) {
    ctg.id <- categories == ctg
    
    # Compute within category distance distribution
    ctg.dist <- dist.mat[ctg.id, ctg.id]
    ctg.dist <- ctg.dist[upper.tri(ctg.dist)]
    
    # Compute category neighborhood distance distribution
    nbhd.dist <- apply(dist.mat[ctg.id, ], MAR=1, function(z) sort(z))
    nbhd.dist <- nbhd.dist[2:sum(ctg.id),]
    
    out.within <- summary_fun(ctg.dist, dmso.dist)
    out.between <- summary_fun(ctg.dist, nbhd.dist)
    return(c(out.within, out.between))
  })
  
  out <- data.frame(t(dist.summary)) %>%
    mutate(Category=colnames(dist.summary)) %>%
    dplyr::rename(DistWithin=X1, DistBetween=X2)
  
  return(out)
}

if (FALSE) {
ecdf <- function(x) {
  n <- length(x)
  return(cbind(sort(x), (1:n) / n))
}

dmso.ecdf <- data.frame(ecdf(dmso.dist), Type='DMSO')
ctg.ecdf <- data.frame(ecdf(ctg.dist), Type='Ctg')
nbhd.ecdf <- data.frame(ecdf(nbhd.dist), Type='Nbhd')

rbind(dmso.ecdf, ctg.ecdf, nbhd.ecdf) %>%
  ggplot(aes(x=X1, y=X2, col=Type)) +
  geom_line()

cvm_stat(ctg.dist, dmso.dist, power=1)
cvm_stat(ctg.dist, nbhd.dist, power=1)
}
################################################################################
# Feature transformation functions
################################################################################
correlation_weight <- function(x, normalize=FALSE) {
  
  # Apply correlation weighting within each cell line
  x <- lapply(unique(x$Cell_Line), function(cl) {
    
    center <- function(z) z - mean(z)
    
    # Compute feature correlation
    xfeat <- dplyr::select(x, matches('(^non|Cell_Line)'))
    xfeat.c <- filter(xfeat, Cell_Line == cl) %>% dplyr::select(-Cell_Line)
    xnfeat.c <- filter(x, Cell_Line == cl) %>% dplyr::select(-matches('(^non)'))
    
    # l2 normalize features
    l2norm <- function(z) z / sqrt(sum(z ^ 2))
    if (normalize)  xfeat.c <- mutate_if(xfeat.c, is.numeric, center)
    
    # Center features
    xfeat.c <- mutate_if(xfeat.c, is.numeric, center)
    
    # Compute feature weight
    xcor <- cor(xfeat.c)
    xcor.wt <- rowSums(1 - abs(xcor))
    xcor.wt <- xcor.wt / max(xcor.wt)
    
    # Weight features
    xfeat.c <- data.frame(t(t(xfeat.c) * xcor.wt))
    
    #configs <- umap.defaults
    #configs$n_neighbors <- 5
    #xfeat.c <- umap(xfeat.c)$layout
    
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

opt_similarity <- function(x, categories, weights, k=1) {
  # Evaluate cell line sets with optimal similarity
  
  # Initialize cell line combinations of size k
  cell.combs <- combn(rownames(x), k, simplify=FALSE)
  
  # Compute similarity within each compound category
  unq.category <- unique(categories)
  unq.category <- unq.category[unq.category != '']
  
  similarity.cat <- lapply(cell.combs, function(cl) {
    cat.sim <- apply(xplot[cl,,drop=FALSE], MAR=2, max)
  })
  
  # Filter DMSO from weight vector
  weights <- weights[names(weights) != '']
  
  # Compute weighted similarity
  similarity.score <- sapply(similarity.cat, function(z) {
    sum(z[names(weights)] * c(weights)) / sum(weights)
  })
  
  names(similarity.score) <- sapply(cell.combs, str_c, collapse=', ')
  opt <- similarity.score[which.max(similarity.score)]
  return(list(opt=opt, score=similarity.score))
}
