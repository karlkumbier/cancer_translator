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

bootstrap_sample <- function(x) {
  # Wrapper for bootstrap sampling maintaining class balance
  xg <- data.frame(Cat=x) %>%
    mutate(Idx=1:n()) %>%
    group_by(Cat) %>%
    summarize(Idx=list(Idx))
  
  # Note: no longer strict bootstrap due to distances for replicates
  out <- lapply(xg$Idx, function(i) sample(i, length(i) * (2 / 3)))
  out <- unique(unlist(out))
  return(out)
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
dmso_distance <- function(x, null_summary=max) {
  # Computes distance for each well to DMSO centroid

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

bioactive_difference <- function(x, cell.line, bootstrap=FALSE) {
  # Compute deviation between DMSO, compound ECDFs for select cell line(s)
  xcl <- filter(x, Cell_Line %in% cell.line)
  if (bootstrap) xcl <- bootstrap_table(xcl)
  
  xdmso <- filter(xcl, Compound_ID == 'DMSO')
  xcpd <- filter(xcl, Compound_ID != 'DMSO') %>%
    group_by(Compound_ID) %>%
    top_n(1, DistNorm)
  
  return(cvm_stat(xdmso$DistNorm, xcpd$DistNorm, power=1))
}

bioactive_difference_ctg <- function(x, cell.lines, categories, bootstrap=FALSE) {
  # Compute deviation between DMSO, category ECDFs for select cell line(s)
  out <- sapply(categories, function(ctg) {
    sapply(cell.lines, function(cl) {
      xctg <- filter(x, Compound_Category %in% c(ctg, 'DMSO'))
      bioactive_difference(xctg, cl, bootstrap)
    }) 
  })
  
  return(out)
}

score_bioactive <- function(x, weights) {
  # Evaluate weighted bioactivity score by cell line set
  x <- x[,names(weights)]
  stopifnot(ncol(x) == length(weights))
  score <- colSums(t(x) * weights) / sum(weights)
  return(return(score))
}

bootstrap_table <- function(x) {
  # Bootstrap sample table maintaining DMSO/cpd balance
  #xdmso <- filter(x, Compound_ID == 'DMSO')
  #xdmso <- slice_sample(xdmso, n=nrow(xdmso), replace=TRUE)
  #xcpd <- filter(x, Compound_ID != 'DMSO')
  #xcpd <- slice_sample(xcpd, n=nrow(xcpd), replace=TRUE)
  
  # Add epsilon noise to deal with duplicate samples
  eps <- rnorm(nrow(x), sd=1e-10)
  id.bootstrap <- bootstrap_sample(x$Category)
  out <- x[id.bootstrap,]# %>% mutate(DistNorm=DistNorm + eps)
    #slice_sample(x, n=nrow(x), replace=TRUE) %>% 
    #mutate(DistNorm=DistNorm + eps)
  
  return(out)
}

################################################################################
# Functions for evaluating compound similarity
################################################################################
category_distance <- function(x, categories, summary_fun, dist.mat=NULL, bootstrap=FALSE) {
  
  # Compute pairwise distances
  if (is.null(dist.mat)) dist.mat <- dist(x) %>% as.matrix
  
  # Take bootstrap sample
  if (bootstrap) id.bootstrap <- bootstrap_sample(categories)
  
  # Initialize DMSO median distance
  id.dmso <- which(categories == 'DMSO')
  if (bootstrap) id.dmso <- id.bootstrap[id.bootstrap %in% id.dmso]
  dmso.dist <- dist.mat[id.dmso, id.dmso]
  dmso.dist <- dmso.dist[upper.tri(dmso.dist)]
  
  # Add epsilon noise if bootstrapping
  #if (bootstrap) {
  #  eps <- rnorm(length(dmso.dist), sd=1e-10)
  #  dmso.dist <- dmso.dist + eps
  #}

  # Compute within category distances
  dist.summary <- sapply(unique(categories), function(ctg) {
    # Compute within category distance distribution
    id.ctg <- which(categories == ctg)
    if (bootstrap) id.ctg <- id.bootstrap[id.bootstrap %in% id.ctg]
    ctg.dist <- dist.mat[id.ctg, id.ctg]
    ctg.dist <- ctg.dist[upper.tri(ctg.dist)]
    
    # Add epsilon noise if bootstrapping
    #if (bootstrap) {
    #  eps <- rnorm(length(ctg.dist), sd=1e-10)
    #  ctg.dist <- ctg.dist + eps
    #}
    
    # Compute category neighborhood distance distribution
    if (bootstrap) dist.mat <- dist.mat[,unique(id.bootstrap)]
    nbhd.dist <- apply(dist.mat[id.ctg, ], MAR=1, function(z) sort(z))
    nbhd.dist <- nbhd.dist[2:length(unique(id.ctg)),]
    
    out.within <- summary_fun(ctg.dist, dmso.dist)
    out.between <- summary_fun(ctg.dist, nbhd.dist)
    return(c(out.within, out.between))
  })
  
  out <- data.frame(t(dist.summary)) %>%
    mutate(Category=colnames(dist.summary)) %>%
    dplyr::rename(DistWithin=X1, DistBetween=X2)
  
  return(out)
}

aggregate_similarity <- function(x, cell.lines) {
  # Aggregate as max similarity across cell lines
  out <- filter(x, Cell_Line %in% cell.lines) %>%
    group_by(Category) %>%
    summarize(DistWithin=max(DistWithin), DistBetween=max(DistBetween)) %>%
    mutate(Cell_Set=str_c(cell.lines, collapse=', '))
  
  return(out)
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
