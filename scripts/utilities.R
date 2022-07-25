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
normalize_features <- function(x) {
  # Wrapper function to normalize features
  features <- colnames(x) %>% str_subset('^nonborder')
  xpredict <- dplyr::select(x, ColID, RowID, Compound_ID)
  
  for (f in features) {
    xpredict.f <- mutate(xpredict, Y=x[[f]])
    x[[f]] <- normalize_features_(xpredict.f)
  }
  
  return(x)
}

normalize_features_ <- function(x) {
  # Normalize intensity as residuals from lm using column as predictor
  
  # Initialize DMSO data set for model fitting
  x.fit <- filter(x, Compound_ID == 'DMSO') %>% dplyr::select(-Compound_ID)

  # Initialize outlier samples to drop from model fitting
  outlier.thresh <- 1.5 * IQR(x.fit$Y)
  outlier.med <- median(x.fit$Y)
  id.low <- x.fit$Y < outlier.med - outlier.thresh
  id.high <- x.fit$Y > outlier.med + outlier.thresh
  id.drop <- id.low | id.high
  
  # Fit model and compute residuals
  if (length(unique(x$ColID)) > 2) {
    fit <- lm(Y ~ ColID + RowID, data=x.fit[!id.drop,])
    x.full <- dplyr::select(x, -Compound_ID)
    ypred <- predict(fit, x.full)
    out <- x.full$Y - ypred
  } else {
    fit <- lm(Y ~ RowID, data=x.fit[!id.drop,])
    x.full <- dplyr::select(x, -Compound_ID)
    ypred <- predict(fit, x.full)
    out <- x.full$Y - ypred  
  }
  
  # DMSO-center feature
  centroid <- median(out[x$Compound_ID == 'DMSO'])
  return(out - centroid)
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

phenoactivity <- function(x, cell.line, categories, bootstrap=FALSE) {
  # Compute deviation between DMSO, compound ECDFs for select cell line(s)
  xcl <- filter(x, Cell_Line %in% cell.line)
  if (bootstrap) xcl <- bootstrap_table(xcl)
  
  out <- sapply(unique(categories), function(ctg) {
    xdmso <- filter(xcl, Category == 'DMSO')
    xcpd <- filter(xcl, Category == ctg)
      
    return(cvm_stat(xdmso$DistNorm, xcpd$DistNorm, power=1))
  })
  
  return(out)
}

score_bioactive <- function(x, weights) {
  # Evaluate weighted bioactivity score by cell line set
  x <- x[names(weights)]
  stopifnot(length(x) == length(weights))
  score <- sum(x * weights) / sum(weights)
  return(score)
}

bootstrap_table <- function(x) {
  # Bootstrap sample table maintaining DMSO/cpd balance

  # Add epsilon noise to deal with duplicate samples
  id.bootstrap <- bootstrap_sample(x$Category)
  out <- x[id.bootstrap,]
  
  return(out)
}

################################################################################
# Functions for evaluating compound similarity
################################################################################
phenosimilarity <- function(x, categories, summary_fun, dist.mat=NULL, bootstrap=FALSE) {
  
  # Compute pairwise distances
  if (is.null(dist.mat)) dist.mat <- dist(x) %>% as.matrix
  
  # Take bootstrap sample
  if (bootstrap) id.bootstrap <- bootstrap_sample(categories)
  
  # Compute within category distances
  dist.summary <- sapply(unique(categories), function(ctg) {
    # Compute within category distance distribution
    id.ctg <- which(categories == ctg)
    if (bootstrap) id.ctg <- id.bootstrap[id.bootstrap %in% id.ctg]
    ctg.dist <- dist.mat[id.ctg, id.ctg]
    ctg.dist <- ctg.dist[upper.tri(ctg.dist)]
    
    # Compute category neighborhood distance distribution
    if (bootstrap) dist.mat <- dist.mat[,unique(id.bootstrap)]
    nbhd.dist <- apply(dist.mat[id.ctg, ], MAR=1, function(z) sort(z))
    nbhd.dist <- nbhd.dist[2:length(unique(id.ctg)),]
    
    # Compute proportion of nearest neighbors in category
    if (!bootstrap) {
      nbhd.ovlp <- apply(dist.mat[id.ctg, ], MAR=1, function(z) order(z))
      nbhd.ovlp <- nbhd.ovlp[2:length(unique(id.ctg)),]
      nbhd.prop <- apply(nbhd.ovlp, MAR=2, function(z) length(intersect(z, id.ctg)))
      nbhd.prop <- mean(nbhd.prop / (length(id.ctg) - 1))
    } else {
      nbhd.prop <- 0
    }
    
    
    out.between <- summary_fun(ctg.dist, nbhd.dist)
    return(c(out.between, nbhd.prop))
  })
  
  out <- data.frame(Dist=dist.summary[1,]) %>%
    mutate(Prop=dist.summary[2,]) %>%
    mutate(Category=colnames(dist.summary)) 
    
  return(out)
}

phenosimilarity_distn <- function(x, categories, dist.mat=NULL) {
  
  # Compute pairwise distances
  if (is.null(dist.mat)) dist.mat <- dist(x) %>% as.matrix
  
  # Compute within category distances
  distn <- lapply(unique(categories), function(ctg) {
    # Compute within category distance distribution
    id.ctg <- which(categories == ctg)
    ctg.dist <- dist.mat[id.ctg, id.ctg]
    ctg.dist <- ctg.dist[upper.tri(ctg.dist)]
    
    # Compute category neighborhood distance distribution
    nbhd.dist <- apply(dist.mat[id.ctg, ], MAR=1, function(z) sort(z))
    nbhd.dist <- c(nbhd.dist[2:length(unique(id.ctg)),])
    
    out <- list(ctg.dist, nbhd.dist)
    names(out) <- c('Query', 'Reference')
    return(out)
  })
  
  names(distn) <- unique(categories)
  return(distn)
}


aggregate_similarity <- function(x, cell.lines) {
  # Aggregate as max similarity across cell lines
  out <- filter(x, Cell_Line %in% cell.lines) %>%
    group_by(Category) %>%
    summarize(Dist=max(Dist)) %>%
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

phenosimilarity_cpd <- function(x, categories, compounds, dist.mat=NULL) {
  
  # Compute pairwise distances
  if (is.null(dist.mat)) dist.mat <- dist(x) %>% as.matrix
  
  # Compute within category distances
  dist.summary <- lapply(unique(categories), function(ctg) {
    # Compute within category distance distribution
    id.ctg <- which(categories == ctg)
    cpd.ctg <- compounds[id.ctg]
    
    # Compute average distance between compounds within category
    ctg.dist <- dist.mat[id.ctg, id.ctg]
    diag(ctg.dist) <- NA
    avg.dist <- rowMeans(ctg.dist, na.rm=TRUE)
    
    # Compute average distance to nearest neighborhood
    nbhd.dist <- apply(dist.mat[id.ctg, ], MAR=1, function(z) sort(z))
    nbhd.dist <- nbhd.dist[2:length(unique(id.ctg)),]
    avg.nbhd.dist <- colMeans(nbhd.dist)
    
    out <- data.frame(Compound_ID=cpd.ctg, Category=ctg, Dist=avg.nbhd.dist/avg.dist)
    return(out)
  })
  
    
  return(rbindlist(dist.summary))
}
