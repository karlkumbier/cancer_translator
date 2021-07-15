################################################################################
# Modeling functions
################################################################################
fit_cell_line <- function(x, 
                          cell.line, 
                          model, 
                          model_predict, 
                          holdout='random', 
                          prop=0.9,
                          reps=50) {
  
  xc <- filter(x, Cell_Line == cell.line)
  xx <- select(xc, matches('^nonborder'))
  yy <- as.factor(xc$Compound_Category)
  
  # Set training set
  if (holdout == 'random') {
    # Randomly sample training set
    id.train <- createDataPartition(yy, times=reps, p=prop)
    
  } else {
    # Randomly hold-out single compound for testing
    cpd.table <- select(xc, Compound_Category, Compound_ID) %>%
      group_by(Compound_Category) %>%
      summarize(Compound_ID=list(unique(Compound_ID)))
    
    cpd.train <- lapply(1:reps, function(i) {
      sapply(cpd.table$Compound_ID, sample, size=length(cpd.table$Compound_ID) - 1)
    })
    
    id.train <- lapply(cpd.train, function(i) {
      which(xc$Compound_ID %in% i)
    })
  }
  
  # Fit models to each training set
  out <- lapply(id.train, fit_cell_line_, 
                x=xx, 
                y=yy, 
                model=model, 
                model_predict=model_predict
  )
                
  # Aggregate model predictions
  xc.select <- select(xc, !matches('^nonborder')) %>% mutate(ID=1:n())
  out <- rbindlist(out) %>%
    left_join(xc.select, by='ID') %>%
    select(-ID)

  return(out)
}

fit_cell_line_ <- function(x, y, id.train, model, model_predict) {
  
  # Upsample for class balance
  id.test <- setdiff(1:nrow(x), id.train)
  xup <- upSample(x[id.train,], y[id.train])
  ytrain <- xup$Class
  xtrain <- select(xup, -Class)
  
  # Fit model
  fit <- model(xtrain, ytrain)
  
  # Evaluate model class predictions
  yy <- as.numeric(y) - 1
  ypred.class <- model_predict(fit, x[id.test,], 'class')
  ypred <- model_predict(fit, x[id.test,], nclass=length(unique(yy)))
  
  out <- data.table(ID=id.test, YpredCl=ypred.class) %>%
    mutate(Ypred=ypred) %>%
    mutate(Ytrue=yy[id.test])
  return(out)
}

irf <- function(x, y) {
  fit <- iRF(x=x, y=y, type='ranger', n.iter=1, n.core=6)
  return(fit)
}

irf_predict <- function(fit, x, type='response', nclass=NULL) {
  
  if (type == 'class') {
    ypred <- predict(fit$rf.list, x)$predictions
  } else {
    ypred <- predict(fit$rf.list, x, predict.all=TRUE)$predictions
    ypred <- sapply(0:(nclass - 1), function(cl) rowMeans(ypred == cl))
    ypred <- lapply(1:nrow(ypred), function(k) ypred[k,])
  }
  
  return(ypred)
}

################################################################################
# Helper functions
################################################################################
summarize_cell_line <- function(x, cell.line) {
  # Compute average bioactivity score across multiple doses of compound
  x <- filter(x, Cell_Line == cell.line) %>%
    group_by(Compound_ID, Compound_Category, Target) %>%
    summarize(DistNorm=mean(DistNorm))
  return(x)
}

expand_category <- function(x, category) {
  # Expand dataset to include DMSO with each compound category
  xcomp <- filter(x, Compound_Category == category)
  
  # Select dmso cell lines
  plates <- xcomp$PlateID
  xdmso <- filter(x, PlateID %in% plates & Compound_ID == 'DMSO')
  
  out <- rbind(xcomp, xdmso) %>%
    mutate(Category=category) %>%
    mutate(Group=ifelse(Compound_ID == 'DMSO', 'DMSO', 'Compound'))
  return(out)
}

bag_predictions <- function(ypred) {
  # Average compound predictions across cell lines
  ypred.bag <- lapply(unique(ypred$Compound_Dose), function(g) {
    out <- filter(ypred, Compound_Dose == g) %>% 
      select(Ypred, Compound_Category, Ytrue)
    
    ypred <- colMeans(do.call(rbind, out$Ypred))
    out <- data.table(
      Compound_Dose=g, 
      YpredBag=(which.max(ypred) - 1),
      Compound_Category=out$Compound_Category[1],
      Ytrue=out$Ytrue[1]
    )
    
    return(out)
  })
  
  return(rbindlist(ypred.bag))
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
  fit <- lm(y ~ as.numeric(col))
  return(fit$residuals)
}


################################################################################
# Functions for evaluating bioactivity of compounds
################################################################################
bioactivity <- function(x) {
  # Wrapper function to compute distance for each well to DMSO

  # Compute center of DMSO point cloud
  features <- '(^PC|^nonborder)'
  xctl <- filter(x, Compound_ID == 'DMSO') %>% select(matches(features))
  xctl <- colMeans(xctl)
  
  # Compoute distance from DMSO center for each well
  xwell <- select(x, matches(features))
  well.dist <- sqrt(colMeans((t(xwell) - xctl) ^ 2))

  # Generate null distribution of maximum distance DMSO to DMSO center
  null <- generate_null(x, n.core=n.core)
  
  # Normalize relative to control
  xdist <- data.frame(Dist=well.dist) %>%
    mutate(DistNorm=Dist / mean(null)) %>%
    mutate(pval=sapply(Dist, function(d) mean(null > d)))
  
  xmeta <- dplyr::select(x, !matches('(^PC|^nonborder)'))
  return(cbind(xdist, xmeta))
}

null_dist <- function(x) {
  # Compute maximum distance between DMSO wells and point DMSO center
  xctl <- filter(x, Compound_ID == 'DMSO') %>% select(matches('(^PC|^nonborder)'))
  xctl.mean <- colMeans(xctl)
  
  # Compute maximum distance
  xdist <- sqrt(colMeans((t(xctl) - xctl.mean) ^ 2))
  return(max(xdist))
}

generate_null <- function(x, n.samples=100, n.core=1) {
  # Wrapper function for generating null distributions through subsampling
  null.dist <- mclapply(1:n.samples, function(i) {
    set.seed(i)
    return(generate_null_(x))
  }, mc.cores=n.core)
  
  return(unlist(null.dist))
}

generate_null_ <- function(x) {
  # Generate single instance of DMSO self-distance from subsample
  x <- sample_frac(x, 2/3)
  null.dist <- null_dist(x)
  return(null.dist)
}

################################################################################
# Plotting utilities
################################################################################
make_segments <- function(z, x='PC1', y='PC2') {
  # Generate data frame for segment plots
  cell.lines <- unique(z$Cell_Line)
  out <- lapply(cell.lines, function(ii) {
    zi <- filter(z, Cell_Line == ii)
    n <- nrow(zi) - 1
    
    out.i <- data.frame(x=head(zi[[x]], n), 
                        y=head(zi[[y]], n), 
                        xend=tail(zi[[x]], n),
                        yend=tail(zi[[y]], n),
                        Cell_Line=ii
    )
    return(out.i)
  })
  
  return(rbindlist(out))
}

plot_compound <- function(x, compound) {
  # Plot compound in PCs 1, 2
  x.target <- filter(x, Compound_ID == compound) %>%
    group_by(Compound_ID, Dose_Category, Cell_Line) %>%
    mutate_if(is.numeric, mean) %>%
    select(PC1, PC2, Compound_ID, Dose_Category, Cell_Line) %>%
    distinct()
  
  plot.seg <- length(unique(x.target$Dose_Category)) > 1
  if (plot.seg) {
    x.segment <- make_segments(x.target) %>% 
      mutate(Compound_ID=unique(x.target$Compound_ID))
  }
  
  x.dmso <- filter(x, Compound_ID == 'DMSO') %>%
    select(PC1, PC2, Compound_ID, Dose_Category, Cell_Line)
  
  p <- rbind(x.target, x.dmso) %>%
    mutate(Dose_Category=as.numeric(Dose_Category)) %>%
    mutate(Dose_Category=6 - Dose_Category) %>%
    mutate(Dose_Category=ifelse(Compound_ID == 'DMSO', 1, Dose_Category)) %>%
    mutate(Compound_ID=factor(Compound_ID, levels=c('DMSO', compound))) %>%
    ggplot() +
    geom_point(aes(x=PC1, y=PC2, col=Compound_ID, size=Dose_Category)) +
    facet_wrap(~Cell_Line) +
    theme_bw()
  
  if (plot.seg) {
    p <- p + geom_segment(data=x.segment, aes(x=x, y=y, xend=xend, yend=yend, col=Compound_ID))
  }
  
  return(p)
}

