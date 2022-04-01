median_impute <- function(x) {
  x[is.na(x)] <- median(x, na.rm=TRUE)
  return(x)
}

################################################################################
# Modeling functions
################################################################################
fit_marker <- function(x, 
                       markers, 
                       model, 
                       model_predict, 
                       marker.features,
                       cpd.train=NULL,
                       trt.train=NULL,
                       holdout='random', 
                       prop=0.9,
                       reps=50) {
  
  features <- unlist(marker.features[markers])
  markers.c <- str_c(features, collapse='|')
  marker.re <- str_c('(', markers.c, ')')
  
  xx <- select(x, matches('^nonborder')) %>% select(matches(marker.re))
  yy <- as.factor(x$Compound_Category)
  
  # Set training set
  if (holdout == 'random') {
    if (is.null(trt.train)) {
      # Randomly hold-out treatments for testing
      cpd.table <- select(x, Compound_Category, Treatment) %>%
        group_by(Compound_Category) %>%
        summarize(Treatment=list(unique(Treatment)))
      
      trt.train <- lapply(1:reps, function(i) {
        out <- sapply(cpd.table$Treatment, function(z) {
          sample(z, length(z) * prop)
        })
        
        return(unlist(out))
      })
    }
    
    # Intiailize training set indices
    id.train <- lapply(trt.train, function(i) {
      which(x$Treatment %in% i)
    })
    
  } else {
    if (is.null(cpd.train)) {
      # Randomly hold-out compounds for testing
      cpd.table <- select(x, Compound_Category, Compound_ID) %>%
        group_by(Compound_Category) %>%
        summarize(Compound_ID=list(unique(Compound_ID)))
      
      cpd.train <- lapply(1:reps, function(i) {
        sapply(cpd.table$Compound_ID, sample, size=length(cpd.table$Compound_ID) - 1)
      })
    }
    
    id.train <- lapply(cpd.train, function(i) {
      which(x$Compound_ID %in% i)
    })
  }
  
  # Fit models to each training set
  out <- lapply(1:length(id.train), function(i) {
    out <- fit_wrap(xx, yy, id.train[[i]], model, model_predict)
    out <- mutate(out, Rep=i)
    return(out)
  })
  
  # Aggregate model predictions
  xx.select <- select(x, !matches('^nonborder')) %>% 
    mutate(ID=1:n())
  
  out <- rbindlist(out) %>%
    left_join(xx.select, by='ID') %>%
    select(-ID) %>%
    mutate(Markerset=str_c(markers, collapse='_'))
  
  return(out)
}

fit_cell_line <- function(x, 
                          cell.line, 
                          model, 
                          model_predict, 
                          trt.train=NULL,
                          cpd.train=NULL,
                          holdout='random', 
                          prop=0.9,
                          reps=50,
                          n.core=1,
                          full.samples=6) {
  
  if (cell.line != 'Full') {
    xc <- filter(x, Cell_Line == cell.line) %>%
      select(matches('(^nonborder|Compound_ID|Cell_Line|Compound_Category|Dose_Category|Treatment)'))
    
    xx <- select(xc, matches('^nonborder'))
    yy <- as.factor(xc$Compound_Category)
  } else {
    xc <- x %>%
      select(matches('(^nonborder|Compound_ID|Cell_Line|Compound_Category|Dose_Category|Treatment)')) 
    
    xdata <- expand_cell_data(xc, full.samples=full.samples, n.core=n.core)
    xx <- xdata$x
    yy <- xdata$y
    
    xc <- data.frame(
      Compound_Category=yy, 
      Cell_Line='Full', 
      Dose_Category=xdata$dose,
      Compound_ID=xdata$compound
    )
  }
  
  # Set training set
  if (holdout == 'random') {
    # Randomly hold out treatment compound for testing
    if (is.null(trt.train)) {
      cpd.table <- select(xc, Compound_Category, Treatment) %>%
        group_by(Compound_Category) %>%
        summarize(Treatment=list(unique(Treatment)))
      
      trt.train <- lapply(1:reps, function(i) {
        out <- sapply(cpd.table$Treatment, function(z) {
          sample(z, length(z) * prop)
        })
        
        return(unlist(out))
      })
    }
    
    # Initialize indices of training samples
    id.train <- lapply(trt.train, function(i) {
      which(xc$Treatment %in% i)
    })
    
  } else {
    # Randomly hold out compound for testing
    if (is.null(cpd.train)) {
      cpd.table <- select(xc, Compound_Category, Compound_ID) %>%
        group_by(Compound_Category) %>%
        summarize(Compound_ID=list(unique(Compound_ID)))
      
      cpd.train <- lapply(1:reps, function(i) {
        sapply(cpd.table$Compound_ID, function(z) {
          sample(z, length(z) - 1)
        })
      })
    }
    
    # Initialize indices of training samples
    id.train <- lapply(cpd.train, function(i) {
      which(xc$Compound_ID %in% i)
    })
  }
  
  # Fit models to each training set
  out <- lapply(1:length(id.train), function(i) {
    out <- fit_wrap(xx, yy, id.train[[i]], model, model_predict)
    out <- mutate(out, Rep=i)
    return(out)
  })
                
  # Aggregate model predictions
  xc.select <- select(xc, !matches('^nonborder')) %>% mutate(ID=1:n())
  out <- rbindlist(out) %>%
    left_join(xc.select, by='ID') %>%
    select(-ID)

  return(out)
}

expand_cell_data <- function(x, full.samples=6, n.core=1) {
  # Generate feature set from multiple cell lines for each compound/dose pair
  groups <- str_c(x$Compound_ID, ', ', x$Dose_Category)
  
  # Generate feature matrix
  groups.gen <- rep(unique(groups), full.samples)
  out.x <- mclapply(groups.gen, function(g) {
    xx <- filter(x, groups == g) %>% 
      group_by(Cell_Line) %>%
      sample_n(1) %>%
      ungroup() %>%
      select(matches('^nonborder'))
    return(unlist(c(xx)))
  }, mc.cores=n.core)
  
  out.x <- do.call(rbind, out.x)
  
  out.y <- mclapply(groups.gen, function(g) {
    yy <- unique(filter(x, groups == g)$Compound_Category)
  }, mc.cores=n.core)
  
  out.y <- as.factor(unlist(out.y))
  dose <- str_remove_all(groups.gen, '^.*, ')
  compound <- str_remove_all(groups.gen, ',.*$')
  return(list(x=out.x, y=out.y, dose=dose, compound=compound))
}

fit_wrap <- function(x, y, id.train, model, model_predict) {
  
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
  fit <- iRF(x=x, y=y, type='ranger', n.iter=1, n.core=1, verbose=FALSE)
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

bag_predictions <- function(ypred, cell.lines=NULL) {
  
  # Filter to selected cell lines
  if (!is.null(cell.lines)) ypred <- filter(ypred, Cell_Line %in% cell.lines)
  
  # Average compound predictions across cell lines
  ypred.bag <- lapply(unique(ypred$Treatment), function(g) {
    out <- filter(ypred, Treatment == g) %>% 
      select(Ypred, Compound_Category, Ytrue, Cell_Line, Rep)
    
    out <- lapply(unique(out$Rep), function(rep) {
      
      out.r <- filter(out, Rep == rep)
      ypred.bag <- colMeans(do.call(rbind, out.r$Ypred))

      out.r <- data.table(
        Treatment=g, 
        YpredBag=(which.max(ypred.bag) - 1),
        Compound_Category=out$Compound_Category[1],
        Ytrue=out$Ytrue[1],
        Rep=rep
      )
    })
    
    return(rbindlist(out))
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
  out <- rep(NA, length(y))
  names(out) <- 1:length(y)
  
  fit <- lm(y ~ as.numeric(col))
  out[names(fit$residuals)] <- fit$residuals
  return(out)
}


################################################################################
# Functions for evaluating bioactivity of compounds
################################################################################
bioactivity <- function(x, xcov.inv.sqrt, null_summary=max) {
  # Wrapper function to compute distance for each well to DMSO

  # Feature-center data
  center <- function(z) z - mean(z)
  features <- '(^PC|^nonborder)'
  xwell <- select(x, matches(features))
  
  # Transform through inverse sqrt covariance
  qwell <- as.matrix(xwell) %*% xcov.inv.sqrt
  
  # Initialize DMSO point cloud center
  xctl <- apply(filter(xwell, x$Compound_ID == 'DMSO'), MAR=2, median)
  qctl <- c(xcov.inv.sqrt %*% xctl)
  
  # Compute mahalanobis distance
  well.dist <- sqrt(colMeans((t(qwell) - qctl) ^ 2))
  
  # Generate null distribution of maximum distance DMSO to DMSO center
  null <- well.dist[x$Compound_ID == 'DMSO']
  null.thresh <- null_summary(null)
  
  xdist <- data.frame(Dist=well.dist) %>% mutate(DistNorm=Dist / null.thresh)
  xmeta <- dplyr::select(x, !matches('(^PC|^nonborder)'))
  
  return(cbind(xdist, xmeta))
}


null_dist <- function(x, xcov.inv.sqrt) {
  # Wrapper function to compute null distribution of DMSO to centroid distances
  
  # Feature-center data
  center <- function(z) z - mean(z)
  features <- '(^PC|^nonborder)'
  xwell <- select(x, matches(features))
  
  # Transform through inverse sqrt covariance
  qwell <- as.matrix(xwell) %*% xcov.inv.sqrt
  
  # Initialize DMSO point cloud center
  xctl <- apply(filter(xwell, x$Compound_ID == 'DMSO'), MAR=2, median)
  qctl <- c(xcov.inv.sqrt %*% xctl)
  
  # Compute mahalanobis distance
  well.dist <- sqrt(colMeans((t(qwell) - qctl) ^ 2))
  return(well.dist[x$Compound_ID == 'DMSO'])
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

