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

