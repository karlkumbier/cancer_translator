library(data.table)
library(tidyverse)
library(rprofiler)
library(superheat)
library(tidytext)

################################################################################
# Helper functions
################################################################################
summarize_cell_line <- function(x, cell.line) {
  x <- filter(x, Cell_Line == cell.line) %>%
    group_by(Compound_ID, Compound_Category) %>%
    summarize(DistNorm=mean(DistNorm))
  return(x)
}

intensity_normalize <- function(x) {
  # Wrapper function to normalize all ntensity features
  col <- as.numeric(x$Col)
  #intensity.features <- colnames(x) %>% str_subset('Intensity')
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


compound_dist <- function(x) {
  # Wrapper function to compute distance for each well to DMSO
  wells <- x$WellID
  well.dist <- lapply(wells, well_dist, x=x)
  
  # Generate null distribution
  null <- generate_null(x, n.core=n.core)
  
  # Normalize relative to control
  xdist <- data.frame(Dist=unlist(well.dist)) %>%
    mutate(DistNorm=Dist / mean(null)) %>%
    mutate(pval=sapply(Dist, function(d) mean(null > d)))
    
  xmeta <- dplyr::select(x, !matches('(^PC|^nonborder)'))
  return(cbind(xdist, xmeta))
}

well_dist <- function(x, well.id) {
  # Compute average distance between compounds of same type and to control
  # TODO: vectorize
  xf <- filter(x, WellID == well.id) %>% select(matches('(^PC|^nonborder)'))
  xctl <- filter(x, Compound_ID == 'DMSO') %>% select(matches('(^PC|^nonborder)'))
  xctl.mean <- colMeans(xctl)
  
  xdist <- colMeans((t(as.matrix(xf)) - xctl.mean) ^ 2)
  xdist <- mean(xdist)
  
  return(xdist)
}

compound_dist_ <- function(x, cpd) {
  # Compute average distance between compounds of same type and to control
  xf <- filter(x, Compound_ID == cpd) %>% select(matches('(^PC|^nonborder)'))
  xctl <- filter(x, Compound_ID == 'DMSO') %>% select(matches('(^PC|^nonborder)'))
  xctl.mean <- colMeans(xctl)
  
  if (cpd == 'DMSO') {
    xdist <- colMeans((t(as.matrix(xf)) - xctl.mean) ^ 2)
    xdist <- max(xdist)  
  } else {
    xdist <- colMeans((t(as.matrix(xf)) - xctl.mean) ^ 2)
    xdist <- mean(xdist)
  }
  
  return(data.frame(Compound_ID=cpd, CpdDist=xdist))
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
  dist <- compound_dist_(x, 'DMSO')
  return(dist$CpdDist)
}

################################################################################
# Setup
################################################################################
col.pal <- RColorBrewer::brewer.pal(11, 'RdYlBu')
col.pal[6] <- '#FFFFFF'
intensity.normalize <- TRUE
n.core <- 16

setwd('~/github/cancer_translator/')
data.dir <- 'data/screens/LH_CDC_1/'
load(str_c(data.dir, 'ks_profiles.Rdata'))

xks <- mutate(xks, ID=str_c(PlateID, '_', WellID)) %>%
  dplyr::select(-PlateID, -WellID)

xmeta <- mutate(xmeta, ID=str_c(PlateID, '_', WellID))
x <- left_join(xks, xmeta, by='ID') %>%
  mutate(Control=Compound_ID == 'DMSO') %>%
  distinct() %>%
  filter(!is.na(PlateID)) %>%
  mutate(Col=as.numeric(sapply(str_split(WellID, '-'), tail, 1)))

# TODO: why do we need distinct?
# TODO: what are NA plates?

# Normalize intensity within each plate
if (intensity.normalize) {
  plates <- unique(x$PlateID)
  x <- lapply(plates, function(p) {
    intensity_normalize(filter(x, PlateID == p))
  })
  
  x <- rbindlist(x)
}

################################################################################
# Simple EDA - counts and intensity
################################################################################
# Cell count plots
filter(x, Control) %>%
  ggplot(aes(x=NCells)) +
  geom_histogram() +
  theme_bw() +
  ggtitle('Cell counts, DMSO wells')

filter(x, Control) %>%
  ggplot(aes(x=NCells, fill=Cell_Line)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~PlateID) +
  ggtitle('Cell counts, DMSO wells')


filter(x, Control) %>%
  ggplot(aes(x=NCells)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~Cell_Line) +
  ggtitle('Cell counts, DMSO wells')

# Drop plate 2021018028 - strange cell count distribution
x <- filter(x, PlateID != '2021018028')

# Feature distribution plots
xcontrol <- filter(x, Control)
plate.id <- xcontrol$PlateID
cell.line <- xcontrol$Cell_Line
xcontrol <- dplyr::select(xcontrol, matches('^nonborder'))
colnames(xcontrol) <- str_remove_all(colnames(xcontrol), 'nonborder\\.\\.\\.')

superheat(xcontrol,
          membership.rows=cell.line,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=col.pal,
          heat.pal.values=seq(0, 1, by=0.1),
          bottom.label.text.angle=90,
          bottom.label.text.size=3,
          bottom.label.size=0.75)

xcontrol[xcontrol < -0.3] <- -0.3
xcontrol[xcontrol > 0.3] <- 0.3
superheat(xcontrol,
          membership.rows=cell.line,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=col.pal,
          heat.pal.values=seq(0, 1, by=0.1),
          bottom.label.text.angle=90,
          bottom.label.text.size=3,
          bottom.label.size=0.75)

# Intensity distribution plots
filter(x, Control) %>%
  dplyr::select(matches('(Intensity.*Mean|PlateID|WellID)')) %>%
  mutate(PlateID=as.factor(PlateID)) %>%
  mutate(Col=sapply(str_split(WellID, '-'), tail, 1)) %>%
  reshape2::melt() %>%
  mutate(variable=str_remove_all(variable, 'nonborder\\.\\.\\.')) %>%
  ggplot(aes(x=PlateID, y=value, fill=Col)) +
  geom_boxplot() +
  facet_wrap(~variable) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  ylab('KS')


################################################################################
# Full dataset PCA
################################################################################
xfeat <- dplyr::select(x, matches('nonborder'))
xpca <- prcomp(xfeat)
x <- cbind(xpca$x, x)

data.frame(PctVar=cumsum(xpca$sdev ^ 2) / sum(xpca$sdev ^ 2)) %>%
  mutate(NPC=1:n()) %>%
  ggplot(aes(x=NPC, y=PctVar)) +
  geom_point() +
  geom_line() +
  theme_bw()

pc.list <- list(c('PC1', 'PC2'), c('PC3', 'PC4'), c('PC5', 'PC6'))
for (pcs in pc.list) {
  p <- ggplot(x, aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=Control), alpha=0.5) +
    facet_wrap(~PlateID) +
    theme_bw()
  plot(p)
  
  p <- ggplot(x, aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=log(NCells), shape=Control), alpha=0.5) +
    facet_wrap(~PlateID) +
    theme_bw() +
    scale_color_gradientn(colors=viridis::inferno(10)[1:9])
  plot(p)
  
  p <- filter(x, Control) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=Cell_Line), alpha=0.5) +
    theme_bw +
    ggtitle('DMSO wells by cell line')
  plot(p)
  
  p <- sample_frac(x, 0.05) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=Compound_Usage), alpha=0.5) +
    theme_bw() +
    ggtitle('Compound usage, randomly sampled wells (0.1)')
  plot(p)
}

################################################################################
# Bioactivity
################################################################################

# Normalize data matrix
l2norm <- function(z) z / sqrt(sum(z ^ 2))
x <- mutate_if(x, is.numeric, l2norm)

xdist <- mclapply(unique(x$PlateID), function(p) {
    out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
    return(compound_dist(out))
}, mc.cores=n.core)

# Note: not taking into accoutn dose here
xdist <- rbindlist(xdist)
group_by(xdist, Compound_Category, Cell_Line) %>%
  summarize(pval=mean(pval), DistNorm=mean(DistNorm)) %>%
  filter(!is.na(Compound_Category)) %>%
  ggplot(aes(x=reorder_within(Compound_Category, DistNorm, Cell_Line), y=DistNorm)) +
  geom_bar(stat='identity', position='dodge', aes(fill=pval)) +
  theme_bw() +
  facet_wrap(Cell_Line~., scales='free') +
  coord_flip()

xdist.786 <- summarize_cell_line(xdist, '786-0')
xdist.a549 <- summarize_cell_line(xdist, 'A549')
xdist.ovcar4 <- summarize_cell_line(xdist, 'OVCAR4')

left_join(xdist.786, xdist.a549, by='Compound_ID') %>%
  filter(!is.na(Compound_Category.y)) %>%
  ggplot(aes(x=log(DistNorm.x), y=log(DistNorm.y), col=Compound_Category.y)) +
  geom_point() +
  theme_bw() +
  xlab('786-0 log l2 distance from DMSO') +
  ylab('A549 log l2 distance from DMSO') +
  geom_hline(yintercept=0, col='grey') +
  geom_vline(xintercept=0, col='grey')

left_join(xdist.ovcar4, xdist.a549, by='Compound_ID') %>%
  filter(!is.na(Compound_Category.y)) %>%
  ggplot(aes(x=log(DistNorm.x), y=log(DistNorm.y), col=Compound_Category.y)) +
  geom_point() +
  theme_bw() +
  xlab('OVCAR4 log l2 distance from DMSO') +
  ylab('A549 log l2 distance from DMSO') +
  geom_hline(yintercept=0, col='grey') +
  geom_vline(xintercept=0, col='grey')

left_join(xdist.ovcar4, xdist.786, by='Compound_ID') %>%
  filter(!is.na(Compound_Category.y)) %>%
  ggplot(aes(x=log(DistNorm.x), y=log(DistNorm.y), col=Compound_Category.y)) +
  geom_point() +
  theme_bw() +
  xlab('OVCAR4 log l2 distance from DMSO') +
  ylab('786-0 log l2 distance from DMSO') +
  geom_hline(yintercept=0, col='grey') +
  geom_vline(xintercept=0, col='grey')
  
  

  
