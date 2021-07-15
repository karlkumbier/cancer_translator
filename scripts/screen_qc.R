#+ setup, echo=FALSE, message=FALSE, warning=FALSE
library(data.table)
library(tidyverse)
library(rprofiler)
library(superheat)
library(tidytext)
library(parallel)
library(iRF)
library(caret)

################################################################################
# Setup
################################################################################
col.pal <- RColorBrewer::brewer.pal(11, 'RdYlBu')
col.pal[6] <- '#FFFFFF'
intensity.normalize <- TRUE
n.core <- 16

setwd('~/github/cancer_translator/')
source('scripts/utilities.R')
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

#' #QC
#' The figures below report the distribution of cell counts by plate/cell line.
#' We drop plate `2021018028` which appears to have a strange count distribution 
#' relative to other A549 plates.
#+ eda, fig.height=8, fig.width=12
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

#' Preliminary EDA indicated strong "column" effects, particularly for intensity 
#' features. We normalize all features relative to well position by regressing
#' features on column ID — i.e. we remove the portion of a feature that can be
#' explained by a wells column position. Plots below report the distribution of 
#' features for DMSO wells after normalization.
#+ superheat, fig.height=12, fig.width=18
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

filter(x, Control) %>%
  dplyr::select(matches('(Intensity.*Mean|PlateID|WellID|Cell_Line)')) %>%
  mutate(PlateID=as.factor(PlateID)) %>%
  mutate(Col=sapply(str_split(WellID, '-'), tail, 1)) %>%
  reshape2::melt() %>%
  mutate(variable=str_remove_all(variable, 'nonborder\\.\\.\\.')) %>%
  ggplot(aes(x=Cell_Line, y=value, fill=Col)) +
  geom_boxplot() +
  facet_wrap(~variable) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  ylab('KS')

#' # PCA analysis
#' To qualitatively assess the feature distribution of different compound 
#' classes we plot PCA projections grouped by various metadata features. Broadly
#' speaking, the primary sources of variation in the data correspond to compound
#' categories.
#+ pca, fig.height=6, fig.width=12
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

pc.list <- list(c('PC1', 'PC2'), c('PC3', 'PC4'))
for (pcs in pc.list) {
  
  p <- filter(x, !str_detect(Compound_Usage, '(query|reference)')) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=Compound_Usage), alpha=0.5) +
    theme_bw() +
    ggtitle('Positive v. negative controls') +
    facet_wrap(~Cell_Line)
  plot(p)
  
  p <- filter(x, Compound_Usage == 'reference_cpd' | Compound_ID == 'DMSO') %>%
    mutate(DMSO=Compound_ID == 'DMSO') %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=Compound_Category, shape=DMSO, alpha=DMSO)) +
    theme_bw() +
    scale_shape_manual(values=c(19, 1)) +
    scale_alpha_manual(values=c(0.6, 0.3)) +
    ggtitle('Compound category') +
    facet_wrap(~Cell_Line)
  plot(p)
  
  p <- filter(x, Compound_Usage == 'reference_cpd' | Compound_ID == 'DMSO') %>%
    mutate(DMSO=Compound_ID == 'DMSO') %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=Pathway, shape=DMSO, alpha=DMSO)) +
    theme_bw() +
    scale_shape_manual(values=c(19, 1)) +
    scale_alpha_manual(values=c(0.6, 0.3)) +
    ggtitle('Pathway') +
    facet_wrap(~Cell_Line)
  plot(p)
  
  p <- filter(x, Compound_Usage == 'reference_cpd' | Compound_ID == 'DMSO') %>%
    mutate(DMSO=Compound_ID == 'DMSO') %>%
    mutate(Target=ifelse(Compound_ID == 'DMSO', 'DMSO', Target)) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=Target, shape=DMSO, alpha=DMSO)) +
    theme_bw() +
    scale_shape_manual(values=c(19, 1)) +
    scale_alpha_manual(values=c(0.6, 0.3)) +
    ggtitle('Target') +
    facet_wrap(~Cell_Line)
  plot(p)
}

#' # Bioactivity filtering
#' Before training compound classifiers, we filter out compounds/doses that
#' cannot be distinguished from DMSO. Specifically, we (i) take a subsample
#' of DMSO wells (ii) compute the center of the DMSO point cloud by taking
#' average feature values across each well in the subsample (iii) compute the
#' l2 distance between each DMSO well in the subsample and the DMSO point cloud 
#' center (iv) define the DMSO maximal distance as the maximum distance (from 
#' iii) over all subsampled wells.
#' 
#' We then ask whether a given well (any compound) is further from the DMSO
#' point cloud center than the maximal DMSO distance. Repeating this process
#' across many subsamples allows us to generate a bioactivity p-value.
#+ bioactivity
################################################################################
# Bioactivity
################################################################################
# Normalize data matrix
l2norm <- function(z) z / sqrt(sum(z ^ 2))
x <- mutate_if(x, is.numeric, l2norm)

# Compute bioactivity scores for each well
xdist <- mclapply(unique(x$PlateID), function(p) {
    out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
    return(bioactivity(out))
}, mc.cores=n.core)

xdist <- rbindlist(xdist)

# Filter to reference compound doses that are bioactive
xdist.select <- filter(xdist, Compound_Usage == 'reference_cpd') %>%
  filter(pval <= 0.05) %>%
  mutate(Compound_Dose=str_c(Compound_ID, ', ', Dose_Category))

#' # Modeling
#' For each cell line, we train classifiers to predict compound category from
#' phenotypic profiling features. Compounds/doses are filtered to include only 
#' those that are bioactive in at least one cell line. The test set is generated
#' by randomly sampling compounds (and taking all bioactive doses).
#+ modeling, fig.height=8, fig.width=15, message=FALSE, warnings=FALSE, echo=FALSE
################################################################################
# Modeling
################################################################################
xf <- mutate(x, Compound_Dose=str_c(Compound_ID, ', ', Dose_Category)) %>%
  filter(Compound_Dose %in% xdist.select$Compound_Dose) %>%
  select(!matches('^PC')) %>%
  filter(Compound_Usage == 'reference_cpd')

# Fit models for each cell line
cell.lines <- unique(xf$Cell_Line)
models <- lapply(cell.lines, function(cell.line) {
  xc <- filter(xf, Cell_Line == cell.line)
  
  
  # Fit model
  xx <- select(xc, matches('^nonborder'))
  yy <- as.factor(xc$Compound_Category)
  
  # Set training set
  id.train <- createDataPartition(yy, times=25, p=0.9)
  
  out <- lapply(1:length(id.train), function(i) {
    
    ii <- id.train[[i]]
    xup <- upSample(xx[ii,], yy[ii])
    ytrain <- xup$Class
    xtrain <- select(xup, -Class)
    
    fit <- iRF(x=xtrain, 
               y=ytrain, 
               n.iter=1, 
               type='ranger', 
               n.core=n.core)
  
    # Evaluate model class predictions
    yy <- as.numeric(yy) - 1
    ypred.class <- predict(fit$rf.list, xx[-ii,])$predictions
    
    # Evaluate model probability predictions
    ypred <- predict(fit$rf.list, xx[-ii,], predict.all=TRUE)$predictions
    ypred <- sapply(sort(unique(yy)), function(cl) rowMeans(ypred == cl))
    ypred <- lapply(1:nrow(ypred), function(k) ypred[k,])
  
    out <- select(xc, -matches('^nonborder'))[-ii,] %>%
      mutate(YpredCl=ypred.class) %>%
      mutate(Ypred=ypred) %>%
      mutate(Ytrue=yy[-ii]) %>%
      mutate(Rep=i)
    
    return(out)
  })
  
  out <- rbindlist(out)
})

models <- rbindlist(models)

# Group model predictions across cell lines
model.bag <- lapply(unique(models$Compound_Dose), function(g) {
  out <- filter(models, Compound_Dose == g) %>% 
    select(Compound_ID, Ypred, Ytrue, Compound_Category)
  
  ypred <- colMeans(do.call(rbind, out$Ypred))
  out <- data.table(Compound_ID=out$Compound_ID[1]) %>%
    mutate(Ytrue=out$Ytrue[1]) %>%
    mutate(Compound_Dose=g) %>%
    mutate(Ypred=list(ypred)) %>%
    mutate(Compound_Category=out$Compound_Category[1]) %>%
    mutate(YpredCl=(which.max(ypred) - 1))
  return(out)
})

model.bag <- rbindlist(model.bag)
bagg.acc <- mean(model.bag$YpredCl == model.bag$Ytrue)

# Accuracy by cell line
group_by(models, Cell_Line) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  ggplot(aes(x=Cell_Line, y=Accuracy)) +
  geom_bar(stat='identity', fill='#0088D1') +
  geom_text(aes(x=Cell_Line, y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  geom_hline(yintercept=bagg.acc, col='red') +
  ylim(c(0, 1))

# Accuracy by compound category
xplot.cell <- group_by(models, Cell_Line, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue))

xplot.bag <- group_by(model.bag, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Cell_Line='Aggregate')

rbind(xplot.cell, xplot.bag) %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=c('Aggregate', 'A549', '786-0', 'OVCAR4'))) %>%
  ggplot(aes(x=Cell_Line, y=Accuracy, fill=Compound_Category)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label=Accuracy), position=position_dodge(width=0.9), vjust=-0.02) +
  theme_bw() +
  ylim(0:1)

