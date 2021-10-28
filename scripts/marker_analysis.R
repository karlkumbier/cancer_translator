#+ setup, echo=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(iRF)
library(caret)

col.pal <- RColorBrewer::brewer.pal(11, 'RdYlBu')
col.pal[6] <- '#FFFFFF'
intensity.normalize <- TRUE
n.core <- 6

setwd('~/github/cancer_translator/')
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

intensity.normalize <- TRUE
fin <- str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata')
load(fin)

# Filter to A549
x <- filter(x, Cell_Line == 'A549')

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
setwd('~/github/cancer_translator/')

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
  filter(pval == 0) %>%
  mutate(Compound_Dose=str_c(Compound_ID, ', ', Dose_Category))

#' # Modeling
#' For each cell line, we train classifiers to predict compound category from
#' phenotypic profiling features. Compounds/doses are filtered to include only 
#' those that are bioactive in at least one cell line.
#+ modeling, fig.height=8, fig.width=18, message=FALSE, warnings=FALSE, echo=FALSE
################################################################################
# Modeling
################################################################################
xf <- mutate(x, Compound_Dose=str_c(Compound_ID, ', ', Dose_Category)) %>%
  filter(Compound_Dose %in% xdist.select$Compound_Dose) %>%
  select(!matches('^PC')) %>%
  filter(Compound_Usage == 'reference_cpd')

# Table of unique compounds by category
compound.table <- group_by(xf, Compound_Category) %>%
  summarize(Ncpd=length(unique(Compound_ID)))

# Filter out compounds with low prevalence
cpd.keep <- filter(compound.table, Ncpd == 5)
xf <- filter(xf, Compound_Category %in% cpd.keep$Compound_Category)

#' First, we assess performance relative to a random holdout set. That is, 80%
#' of wells are randomly sampled to train models and the remaining 20% are used 
#' to assess accuracy. Note: a compound can appear in both the training and test 
#' sets but at different doses.
#+ random_holdout, fig.height=8, fig.width=18, message=FALSE, echo=FALSE
################################################################################
# Random holdout predictions
################################################################################
markers <- c('HOECHST', 'Mitotracker', 'SYTO14', 'Alexa')
markers.re <- str_c('(', str_c(markers, collapse='|'), ')')

id.morphology <- !str_detect(colnames(xf), markers.re  )
id.feat <- str_detect(colnames(xf), 'nonborder')
id <- id.feat & id.morphology
colnames(xf)[id] <- str_c(colnames(xf)[id], 'morphology')

markers1 <- lapply(markers, function(m) c(m, 'morphology'))
markers2 <- combn(markers, 2, simplify=FALSE)
markers2 <- lapply(markers2, function(m) c(m, 'morphology'))

markers4 <- list(markers)
markers11 <- lapply(markers, function(m) m)

# Fit models for each cell line
markers <- c(markers1, markers2, markers11, markers4)
ypred <- lapply(markers, 
                fit_marker,
                x=xf,
                model=irf,
                model_predict=irf_predict,
                prop=0.8
)

ypred <- rbindlist(ypred)

xplot.marker <- group_by(ypred, Marker) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.marker %>%
  mutate(Marker=str_replace_all(Marker, '\\|', ', ')) %>%
  ggplot(aes(x=reorder(Marker, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', fill='#0088D1') +
  geom_text(aes(x=Marker, y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  ylim(c(0, 1)) +
  theme(axis.text.x=element_text(angle=90))

# Accuracy by compound category
xplot.marker <- group_by(ypred, Marker, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue))

xplot.marker %>%
  mutate(Marker=str_replace_all(Marker, '\\|', ', ')) %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  ggplot(aes(x=reorder(Marker, Accuracy), y=Accuracy, fill=Compound_Category)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label=Accuracy), position=position_dodge(width=0.9), vjust=-0.02, size=3) +
  theme_bw() +
  ylim(0:1) +
  theme(axis.text.x=element_text(angle=90))

#' Next, we assess performance relative to a randomly held out compounds. That 
#' is, 4 compounds from each category randomly sampled to train models and the 
#' remaining compound from each category is used to assess accuracy. Note: 
#' models here are evaluated on compounds they have never seen
#+ compound_holdout, fig.height=8, fig.width=18, message=FALSE, echo=FALSE
################################################################################
# Compound holdout predictions
################################################################################
# Fit models for each cell line
ypred <- lapply(markers, 
                fit_marker,
                x=xf,
                model=irf,
                model_predict=irf_predict,
                holdout='compound'
)

ypred <- rbindlist(ypred)

xplot.marker <- group_by(ypred, Marker) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.marker %>%
  mutate(Marker=str_replace_all(Marker, '\\|', ', ')) %>%
  ggplot(aes(x=reorder(Marker, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', fill='#0088D1') +
  geom_text(aes(x=Marker, y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  ylim(c(0, 1)) +
  theme(axis.text.x=element_text(angle=90))

# Accuracy by compound category
xplot.marker <- group_by(ypred, Marker, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue))

xplot.marker %>%
  mutate(Marker=str_replace_all(Marker, '\\|', ', ')) %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  ggplot(aes(x=reorder(Marker, Accuracy), y=Accuracy, fill=Compound_Category)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label=Accuracy), position=position_dodge(width=0.9), vjust=-0.02, size=3) +
  theme_bw() +
  ylim(0:1) +
  theme(axis.text.x=element_text(angle=90))
