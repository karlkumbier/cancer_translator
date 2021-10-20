#+ setup, echo=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(iRF)
library(caret)
library(superheat)

col.pal <- RColorBrewer::brewer.pal(11, 'RdYlBu')
col.pal[6] <- '#FFFFFF'
intensity.normalize <- TRUE
n.core <- 12
write.tables <- FALSE

setwd('~/github/cancer_translator/')
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

intensity.normalize <- TRUE
fin <- str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata')
load(fin)

# Filter to A549
x <- filter(x, Cell_Line == 'A549')

options(knitr.table.format = function() {
  if (knitr::is_latex_output())
    "latex" else "pipe"
})

#' # Overview
#' In this notebook, we consider the utility of different marker sets for 
#' compound target prediction. All analyses are restricted to A549 cell line.
#' 
#' # Bioactivity filtering
#' Before training compound classifiers, we filter out compounds/dose pairs that
#' cannot be distinguished from DMSO. Bioactivity relative to DMSO is evaluated 
#' as follows:
#' 
#' 1. Subsample DMSO wells
#' 2. Compute the center of the DMSO point cloud by taking average feature 
#' values across each well in the subsample.
#' 3. Compute the l2 distance between each DMSO well in the subsample and the 
#' DMSO point cloud center.
#' 4. Define the DMSO maximal distance as the maximum distance (from 3) over all 
#' subsampled wells.
#' 5. Any well (i.e. compound) that is further from the DMSO point cloud center 
#' than the maximal DMSO distance is called as bioactive in the given subsample. 
#'
#' Repeating this process across many subsamples allows us to generate a 
#' bioactivity p-value. Tablesbelow summarize bioactivity by plate, compound 
#' category (in reference set), and cell line.
#+ bioactivity
################################################################################
# Bioactivity
################################################################################
setwd('~/github/cancer_translator/')

# Normalize data matrix
l2norm <- function(z) z / sqrt(sum(z ^ 2))
x <- mutate_if(x, is.numeric, l2norm)

# Compute bioactivity scores for each well by plate
xdist <- mclapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(bioactivity(out))
}, mc.cores=n.core)

xdist <- rbindlist(xdist)

# Bioactivity by plate
bioactive.plate <- group_by(xdist, PlateID) %>% 
  summarize(ProportionBioactive=mean(pval == 0))

knitr::kable(bioactive.plate)

if (write.tables) {
  write.csv(file='results/bioactivity_plate.csv', bioactive.plate, quote=FALSE)
}

# Bioactivity by compound category
bioactive.category <- filter(xdist, Compound_Usage == 'reference_cpd') %>%
  group_by(Compound_Category) %>% 
  summarize(ProportionBioactive=mean(pval == 0))

knitr::kable(bioactive.category)

if (write.tables) {
write.csv(file='results/bioactivity_category.csv', bioactive.category, quote=FALSE)
}

# Filter to reference compound doses that are bioactive
xdist.select <- filter(xdist, Compound_Usage == 'reference_cpd') %>%
  filter(pval == 0) %>%
  mutate(Compound_Dose=str_c(Compound_ID, ', ', Dose_Category))

#' # Modeling
#' For each marker set, we train classifiers to predict compound category from
#' phenotypic profiling features. Compounds/doses are filtered to include only 
#' those that are bioactive, while categories are filtered to include only 
#' reference compound categories with >= 5 compounds.
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
cpd.keep <- filter(compound.table, Ncpd >= 5)
xf <- filter(xf, Compound_Category %in% cpd.keep$Compound_Category)

#' Below, we assess the contribution of different channels to target prediction
#' accuracy. Channels are defined as one of: `HOECHST`, `Mitotracker`, `SYTO14`,
#' `Alexa`, or `morphology` — where `morphology` is any feature not specific to 
#' one of the first four channels. We train each model using the subset of 
#' features indicated and evaluate performance on a hold-out test set.
#' 
#' First, we assess performance relative to a random holdout set. That is, 80%
#' of wells are randomly sampled to train models and the remaining 20% are used 
#' to assess accuracy. Note: a compound can appear in both the training and test 
#' sets but at different doses.
#+ random_holdout, fig.height=12, fig.width=18, message=FALSE, echo=FALSE
################################################################################
# Random holdout predictions
################################################################################
markers <- c('HOECHST', 'Mitotracker', 'SYTO14', 'Alexa')
markers.re <- str_c('(', str_c(markers, collapse='|'), ')')

id.morphology <- !str_detect(colnames(xf), markers.re )
id.feat <- str_detect(colnames(xf), 'nonborder')
id <- id.feat & id.morphology
colnames(xf)[id] <- str_c(colnames(xf)[id], 'morphology')

# Generate marker set combinations to be evaluated
markers1 <- markers
markers1.morph <- lapply(markers, function(m) c(m, 'morphology'))

markers2 <- combn(markers, 2, simplify=FALSE)
markers2.morph <- lapply(markers2, function(m) c(m, 'morphology'))

markers3 <- combn(markers, 3, simplify=FALSE)
markers3.morph <- lapply(markers3, function(m) c(m, 'morphology'))

markers4 <- list(markers)
markers4.morph <- lapply(markers4, function(m) c(m, 'morphology'))

markers <- c(
  'morphology',
  markers1, 
  markers1.morph, 
  markers2, 
  markers2.morph, 
  markers3, 
  markers3.morph, 
  markers4, 
  markers4.morph
)

group <- c(
  'morphology',
  rep(c('1 marker', '1 marker, morphology'), each=length(markers1)),
  rep(c('2 marker', '2 marker, morphology'), each=length(markers2)),
  rep(c('3 marker', '3 marker, morphology'), each=length(markers3)),
  rep(c('4 marker', '4 marker, morphology'), each=length(markers4))
)

group.key <- data.table(Marker=markers, Group=group) %>%
  mutate(Marker=sapply(Marker, function(m) str_c(m, collapse='|')))

# Fit models for each cell line
ypred <- mclapply(markers, 
                fit_marker,
                x=xf,
                model=irf,
                model_predict=irf_predict,
                prop=0.8,
                mc.cores=n.core
)

ypred <- rbindlist(ypred)

xplot.marker <- group_by(ypred, Marker) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  left_join(group.key, by='Marker')

xplot.marker %>%
  mutate(Marker=str_replace_all(Marker, '\\|', ', ')) %>%
  ggplot(aes(x=reorder(Marker, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Group)) +
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
  ggplot(aes(x=reorder(Compound_Category, Accuracy), y=Accuracy, fill=Compound_Category)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label=Accuracy), position=position_dodge(width=0.9), vjust=-0.02) +
  theme_bw() +
  theme(legend.position='none') +
  ylim(0:1) +
  facet_wrap(~Marker) +
  theme(axis.text.x=element_text(angle=90))


# Heatmap visualization by category
ncpd <- length(unique(xplot.marker$Compound_Category))
xplot <- matrix(xplot.marker$Accuracy, nrow=ncpd)

rownames(xplot) <- unique(xplot.marker$Compound_Category)
colnames(xplot) <- unique(xplot.marker$Marker)
group.id <- group.key$Group[match(colnames(xplot), group.key$Marker)]

# Threshold at minimum 0.7 for visualization
colnames(xplot) <- str_replace_all(colnames(xplot), '\\|', ', ')
xplot[xplot < 0.7] <- 0.7

superheat(xplot, 
          membership.cols=group.id,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=viridis::inferno(10),
          heat.pal.values=seq(0, 1, by=0.1),
          bottom.label='variable',
          bottom.label.text.angle=90,
          bottom.label.size=1)

#' Next, we assess performance relative to a randomly held out compounds. That 
#' is, 4 compounds from each category randomly sampled to train models and the 
#' remaining compound from each category is used to assess accuracy. Note: 
#' models here are evaluated on compounds they have never seen
#+ compound_holdout, fig.height=12, fig.width=18, message=FALSE, echo=FALSE
################################################################################
# Compound holdout predictions
################################################################################
# Fit models for each cell line
ypred <- mclapply(markers, 
                fit_marker,
                x=xf,
                model=irf,
                model_predict=irf_predict,
                holdout='compound',
                mc.cores=n.core
)

ypred <- rbindlist(ypred)

xplot.marker <- group_by(ypred, Marker) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  left_join(group.key, by='Marker')

xplot.marker %>%
  mutate(Marker=str_replace_all(Marker, '\\|', ', ')) %>%
  ggplot(aes(x=reorder(Marker, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Group)) +
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
  ggplot(aes(x=reorder(Compound_Category, Accuracy), y=Accuracy, fill=Compound_Category)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label=Accuracy), position=position_dodge(width=0.9), vjust=-0.02) +
  theme_bw() +
  theme(legend.position='none') +
  ylim(0:1) +
  facet_wrap(~Marker) +
  theme(axis.text.x=element_text(angle=90))


# Heatmap visualization by category
ncpd <- length(unique(xplot.marker$Compound_Category))
xplot <- matrix(xplot.marker$Accuracy, nrow=ncpd)

rownames(xplot) <- unique(xplot.marker$Compound_Category)
colnames(xplot) <- unique(xplot.marker$Marker)
group.id <- group.key$Group[match(colnames(xplot), group.key$Marker)]

# Threshold at minimum 0.5 for visualization
colnames(xplot) <- str_replace_all(colnames(xplot), '\\|', ', ')
xplot[xplot < 0.5] <- 0.5

superheat(xplot, 
          membership.cols=group.id,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=viridis::inferno(10),
          heat.pal.values=seq(0, 1, by=0.1),
          bottom.label='variable',
          bottom.label.text.angle=90,
          bottom.label.size=1)
