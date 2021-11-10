#+ setup, echo=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(iRF)
library(caret)
library(superheat)
library(ggsci)

setwd('~/github/cancer_translator/')
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

################################################################################
# Initialize analysis parameters
################################################################################
intensity.normalize <- TRUE
n.core <- 6
cell.line <- 'A549'
prop.train <- 0.8

compound.set <- c(
  'HSP90',
  'MT',
  'mTOR',
  'DNA',
  'Proteasome',
  'HDAC'
)

fin <- str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata')
load(fin)

# Initialize output files
output.dir <- 'results/marker/'
dir.create(output.dir, showWarnings=TRUE)
bioactivity.file <- str_c(output.dir, 'bioactivity_', cell.line, '.Rdata') 
output.file <- str_c(output.dir, 'classification_', cell.line, '.Rdata') 

# Initialize color palettes
heat.pal <- c('#FFFFFF', pal_material('light-green')(10))
group.pal <- pal_jco()(10)[c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10)]

# Filter to selected cell line
x <- filter(x, Cell_Line == cell.line)

# Clean compound category names
select_category <- function(x) na.omit(x)[1]

xcat.key <- select(x, Compound_ID, Compound_Category) %>%
  group_by(Compound_ID) %>%
  summarize(Compound_Category=select_category(Compound_Category))

################################################################################
# Initialize marer sets evaluated in analysis
################################################################################
# Initialize feature sets for each marker
colnames(x) <- str_replace_all(colnames(x), 'MitoTracker', 'Mitotracker')

markers <- c('HOECHST', 'Mitotracker', 'SYTO14', 'Alexa')
markers.re <- str_c('(', str_c(markers, collapse='|'), ')')
features.marker <- lapply(markers, function(m) str_subset(colnames(x), m))

id.morphology <- !str_detect(colnames(x), markers.re)
features.morphology <- colnames(x)[id.morphology] %>% str_subset('nonborder')

features <- c(features.marker, list(features.morphology))
names(features) <- c(markers, 'morphology')

# Initialize marker set combinations
n.marker <- length(markers)
marker.comb <- lapply(1:n.marker, function(k) combn(markers, k, simplify=FALSE))
marker.comb <- unlist(marker.comb, recursive=FALSE)
n.markers <- sapply(marker.comb, length)

marker.comb.m <- lapply(marker.comb, function(m) c(m, 'morphology'))
marker.comb <- c(marker.comb, marker.comb.m, 'morphology')

# Initialize marker set groups
groups <- str_c(n.markers, ' markers')
groups <- c(groups, str_c(n.markers, ' markers, morphology'), 'morphology')
names(groups) <- sapply(marker.comb, str_c, collapse='_')

#' # Overview
#' This notebook considers the problems of classifying compound category based 
#' on a the cell painting markerset. In addition to models based on the full 
#' markerset, we consider how prediction accuracy varies across subsets of 
#' markers.
#'
#' ## Bioactivity filtering
#' Before training compound classifiers, we filter out treatments (i.e. 
#' compounds/dose pairs) that cannot be distinguished from DMSO according to 
#' bioactivity scores (see bioactivity analysis).
#' 
#' **Note:** bioactivity filtering is based on all markers.
#' 
#' ## Compound filtering
#' The following analyses evaluate compound classification relative to the 
#' reference compound categories used in ORACL.
#+ load_bioactivity
# Load bioactivity table
setwd('~/github/cancer_translator/')
load(bioactivity.file)

# Filter to reference compound doses that are bioactive in the full markerset
xdist.select <- filter(xdist, Compound_Usage == 'reference_cpd') %>%
  filter(Markerset == 'HOECHST_Mitotracker_SYTO14_Alexa_morphology') %>%
  filter(pval == 0) %>%
  mutate(Treatment=str_c(Compound_ID, ', ', Dose_Category))

xf <- mutate(x, Treatment=str_c(Compound_ID, ', ', Dose_Category)) %>%
  filter(Treatment %in% xdist.select$Treatment) %>%
  select(!matches('^PC')) %>%
  filter(Compound_Usage == 'reference_cpd') %>%
  filter(Compound_Category %in% compound.set)

# Table of unique compounds by category
compound.table <- group_by(xf, Compound_Category) %>%
  summarize(Ncpd=length(unique(Compound_ID)))

print(compound.table)

#' # Modeling
#' For each cell line, we train classifiers to predict compound category from
#' phenotypic profiling features.
#' 
#' First, we assess performance relative to a random holdout set. `r prop.train`
#' of wells are randomly sampled to train models and the remaining wells are 
#' used to assess accuracy. **Note:** based on this formulation, a compound can 
#' appear in both the training and test sets but as a different (well) replicate.
#+ random_holdout, fig.height=8, fig.width=15, message=FALSE, echo=FALSE
################################################################################
# Random holdout predictions
################################################################################
# Fit models for each cell line
ypred <- mclapply(
  marker.comb, 
  fit_marker,
  x=xf,
  model=irf,
  model_predict=irf_predict,
  prop=prop.train,
  mc.cores=n.core
)

ypred <- rbindlist(ypred) %>%
  mutate(Treatment=str_c(Compound_ID, ', ', Dose_Category))

#+ random_holdout_acc, fig.height=8, fig.width=18, message=FALSE
# Aggregate model predictions for visualization
xplot.group <- group_by(ypred, Markerset) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  mutate(Group=groups[Markerset])

xplot.group %>%
  ggplot(aes(x=reorder(Markerset, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Group)) +
  geom_text(aes(label=Accuracy), size=2, nudge_y=0.02) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  ylim(c(0, 1.05))

#+ random_holdout_acc_dose, fig.height=12, fig.width=18, message=FALSE
# Accuracy by cell line, dose
xplot <- ypred %>%
  mutate(Dose=str_remove_all(Treatment, '^.*, ')) %>%
  group_by(Markerset, Dose) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Group=groups[Markerset]) %>%
  arrange(Dose, Markerset)

# Accuracy by compound category heat map
compounds <- unique(xplot$Compound_Category)
markersets <- unique(xplot$Markerset)
xplot.hm <- matrix(xplot$Accuracy, nrow=length(markersets))

colnames(xplot.hm) <- compounds
rownames(xplot.hm) <- markersets

# Initialize row attributes
group <- groups[markersets]
row.pal <- group.pal[as.factor(group)]

superheat(
  xplot.hm,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  membership.rows=group,
  left.label.col=row.pal,
  left.label='variable',
  left.label.text.size=3,
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label.text.angle=90
)


#+ random_holdout_acc_hm, fig.height=12, fig.width=24, message=FALSE
# Accuracy by compound category
xplot <- group_by(ypred, Markerset, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Group=groups[Markerset]) %>%
  arrange(Compound_Category, Markerset)

# Accuracy by compound category heat map
compounds <- unique(xplot$Compound_Category)
markersets <- unique(xplot$Markerset)
xplot.hm <- matrix(xplot$Accuracy, nrow=length(markersets))

colnames(xplot.hm) <- compounds
rownames(xplot.hm) <- markersets

# Initialize row attributes
group <- groups[markersets]
row.pal <- group.pal[as.factor(group)]

superheat(
  xplot.hm,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  membership.rows=group,
  left.label.col=row.pal,
  left.label='variable',
  left.label.text.size=3,
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label.text.angle=90
)

ypred.random <- ypred

#' #### Compound holdout
#' Next, we assess performance relative to a randomly held out compounds. That 
#' is, 4 compounds from each category randomly sampled to train models and the 
#' remaining compound from each category is used to assess accuracy. **Note:** 
#' models here are evaluated on compounds they have never seen.
#' 
#' Plots below are the same as above but using hold-out compounds.
#+ compound_holdout, fig.height=8, fig.width=18, message=FALSE, echo=FALSE
################################################################################
# Compound holdout predictions
################################################################################
# Fit models for each cell line
ypred <- mclapply(
  marker.comb, 
  fit_marker,
  x=xf,
  model=irf,
  model_predict=irf_predict,
  prop=prop.train,
  holdout='compound',
  mc.cores=n.core
)

ypred <- rbindlist(ypred) %>%
  mutate(ypred, Treatment=str_c(Compound_ID, ', ', Dose_Category))

#+ cpd_holdout_acc, fig.height=8, fig.width=18, message=FALSE
# Aggregate model predictions for visualization
xplot.group <- group_by(ypred, Markerset) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  mutate(Group=groups[Markerset])

xplot.group %>%
  ggplot(aes(x=reorder(Markerset, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Group)) +
  geom_text(aes(label=Accuracy), size=2, nudge_y=0.02) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  ylim(c(0, 1.05))

#+ cpd_holdout_acc_dose, fig.height=12, fig.width=18, message=FALSE
# Accuracy by cell line, dose
# Accuracy by cell line, dose
xplot <- ypred %>%
  mutate(Dose=str_remove_all(Treatment, '^.*, ')) %>%
  group_by(Markerset, Dose) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Group=groups[Markerset]) %>%
  arrange(Dose, Markerset)

# Accuracy by compound category heat map
compounds <- unique(xplot$Compound_Category)
markersets <- unique(xplot$Markerset)
xplot.hm <- matrix(xplot$Accuracy, nrow=length(markersets))

colnames(xplot.hm) <- compounds
rownames(xplot.hm) <- markersets

# Initialize row attributes
group <- groups[markersets]
row.pal <- group.pal[as.factor(group)]

superheat(
  xplot.hm,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  membership.rows=group,
  left.label.col=row.pal,
  left.label='variable',
  left.label.text.size=3,
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label.text.angle=90
)

#+ cpd_holdout_acc_hm, fig.height=12, fig.width=24, message=FALSE
# Accuracy by compound category
xplot <- group_by(ypred, Markerset, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Group=groups[Markerset]) %>%
  arrange(Compound_Category, Markerset)

# Accuracy by compound category heat map
compounds <- unique(xplot$Compound_Category)
markersets <- unique(xplot$Markerset)
xplot.hm <- matrix(xplot$Accuracy, nrow=length(markersets))

colnames(xplot.hm) <- compounds
rownames(xplot.hm) <- markersets

# Initialize row attributes
group <- groups[markersets]
row.pal <- group.pal[as.factor(group)]

superheat(
  xplot.hm,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  membership.rows=group,
  left.label.col=row.pal,
  left.label='variable',
  left.label.text.size=3,
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label.text.angle=90
)

ypred.compound <- ypred