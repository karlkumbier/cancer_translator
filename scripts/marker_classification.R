#+ setup, echo=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(iRF)
library(caret)
library(superheat)
library(ggsci)

analysis.dir <- '~/github/cancer_translator/'
setwd(analysis.dir)
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

################################################################################
# Initialize analysis parameters
################################################################################
intensity.normalize <- TRUE
n.core <- 16
cell.line <- 'A549'
prop.train <- 0.8
reps <- 50

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
dir.create(output.dir, showWarnings=FALSE)
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
#' markers, compounds doses, and compound categories.
#' 
#' # Key takeaways
#' 
#' 1. Marker-specific features profile marginal improvement in classification 
#' accuracy relative to morphology features.
#' 
#' 2. Marker-specific features alone show a notable drop in classification 
#' accuracy without morphology features, except when considering all markers 
#' (i.e. marker-specific features from all 4 channels).
#' 
#' 3. Trends 1, 2 hold across compound categories and doses, though we observe 
#' higher classification accuracy in certain categories (e.g. HDAC, mTOR).
#'
#' ## Bioactivity filtering
#' Before training compound classifiers, we filter out treatments (i.e. 
#' compounds/dose pairs) that cannot be distinguished from DMSO according to 
#' bioactivity scores (see bioactivity analysis).
#' 
#' **Note:** bioactivity filtering is based on all markers. We set bioactivity 
#' filter as p-value = 0, which results in a larger treatment set than distance 
#' > 2 filtering (see bioactivity notebook).
#' 
#' ## Compound filtering
#' The following analyses evaluate compound classification relative to the 
#' reference compound categories used in ORACL.
#+ load_bioactivity, message=FALSE, echo=FALSE, warning=FALSE
# Load bioactivity table
setwd(analysis.dir)
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
#' For each marker, we train classifiers to predict compound category from
#' phenotypic profiling features.
#' 
#' First, we assess performance relative to a random holdout set. `r prop.train`
#' of wells are randomly sampled to train models and the remaining wells are 
#' used to assess accuracy. **Note:** based on this formulation, a compound can 
#' appear in both the training and test sets but as a different (well) replicate.
#+ random_holdout, fig.height=8, fig.width=15, message=FALSE, echo=FALSE, warning=FALSE
################################################################################
# Random holdout predictions
################################################################################
# Initialize treatments used in training set
cpd.table <- select(xf, Compound_Category, Treatment) %>%
  group_by(Compound_Category) %>%
  summarize(Treatment=list(unique(Treatment)))

trt.train <- lapply(1:reps, function(i) {
  sapply(cpd.table$Treatment, function(z) {
    sample(z, length(z) * prop.train)
  })
})

trt.train <- lapply(trt.train, unlist)

# Fit models for each marker
ypred <- mclapply(
  marker.comb, 
  fit_marker,
  x=xf,
  model=irf,
  model_predict=irf_predict,
  marker.features=features,
  trt.train=trt.train,
  mc.cores=n.core
)

ypred <- rbindlist(ypred) %>%
  mutate(Treatment=str_c(Compound_ID, ', ', Dose_Category))

#+ random_holdout_acc, fig.height=8, fig.width=18, message=FALSE, echo=FALSE, warning=FALSE
# Aggregate model predictions for visualization
xplot.group <- group_by(ypred, Markerset, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=Accuracy) %>%
  mutate(Group=groups[Markerset])

xtext.group <- group_by(xplot.group, Markerset, Group) %>%
  summarize(Accuracy=round(mean(Accuracy), 3))

xplot.group %>%
  ggplot(aes(x=reorder(Markerset, Accuracy), y=Accuracy, col=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.7) +
  geom_text(data=xtext.group, aes(y=1, label=Accuracy), size=2) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by markerset') +
  ylim(c(0, 1.05))

group_by(xplot.group, Rep) %>%
  mutate(Accuracy=Accuracy / Accuracy[Markerset == 'morphology']) %>%
  ggplot(aes(x=reorder(Markerset, Accuracy), col=Group)) +
  geom_boxplot(aes(fill=Group, y=Accuracy), alpha=0.7, outlier.shape=NA) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by markerset', 'relative to morphology') +
  geom_hline(yintercept=1, col='grey', linetype=2)

group_by(xplot.group, Rep) %>%
  mutate(Accuracy=Accuracy / Accuracy[Markerset == 'morphology']) %>%
  ggplot(aes(x=reorder(Group, Accuracy), col=Group)) +
  geom_boxplot(aes(fill=Group, y=Accuracy), alpha=0.7, outlier.shape=NA) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  theme(legend.position='none') +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by markerset group', 'relative to morphology') +
  geom_hline(yintercept=1, col='grey', linetype=2)

#+ random_holdout_acc_dose, fig.height=12, fig.width=18, message=FALSE, echo=FALSE, warning=FALSE
# Accuracy by cell line, dose
xplot <- ypred %>%
  mutate(Dose=str_remove_all(Treatment, '^.*, ')) %>%
  group_by(Markerset, Dose, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Group=groups[Markerset]) %>%
  arrange(Dose, Markerset)

xplot %>%
  ggplot(aes(x=Dose, y=Accuracy, col=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by markerset group, dose') +
  ylim(c(0, 1.05))

group_by(xplot, Rep, Dose) %>%
  mutate(Accuracy=Accuracy / Accuracy[Markerset == 'morphology']) %>%
  ggplot(aes(x=Dose, y=Accuracy, col=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by markerset group, dose', 'relative to morphology') +
  geom_hline(yintercept=1, col='grey', linetype=2)

#+ random_holdout_acc_hm, fig.height=12, fig.width=24, message=FALSE
# Accuracy by compound category
xplot <- group_by(ypred, Markerset, Compound_Category, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Group=groups[Markerset]) %>%
  arrange(Compound_Category, Markerset)

xplot %>%
  ggplot(aes(x=Compound_Category, y=Accuracy, col=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ylim(c(0, 1.05)) +
  ggtitle('Accuracy over hold-out replicates by compound category')

group_by(xplot, Rep, Compound_Category) %>%
  mutate(Accuracy=Accuracy / Accuracy[Markerset == 'morphology']) %>%
  ggplot(aes(x=Compound_Category, y=Accuracy, col=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by compound category', 'relative to morphology') +
  geom_hline(yintercept=1, col='grey', linetype=2)

#' #### Compound holdout
#' Next, we assess performance relative to a randomly held out compounds. That 
#' is, 4 compounds from each category randomly sampled to train models and the 
#' remaining compound from each category is used to assess accuracy. **Note:** 
#' models here are evaluated on compounds they have never seen.
#' 
#' Plots below are the same as above but using hold-out compounds.
#+ compound_holdout, fig.height=8, fig.width=18, message=FALSE, echo=FALSE, warning=FALSE
################################################################################
# Compound holdout predictions
################################################################################
# Intialize compounds for training
cpd.table <- select(xf, Compound_Category, Compound_ID) %>%
  group_by(Compound_Category) %>%
  summarize(Compound_ID=list(unique(Compound_ID)))

cpd.train <- lapply(1:reps, function(i) {
  sapply(cpd.table$Compound_ID, function(z) {
    sample(z, length(z) - 1)
  })
})

cpd.train <- lapply(cpd.train, unlist)

# Fit models for each cell line
ypred <- mclapply(
  marker.comb, 
  fit_marker,
  x=xf,
  model=irf,
  model_predict=irf_predict,
  marker.features=features,
  cpd.train=cpd.train,
  holdout='compound',
  mc.cores=n.core
)

ypred <- rbindlist(ypred) %>%
  mutate(Treatment=str_c(Compound_ID, ', ', Dose_Category))

#+ cpd_holdout_acc, fig.height=8, fig.width=18, message=FALSE, echo=FALSE, warning=FALSE
# Aggregate model predictions for visualization
xplot.group <- group_by(ypred, Markerset, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=Accuracy) %>%
  mutate(Group=groups[Markerset])

xtext.group <- group_by(xplot.group, Markerset, Group) %>%
  summarize(Accuracy=round(mean(Accuracy), 3))

xplot.group %>%
  ggplot(aes(x=reorder(Markerset, Accuracy), y=Accuracy, col=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.7) +
  geom_text(data=xtext.group, aes(y=1, label=Accuracy), size=2) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by markerset') +
  ylim(c(0, 1.05))

group_by(xplot.group, Rep) %>%
  mutate(Accuracy=Accuracy / Accuracy[Markerset == 'morphology']) %>%
  ggplot(aes(x=reorder(Markerset, Accuracy), col=Group)) +
  geom_boxplot(aes(fill=Group, y=Accuracy), alpha=0.7, outlier.shape=NA) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by markerset', 'relative to morphology') +
  geom_hline(yintercept=1, col='grey', linetype=2)

group_by(xplot.group, Rep) %>%
  mutate(Accuracy=Accuracy / Accuracy[Markerset == 'morphology']) %>%
  ggplot(aes(x=reorder(Group, Accuracy), col=Group)) +
  geom_boxplot(aes(fill=Group, y=Accuracy), alpha=0.7, outlier.shape=NA) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  theme(legend.position='none') +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by markerset group', 'relative to morphology') +
  geom_hline(yintercept=1, col='grey', linetype=2)

#+ cpd_holdout_acc_dose, fig.height=12, fig.width=18, message=FALSE, echo=FALSE, warning=FALSE
# Accuracy by cell line, dose
xplot <- ypred %>%
  mutate(Dose=str_remove_all(Treatment, '^.*, ')) %>%
  group_by(Markerset, Dose, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Group=groups[Markerset]) %>%
  arrange(Dose, Markerset)

xplot %>%
  ggplot(aes(x=Dose, y=Accuracy, col=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by markerset group, dose') +
  ylim(c(0, 1.05))

group_by(xplot, Rep, Dose) %>%
  mutate(Accuracy=Accuracy / Accuracy[Markerset == 'morphology']) %>%
  ggplot(aes(x=Dose, y=Accuracy, col=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by markerset group, dose', 'relative to morphology') +
  geom_hline(yintercept=1, col='grey', linetype=2)

#+ cpd_holdout_acc_hm, fig.height=12, fig.width=24, message=FALSE
# Accuracy by compound category
xplot <- group_by(ypred, Markerset, Compound_Category, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Group=groups[Markerset]) %>%
  arrange(Compound_Category, Markerset)

xplot %>%
  ggplot(aes(x=Compound_Category, y=Accuracy, col=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ylim(c(0, 1.05))

group_by(xplot, Rep, Compound_Category) %>%
  mutate(Accuracy=Accuracy / Accuracy[Markerset == 'morphology']) %>%
  ggplot(aes(x=Compound_Category, y=Accuracy, col=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by markerset, dose', 'relative to morphology') +
  geom_hline(yintercept=1, col='grey', linetype=2)

xplot <- group_by(ypred, Markerset, Compound_Category, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Group=groups[Markerset]) %>%
  arrange(Compound_Category, Markerset)

xplot %>%
  ggplot(aes(x=Compound_Category, y=Accuracy, col=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ylim(c(0, 1.05)) +
  ggtitle('Accuracy over hold-out replicates by compound category')

group_by(xplot, Rep, Compound_Category) %>%
  mutate(Accuracy=Accuracy / Accuracy[Markerset == 'morphology']) %>%
  ggplot(aes(x=Compound_Category, y=Accuracy, col=Group)) +
  geom_boxplot(aes(fill=Group), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  scale_color_manual(values=group.pal) +
  ggtitle('Accuracy over hold-out replicates by compound category', 'relative to morphology') +
  geom_hline(yintercept=1, col='grey', linetype=2)