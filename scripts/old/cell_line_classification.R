#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(iRF)
library(caret)
library(superheat)
library(ggsci)

intensity.normalize <- TRUE
n.core <- 16
prop.train <- 0.8
reps <- 100

compound.set <- c(
  'HSP90',
  'MT',
  'mTOR',
  'DNA',
  'Proteasome',
  'HDAC'
)

analysis.dir <- '~/github/cancer_translator/'
setwd(analysis.dir)
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

output.dir <- 'results/cell_line/'
dir.create(output.dir, showWarnings=FALSE)
bioactivity.file <- str_c(output.dir, 'bioactivity.Rdata')
output.file <- str_c(output.dir, 'classification.Rdata') 

intensity.normalize <- TRUE
fin <- str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata')
load(fin)

# Clean compound category names
select_category <- function(x) na.omit(x)[1]

xcat.key <- select(x, Compound_ID, Compound_Category) %>%
  group_by(Compound_ID) %>%
  summarize(Compound_Category=select_category(Compound_Category))

# Initialize color palettes
heat.pal <- c('#FFFFFF', pal_material("light-green")(10))
ncell.pal <- pal_material("purple")(10)
compound.pal <- pal_jco()(10)
cell.pal <- pal_nejm()(8)

#' # Overview
#' This notebook considers the problems of classifying compound category based 
#' on a panel of 6 cell lines: OVCAR4, A549, DU145, ALS-WT, HEPG2, and 786-0. In 
#' addition to models based on individual cell lines, we consider whether 
#' classification can be improved by aggregating information across cell lines.
#'
#' # Key takeaways
#' 
#' 1. Increasing from 1 to 2 cell lines improves compound category 
#' classification accuracy for reference set compounds. Additional cell lines
#' (beyond 2) do not improve accuracy. The best performing pairs include A549 
#' with another cell line (e.g. A549, HEPG2).
#' 
#' 2. Classification accuracy is tied to dose, with lower doses being more 
#' difficult to classify.
#' 
#' 3. The largest improvements in classification accuracy (from 1 to 2 cell 
#' lines) are observed in HDAC compounds.
#' 
#' 
#'
#' ## Bioactivity filtering
#' Before training compound classifiers, we filter out treatments (i.e. 
#' compounds/dose pairs) that cannot be distinguished from DMSO according to 
#' bioactivity scores (see bioactivity analysis).
#' 
#' **Note:** bioactivity filtering is based on all cell lines. All 
#' treatments deemed bioactive wrt at least one cell line appear in our 
#' classification analysis. This allows us to penalize cell line models (wrt 
#' classification performances) that fail to detect bioactive treatments that 
#' other cell line models capture. We set bioactivity filter as p-value = 0, 
#' which results in a larger treatment set than distance > 2 filtering (see 
#' bioactivity notebook).
#' 
#' The following analyses evaluate compound classification relative to the 
#' reference compound categories used in ORACL.
#+ load_bioactivity
# Load bioactivity table
setwd('~/github/cancer_translator/')
load(bioactivity.file)

# Filter to reference compound doses that are bioactive in > 1 cell line
xdist.select <- filter(xdist, Compound_Usage == 'reference_cpd') %>%
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
#' used to assess accuracy. 
#' 
#' **Note:** based on this formulation, a compound can appear in both the 
#' training and test sets but as a different (well) replicate.
#+ random_holdout, message=FALSE, echo=FALSE, warning=FALSE
################################################################################
# Random holdout predictions
################################################################################
# Fit models for each cell line
cell.lines <- unique(xf$Cell_Line)

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

# Fit models
ypred <- mclapply(
  cell.lines, 
  fit_cell_line,
  x=xf,
  model=irf,
  model_predict=irf_predict,
  trt.train=trt.train,
  mc.cores=n.core
)

ypred <- rbindlist(ypred, use.names=TRUE) %>%
  mutate(Treatment=str_c(Compound_ID, ', ', Dose_Category))

#' In addition to predictions based on individual cell lines, we consider an
#' aggregate model that combines information/predictions across all cell 
#' lines by bagging (i.e. average model predictions).
#' 
#' **Note:** We previously considered a full model based on concatenated profiles. 
#' Performance was comparable to bagging and training takes considerably longer
#' for combinatorial cell line sets, so we do not report results for the full
#' model here.
#+ random_holdout_acc, fig.height=12, fig.width=24, message=FALSE, echo=FALSE, warning=FALSE
# Initialize cell line sets
n.cell.line <- length(cell.lines)

cell.sets <- lapply(2:n.cell.line, function(k) {
  combn(cell.lines, k, simplify=FALSE)
})

cell.sets <- unlist(cell.sets, recursive=FALSE)

# Bag predictions over cell line sets
ypred.bag <- mclapply(cell.sets, function(s) {
  out <- bag_predictions(ypred, s) %>%
    mutate(Cell_Line=str_c(s, collapse=', '))
  return(out)
}, mc.cores=n.core)

# Aggregate model predictions for visualization
xplot.bag <- rbindlist(ypred.bag) %>%
  group_by(Cell_Line, Rep) %>%
  summarize(Accuracy=mean(Ytrue == YpredBag), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  select(Cell_Line, Rep, Accuracy)

xplot.cell <- group_by(ypred, Cell_Line, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Ncells=str_count(Cell_Line, ',') + 1)

xplot.text <- group_by(xplot.group, Cell_Line, Ncells) %>%
  summarize(Accuracy=round(mean(Accuracy), 3))

#' The figures below report the distribtuion of prediction accuracy across 
#' hold-out replicates. I.e. each point represents a single replicate for which
#' a randomly sampled collection of treatments are held-out for testing. Since
#' we observe a high degree of variability across hold-out sets, we also report 
#' accuracy relative to A549 classification for each hold-out set.
#+ random_holdout_acc_plots, fig.height=8, fig.width=18, message=FALSE, echo=FALSE, warning=FALSE
xplot.group %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), col=Ncells)) +
  geom_boxplot(aes(fill=Ncells, y=Accuracy), alpha=0.7, outlier.shape=NA) +
  geom_jitter(aes(y=Accuracy), alpha=0.9, width=0.15) +
  geom_text(data=xplot.text, aes(y=1, label=Accuracy), size=2.25, nudge_y=0.02) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  scale_color_gradientn(colors=ncell.pal[-1]) +
  ylim(c(0, 1.05)) +
  theme(legend.position='none') +
  xlab('Number of cell lines used') +
  ggtitle('Accuracy over hold-out replicates by cell line')

group_by(xplot.group, Rep) %>%
  mutate(Accuracy=Accuracy / Accuracy[Cell_Line == 'A549']) %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), col=Ncells)) +
  geom_boxplot(aes(fill=Ncells, y=Accuracy), alpha=0.7, outlier.shape=NA) +
  geom_jitter(aes(y=Accuracy), alpha=0.9, width=0.15) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  scale_color_gradientn(colors=ncell.pal[-1]) +
  theme(legend.position='none') +
  xlab('Number of cell lines used') +
  ylab('Relative accuracy') +
  ggtitle('Accuracy over hold-out replicates by cell line', 'relative to A549') +
  geom_hline(yintercept=1, col='grey', linetype=2)

group_by(xplot.group, Rep) %>%
  mutate(Ncells=as.factor(Ncells)) %>%
  mutate(Accuracy=Accuracy / Accuracy[Cell_Line == 'A549']) %>%
  ggplot(aes(x=Ncells, col=Ncells)) +
  geom_boxplot(aes(fill=Ncells, y=Accuracy), alpha=0.7, outlier.shape=NA) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=ncell.pal[-1]) +
  scale_color_manual(values=ncell.pal[-1]) +
  theme(legend.position='none') +
  xlab('Number of cell lines used') +
  ggtitle('Accuracy over hold-out replicates by cell line', 'relative to A549') +
  geom_hline(yintercept=1, col='grey', linetype=2)

#+ random_holdout_acc_dose, fig.height=8, fig.width=18, message=FALSE
# Accuracy by cell line, dose
xplot.bag <- rbindlist(ypred.bag) %>%
  mutate(Dose=str_remove_all(Treatment, '^.*, ')) %>%
  ungroup() %>%
  group_by(Dose, Cell_Line, Rep) %>%
  summarize(Accuracy=mean(Ytrue == YpredBag), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  rename(Dose_Category=Dose) %>%
  select(Cell_Line, Dose_Category, Rep, Accuracy)

xplot.cell <- group_by(ypred, Cell_Line, Dose_Category, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  select(Cell_Line, Dose_Category, Rep, Accuracy)

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Ncells=str_count(Cell_Line, ',') + 1)

xplot.group %>%
  mutate(Ncells=as.factor(Ncells)) %>%
  ggplot(aes(x=Dose_Category, y=Accuracy)) +
  geom_boxplot(aes(fill=Ncells, col=Ncells), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=ncell.pal[-1]) +
  scale_color_manual(values=ncell.pal[-1]) +
  ylim(c(0, 1.05)) +
  ggtitle('Accuracy by number cells, dose')

group_by(xplot.group, Rep, Dose_Category) %>%
  mutate(Accuracy=Accuracy / Accuracy[Cell_Line == 'A549']) %>%
  mutate(Ncells=as.factor(Ncells)) %>%
  ggplot(aes(x=Dose_Category, y=Accuracy)) +
  geom_boxplot(aes(fill=Ncells, col=Ncells), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=ncell.pal[-1]) +
  scale_color_manual(values=ncell.pal[-1]) +
  ggtitle('Accuracy by number cells, dose', 'relative to A549') +
  geom_hline(yintercept=1, col='grey', linetype=2)


#+ random_holdout_acc_hm, fig.height=8, fig.width=24, message=FALSE
# Accuracy by compound category
xplot.cell <- group_by(ypred, Cell_Line, Compound_Category, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop')

xplot.bag <- rbindlist(ypred.bag) %>%
  group_by(Cell_Line, Compound_Category, Rep) %>%
  summarize(Accuracy=mean(YpredBag == Ytrue), .groups='drop')

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Ncells=str_count(Cell_Line, ',') + 1) %>%
  arrange( Compound_Category, Cell_Line)

xplot.group %>%
  mutate(Ncells=as.factor(Ncells)) %>%
  ggplot(aes(x=Compound_Category, y=Accuracy)) +
  geom_boxplot(aes(fill=Ncells, col=Ncells), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=ncell.pal[-1]) +
  scale_color_manual(values=ncell.pal[-1]) +
  ylim(c(0, 1.05)) +
  ggtitle('Accuracy by number cells, compound category')

group_by(xplot.group, Rep, Compound_Category) %>%
  mutate(Accuracy=Accuracy / Accuracy[Cell_Line == 'A549']) %>%
  mutate(Ncells=as.factor(Ncells)) %>%
  ggplot(aes(x=Compound_Category, y=Accuracy)) +
  geom_boxplot(aes(fill=Ncells, col=Ncells), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=ncell.pal[-1]) +
  scale_color_manual(values=ncell.pal[-1]) +
  ggtitle('Accuracy by number cells, compound category', 'relative to A549') +
  geom_hline(yintercept=1, col='grey', linetype=2)


#' #### Compound holdout
#' Next, we assess performance relative to a randomly held out compounds. That 
#' is, 4 compounds from each category randomly sampled to train models and the 
#' remaining compound from each category is used to assess accuracy. **Note:** 
#' models here are evaluated on compounds they have never seen.
#' 
#' Plots below are the same as above but using hold-out compounds.
#+ compound_holdout, message=FALSE, echo=FALSE, warning=FALSE
################################################################################
# Compound holdout predictions
################################################################################
# Initialize cell lines for modeling
cell.lines <- c(unique(xf$Cell_Line))

# Intialize compounds for training
cpd.table <- select(xf, Compound_Category, Compound_ID) %>%
  group_by(Compound_Category) %>%
  summarize(Compound_ID=list(unique(Compound_ID)))

cpd.train <- lapply(1:reps, function(i) {
  sapply(cpd.table$Compound_ID, function(z) {
    sample(z, length(z) - 1)
  })
})

cpd.train <- lapply(cpd.train, c)

# Fit model
ypred <- mclapply(
  cell.lines, 
  fit_cell_line,
  x=xf,
  model=irf,
  cpd.train=cpd.train,
  model_predict=irf_predict,
  holdout='compound',
  mc.cores=n.core
)

ypred <- rbindlist(ypred, use.names=TRUE) %>%
  mutate(Treatment=str_c(Compound_ID, ', ', Dose_Category))

n.cell.line <- length(cell.lines)

cell.sets <- lapply(2:n.cell.line, function(k) {
  combn(cell.lines, k, simplify=FALSE)
})

cell.sets <- unlist(cell.sets, recursive=FALSE)

# Bag predictions over cell line sets
ypred.bag <- mclapply(cell.sets, function(s) {
  out <- bag_predictions(ypred, s) %>%
    mutate(Cell_Line=str_c(s, collapse=', '))
  return(out)
}, mc.cores=n.core)

# Aggregate model predictions for visualization
xplot.bag <- rbindlist(ypred.bag) %>%
  group_by(Cell_Line, Rep) %>%
  summarize(Accuracy=mean(Ytrue == YpredBag), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  select(Cell_Line, Rep, Accuracy)

xplot.cell <- group_by(ypred, Cell_Line, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Ncells=str_count(Cell_Line, ',') + 1)

xplot.text <- group_by(xplot.group, Cell_Line, Ncells) %>%
  summarize(Accuracy=round(mean(Accuracy), 3))

#' The figures below report the distribtuion of prediction accuracy across 
#' hold-out replicates. I.e. each point represents a single replicate for which
#' andomly sampled compounds are held-out for testing. Since
#' we observe a high degree of variability across hold-out sets, we also report 
#' accuracy relative to A549 classification for each hold-out set.
#+ cpd_holdout_acc_plots, fig.height=12, fig.width=24, message=FALSE, echo=FALSE, warning=FALSE
xplot.group %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), col=Ncells)) +
  geom_boxplot(aes(fill=Ncells, y=Accuracy), alpha=0.7, outlier.shape=NA) +
  geom_jitter(aes(y=Accuracy), alpha=0.9, width=0.15) +
  geom_text(data=xplot.text, aes(y=1, label=Accuracy), size=2.25, nudge_y=0.02) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  scale_color_gradientn(colors=ncell.pal[-1]) +
  ylim(c(0, 1.05)) +
  theme(legend.position='none') +
  xlab('Number of cell lines used') +
  ggtitle('Accuracy over hold-out replicates by cell line')

group_by(xplot.group, Rep) %>%
  mutate(Accuracy=Accuracy / Accuracy[Cell_Line == 'A549']) %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), col=Ncells)) +
  geom_boxplot(aes(fill=Ncells, y=Accuracy), alpha=0.7, outlier.shape=NA) +
  geom_jitter(aes(y=Accuracy), alpha=0.9, width=0.15) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  scale_color_gradientn(colors=ncell.pal[-1]) +
  theme(legend.position='none') +
  xlab('Number of cell lines used') +
  ggtitle('Accuracy over hold-out replicates by cell line', 'relative to A549') +
  geom_hline(yintercept=1, col='grey', linetype=2)


group_by(xplot.group, Rep) %>%
  mutate(Ncells=as.factor(Ncells)) %>%
  mutate(Accuracy=Accuracy / Accuracy[Cell_Line == 'A549']) %>%
  ggplot(aes(x=Ncells, col=Ncells)) +
  geom_boxplot(aes(fill=Ncells, y=Accuracy), alpha=0.7, outlier.shape=NA) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=ncell.pal[-1]) +
  scale_color_manual(values=ncell.pal[-1]) +
  theme(legend.position='none') +
  xlab('Number of cell lines used') +
  ggtitle('Accuracy over hold-out replicates by cell line', 'relative to A549') +
  geom_hline(yintercept=1, col='grey', linetype=2)


#+ cpd_holdout_acc_dose, fig.height=8, fig.width=18, message=FALSE, echo=FALSE, warning=FALSE
# Accuracy by cell line, dose
xplot.bag <- rbindlist(ypred.bag) %>%
  mutate(Dose=str_remove_all(Treatment, '^.*, ')) %>%
  ungroup() %>%
  group_by(Dose, Cell_Line, Rep) %>%
  summarize(Accuracy=mean(Ytrue == YpredBag), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  rename(Dose_Category=Dose) %>%
  select(Cell_Line, Dose_Category, Rep, Accuracy)

xplot.cell <- group_by(ypred, Cell_Line, Dose_Category, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  select(Cell_Line, Dose_Category, Rep, Accuracy)

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Ncells=str_count(Cell_Line, ',') + 1)

xplot.group %>%
  mutate(Ncells=as.factor(Ncells)) %>%
  ggplot(aes(x=Dose_Category, y=Accuracy)) +
  geom_boxplot(aes(fill=Ncells, col=Ncells), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=ncell.pal[-1]) +
  scale_color_manual(values=ncell.pal[-1]) +
  ylim(c(0, 1.05)) +
  ggtitle('Accuracy by number cells, dose')

group_by(xplot.group, Rep, Dose_Category) %>%
  mutate(Accuracy=Accuracy / Accuracy[Cell_Line == 'A549']) %>%
  mutate(Ncells=as.factor(Ncells)) %>%
  ggplot(aes(x=Dose_Category, y=Accuracy)) +
  geom_boxplot(aes(fill=Ncells, col=Ncells), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=ncell.pal[-1]) +
  scale_color_manual(values=ncell.pal[-1]) +
  ggtitle('Accuracy by number cells, dose', 'relative to A549') +
  geom_hline(yintercept=1, col='grey', linetype=2)


#+ cpd_holdout_acc_hm, fig.height=8, fig.width=24, message=FALSE
# Accuracy by compound category
xplot.cell <- group_by(ypred, Cell_Line, Compound_Category, Rep) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop')

xplot.bag <- rbindlist(ypred.bag) %>%
  group_by(Cell_Line, Compound_Category, Rep) %>%
  summarize(Accuracy=mean(YpredBag == Ytrue), .groups='drop')

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Ncells=str_count(Cell_Line, ',') + 1) %>%
  arrange( Compound_Category, Cell_Line)

xplot.group %>%
  mutate(Ncells=as.factor(Ncells)) %>%
  ggplot(aes(x=Compound_Category, y=Accuracy)) +
  geom_boxplot(aes(fill=Ncells, col=Ncells), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=ncell.pal[-1]) +
  scale_color_manual(values=ncell.pal[-1]) +
  ylim(c(0, 1.05)) +
  ggtitle('Accuracy by number cells, compound category')

group_by(xplot.group, Rep, Compound_Category) %>%
  mutate(Accuracy=Accuracy / Accuracy[Cell_Line == 'A549']) %>%
  mutate(Ncells=as.factor(Ncells)) %>%
  ggplot(aes(x=Compound_Category, y=Accuracy)) +
  geom_boxplot(aes(fill=Ncells, col=Ncells), alpha=0.7) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=ncell.pal[-1]) +
  scale_color_manual(values=ncell.pal[-1]) +
  ggtitle('Accuracy by number cells, compound category', 'relative to A549') +
  geom_hline(yintercept=1, col='grey', linetype=2)
