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

setwd('~/github/cancer_translator/')
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

output.dir <- 'results/cell_line/'
dir.create(output.dir, showWarnings=TRUE)
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
heat.pal <- c('#FFFFFF', pal_material("light-blue")(10))
ncell.pal <- pal_material("light-green")(10)
compound.pal <- pal_jco()(10)
cell.pal <- pal_nejm()(8)

#' # Overview
#' This notebook considers the problems of classifying compound category based 
#' on a panel of 6 cell lines: OVCAR4, A549, DU145, ALS-WT, HEPG2, and 786-0. In 
#' addition to models based on individual cell lines, we consider whether 
#' classification can be improved by aggregating information across cell lines.
#'
#' ## Bioactivity filtering
#' Before training compound classifiers, we filter out treatments (i.e. 
#' compounds/dose pairs) that cannot be distinguished from DMSO according to 
#' bioactivity scores (see bioactivity analysis).
#' 
#' **Note:** bioactivity filtering is based on all cell lines. I.e. all 
#' treatments deemed bioactive wrt at least one cell line appear in our 
#' classification analysis. This allows us to penalize cell line models (wrt 
#' classification performances) that fail to detect bioactive treatments that 
#' other cell line models capture.
#' 
#' ## Compound filtering
#' The following analyses evaluate compound classification relative to the 
#' reference compound categories used in ORACL.
#+ load_bioactivity
# Load bioactivity table
load(bioactivity.file)

# Filter to reference compound doses that are bioactive in > 1 cell line
xdist.select <- filter(xdist, Compound_Usage == 'reference_cpd') %>%
  filter(pval == 0) %>%
  mutate(Compound_Dose=str_c(Compound_ID, ', ', Dose_Category))


xf <- mutate(x, Compound_Dose=str_c(Compound_ID, ', ', Dose_Category)) %>%
  filter(Compound_Dose %in% xdist.select$Compound_Dose) %>%
  select(!matches('^PC')) %>%
  filter(Compound_Usage == 'reference_cpd')

# Table of unique compounds by category
compound.table <- group_by(xf, Compound_Category) %>%
  summarize(Ncpd=length(unique(Compound_ID)))

# Filter to ORACL compounds
cpd.keep <- filter(compound.table, Ncpd == 5)
xf <- filter(xf, Compound_Category %in% cpd.keep$Compound_Category)

#' # Modeling
#' For each cell line, we train classifiers to predict compound category from
#' phenotypic profiling features.
#' 
#' First, we assess performance relative to a random holdout set. `r prop.train`
#' of wells are randomly sampled to train models and the remaining wells are 
#' used to assess accuracy. **Note:** based on this formulation, a compound can 
#' appear in both the training and test sets but at different doses.
#+ random_holdout, fig.height=8, fig.width=15, message=FALSE, echo=FALSE
################################################################################
# Random holdout predictions
################################################################################
# Fit models for each cell line
cell.lines <- c(unique(xf$Cell_Line), 'Full')
ypred <- mclapply(
  cell.lines, 
  fit_cell_line,
  x=xf,
  model=irf,
  model_predict=irf_predict,
  prop=prop.train,
  mc.cores=n.core
)

ypred <- rbindlist(ypred, use.names=TRUE)

#' In addition to predictions based on individual cell lines, we consider an
#' 2 aggregate model that combines information/predictions across all cell 
#' lines. First, we aggregate predictions by bagging (i.e. average model 
#' predictions) across cell lines. Second, we consider "full" models trained on 
#' data from all cell lines by concatenating features across cell lines.
#+ random_holdout_acc, fig.height=8, fig.width=15, message=FALSE
# Accuracy by cell line
models <- c(unique(ypred$Cell_Line), 'Bagged')
ypred.bag <- bag_predictions(filter(ypred, Cell_Line != 'Full'))

xplot.bag <- data.frame(Cell_Line='Bagged') %>%
  mutate(Accuracy=mean(ypred.bag$Ytrue == ypred.bag$YpredBag)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.cell <- group_by(ypred, Cell_Line) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=models))

xplot.group %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Cell_Line)) +
  geom_text(aes(x=Cell_Line, y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  theme(legend.position='none') +
  scale_fill_nejm() +
  ylim(c(0, 1.05))

# Accuracy by cell line, dose
xplot.bag <- mutate(ypred.bag, Dose=str_remove_all(Compound_Dose, '^.*, ')) %>%
  ungroup() %>%
  group_by(Dose) %>%
  summarize(Accuracy=mean(Ytrue == YpredBag)) %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  mutate(Cell_Line='Bagged') %>%
  rename(Dose_Category=Dose) %>%
  select(Cell_Line, Dose_Category, Accuracy)

xplot.cell <- group_by(ypred, Cell_Line, Dose_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=models))

xplot.group %>%
  ggplot(aes(x=reorder_within(Cell_Line, Accuracy, Dose_Category), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Cell_Line)) +
  geom_text(aes(x=reorder_within(Cell_Line, Accuracy, Dose_Category), y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  facet_wrap(~Dose_Category, scales='free_x') +
  scale_fill_nejm() +
  ylim(c(0, 1.05)) +
  ggtitle('Classification accuracy by cell line, dose')

xplot.group %>%
  ggplot(aes(x=reorder_within(Dose_Category, Accuracy, Cell_Line), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Cell_Line)) +
  geom_text(aes(x=reorder_within(Dose_Category, Accuracy, Cell_Line), y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  facet_wrap(~Cell_Line, scales='free_x') +
  scale_fill_nejm() +
  ylim(c(0, 1.05)) +
  ggtitle('Classification accuracy by cell line, dose')

# Accuracy by compound category
xplot.cell <- group_by(ypred, Cell_Line, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop')

xplot.bag <- group_by(ypred.bag, Compound_Category) %>%
  summarize(Accuracy=mean(YpredBag == Ytrue), .groups='drop') %>%
  mutate(Cell_Line='Bagged')

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=models))

xplot.group %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy, fill=Compound_Category)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label=Accuracy), position=position_dodge(width=0.9), vjust=-0.02) +
  theme_bw() +
  scale_fill_jco() +
  ylim(0:1)

# Accuracy by compound category heat map
compounds <- unique(xplot.group$Compound_Category)
cell.lines <- unique(xplot.group$Cell_Line)
xplot.group.hm <- matrix(xplot.group$Accuracy, nrow=length(compounds))

rownames(xplot.group.hm) <- compounds
colnames(xplot.group.hm) <- cell.lines
superheat(
  xplot.group.hm,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  heat.pal=viridis::inferno(10),
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label.text.angle=90
)

#' To assess model performance wrt sampling variability, we compute the average
#' accuracy for each compound/dose/cell line category. 
#+ pred_variaiblity_rand, fig.height=8, fig.width=12
xplot.cell <- mutate(ypred, Group=str_c(Compound_ID, ', ', Dose_Category)) %>%
  group_by(Cell_Line, Group, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.bag <- rename(ypred.bag, Group=Compound_Dose) %>%
  mutate(Cell_Line='Bagged') %>%
  select(Cell_Line, Group, Compound_Category, Accuracy)

xplot.group <- rbind(xplot.cell, xplot.bag) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=models))

xplot.group %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy, fill=Cell_Line)) +
  geom_boxplot(aes(fill=Cell_Line, col=Cell_Line), alpha=0.7, outlier.shape=NA) +
  geom_jitter(aes(fill=Cell_Line, col=Cell_Line), height=0, width=0.15) +
  theme_bw() +
  scale_fill_nejm() +
  scale_color_nejm() +
  theme(legend.position='none') +
  ylim(0:1)

xplot.group %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  mutate(Dose=str_remove_all(Group, '^.*, ')) %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy, fill=Cell_Line)) +
  geom_boxplot(aes(fill=Cell_Line, col=Cell_Line), alpha=0.7, outlier.shape=NA) +
  geom_jitter(aes(fill=Cell_Line, col=Cell_Line), height=0, width=0.15) +
  theme_bw() +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~Dose) +
  theme(legend.position='none') +
  ylim(0:1) +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Classification accuracy for compound/dose pairs by cell line, dose category')

ypred.random <- ypred
ypred.bag.random <- ypred.bag

#' #### Random holdout
#' Next, we assess performance relative to a randomly held out compounds. That 
#' is, 4 compounds from each category randomly sampled to train models and the 
#' remaining compound from each category is used to assess accuracy. **Note:** 
#' models here are evaluated on compounds they have never seen.
#' 
#' Plots below are the same as above but using hold-out compounds.
#+ compound_holdout, fig.height=8, fig.width=15, message=FALSE, echo=FALSE
################################################################################
# Compound holdout predictions
################################################################################
# Fit models for each cell line
cell.lines <- c(unique(xf$Cell_Line), 'Full')
ypred <- mclapply(
  cell.lines, 
  fit_cell_line,
  x=xf,
  model=irf,
  model_predict=irf_predict,
  holdout='compound',
  mc.cores=n.core
)

ypred <- rbindlist(ypred, use.names=TRUE)

# Accuracy by cell line
ypred.bag <- bag_predictions(filter(ypred, Cell_Line != 'Full'))
xplot.bag <- data.frame(Cell_Line='Bagged') %>%
  mutate(Accuracy=mean(ypred.bag$Ytrue == ypred.bag$YpredBag)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.cell <- group_by(ypred, Cell_Line) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=models))
  
xplot.group %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Cell_Line)) +
  geom_text(aes(x=Cell_Line, y=Accuracy + 0.02, label=Accuracy)) +
  scale_fill_manual(values=cell.pal) +
  theme_bw() +
  theme(legend.position='none') +
  ylim(c(0, 1))

# Accuracy by cell line, dose
xplot.bag <- mutate(ypred.bag, Dose=str_remove_all(Compound_Dose, '^.*, ')) %>%
  ungroup() %>%
  group_by(Dose) %>%
  summarize(Accuracy=mean(Ytrue == YpredBag)) %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  mutate(Cell_Line='Bagged') %>%
  rename(Dose_Category=Dose) %>%
  select(Cell_Line, Dose_Category, Accuracy)

xplot.cell <- group_by(ypred, Cell_Line, Dose_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=models))

xplot.group %>%
  ggplot(aes(x=reorder_within(Cell_Line, Accuracy, Dose_Category), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Cell_Line)) +
  geom_text(aes(x=reorder_within(Cell_Line, Accuracy, Dose_Category), y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  facet_wrap(~Dose_Category, scales='free_x') +
  scale_fill_manual(values=cell.pal) +
  ylim(c(0, 1.05)) +
  ggtitle('Classification accuracy by cell line, dose')

xplot.group %>%
  ggplot(aes(x=reorder_within(Dose_Category, Accuracy, Cell_Line), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Cell_Line)) +
  geom_text(aes(x=reorder_within(Dose_Category, Accuracy, Cell_Line), y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  facet_wrap(~Cell_Line, scales='free_x') +
  scale_fill_manual(values=cell.pal) +
  ylim(c(0, 1.05)) +
  ggtitle('Classification accuracy by cell line, dose')

# Accuracy by compound category
xplot.cell <- group_by(ypred, Cell_Line, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue))

xplot.bag <- group_by(ypred.bag, Compound_Category) %>%
  summarize(Accuracy=mean(YpredBag == Ytrue), .groups='drop') %>%
  mutate(Cell_Line='Bagged')

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=models))

xplot.group %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=unique(xplot.group$Cell_Line))) %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy, fill=Compound_Category)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label=Accuracy), position=position_dodge(width=0.9), vjust=-0.02) +
  theme_bw() +
  scale_fill_manual(values=compound.pal) +
  ylim(0:1)

# Accuracy by compound category heat map
compounds <- unique(xplot.group$Compound_Category)
cell.lines <- unique(xplot.group$Cell_Line)
xplot.group.hm <- matrix(xplot.group$Accuracy, nrow=length(compounds))

rownames(xplot.group.hm) <- compounds
colnames(xplot.group.hm) <- cell.lines
xplot.group.hm[xplot.group.hm < 0.5] <- 0.5

superheat(
  xplot.group.hm,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  heat.pal=viridis::inferno(10),
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label.text.angle=90
)

#' To assess model performance wrt sampling variability, we compute the average
#' accuracy for each compound/dose/cell line category. 
#+ pred_variaiblity_cpd, fig.height=8, fig.width=12
xplot.cell <- mutate(ypred, Group=str_c(Compound_ID, ', ', Dose_Category)) %>%
  group_by(Cell_Line, Group, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.bag <- rename(ypred.bag, Group=Compound_Dose) %>%
  mutate(Cell_Line='Bagged') %>%
  select(Cell_Line, Group, Compound_Category, Accuracy)

xplot.group <- rbind(xplot.cell, xplot.bag) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=models))

xplot.group %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy, fill=Cell_Line)) +
  geom_boxplot(aes(fill=Cell_Line, col=Cell_Line), alpha=0.7, outlier.shape=NA) +
  geom_jitter(aes(fill=Cell_Line, col=Cell_Line), height=0, width=0.15) +
  theme_bw() +
  scale_fill_nejm() +
  scale_color_nejm() +
  theme(legend.position='none') +
  ylim(0:1)

xplot.group %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  mutate(Dose=str_remove_all(Group, '^.*, ')) %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy, fill=Cell_Line)) +
  geom_boxplot(aes(fill=Cell_Line, col=Cell_Line), alpha=0.7, outlier.shape=NA) +
  geom_jitter(aes(fill=Cell_Line, col=Cell_Line), height=0, width=0.15) +
  theme_bw() +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~Dose) +
  ylim(0:1) +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Classification accuracy for compound/dose pairs by cell line, dose category')

save(file=output.file, xdist, bioactive.cat, ypred.random, ypred.bag.random)