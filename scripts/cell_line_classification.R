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

compound.set <- c(
  'HSP90',
  'MT',
  'mTOR',
  'DNA',
  'Proteasome',
  'HDAC'
)

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
#' used to assess accuracy. **Note:** based on this formulation, a compound can 
#' appear in both the training and test sets but as a different (well) replicate.
#+ random_holdout, fig.height=8, fig.width=15, message=FALSE, echo=FALSE
################################################################################
# Random holdout predictions
################################################################################
# Fit models for each cell line
cell.lines <- unique(xf$Cell_Line)

ypred <- mclapply(
  cell.lines, 
  fit_cell_line,
  x=xf,
  model=irf,
  model_predict=irf_predict,
  prop=prop.train,
  mc.cores=n.core
)

ypred <- rbindlist(ypred, use.names=TRUE) %>%
  mutate(Treatment=str_c(Compound_ID, ', ', Dose_Category))

#' In addition to predictions based on individual cell lines, we consider an
#' aggregate model that combines information/predictions across all cell 
#' lines by bagging (i.e. average model predictions).
#' 
#' **Note:** We also considered a full model based on concatenated profiles. 
#' Performance was comparable to bagging and training takes considerably longer
#' for combinatorial cell line sets, so we do not report results for the full
#' model here.
#+ random_holdout_acc, fig.height=8, fig.width=15, message=FALSE
# Accuracy by cell line
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
  group_by(Cell_Line) %>%
  summarize(Accuracy=mean(Ytrue == YpredBag), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  select(Cell_Line, Accuracy)

xplot.cell <- group_by(ypred, Cell_Line) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Ncells=str_count(Cell_Line, ',') + 1)

xplot.group %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Ncells)) +
  geom_text(aes(x=Cell_Line, y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  ylim(c(0, 1.05))

# Accuracy by cell line, dose
xplot.bag <- rbindlist(ypred.bag) %>%
  mutate(Dose=str_remove_all(Treatment, '^.*, ')) %>%
  ungroup() %>%
  group_by(Dose, Cell_Line) %>%
  summarize(Accuracy=mean(Ytrue == YpredBag), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  rename(Dose_Category=Dose) %>%
  select(Cell_Line, Dose_Category, Accuracy)

xplot.cell <- group_by(ypred, Cell_Line, Dose_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  select(Cell_Line, Dose_Category, Accuracy)

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Ncells=str_count(Cell_Line, ',') + 1)

xplot.group %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Ncells)) +
  geom_text(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  facet_grid(Dose_Category~.) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  ylim(c(0, 1.05)) +
  ggtitle('Classification accuracy by cell line, dose')

# Accuracy by compound category
xplot.cell <- group_by(ypred, Cell_Line, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop')

xplot.bag <- rbindlist(ypred.bag) %>%
  group_by(Cell_Line, Compound_Category) %>%
  summarize(Accuracy=mean(YpredBag == Ytrue), .groups='drop')

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Ncells=str_count(Cell_Line, ',') + 1) %>%
  arrange( Compound_Category, Cell_Line)

# Accuracy by compound category heat map
compounds <- unique(xplot.group$Compound_Category)
cell.lines <- unique(xplot.group$Cell_Line)
xplot.group.hm <- matrix(xplot.group$Accuracy, nrow=length(cell.lines))

colnames(xplot.group.hm) <- compounds
rownames(xplot.group.hm) <- cell.lines

# Initialize row attributes
ncell <- select(xplot.group, Cell_Line, Ncells) %>% distinct()
ncell.col <- ncell.pal[-1][floor(ncell$Ncells * 1.5)]

superheat(
  xplot.group.hm,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  membership.rows=ncell$Ncells,
  left.label.col=ncell.col,
  left.label='variable',
  left.label.text.size=3,
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label.text.angle=90
)

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
cell.lines <- c(unique(xf$Cell_Line))

ypred <- mclapply(
  cell.lines, 
  fit_cell_line,
  x=xf,
  model=irf,
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
  group_by(Cell_Line) %>%
  summarize(Accuracy=mean(Ytrue == YpredBag), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  select(Cell_Line, Accuracy)

xplot.cell <- group_by(ypred, Cell_Line) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Ncells=str_count(Cell_Line, ',') + 1)

xplot.group %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Ncells)) +
  geom_text(aes(x=Cell_Line, y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  ylim(c(0, 1.05))

# Accuracy by cell line, dose
xplot.bag <- rbindlist(ypred.bag) %>%
  mutate(Dose=str_remove_all(Treatment, '^.*, ')) %>%
  ungroup() %>%
  group_by(Dose, Cell_Line) %>%
  summarize(Accuracy=mean(Ytrue == YpredBag), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  rename(Dose_Category=Dose) %>%
  select(Cell_Line, Dose_Category, Accuracy)

xplot.cell <- group_by(ypred, Cell_Line, Dose_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop') %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  select(Cell_Line, Dose_Category, Accuracy)

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Ncells=str_count(Cell_Line, ',') + 1)

xplot.group %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', aes(fill=Ncells)) +
  geom_text(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  facet_grid(Dose_Category~.) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  ylim(c(0, 1.05)) +
  ggtitle('Classification accuracy by cell line, dose')

# Accuracy by compound category
xplot.cell <- group_by(ypred, Cell_Line, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue), .groups='drop')

xplot.bag <- rbindlist(ypred.bag) %>%
  group_by(Cell_Line, Compound_Category) %>%
  summarize(Accuracy=mean(YpredBag == Ytrue), .groups='drop')

xplot.group <- rbind(xplot.bag, xplot.cell) %>%
  mutate(Ncells=str_count(Cell_Line, ',') + 1) %>%
  arrange( Compound_Category, Cell_Line)

# Accuracy by compound category heat map
compounds <- unique(xplot.group$Compound_Category)
cell.lines <- unique(xplot.group$Cell_Line)
xplot.group.hm <- matrix(xplot.group$Accuracy, nrow=length(cell.lines))

colnames(xplot.group.hm) <- compounds
rownames(xplot.group.hm) <- cell.lines

# Initialize row attributes
ncell <- select(xplot.group, Cell_Line, Ncells) %>% distinct()
ncell.col <- ncell.pal[-1][floor(ncell$Ncells * 1.5)]

superheat(
  xplot.group.hm,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  membership.rows=ncell$Ncells,
  left.label.col=ncell.col,
  left.label='variable',
  left.label.text.size=3,
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label.text.angle=90
)