#+ setup, echo=FALSE
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
n.core <- 6

setwd('~/github/cancer_translator/')
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

intensity.normalize <- TRUE
fin <- str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata')
load(fin)

# TODO: clean compound category names
select_category <- function(x) na.omit(x)[1]

xcat.key <- select(x, Compound_ID, Compound_Category) %>%
  group_by(Compound_ID) %>%
  summarize(Compound_Category=select_category(Compound_Category))

#' # Bioactivity filtering
#' Before training compound classifiers, we filter out compounds/doses that
#' cannot be distinguished from DMSO. Specifically, we (i) take a subsample
#' of DMSO wells (ii) compute the center of the DMSO point cloud (i.e.
#' average feature values across each well in the subsample) (iii) compute the
#' l2 distance between each DMSO well in the subsample and the DMSO point cloud 
#' center (iv) define the DMSO maximal distance as the maximum distance (from 
#' iii) over all subsampled wells.
#' 
#' We then ask whether a given well compound/dose is further from the DMSO
#' point cloud center than the maximal DMSO distance. Repeating this process
#' across many subsamples allows us to generate a bioactivity p-value. Tables 
#' below summarize bioactivity by plate, compound category (in reference set),
#' and cell line.
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

# Group bioactivity scores by cell line/compound
n.cell.line <- length(unique(x$Cell_Line))
xgroup <- group_by(xdist, Cell_Line, Compound_ID, Dose_Category) %>% 
  summarize(DistNorm=mean(DistNorm), .groups='drop') %>%
  group_by(Compound_ID) %>%
  mutate(Count=n()) %>%
  filter(Count == n.cell.line) %>%
  ungroup() %>%
  arrange(Compound_ID, Cell_Line) %>%
  left_join(xcat.key, by='Compound_ID')

# Format cell x compound bioactivity for visualization
xplot <- matrix(xgroup$DistNorm, nrow=n.cell.line)
rownames(xplot) <- unique(xgroup$Cell_Line)
colnames(xplot) <- unique(xgroup$Compound_ID)
category <- matrix(xgroup$Compound_Category, nrow=n.cell.line)[1,]

# Plot bioactivity by cell line, compound
xplot.t <- xplot
xplot.t[xplot.t < 1] <- 0
superheat(log(xplot.t + 1, base=2),
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=viridis::inferno(10), 
          heat.pal.values=seq(0, 1, by=0.1),
          title='Bioactivity by cell line, compound (log scale)')

# Plot bioactivity by cell line, coompound category
cat.count <- table(category)
cat.keep <- setdiff(names(cat.count[cat.count > 20]), 'Others')
superheat(log(xplot.t + 1, base=2)[,category %in% cat.keep],
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          membership.cols=category[category %in% cat.keep],
          bottom.label.text.angle=90,
          bottom.label.size=0.75,
          heat.pal=viridis::inferno(10), 
          heat.pal.values=seq(0, 1, by=0.1),
          title='Bioactivity by cell line, compound (log scale)\nprevalent categories')

# Generate table of bioactivity by cell line/category
bioactive.cell.cat <- filter(xdist, Compound_Usage == 'reference_cpd') %>%
  group_by(Compound_Category, Cell_Line) %>% 
  summarize(PropBioactive=mean(pval == 0))

bioactive.cell.cat %>%
  ggplot(aes(x=reorder_within(Compound_Category, PropBioactive, Cell_Line), y=PropBioactive)) +
  geom_bar(stat='identity', aes(fill=Compound_Category)) +
  theme_bw() +
  facet_wrap(~Cell_Line, scales='free_x') +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_hue(l=60) +
  ggtitle('Bioactivity by category/cell line')
  
write.csv(file='results/bioactivity_cell_category.csv', bioactive.cell.cat, quote=FALSE)


# Generate table of bioactivity by plate
bioactive.plate <- group_by(xdist, PlateID) %>% 
  summarize(PropBioactive=mean(pval == 0))

bioactive.plate %>%
  ggplot(aes(x=reorder(PlateID, PropBioactive), y=PropBioactive)) +
  geom_bar(stat='identity') +
  theme_bw() +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_hue(l=60) +
  ggtitle('Bioactivity by plate')

write.csv(file='results/bioactivity_plate.csv', bioactive.plate, quote=FALSE)

# Filter to reference compound doses that are bioactive in > 1 cell line
xdist.select <- filter(xdist, Compound_Usage == 'reference_cpd') %>%
  filter(pval == 0) %>%
  mutate(Compound_Dose=str_c(Compound_ID, ', ', Dose_Category))

#' # Modeling
#' For each cell line, we train classifiers to predict compound category from
#' phenotypic profiling features. Compounds/doses are filtered to include only 
#' those that are bioactive in at least one cell line.
#+ modeling, fig.height=8, fig.width=15, message=FALSE, warnings=FALSE, echo=FALSE
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
#+ random_holdout, fig.height=8, fig.width=15, message=FALSE, echo=FALSE
################################################################################
# Random holdout predictions
################################################################################
# Fit models for each cell line
ypred <- lapply(unique(xf$Cell_Line), 
                fit_cell_line,
                x=xf,
                model=irf,
                model_predict=irf_predict,
                prop=0.8
)

ypred <- rbindlist(ypred)
ypred.bag <- bag_predictions(ypred)

# Accuracy by cell line
xplot.bag <- data.frame(Cell_Line='Aggregate') %>%
  mutate(Accuracy=mean(ypred.bag$Ytrue == ypred.bag$YpredBag)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.cell <- group_by(ypred, Cell_Line) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.group <- rbind(xplot.bag, xplot.cell)

xplot.group %>%
  mutate(Cell_Line=factor(Cell_Line, levels=xplot.group$Cell_Line)) %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', fill='#0088D1') +
  geom_text(aes(x=Cell_Line, y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  ylim(c(0, 1))


# Accuracy by compound category
xplot.cell <- group_by(ypred, Cell_Line, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue))

xplot.bag <- group_by(ypred.bag, Compound_Category) %>%
  summarize(Accuracy=mean(YpredBag == Ytrue)) %>%
  mutate(Cell_Line='Aggregate')

xplot.group <- rbind(xplot.bag, xplot.cell)

xplot.group %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=unique(xplot.group$Cell_Line))) %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy, fill=Compound_Category)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label=Accuracy), position=position_dodge(width=0.9), vjust=-0.02) +
  theme_bw() +
  ylim(0:1)

#' Next, we assess performance relative to a randomly held out compounds. That 
#' is, 4 compounds from each category randomly sampled to train models and the 
#' remaining compound from each category is used to assess accuracy. Note: 
#' models here are evaluated on compounds they have never seen
#+ compound_holdout, fig.height=8, fig.width=15, message=FALSE, echo=FALSE
################################################################################
# Compound holdout predictions
################################################################################
# Fit models for each cell line
ypred <- lapply(unique(xf$Cell_Line), fit_cell_line,
                x=xf,
                model=irf,
                model_predict=irf_predict,
                holdout='compound'
)

ypred <- rbindlist(ypred)
ypred.bag <- bag_predictions(ypred)

# Accuracy by cell line
xplot.bag <- data.frame(Cell_Line='Aggregate') %>%
  mutate(Accuracy=mean(ypred.bag$Ytrue == ypred.bag$YpredBag)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.cell <- group_by(ypred, Cell_Line) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue)) %>%
  mutate(Accuracy=round(Accuracy, 2))

xplot.group <- rbind(xplot.bag, xplot.cell)

rbind(xplot.bag, xplot.cell) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=unique(xplot.group$Cell_Line))) %>%
  ggplot(aes(x=reorder(Cell_Line, Accuracy), y=Accuracy)) +
  geom_bar(stat='identity', fill='#0088D1') +
  geom_text(aes(x=Cell_Line, y=Accuracy + 0.02, label=Accuracy)) +
  theme_bw() +
  ylim(c(0, 1))

# Accuracy by compound category
xplot.cell <- group_by(ypred, Cell_Line, Compound_Category) %>%
  summarize(Accuracy=mean(YpredCl == Ytrue))

xplot.bag <- group_by(ypred.bag, Compound_Category) %>%
  summarize(Accuracy=mean(YpredBag == Ytrue)) %>%
  mutate(Cell_Line='Aggregate')

xplot.group <- rbind(xplot.bag, xplot.cell)

rbind(xplot.cell, xplot.bag) %>%
  mutate(Accuracy=round(Accuracy, 2)) %>%
  mutate(Cell_Line=factor(Cell_Line, levels=unique(xplot.group$Cell_Line))) %>%
  ggplot(aes(x=Cell_Line, y=Accuracy, fill=Compound_Category)) +
  geom_bar(stat='identity', position='dodge') +
  geom_text(aes(label=Accuracy), position=position_dodge(width=0.9), vjust=-0.02) +
  theme_bw() +
  ylim(0:1)
