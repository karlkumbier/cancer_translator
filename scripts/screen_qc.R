#+ setup, echo=FALSE, message=FALSE, warning=FALSE
library(data.table)
library(tidyverse)
library(rprofiler)
library(superheat)
library(tidytext)
library(parallel)
library(iRF)
library(caret)
library(ggsci)

################################################################################
# Setup
################################################################################
# Initialize color palettes
col.pal <- RColorBrewer::brewer.pal(11, 'RdYlBu')
col.pal[6] <- '#FFFFFF'

intensity.normalize <- TRUE
n.core <- 6

setwd('~/github/cancer_translator/')
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'
load(str_c(data.dir, 'ks_profiles.Rdata'))

# TODO: check plate loading and clean plateID in data
xks <- select(xks, -PlateID, -WellID)
x <- cbind(xks, xmeta) %>%
  mutate(Control=Compound_ID == 'DMSO') %>%
  distinct() %>%
  filter(!is.na(PlateID)) %>%
  mutate(Col=as.numeric(sapply(str_split(WellID, '-'), tail, 1))) %>%
  mutate_if(is.numeric, median_impute)

#' #QC
#' The figures below report the distribution of cell counts by plate/cell line.
#' We drop plate `2021018028` which appears to have a strange count distribution 
#' relative to other A549 plates.
#+ eda, fig.height=8, fig.width=15
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
  ggplot(aes(x=reorder(Cell_Line, NCells), y=NCells)) +
  geom_boxplot(alpha=0.7, aes(fill=Cell_Line, col=Cell_Line)) +
  theme_bw() +
  scale_fill_nejm() +
  scale_color_nejm() +
  theme(legend.position='none') +
  ggtitle('Cell counts, DMSO wells')

filter(x, Control) %>%
  group_by(Cell_Line, PlateID) %>%
  summarize(NCells=mean(NCells)) %>%
  ggplot(aes(x=reorder(Cell_Line, NCells), y=NCells)) +
  geom_boxplot(alpha=0.7, aes(fill=Cell_Line, col=Cell_Line)) +
  theme_bw() +
  scale_fill_nejm() +
  scale_color_nejm() +  
  theme(legend.position='none') +
  ggtitle('Plate average cell counts, DMSO wells')

#' # Well position effects
#' To assess the effect of well position on feature distribution, we plot
#' feature distribution (KS statistic) for DMSO wells against row/column. We 
#' threshold KS statistics at +/- `r ks.thresh` for visualization.
#+ superheat, fig.height=15, fig.width=21
# Feature distribution plots
ks.thresh <- 0.3

xcontrol <- filter(x, Control)
plate.id <- xcontrol$PlateID
cell.line <- xcontrol$Cell_Line
column <- xcontrol$Col

xcontrol <- dplyr::select(xcontrol, matches('^nonborder'))
colnames(xcontrol) <- str_remove_all(colnames(xcontrol), 'nonborder\\.\\.\\.')

# threshold KS values for visualization
xcontrol[xcontrol < -ks.thresh] <- -ks.thresh
xcontrol[xcontrol > ks.thresh] <- ks.thresh

superheat(xcontrol,
          membership.rows=column,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=col.pal,
          heat.pal.values=seq(0, 1, by=0.1),
          bottom.label.text.angle=90,
          bottom.label.text.size=3,
          bottom.label.size=0.75,
          title='Feature distribution by column')

data.frame(xcontrol) %>%
  select(matches('Intensity')) %>%
  mutate(PlateID=as.factor(plate.id)) %>%
  mutate(Column=as.factor(column)) %>%
  mutate(CellLine=cell.line) %>%
  mutate(Group=str_c(CellLine, ', ', PlateID)) %>%
  reshape2::melt() %>%
  rename(KS=value) %>%
  ggplot(aes(x=variable, y=Column, fill=KS)) +
  geom_tile() +
  facet_wrap(Group~.) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=col.pal)


#' # Normalization
#' To control for well position effects observed above, we normalize all 
#' features by regressing KS value on column ID — i.e. we remove the portion of 
#' a feature that can be explained by a wells column position. Plots below are
#' the same as above but after normalization
#+ normalize, fig.height=15, fig.width=21
# Normalize intensity within each plate
if (intensity.normalize) {
  plates <- unique(x$PlateID)
  x <- lapply(plates, function(p) {
    intensity_normalize(filter(x, PlateID == p))
  })
  
  x <- rbindlist(x)
}

xcontrol <- filter(x, Control)
plate.id <- xcontrol$PlateID
cell.line <- xcontrol$Cell_Line
column <- xcontrol$Col

xcontrol <- dplyr::select(xcontrol, matches('^nonborder'))
colnames(xcontrol) <- str_remove_all(colnames(xcontrol), 'nonborder\\.\\.\\.')

# threshold KS values for visualization
xcontrol[xcontrol < -ks.thresh] <- -ks.thresh
xcontrol[xcontrol > ks.thresh] <- ks.thresh

superheat(xcontrol,
          membership.rows=column,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=col.pal,
          heat.pal.values=seq(0, 1, by=0.1),
          bottom.label.text.angle=90,
          bottom.label.text.size=3,
          bottom.label.size=0.75,
          title='Feature distribution by column')

data.frame(xcontrol) %>%
  select(matches('Intensity')) %>%
  mutate(PlateID=as.factor(plate.id)) %>%
  mutate(Column=as.factor(column)) %>%
  mutate(CellLine=cell.line) %>%
  mutate(Group=str_c(CellLine, ', ', PlateID)) %>%
  reshape2::melt() %>%
  rename(KS=value) %>%
  ggplot(aes(x=variable, y=Column, fill=KS)) +
  geom_tile() +
  facet_wrap(Group~.) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=col.pal)


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
    facet_wrap(~Cell_Line) +
    scale_color_hue(l=50)
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

fout <- str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata')
save(file=fout, x)
