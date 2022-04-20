#+ setup, echo=FALSE, message=FALSE, warning=FALSE
library(data.table)
library(tidyverse)
library(superheat)
library(tidytext)
library(parallel)
library(iRF)
library(caret)
library(ggsci)
library(hrbrthemes)

theme_set(theme_ipsum(base_family='sans'))


################################################################################
# Setup
################################################################################
# Initialize color palettes
col.pal <- RColorBrewer::brewer.pal(11, 'RdYlBu')
col.pal[6] <- '#FFFFFF'

normalize <- TRUE
n.core <- 12
ks.thresh <- 0.3

setwd('~/github/cancer_translator/')
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'
load(str_c(data.dir, 'ks_profiles.Rdata'))

xks <- mutate(xks, ID=str_c(PlateID, '_', WellID)) %>%
  select(-PlateID, -WellID) %>%
  mutate_if(is.numeric, median_impute)

x <- mutate(xmeta, ID=str_c(PlateID, '_', WellID)) %>%
  mutate(Control=Compound_ID == 'DMSO') %>%
  distinct() %>%
  filter(!is.na(PlateID)) %>%
  mutate(Col=as.numeric(sapply(str_split(WellID, '-'), tail, 1))) %>%
  left_join(xks, by='ID')

n.cell.line <- length(unique(x$Cell_Line))

#' # QC
#' The figures below report the distribution of cell counts by plate/cell line.
#+ eda, fig.height=8, fig.width=15
################################################################################
# Simple EDA - counts and intensity
################################################################################
# Cell count plots
filter(x, Control) %>%
  ggplot(aes(x=NCells)) +
  geom_histogram() +
  ggtitle('Cell counts, DMSO wells')

filter(x, Control) %>%
  ggplot(aes(x=reorder(Cell_Line, NCells), y=NCells)) +
  geom_boxplot(alpha=0.7, aes(fill=Cell_Line, col=Cell_Line)) +
  scale_fill_nejm() +
  scale_color_nejm() +
  theme(legend.position='none') +
  ggtitle('Cell counts, DMSO wells')

filter(x, Control) %>%
  ggplot(aes(x=reorder(PlateID, NCells), y=NCells)) +
  geom_boxplot(alpha=0.7, aes(fill=Cell_Line, col=Cell_Line)) +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~Cell_Line, scales='free') +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Cell counts, DMSO wells')


#' # Well position effects
#' To assess the effect of well position on feature distribution, we plot
#' feature distribution (KS statistic) for DMSO wells against row/column. We 
#' threshold KS statistics at +/- `r ks.thresh` for visualization.
#+ superheat, fig.height=20, fig.width=30
# Feature distribution plots
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
  facet_wrap(Group~., nrow=n.cell.line) +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=col.pal)


#' # Normalization
#' To control for well position effects observed above, we normalize all 
#' features by regressing KS value on row, column, and plate ID — i.e. we remove 
#' the portion of a feature that can be explained by plate + plate position. 
#' Plots below are the same as above but after normalization
#+ normalize, fig.height=15, fig.width=21
# Save original data for comparison
xfull <- x

# Normalize intensity within each plate
if (normalize) {
  x <- mclapply(unique(x$PlateID), function(p) {
    normalize_features(filter(x, PlateID == p))
  }, mc.cores=n.core)
  
  x <- rbindlist(x)
}

f <- 'nonborder...Intensity.Nucleus.HOECHST.33342.Mean'
filter(xfull, Compound_ID == 'DMSO') %>%
  mutate(ColID=as.factor(ColID)) %>%
  ggplot(aes_string(x='PlateID', y=f)) +
  geom_boxplot(aes(fill=ColID, col=ColID), alpha=0.7) +
  facet_wrap(~Cell_Line, scales='free_x') +
  theme(axis.text.x=element_blank()) +
  scale_fill_nejm() +
  scale_color_nejm()

filter(x, Compound_ID == 'DMSO') %>%
  mutate(ColID=as.factor(ColID)) %>%
  ggplot(aes_string(x='PlateID', y=f)) +
  geom_boxplot(aes(fill=ColID, col=ColID), alpha=0.7) +
  facet_wrap(~Cell_Line, scales='free_x') +
  theme(axis.text.x=element_blank()) +
  scale_fill_nejm() +
  scale_color_nejm()



#' # PCA analysis
#' To qualitatively assess the feature distribution of different compound 
#' classes we plot PCA projections grouped by various metadata features. Broadly
#' speaking, the primary sources of variation in the data correspond to compound
#' categories.
#+ pca, fig.height=6, fig.width=12
################################################################################
# Full dataset PCA
################################################################################
setwd('~/github/cancer_translator/')

xc <- filter(x, Compound_ID == 'DMSO')
xc.full <- filter(xfull, Compound_ID == 'DMSO')

# Compute PCA on normalized data
xfeat <- dplyr::select(xc, matches('nonborder')) 
xpca <- prcomp(xfeat)
xc <- cbind(xpca$x, xc)

data.frame(PctVar=cumsum(xpca$sdev ^ 2) / sum(xpca$sdev ^ 2)) %>%
  mutate(NPC=1:n()) %>%
  ggplot(aes(x=NPC, y=PctVar)) +
  geom_point() +
  geom_line()

# Compute PCA on unnormalized data
xfeat <- dplyr::select(xc.full, matches('nonborder'))
xpca <- prcomp(xfeat)
xc.full <- cbind(xpca$x, xc.full)

pc.list <- list(c('PC1', 'PC2'), c('PC3', 'PC4'))
for (pcs in pc.list) {
  
  p <- group_by(xc, Cell_Line) %>%
    sample_frac(0.5) %>%
    mutate(ColID=as.factor(ColID)) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=ColID), alpha=0.5) +
    ggtitle('Plate ID', 'Normalized') +
    facet_wrap(~Cell_Line) +
    scale_color_hue(l=50) +
    theme(legend.position='none')
  plot(p)
  
  p <- group_by(xc.full, Cell_Line) %>%
    sample_frac(0.5) %>%
    mutate(ColID=as.factor(ColID)) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=ColID), alpha=0.5) +
    ggtitle('Plate ID', 'Unnormalized') +
    facet_wrap(~Cell_Line) +
    scale_color_hue(l=50) +
    theme(legend.position='none')
  plot(p)
  
  p <- group_by(xc, Cell_Line) %>%
    mutate(PlateID=as.numeric(as.factor(PlateID))) %>%
    mutate(PlateID=as.factor(PlateID)) %>%
    sample_frac(0.5) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=PlateID), alpha=0.5) +
    ggtitle('Plate ID', 'Normalized') +
    facet_wrap(~Cell_Line) +
    scale_color_hue(l=60) +
    theme(legend.position='none')
  plot(p)
  
  p <- group_by(xc.full, Cell_Line) %>%
    mutate(PlateID=as.numeric(as.factor(PlateID))) %>%
    mutate(PlateID=as.factor(PlateID)) %>%
    sample_frac(0.5) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=PlateID), alpha=0.5) +
    ggtitle('Plate ID', 'Unnormalized') +
    facet_wrap(~Cell_Line) +
    scale_color_hue(l=60) +
    theme(legend.position='none')
  plot(p)
}

normalize <- ifelse(normalize, 'normalized', 'unnormalized')
fout <- str_c(data.dir, 'profiles_', normalize, '.Rdata')
save(file=fout, x)
