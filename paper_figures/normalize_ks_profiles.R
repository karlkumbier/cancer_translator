#' This script runs plate position normalization on raw KS profiles. For each
#' feature and plate, we fit a normalization model on DMSO-treated wells by 
#' regressing the feature value (KS profile) on row and column position. The
#' resulting normalization model then run over all wells in the plate, and model
#' predicted feature values are subtracted from observed feature values—i.e.
#' we remove the portion of the feature value that can be explained by plate 
#' position.
#' 
#' Last updated 10/5/2022, Karl Kumbier
library(data.table)
library(tidyverse)
library(superheat)
library(tidytext)
library(parallel)
library(iRF)
library(caret)
library(ggsci)
library(hrbrthemes)

theme_set(
  theme_ipsum(
    plot_title_size=28,
    axis_title_size=24, 
    strip_text_size=28, 
    axis_text_size=18,
    base_size=18, 
    base_family='sans'
  )
)

set.seed(47)

# Initialize color palettes
col.pal <- RColorBrewer::brewer.pal(11, 'RdYlBu')
col.pal[6] <- '#FFFFFF'

# Initialize parameters for analysis
normalize <- TRUE
n.core <- 12
ks.thresh <- 0.3
save.fig <- TRUE
normalize.s <- ifelse(normalize, 'normalized', 'unnormalized')

# Initialize analysis directories, load KS profiles
ccl.dir <- Sys.getenv('CCL_DIR')
data.dir <- str_c(ccl.dir, 'data/screens/LH_CDC_1/')
fig.dir <- str_c(ccl.dir, 'paper_figures/Supplemental/figsupp/')
output.file <- str_c(data.dir, 'profiles_', normalize.s, '.Rdata')
dir.create(fig.dir, recursive=TRUE, showWarnings=FALSE)

load(str_c(data.dir, 'ks_profiles.Rdata'))

source(str_c(ccl.dir, 'paper_figures/utilities.R'))

# Median impute missing values - 2 wells no nearest neighbors
xks <- mutate(xks, ID=str_c(PlateID, '_', WellID)) %>%
  select(-PlateID, -WellID) %>%
  mutate_if(is.numeric, median_impute)

# Merge KS profiles with metadata, initialize DMSO as normalization controls
x <- mutate(xmeta, ID=str_c(PlateID, '_', WellID)) %>%
  mutate(Control=Compound_ID == 'DMSO') %>%
  distinct() %>%
  filter(!is.na(PlateID)) %>%
  mutate(Col=as.numeric(sapply(str_split(WellID, '-'), tail, 1))) %>%
  left_join(xks, by='ID')

n.cell.line <- length(unique(x$Cell_Line))

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
row <- xcontrol$Row

xcontrol <- dplyr::select(xcontrol, matches('^nonborder'))
colnames(xcontrol) <- str_remove_all(colnames(xcontrol), 'nonborder\\.\\.\\.')

# threshold KS values for visualization
xcontrol[xcontrol < -ks.thresh] <- -ks.thresh
xcontrol[xcontrol > ks.thresh] <- ks.thresh

# score positional effect for column ordering
score_column <- function(z) {
  mean(z[column == max(column)]) - mean(z[column == min(column)])
}

score <- apply(xcontrol, MAR=2, score_column)
column <- str_c('Column ', column)

fout <- str_c(fig.dir, 'fig_supp_col_unnormalized_full.png')
if (save.fig) png(fout, height=8, width=16, units='in', res=300)

superheat(
  as.matrix(xcontrol)[,order(score)],
  membership.rows=column,
  pretty.order.rows=TRUE,
  heat.pal=col.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label='none',
  title='Feature distribution by column'
)

if (save.fig) dev.off()

# score positional effect for row ordering
score_row <- function(z) {
  mean(z[row == max(row)]) - mean(z[row == min(row)])
}

score <- apply(xcontrol, MAR=2, score_row)
row.levels <- str_c('Row ', 1:max(row))
row <- factor(str_c('Row ', row), levels=row.levels)

fout <- str_c(fig.dir, 'fig_supp_row_unnormalized_full.png')
if (save.fig) png(fout, height=8, width=16, units='in', res=300)

superheat(
  as.matrix(xcontrol)[,order(score)],
  membership.rows=row,
  pretty.order.rows=TRUE,
  heat.pal=col.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label='none',
  title='Feature distribution by row'
)

if (save.fig) dev.off()

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

fout <- str_c(fig.dir, 'fig_supp_unnormalized_hoechst.pdf')
if (save.fig) pdf(fout, height=8, width=16)

filter(xfull, Compound_ID == 'DMSO') %>%
  mutate(ColID=as.factor(ColID)) %>%
  ggplot(aes_string(x='PlateID', y=f)) +
  geom_boxplot(aes(fill=ColID, col=ColID), alpha=0.7, outlier.size=0.5) +
  facet_wrap(~Cell_Line, scales='free_x') +
  theme(axis.text.x=element_blank()) +
  scale_fill_jama() +
  scale_color_jama() +
  ylab('Hoechst, mean nuclear intensity (KS)') +
  xlab('Plate') +
  theme(legend.position='none') +
  ylim(c(-1.1, 0.3))

if (save.fig) dev.off()


fout <- str_c(fig.dir, 'fig_supp_normalized_hoechst.pdf')
if (save.fig) pdf(fout, height=8, width=16)

filter(x, Compound_ID == 'DMSO') %>%
  mutate(ColID=as.factor(ColID)) %>%
  ggplot(aes_string(x='PlateID', y=f)) +
  geom_boxplot(aes(fill=ColID, col=ColID), alpha=0.7, outlier.size=0.5) +
  facet_wrap(~Cell_Line, scales='free_x') +
  theme(axis.text.x=element_blank()) +
  scale_fill_jama() +
  scale_color_jama() +
  ylab('Hoechst, mean nuclear intensity\n(normalized KS)') +
  xlab('Plate') +
  theme(legend.position='none') +
  ylim(c(-1.1, 0.3))

if (save.fig) dev.off()

#' # Well position effects - post normalization
#' To assess the effect of position normalization over all features, we plot
#' feature distribution (KS statistic) for DMSO wells against row/column. We 
#' threshold KS statistics at +/- `r ks.thresh` for visualization.
#+ superheat_normalized, fig.height=20, fig.width=30
# Feature distribution plots
xcontrol <- filter(x, Control)
plate.id <- xcontrol$PlateID
cell.line <- xcontrol$Cell_Line
column <- xcontrol$Col
row <- xcontrol$Row

xcontrol <- dplyr::select(xcontrol, matches('^nonborder'))
colnames(xcontrol) <- str_remove_all(colnames(xcontrol), 'nonborder\\.\\.\\.')

# threshold KS values for visualization
xcontrol[xcontrol < -ks.thresh] <- -ks.thresh
xcontrol[xcontrol > ks.thresh] <- ks.thresh

# score positional effect for column ordering
score_column <- function(z) {
  mean(z[column == max(column)]) - mean(z[column == min(column)])
}

score <- apply(xcontrol, MAR=2, score_column)
column <- str_c('Column ', column)

fout <- str_c(fig.dir, 'fig_supp_col_normalized_full.png')
if (save.fig) png(fout, height=8, width=16, units='in', res=300)

superheat(
  as.matrix(xcontrol)[,order(score)],
  membership.rows=column,
  pretty.order.rows=TRUE,
  heat.pal=col.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label='none',
  title='Feature distribution by column'
)

if (save.fig) dev.off()

# score positional effect for row ordering
score_row <- function(z) {
  mean(z[row == max(row)]) - mean(z[row == min(row)])
}

score <- apply(xcontrol, MAR=2, score_row)
row.levels <- str_c('Row ', 1:max(row))
row <- factor(str_c('Row ', row), levels=row.levels)

fout <- str_c(fig.dir, 'fig_supp_row_normalized_full.png')
if (save.fig) png(fout, height=8, width=16, units='in', res=300)

superheat(
  as.matrix(xcontrol)[,order(score)],
  membership.rows=row,
  pretty.order.rows=TRUE,
  heat.pal=col.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label='none',
  title='Feature distribution by row'
)

if (save.fig) dev.off()



#' # PCA analysis
#' To qualitatively assess the feature distribution of different compound 
#' classes we plot PCA projections grouped by various metadata features. Broadly
#' speaking, the primary sources of variation in the data correspond to compound
#' categories.
#+ pca, fig.height=6, fig.width=12
################################################################################
# Full dataset PCA
################################################################################
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

pc.list <- list(c('PC1', 'PC2'))
for (pcs in pc.list) {
  
  fout <- str_c(fig.dir, 'fig_supp_col_normalized_pca.pdf')
  if (save.fig) pdf(fout, height=8, width=16)
  p <- group_by(xc, Cell_Line) %>%
    sample_frac(0.5) %>%
    mutate(ColID=as.factor(ColID)) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=ColID), alpha=0.8) +
    ggtitle('Normalized') +
    facet_wrap(~Cell_Line) +
    scale_color_jama() + #viridis_c(limits=c(-10, 30)) +
    theme(legend.position='none')
  plot(p)
  if (save.fig) dev.off()
  
  
  fout <- str_c(fig.dir, 'fig_supp_col_unnormalized_pca.pdf')
  if (save.fig) pdf(fout, height=8, width=16)
  p <- group_by(xc.full, Cell_Line) %>%
    sample_frac(0.5) %>%
    mutate(ColID=as.factor(ColID)) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=ColID), alpha=0.5) +
    ggtitle('Unnormalized') +
    facet_wrap(~Cell_Line) +
    scale_color_jama() + #viridis_c(limits=c(-10, 30)) +
    theme(legend.position='none')
  plot(p)
  if (save.fig) dev.off()
  
  fout <- str_c(fig.dir, 'fig_supp_row_normalized_pca.pdf')
  if (save.fig) pdf(fout, height=8, width=16)
  p <- group_by(xc, Cell_Line) %>%
    sample_frac(0.5) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=RowID), alpha=0.8) +
    ggtitle('Normalized') +
    facet_wrap(~Cell_Line) +
    scale_color_viridis_c() +
    labs(color='Row')
  plot(p)
  if (save.fig) dev.off()
  
  fout <- str_c(fig.dir, 'fig_supp_row_unnormalized_pca.pdf')
  if (save.fig) pdf(fout, height=8, width=16)
  p <- group_by(xc.full, Cell_Line) %>%
    sample_frac(0.5) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=RowID), alpha=0.5) +
    ggtitle('Unnormalized') +
    facet_wrap(~Cell_Line) +
    scale_color_viridis_c() +
    labs(color='Row')
  plot(p)
  if (save.fig) dev.off()
  
}


save(file=output.file, x)
