#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(ggsci)
library(hrbrthemes)
library(readr)


theme_set(
  theme_ipsum(
    axis_title_size=18, 
    strip_text_size=18, 
    axis_text_size=14,
    base_size=18, 
    base_family='sans'
  )
)

################################################################################
# Analysis parameters
################################################################################
# Initialize analysis parameters
n.core <- 16
min.cat <- 5
n.bootstrap <- 50
save.fig <- TRUE

# Initialize color palettes
heat.pal <- viridis::viridis(10)
col.pal <- pal_nejm()(8)
col.pal.2 <- pal_jama()(7)

################################################################################
# Load dataset
################################################################################
setwd('~/github/cancer_translator/')
source('scripts/utilities.R')
source('scripts/load_normalized_data.R')


# Initialize output directories
output.dir <- '~/github/cancer_translator/results/cell_line/'
fig.dir <- '~/github/cancer_translator/results/figures/fig1/'

dir.create(output.dir, showWarnings=FALSE)
output.file <- str_c(output.dir, 'bioactivity.Rdata') 

#' # Overview
#' This notebook considers the problems of bioactivity detection based on 
#' a panel of 6 cell lines: OVCAR4, A549, DU145, ALS-WT, HEPG2, and 786-0. We
#' consider how proportion biactive varies wrt threshold.
#' 
#' **Note:** We assess bioactivity relative to the highest dose level of a given 
#' compound.
#+ bioactivity, fig.height=8, fig.width=14, echo=FALSE

# Mean-center and correlation weight features by cell line
x <- correlation_weight(x, normalize=TRUE)

#' Distribution of distances from DMSO centroids for positive + negative 
#' controls by plate, cell line.
#+ bioactivity_control, fig.height=8, fig.width=14, echo=FALSE
# Compute distances between compounds and DMSO centroids
xdist <- lapply(unique(x$Cell_Line), function(cl) {
  xc <- filter(x, Cell_Line == cl)
  return(dmso_distance(xc, null_summary))
})

# Initialize functions for null distribution summaries
summary_functions <- list()

# Quantile-based summary functions
qtp <- function(x, prob) {
  id.dmso <- x$Category == 'DMSO'
  threshold <- quantile(x$Dist[id.dmso], prob)
  return(mean(x$Dist[!id.dmso] > threshold))
}

summary_functions$qt_90 <- function(x) qtp(x, 0.9)
summary_functions$qt_95 <- function(x) qtp(x, 0.95)
summary_functions$qt_99 <- function(x) qtp(x, 0.99)

# Z-score summary function
summary_functions$z_score <- function(x) {
  id.dmso <- x$Category == 'DMSO'
  threshold <- mean(x$Dist[id.dmso]) + 2 * sd(x$Dist[id.dmso])
  return(mean(x$Dist[!id.dmso] > threshold))
}

# IQR summary function
summary_functions$IQR <- function(x) {
  id.dmso <- x$Category == 'DMSO'
  threshold <- median(x$Dist[id.dmso]) + 1.5 * IQR(x$Dist[id.dmso])
  return(mean(x$Dist[!id.dmso] > threshold))
}

# Max summary function
summary_functions$Max <- function(x) {
  id.dmso <- x$Category == 'DMSO'
  threshold <- max(x$Dist[id.dmso])
  return(mean(x$Dist[!id.dmso] > threshold))
}

out <- sapply(summary_functions, function(f) {
  sapply(xdist, function(z) f(z))
})

rownames(out) <- unique(x$Cell_Line)
reshape2::melt(out) %>%
  ggplot(aes(x=reorder_within(Var1, value, Var2), y=value, fill=Var1)) +
  geom_bar(stat='identity') +
  facet_wrap(~Var2, scales='free_x', ncol=2) +
  scale_fill_jama() +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  ylab('Proportion bioactive') +
  xlab(NULL)
