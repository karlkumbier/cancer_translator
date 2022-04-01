#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)
library(patchwork)
library(hrbrthemes)
library(dbscan)
library(umap)

theme_set(
  theme_ipsum(
    axis_title_size=14, 
    plot_title_size=14, 
    strip_text_size=14
  )
)

################################################################################
# Initialize constant parameters to be used in analysis
################################################################################
# Initialize analysis parameters
intensity.normalize <- FALSE
n.core <- 16
min.cat <- 5
pval.thresh <- 0
ncell.thresh <- 50

# Initialize clustering parameters
#cell.line <- 'OVCAR4'

# Initialize color palettes
heat.pal <- viridis::viridis(10)

################################################################################
# Load dataset
################################################################################
setwd('~/github/cancer_translator/')
source('scripts/utilities.R')

# Load KS profile data
data.dir <- 'data/screens/LH_CDC_1/'
load(str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata'))
x <- dplyr::select(x, -matches('^PC'))

# initialize output directories
output.dir <- '~/github/cancer_translator/results/clustering/'
dir.create(output.dir, showWarnings=FALSE)

# Clean compound categories
x <- group_by(x, Compound_ID) %>%
  mutate(Compound_Category=select_category(Compound_Category)) %>%
  ungroup()

# Filter to compound/dose combinations evaluated in all cell lines
x.treat <- select(x, Cell_Line, Compound_ID, Dose_Category, Compound_Usage) %>%
  distinct() %>%
  group_by(Compound_ID, Dose_Category, Compound_Usage) %>% 
  count() %>%
  filter(n == n.cell.line) %>%
  mutate(ID=str_c(Compound_ID, Dose_Category, Compound_Usage))

x <- mutate(x, ID=str_c(Compound_ID, Dose_Category, Compound_Usage)) %>%
  filter(ID %in% x.treat$ID) %>%
  filter(!is.na(Compound_Usage)) %>%
  dplyr::select(-Compound_Category)

# Load compound target table
x.target <- 'data/Bioactivity_broad_label/' %>%
  str_c('bioactivity_broad_institute_label_20220108.csv') %>%
  fread() %>%
  dplyr::select(Compound_ID, MOA_broad_institute, Target_broad_institute) %>%
  distinct() %>%
  dplyr::rename(Compound_Category=MOA_broad_institute) 

x <- left_join(x, x.target, by='Compound_ID')

xcat.tab <- filter(x, !is.na(Compound_Category)) %>%
  filter(Compound_Category != '') %>%
  dplyr::select(Compound_ID, Compound_Category) %>%
  distinct() %>%
  group_by(Compound_Category) %>%
  summarize(Count=n())

xcat.tab.top <- filter(xcat.tab, Count >= min.cat)
cat.keep <- xcat.tab.top$Compound_Category

# Initialize # of cell lines
n.cell.line <- length(unique(x$Cell_Line))

#' # Overview
#' This notebook considers the problem clustering based on A549 phenotypic 
#' profiles. Compounds are first filtered based on bioactivity — as specified in
#' the `cell_line_broad` notebook.
#+ bioactivity, fig.height=8, fig.width=12, echo=FALSE
# Mean-center features by cell line for subsequent analysis
x <- correlation_weight(x)

#' the `cell_line_broad` notebook. At a high level, bioactivity is based on the 
#' mahalanobis distance (evaluated over phenotypic profiles) between a target
#' well and a collection of DMSO control wells. 
#' 
#' Figures below report mahalanobis distance between a well (i.e. treatment 
#' condition) and the DMSO point cloud center, normalized relative to the 
#' average distance of individual DMSO wells from the DMSO point cloud center, 
#' against cell count. Treatment conditions called as bioactive are highlighted 
#' in orange.
#+ bioactivity, fig.height=8, fig.width=12
################################################################################
# Initialize parameters for bioactivity analysis
################################################################################
# Initialize DMSO distance summary function
null_summary <- function(x) {
  iqr.thresh <- IQR(x)
  return(max(x[x < median(x) + iqr.thresh]))
}

################################################################################
# Compute bioactivity scores
################################################################################
xdist <- lapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(bioactivity(out, null_summary))
})

# Average bioactivity scores within treatment/cell line groups
xgroup <- rbindlist(xdist) %>%
  group_by(Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup() %>%
  group_by(Cell_Line, Compound_ID, Compound_Category) %>% 
  summarize(DistNorm=mean(DistNorm), NCells=mean(NCells), .groups='drop') %>%
  arrange(Compound_ID, Cell_Line) 

mutate(xgroup, Bioactive=DistNorm > 1) %>%
  ggplot(aes(x=NCells, y=DistNorm)) +
  geom_point(aes(col=Bioactive), alpha=0.6) +
  facet_wrap(~Cell_Line) +
  scale_x_log10() +
  scale_color_d3() +
  geom_vline(xintercept=ncell.thresh, lty=2, col='grey')

#' # Clustering analysis
#' The following case study consider the problem of clustering compounds based 
#' on phenotypic profiles. Prior to clustering we filter compounds as follows:
#' 
#' 1. Categories with at least `r min.cat` compounds 
#' 2. Compounds that are called as bioactive (distance > 1)
#' 3. Compounds that do not result in cell death phenotype (Ncells >= `r ncell.thresh`)
#+ clustering_analysis
################################################################################
# Initialize cluster analysis dataset
################################################################################
# Initialize bioactive compound set
compounds.select <- filter(xgroup, DistNorm > 1, NCells >= ncell.thresh)

# Get compound counts for bioactive set and filter based on minimum category size
compound.table <- select(compounds.select, Compound_ID, Compound_Category) %>%
  distinct() %>%
  group_by(Compound_Category) %>%
  filter(Compound_Category %in% cat.keep) %>%
  mutate(Count=n()) %>%
  filter(Count >= min.cat)

x <- filter(x, Compound_ID %in% compound.table$Compound_ID) %>%
  select(-matches('^PC')) %>%
  group_by(Compound_ID, Cell_Line) %>% 
  filter(Dose == max(as.numeric(Dose))) %>%
  group_by(Compound_ID, Compound_Category, Cell_Line) %>%
  summarize_if(is.numeric, mean) %>%
  ungroup()

#' ## Case study for `r cell.line`
#+ a549_analysis, fig.height=8, fig.width=16
# Initialize feature set for select cell line
xc <- filter(x, Cell_Line == cell.line)
xc.feat <- select(xc, matches('^non'))
colnames(xc.feat) <- str_remove_all(colnames(xc.feat), '^nonborder\\.\\.\\.')


#' Below we visualize inferred cluster in in umap space. In the top figure, 
#' points corresponding to wells are colored by cluster. In the bottom figure
#' points are colored by log cell count. Clearly defined clusters are marked by
#' reduced cell counts.
################################################################################
# UMAP visualization
################################################################################
# Initialize umap embedding for visualization
configs <- umap.defaults
configs$n_neighbors <- min.cat
xumap <- umap(xc.feat)

# Initialize data for visualization
xplot <- data.frame(X1=xumap$layout[,1], X2=xumap$layout[,2]) %>%
  mutate(Cluster=as.factor(cl$cluster[1:n()])) %>%
  mutate(Category=xc$Compound_Category) %>%
  mutate(NCells=xc$NCells[1:n()])

# Plot clusters/cell counts in UMAP
mutate(xplot, log_ncells=log(NCells)) %>%
  ggplot(aes(x=X1, y=X2, col=(Category == 'tubulin polymerization inhibitor'))) +
  geom_point(alpha=0.8) +
  ggtitle('UMAP-embedded phenotypic profiles') +
  xlab('UMAP1') +
  ylab('UMAP2') +
  theme(legend.position='none') +
  scale_color_nejm()
