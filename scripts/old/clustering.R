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

#' ### Clustering
#' Below we evaluate the strength of clusters inferred using hdbscan with 
#' minimum cluster size set to the smallest allowable compound category:
#' `r min.cat`.
#+ clustering, fig.height=8, fig.width=12

# Compute hdbscan clusters
cl <- hdbscan(xc.feat, minPts=min.cat)
cluster <- cl$cluster

# Compute clusters + silhouettes for range of k
xc.dist <- dist(as.matrix(xc.feat))

silhouette <- cluster::silhouette(cl$cluster, xc.dist)[,'sil_width']

x.sil <- data.frame(Cluster=as.factor(cl$cluster), Silhouette=silhouette) %>%
  group_by(Cluster) %>%
  arrange(Silhouette) %>%
  mutate(Idx=1:n()) %>%
  ungroup()

filter(x.sil, Cluster != 0) %>%
  ggplot(aes(x=Idx, y=Silhouette, fill=Cluster)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_hue(l=60) +
  facet_wrap(~Cluster, scales='free_y', nrow=3) +
  theme(legend.position='None') +
  theme(axis.text.y=element_blank()) +
  xlab(NULL) +
  ggtitle('Silhouette scores by cluster')

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
  mutate(NCells=xc$NCells[1:n()])

# Initialize cluster to be dropped - negative silhouette
x.sil.avg <- group_by(x.sil, Cluster) %>% summarize(Silhouette=mean(Silhouette))
k.drop <- filter(x.sil.avg, Silhouette <= 0)$Cluster

# Plot clusters/cell counts in UMAP
mutate(xplot, log_ncells=log(NCells)) %>%
  filter(!Cluster %in% k.drop) %>%
  ggplot(aes(x=X1, y=X2, col=Cluster)) +
  geom_point(alpha=0.8) +
  scale_color_hue(l=60) +
  ggtitle('UMAP-embedded phenotypic profiles') +
  xlab('UMAP1') +
  ylab('UMAP2') +
  theme(legend.position='none')

mutate(xplot, log_ncells=log(NCells)) %>%
  filter(!Cluster %in% k.drop) %>%
  ggplot(aes(x=X1, y=X2, col=log_ncells)) +
  geom_point(alpha=0.8) +
  scale_color_viridis_c() +
  ggtitle('UMAP-embedded phenotypic profiles') +
  xlab('UMAP1') +
  ylab('UMAP2')

#' To assess the degree to which inferred clusters represent compound 
#' categories, we compute the proportion of samples in each cluster belonging
#' to a particular compound category. The figure below reports the top 
#' categories (by proportion) represented in each cluster.
#+ category_summary, fig.height=16, fig.width=16
################################################################################
# Cluster category summary
################################################################################
data.frame(Cluster=as.factor(cl$cluster), Category=xc$Compound_Category) %>%
  filter(Cluster != 0) %>%
  group_by(Cluster) %>%
  mutate(ClusterN=str_c(Cluster, ', n = ', n())) %>%
  group_by(Cluster, Category, ClusterN) %>%
  summarize(Count=n(), .groups='drop') %>%
  group_by(Cluster) %>%
  mutate(Proportion=Count / sum(Count), Count=n()) %>%
  top_n(5, Proportion) %>%
  ggplot(aes(x=reorder_within(Category, Proportion, Cluster), y=Proportion)) +
  geom_col(position='dodge', width=0.8, aes(fill=Cluster)) +
  facet_wrap(ClusterN~., scales='free_y', nrow=3) +
  scale_fill_hue(l=60) +
  coord_flip() +
  xlab(NULL) +
  theme(legend.position='none') +
  ggtitle('Compound categories enrichment by cluster')
  

################################################################################
# Objective function:
# - Delta (by category) across clusters (relative to next best)
#   - I.e. is this compound uniquely enriched within cluster
# - Delta (by category) within cluster (relative to average)
#   - I.e. is this cluster enriched for a small number of compounds
# - Sum and weight by silhouette
#   - I.e. does this occur in a "strong" cluster
################################################################################

enrichment_score_ <- function(category, cluster, silhouette, k) {

  # Compute cluster category table 
  category.k <- category[cluster == k]  
  category.table.k <- c(table(category.k))
  prec.k <- max(category.table.k) / sum(category.table.k)
  
  # Compute proportion of compounds identified from max category
  max.cat <- names(category.table.k)[which.max(category.table.k)]
  rec.k <- max(category.table.k) / table(category)[max.cat]
  
  # Initialize silhouette score for cluster k
  sil.k <- silhouette$Silhouette[silhouette$Cluster == k]
  sil.k <- max(0, sil.k)
  
  out <- data.frame(Silhouette=sil.k) %>%
    mutate(Prec=prec.k) %>%
    mutate(Rec=rec.k) %>%
    mutate(Category=max.cat) %>%
    mutate(Cluster=k)
  
  return(out)
}

enrichment_score <- function(category, cluster, silhouette) {
  # Wrapper function for computing enrichment across all clusters
  out <- lapply(sort(unique(cluster)), function(k) {
    enrichment_score_(category, cluster, silhouette, k)
  })
  
  return(rbindlist(out))
}


loss <- enrichment_score(xc$Compound_Category, cl$cluster, x.sil.avg)

# Initialize analysis output
x.cluster.tab <- data.frame(Compound_ID=xc$Compound_ID) %>%
  mutate(Compound_Category=xc$Compound_Category) %>%
  mutate(Cluster=cluster) %>%
  mutate(CellLine=cell.line)

save(file=str_c(output.dir, cell.line, '.Rdata'), x.cluster.tab, loss)
