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

theme_set(theme_ipsum(axis_title_size=14, plot_title_size=14, strip_text_size=14))

select_category <- function(x) {
  # Select majority value from vector of categories
  x <- na.omit(x)
  if (length(x) == 0) return(NA)
  xtab <- table(x)
  return(names(xtab)[which.max(xtab)])
}

################################################################################
# Initialize functional parameters to be used in analysis
################################################################################
# Bioactivity thresholding - beyond 90th percentile of DMSO wells
qt90 <- function(x) quantile(x, 0.90)

# Feature normalization - outlier trim + l2 norm
l2norm <- function(z) z / sqrt(sum(z ^ 2))

normalization <- function(z) {
  z[z > quantile(z, 0.95)] <- quantile(z, 0.95)
  z[z < quantile(z, 0.05)] <- quantile(z, 0.05)
  return(l2norm(z))
}

################################################################################
# Initialize constant parameters to be used in analysis
################################################################################
# Initialize analysis parameters
intensity.normalize <- TRUE
n.core <- 16
min.cat <- 5
pval.thresh <- 0
ncell.thresh <- 50

# Initialize clustering parameters
#cell.line <- 'OVCAR4'
pct.var.thresh <- 0.95

# Initialize color palettes
heat.pal <- viridis::viridis(10)

################################################################################
# Load dataset
################################################################################
# Initialize input/output directories and load KS profile data
analysis.dir <- '~/github/cancer_translator/'
output.dir <- str_c(analysis.dir, 'results/clustering/')
dir.create(output.dir, showWarnings=FALSE)
  
setwd(analysis.dir)
source('scripts/utilities.R')

data.dir <- 'data/screens/LH_CDC_1/'
load(str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata'))

# Clean compound categories
x <- group_by(x, Compound_ID) %>%
  mutate(Compound_Category=select_category(Compound_Category)) %>%
  ungroup()

n.cell.line <- length(unique(x$Cell_Line))

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

################################################################################
# Load compound category labels
################################################################################
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

#' # Overview
#' This notebook considers the problem clustering based on A549 phenotypic 
#' profiles. Compounds are first filtered based on bioactivity — as specified in
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
# Compute bioactivity scores
################################################################################
# l2-normalize data matrix
ncells <- x$NCells
x <- mutate_if(x, is.numeric, l2norm) %>% mutate(NCells=ncells)

# Compute inverse covariance matrix
xfeat <- dplyr::select(x, matches('^nonborder'))
xcov.inv <- solve(cov(xfeat))

# Compute bioactivity scores for each well
xdist <- lapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(bioactivity(out, qt90, xcov.inv=xcov.inv, n.core=n.core))
})

# Average bioactivity scores within treatment/cell line groups
xgroup <- rbindlist(xdist) %>%
  group_by(Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup() %>%
  group_by(Cell_Line, Compound_ID, Compound_Category) %>% 
  summarize(DistNorm=mean(DistNorm), pval=mean(pval), NCells=mean(NCells), .groups='drop') %>%
  arrange(Compound_ID, Cell_Line) 

mutate(xgroup, Bioactive=pval == 0) %>%
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
#' 2. Compounds that are called as bioactive (pval $\le$ `r pval.thresh`)
#' 3. Compounds that do not result in cell death phenotype (Ncells $\ge$ `r ncell.thresh`)
#' 
#' For each compound in the resulting set, we use the maximum dose level for 
#' clustering analysis.
#+ clustering_analysis

# Initialize bioactive, non-death compound set
compounds.select <- filter(xgroup, pval <= pval.thresh, NCells >= ncell.thresh)

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
#' 
#' Analysis summary:
#' 
#' - KS profile features are outlier trimmed @ 5th, 95th percentiles and l2 
#' normalized
#' - KS profiles are projected into PC-space and PCs that explain up to 
#' `r pct.var.thresh` percent of variance are maintained.
#' 
#' These steps are intended to (i) limit the impact of "extreme" wells on 
#' computed PCs and (ii) reduce the influence of highly redundant features on
#' downstream clustering analysis.
#+ clustering_analysis_preprocess, fig.height=8, fig.width=16

# Initialize feature set for select cell line
xc <- filter(x, Cell_Line == cell.line)
xc.feat <- select(xc, matches('^non')) %>% mutate_all(normalization)
colnames(xc.feat) <- str_remove_all(colnames(xc.feat), '^nonborder\\.\\.\\.')

# Compute principle component features for select cell line
xc.pca <- prcomp(xc.feat, scale=FALSE)
var.explained <- xc.pca$sdev ^ 2 / sum(xc.pca$sdev ^ 2)
pct.var <- cumsum(var.explained)
id.select <- pct.var < pct.var.thresh
xc.feat.pca <- xc.pca$x[,id.select]

#' the figure below reports features used in subsequent clustering analysis. 
#' Rows correspond to individual wells (treatment conditions) and columns PC
#' features. Bars at the top of the plot report the percent variance captured by 
#' each PC feature.
#+ pca_heatmap, fig.height=12, fig.width=12
superheat(
  xc.feat.pca,
  row.dendrogram=TRUE,
  pretty.order.rows=TRUE,
  yt=var.explained[id.select],
  yt.plot.type='bar',
  yt.axis.name='Variance explained',
  heat.pal=heat.pal,
  bottom.label.text.angle=90,
  bottom.label.text.size=3,
  heat.pal.values=seq(0, 1, by=0.1),
  title='PCA features'
)  

#' ### PCA clustering
#' Below we evaluate the strength of clusters learned from PCs. Clusters are 
#' inferred using hdbscan, setting the  with minimum cluster size parameter 
#' equal to the minimum category size defined above (`r min.cat`). We assess
#' the strength of inferred clusters based on silhouette scores. The figure
#' below reports silhouette scores by cluster.
#' To assess the strength of inferred clusters, we compute silhouette scores
#+ clustering, fig.height=8, fig.width=12

# Compute clusters
cl <- hdbscan(xc.feat.pca, minPts=min.cat)

# Compute silhouettes for inferred clusters
xc.dist <- dist(xc.feat.pca)
silhouette <- cluster::silhouette(cl$cluster, xc.dist)[,'sil_width']

x.sil <- data.frame(Cluster=as.factor(cl$cluster), Silhouette=silhouette) %>%
  group_by(Cluster) %>%
  arrange(Silhouette) %>%
  mutate(Idx=1:n()) %>%
  ungroup()

ggplot(x.sil, aes(x=Idx, y=Silhouette, fill=Cluster)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_hue(l=60) +
  facet_wrap(~Cluster, scales='free_y') +
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
xumap <- umap(xc.feat.pca)

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
  scale_color_hue(l=60, drop=FALSE) +
  ggtitle('UMAP-embedded phenotypic profiles')

mutate(xplot, log_ncells=log(NCells)) %>%
  filter(!Cluster %in% k.drop) %>%
  ggplot(aes(x=X1, y=X2, col=log_ncells)) +
  geom_point(alpha=0.8) +
  scale_color_viridis_c() +
  ggtitle('UMAP-embedded phenotypic profiles')

#' To assess the degree to which inferred clusters represent compound 
#' categories, we compute the proportion of samples in each cluster belonging
#' to a particular compound category. The figure below reports the top 
#' categories (by proportion) represented in each cluster.
#+ category_summary, fig.height=16, fig.width=16
################################################################################
# Cluster category summary
################################################################################
data.frame(Cluster=as.factor(cl$cluster), Category=xc$Compound_Category) %>%
  group_by(Cluster) %>%
  mutate(ClusterN=str_c(Cluster, ', n = ', n())) %>%
  group_by(Cluster, Category, ClusterN) %>%
  summarize(Count=n(), .groups='drop') %>%
  group_by(Cluster) %>%
  mutate(Proportion=Count / sum(Count), Count=n()) %>%
  top_n(5, Proportion) %>%
  ggplot(aes(x=reorder_within(Category, Proportion, Cluster), y=Proportion)) +
  geom_col(position='dodge', width=0.8, aes(fill=Cluster)) +
  facet_wrap(ClusterN~., scales='free_y') +
  scale_fill_hue(l=60) +
  coord_flip() +
  xlab(NULL) +
  theme(legend.position='none') +
  ggtitle('Compound categories enrichment by cluster')

#' The figure below reports the enrichment of compound categories within each
#' cluster based on based on hypergeometric tests. The two largest clusters are
#' enriched for many categories and therefore unlikely to be useful for compound
#' classification. In contrast, the smaller clusters are targeted towards 
#' specific compound categories.  
#+ cluster_heatmap, fig.height=12, fig.width=16, warning=FALSE
################################################################################
# Compound clustering heatmap
################################################################################
# Compute enrichment of compound categories by cluster
cluster <- cl$cluster

enrichment <- sapply(min(cluster):max(cluster), function(k) {
    hyper.test <- sapply(unique(xc$Compound_Category), function(cat) {
      x <- sum(xc$Compound_Category == cat & cluster == k)
      m <- sum(xc$Compound_Category == cat)
      n <- sum(xc$Compound_Category != cat)
      if (x == 0) return(1)
      return(1 - phyper(x, m, n, sum(cluster == k)))
  })
})

# Threshold enrichment for visualization
xplot <- enrichment
xplot[xplot < 1e-8] <- 1e-8
xplot <- -log(xplot, base=10)
colnames(xplot) <- min(cluster):max(cluster)

# Filter to compounds enriched in > 1 cluster
id.keep <- rowSums(xplot > 2) > 0

superheat(
  xplot[id.keep,],
  pretty.order.cols=TRUE,
  pretty.order.rows=TRUE,
  yt=log(c(table(cluster)), base=2),
  yt.plot.type='bar',
  yt.axis.name='Log # compounds',
  left.label.text.size=3,
  left.label.size=0.3,
  heat.pal.values=seq(0, 1, by=0.1)
)

#' To assess the overall quality of inferred clusters relative to a given cell 
#' line, we define an objective function that captures (i) strength of 
#' clustering (ii) how well clusters recapitulate compound category labels.
#' Specifically, we define:
#' 
#' 1. *Precision:* of a cluster is given by proportion of samples in a given 
#' cluster belonging to the maximally represented category. I.e. how much does 
#' the cluster concentrate around a specific compound category.
#' 2. *Recall:* of a cluster is given by the proportion of samples from the 
#' maximally represented category that fall in a given cluster. I.e. how much
#' of the top category does a given cluster recover?
#' 3. *Silhouette:* is given by the average silhouette score over all samples
#' in a given cluster. I.e. how tightly do samples within a cluster group 
#' relative to samples between clusters.
#' 
#' Our objective function is given by the geometric mean of the three values 
#' above. Below we compute objective funciton values for cell line: `r cell.line`
#+ enrichment
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
