#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)
library(patchwork)
library(hrbrthemes)
theme_set(theme_ipsum())

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

# Cluster enrichment - gap between top 2 -log10(p-value)
cluster_enrichment <- function(z) {
  z[z == 0] <- min(z[z != 0])
  logp <- -log(z, base=10)
  p.gap <- apply(logp, MAR=1, function(zz) tail(diff(sort(zz)), 1))
  return(median(p.gap))
}

# Cluster stability - correlation between clusterings
cluster_stability <- function(x, k) {
  # Compute kmeans cluster stability
  c1 <- kmeans(x, centers=k, nstart=3, iter.max=25)$cluster
  c1 <- sapply(1:k, function(kk) as.numeric(c1 == kk))
  
  c2 <- kmeans(x, centers=k, nstart=3, iter.max=25)$cluster
  c2 <- sapply(1:k, function(kk) as.numeric(c2 == kk))

  csim <- cor(c1, c2)
  row.max <- apply(csim, MAR=2, max)
  col.max <- apply(csim, MAR=2, max)
  return((mean(row.max) + mean(col.max)) / 2)
}

################################################################################
# Initialize constant parameters to be used in analysis
################################################################################
# Initialize analysis parameters
intensity.normalize <- TRUE
n.core <- 16
min.cat <- 5
n.subsample <- 100
pval.thresh <- 0
ncell.thresh <- 500

# Initialize clustering parameters
cell.line <- 'OVCAR4'
kvals <- 2:25
pct.var.thresh <- 0.9
n.stab <- 10

# Initialize color palettes
heat.pal <- viridis::viridis(10)

################################################################################
# Load dataset
################################################################################
# Initialize input/output directories and load KS profile data
analysis.dir <- '~/github/cancer_translator/'
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
#' the `cell_line_broad` notebook.
#+ bioactivity, fig.height=8, fig.width=12, echo=FALSE
################################################################################
# Compute bioactivity scores
################################################################################
# Normalize data matrix
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

# Group bioactivity scores by treatment, cell line
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
  scale_color_d3()

#' # Clustering analysis
#' The following case study consider the problem of clustering compounds based 
#' on phenotypic profiles. Prior to clustering we filter to a compound set of:
#' 
#' 1. Categories with at least `r min.cat` compounds 
#' 2. Compounds that are called as bioactive (pval <= `r pval.thresh`)
#' 3. Compounds that do not result in cell death phenotype (Ncells >= `r ncell.thresh`)
#+ clustering
################################################################################
# Initialize cluster analysis dataset
################################################################################
# Initialize bioactive compound set
compounds.select <- filter(xgroup, pval <= pval.thresh, NCells >= ncell.thresh)

# Get compound counts for bioactive set
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

#' ## Case study for A549
#' 
#' Analysis summary:
#' - Features normalized (rank v. thresholded and l2 normalized)
#' - PCA features computed (unnormalized, normalized)
#' - Feature set = PCs that explain up to `r pct.var.thresh` variance
#+ a549_analysis, fig.height=8, fig.width=16

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

#+ pca_heatmap, fig.height=12, fig.width=12
superheat(
  xc.feat,
  row.dendrogram=TRUE,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  heat.pal=heat.pal,
  bottom.label.text.angle=90,
  bottom.label.text.size=3,
  bottom.label.size=0.75,
  heat.pal.values=seq(0, 1, by=0.1),
  title='Raw features'
)  

superheat(
  xc.feat.pca,
  row.dendrogram=TRUE,
  pretty.order.rows=TRUE,
  yt=var.explained[id.select],
  yt.plot.type='bar',
  yt.axis.name='Variance',
  heat.pal=heat.pal,
  bottom.label.text.angle=90,
  bottom.label.text.size=3,
  heat.pal.values=seq(0, 1, by=0.1),
  title='PCA features'
)  

#+ pca_assesment_heatmap, fig.height=12, fig.width=24
################################################################################
# Top compound categories PC distributions
################################################################################
xpca.cat <- reshape2::melt(xc.feat.pca) %>%
  mutate(Category=xc$Compound_Category[Var1]) %>%
  group_by(Var2, Category) %>%
  summarize(MedCat=median(value), .groups='drop')

xplot <- matrix(xpca.cat$MedCat, ncol=ncol(xc.feat.pca))
rownames(xplot) <- unique(xpca.cat$Category)
colnames(xplot) <- unique(xpca.cat$Var2)

superheat(
  apply(xplot, MAR=2, rank),
  pretty.order.rows=TRUE,
  row.dendrogram=TRUE,
  left.label.text.size=2.5,
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1)
)


#' ### PCA clustering
#' Below we evaluate the strength of clusters learned from PCAs. Clusters are 
#' inferred using hdbscan.
#+ clustering_unnormalized, fig.height=8, fig.width=12
# Cluster compounds for select cell line
library(dbscan)

# HDBSCAN cluster compounds
cl <- hdbscan(xc.feat.pca, minPts=5)
cluster.pal <- RColorBrewer::brewer.pal(9, 'Set3')

# Initialize data for PCA plots
n.group <- floor(ncol(xc.feat.pca) / 2)
xc.pca.group <- lapply(1:n.group, function(g) {
  ii <- ((g - 1) * 2 + 1):(g * 2)
  
  out <- data.frame(X1=xc.feat.pca[,ii[1]], X2=xc.feat.pca[,ii[2]]) %>% 
    mutate(Cluster=as.factor(cl$cluster[1:n()])) %>%
    mutate(Group=g) %>%
    mutate(NCells=xc$NCells[1:n()])
  
  return(out)
})
  
  
# Plot clusters/cell counts in PCA space
rbindlist(xc.pca.group) %>%
  mutate(log_ncells=log(NCells)) %>%
  ggplot(aes(x=X1, y=X2, col=Cluster)) +
  geom_point(alpha=0.6) +
  facet_wrap(~Group) +
  scale_color_d3()

rbindlist(xc.pca.group) %>%
  mutate(log_ncells=log(NCells)) %>%
  ggplot(aes(x=X1, y=X2, col=log_ncells)) +
  geom_point() +
  facet_wrap(~Group) +
  scale_color_material('light-green')


# Compute clusters + silhouettes for range of k
xc.dist <- dist(xc.feat.pca)
silhouette <- cluster::silhouette(cl$cluster, xc.dist)[,'sil_width']

data.frame(Cluster=as.factor(cl$cluster), Silhouette=silhouette) %>%
  group_by(Cluster) %>%
  arrange(Silhouette) %>%
  mutate(Idx=1:n()) %>%
  ungroup() %>%
  ggplot(aes(x=Idx, y=Silhouette, fill=Cluster)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_d3() +
  facet_wrap(~Cluster, scales='free_y') +
  theme(legend.position='None') +
  theme(axis.text.y=element_blank()) +
  xlab(NULL) +
  ggtitle('Silhouette scores by cluster')

data.frame(Cluster=as.factor(cl$cluster), Category=xc$Compound_Category) %>%
  group_by(Cluster, Category) %>%
  summarize(Count=n(), .groups='drop') %>%
  group_by(Cluster) %>%
  mutate(Proportion=Count / sum(Count)) %>%
  top_n(10, Proportion) %>%
  ggplot(aes(x=reorder_within(Category, Proportion, Cluster), y=Proportion, fill=Category)) +
  geom_col(position='dodge', width=0.8) +
  facet_wrap(Cluster~., scales='free_y') +
  scale_fill_hue(l=60) +
  coord_flip() +
  xlab(NULL) +
  theme(legend.position='none')
  
  
#+ cluster_heatmap, fig.height=16, fig.width=16
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
xplot[xplot < 1e-5] <- 1e-5
xplot <- -log(xplot, base=10)

superheat(
  xplot,
  pretty.order.cols=TRUE,
  pretty.order.rows=TRUE,
  yt=log(c(table(cluster)), base=2),
  yt.plot.type='bar',
  yt.axis.name='Log # compounds',
  left.label.text.size=3
)

