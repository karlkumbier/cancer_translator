#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)
library(patchwork)
library(hrbrthemes)

select_category <- function(x) {
  # Select majority value from vector of categories
  x <- na.omit(x)
  if (length(x) == 0) return(NA)
  xtab <- table(x)
  return(names(xtab)[which.max(xtab)])
}

# Initialize normalization and scoring functions for clustering analysis
normalization <- function(z) {
  return(rank(z) / length(z))
  #z[z > quantile(z, 0.95)] <- quantile(z, 0.95)
  #z[z < quantile(z, 0.05)] <- quantile(z, 0.05)
  #return(z / sqrt(sum(z ^ 2)))
}

score_fun <- function(z) {
  # Compute enrichment scores from cluster enrichment matrix
  z[z == 0] <- min(z[z != 0])
  logp <- -log(z, base=10)
  p.gap <- apply(logp, MAR=1, function(zz) tail(diff(sort(zz)), 1))
  return(median(p.gap))
}

cluster_stability <- function(x, k) {
  # Compute kmeans cluster stability
  c1 <- kmeans(x, centers=k, nstart=3, iter.max=25)$cluster
  c1 <- sapply(1:k, function(kk) as.numeric(c1 == kk))
  #c1[c1 == 0] <- -1
  
  c2 <- kmeans(x, centers=k, nstart=3, iter.max=25)$cluster
  c2 <- sapply(1:k, function(kk) as.numeric(c2 == kk))
  #c2[c2 == 0] <- -1
  
  #csim <- (t(c1) %*% c2) / nrow(c1)
  csim <- cor(c1, c2)
  row.max <- apply(csim, MAR=2, max)
  col.max <- apply(csim, MAR=2, max)
  return((mean(row.max) + mean(col.max)) / 2)
}

# Initialize analysis parameters
intensity.normalize <- TRUE
n.core <- 16
min.cat <- 5
n.subsample <- 100

# Initialize color palettes
heat.pal <- viridis::viridis(10)
heat.pal.signed <- RColorBrewer::brewer.pal(9, 'RdBu')

################################################################################
# Load dataset
################################################################################
# Initialize input/output directories and load KS profile data
analysis.dir <- '~/github/cancer_translator/'
setwd(analysis.dir)
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

output.dir <- 'results/cell_line/'
dir.create(output.dir, showWarnings=FALSE)
output.file <- str_c(output.dir, 'bioactivity.Rdata') 

intensity.normalize <- TRUE
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
#' profiles. Compounds are first filtered basedon bioactivity — as specified in
#' the `cell_line_broad` notebook.
#+ bioactivity, fig.height=8, fig.width=14, echo=FALSE
################################################################################
# Compute bioactivity scores
################################################################################
# Normalize data matrix
l2norm <- function(z) z / sqrt(sum(z ^ 2))
qt90 <- function(x) quantile(x, 0.99)
x <- mutate_if(x, is.numeric, l2norm)

# Compute bioactivity scores for each well
xdist <- mclapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(bioactivity(out, qt90))
}, mc.cores=n.core)

# Group bioactivity scores by treatment, cell line
xgroup <- rbindlist(xdist) %>%
  group_by(Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup() %>%
  group_by(Cell_Line, Compound_ID, Compound_Category) %>% 
  summarize(DistNorm=mean(DistNorm), pval=mean(pval), .groups='drop') %>%
  arrange(Compound_ID, Cell_Line) 

#' # Clustering analysis
#' The following case study consider the problem of clustering compounds based 
#' on phenotypic profiles. We filter to categories with at least `R min.cat` 
#' compounds and remove inactive compounds before clustering.
#+ clustering
################################################################################
# Initialize cluster analysis dataset
################################################################################
# Initialize bioactive compound set
compounds.select <- filter(xgroup, pval == 0)

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

################################################################################
# Initialize cluster analysis parameters
################################################################################
# Initialize clustering parameters
kvals <- 2:50
pct.var.thresh <- 0.9
n.stab <- 10

#' ## Case study for A549
#' 
#' Analysis summary:
#' - Features normalized (rank v. thresholded and l2 normalized)
#' - PCA features computed (unnormalized, normalized)
#' - Feature set = PCs that explain up to `r pct.var.thresh` variance
#+ a549_analysis

# Initialize feature set for select cell line
xc <- filter(x, Cell_Line == 'A549')
xc.feat <- select(xc, matches('^non')) %>% mutate_all(normalization)
colnames(xc.feat) <- str_remove_all(colnames(xc.feat), '^nonborder\\.\\.\\.')

# Compute principle component features for select cell line
xc.pca <- prcomp(xc.feat)
var.explained <- xc.pca$sdev ^ 2 / sum(xc.pca$sdev ^ 2)
pct.var <- cumsum(var.explained)
id.select <- pct.var < pct.var.thresh
xc.feat.pca <- xc.pca$x[,id.select]

#+ pca_heatmap, fig.height=12, fig.width=12
superheat(
  xc.feat,
  row.dendrogram=TRUE,
  pretty.order.rows=TRUE,
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

#' ### PC qualitative assessment
#' Below we consider whether PC components that capture low variance in the 
#' data are enriched for compound categories. In other words, whether these 
#' components are biologicially relevant or noise.
#+ pca_assessment, fig.height=6, fig.width=8
################################################################################
# Top compound IDs by PC
################################################################################
pc.l2 <- apply(xc.feat.pca, MAR=2, function(z) sqrt(sum(z ^ 2)))
pc.l2 <- round(pc.l2 / max(pc.l2), 3)

data.frame(PC=names(pc.l2), l2=pc.l2) %>%
  mutate(PC=factor(PC, levels=PC)) %>%
  ggplot(aes(x=PC, y=l2)) +
  geom_bar(stat='identity') +
  theme_ipsum()

#+ pca_assessment_top, fig.height=18, fig.width=24
reshape2::melt(xc.feat.pca) %>%
  group_by(Var2) %>%
  top_n(10, value) %>%
  mutate(Category=xc$Compound_Category[Var1]) %>%
  mutate(ID=xc$Compound_ID[Var1]) %>%
  ggplot(aes(x=reorder_within(ID, value, Var2), y=value, fill=Category)) +
  geom_bar(stat='identity') +
  coord_flip() +
  facet_wrap(~Var2, scales='free') +
  theme_ipsum() +
  theme(legend.position='none') +
  theme(axis.text.x=element_blank())

reshape2::melt(xc.feat.pca) %>%
  group_by(Var2) %>%
  top_n(10, -value) %>%
  mutate(Category=xc$Compound_Category[Var1]) %>%
  mutate(ID=xc$Compound_ID[Var1]) %>%
  ggplot(aes(x=reorder_within(ID, value, Var2), y=value, fill=Category)) +
  geom_bar(stat='identity') +
  coord_flip() +
  facet_wrap(~Var2, scales='free') +
  theme_ipsum() +
  theme(legend.position='none') +
  theme(axis.text.x=element_blank())

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

superheat(
  apply(xplot, MAR=2, rank)[,1:4],
  pretty.order.rows=TRUE,
  row.dendrogram=TRUE,
  left.label.text.size=2.5,
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1)
)

superheat(
  apply(xplot, MAR=2, rank)[,13:16],
  pretty.order.rows=TRUE,
  row.dendrogram=TRUE,
  left.label.text.size=2.5,
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1)
)

#' ### PCA clustering unnormalized
#' Below we evaluate the strength of clusters (wrt silhouette score, enrichment)
#' relative to unnormalized PCs
#+ clustering_unnormalized, fig.height=8, fig.width=12
# Cluster compounds for select cell line
library(dbscan)

cl <- hdbscan(xc.feat.pca, minPts=5)
pairs(xc.feat.pca[,1:2], col=cl$cluster, pch=20)


cluster.analysis <- lapply(kvals, function(k) {
  cluster <- kmeans(xc.feat.pca, centers=k, nstart=10, iter.max=25)$cluster
  
  # Cluster stability analysis
  stability <- replicate(n.stab, cluster_stability(xc.feat.pca, k))
  return(list(cluster=cluster, stability=stability))
})

clusters <- lapply(cluster.analysis, function(z) z$cluster)
stability <- sapply(cluster.analysis, function(z) mean(z$stability))

# Compute clusters + silhouettes for range of k
xc.dist <- dist(xc.feat.pca)

silhouettes <- sapply(clusters, function(cluster) {
  return(cluster::silhouette(cluster, xc.dist)[,'sil_width'])
})

# Compute category enrichment at select resolution
enrichments <- lapply(clusters, function(cluster) {
  
  # Compute enrichment of compound categories by cluster
  enrichment <- sapply(1:max(cluster), function(k) {
    hyper.test <- sapply(unique(xc$Compound_Category), function(cat) {
      x <- sum(xc$Compound_Category == cat & cluster == k)
      m <- sum(xc$Compound_Category == cat)
      n <- sum(xc$Compound_Category != cat)
      return(1 - phyper(x, m, n, sum(cluster == k)))
    })
    
    return(hyper.test)
  })
  
  return(data.frame(enrichment))
})

xplot.unnorm <- data.frame(K=kvals) %>%
  mutate(Silhouette=colMeans(silhouettes)) %>%
  mutate(Stability=stability) %>%
  mutate(PC_Normalize=FALSE)

# Initialize visualization parameters
p1 <- xplot.unnorm %>%
  ggplot(aes(x=K, y=Silhouette, col=PC_Normalize)) +
  geom_line(size=1.1) +
  theme_ipsum() +
  scale_color_d3() +
  theme(legend.position='none')

p2 <- xplot.unnorm %>% 
  ggplot(aes(x=K, y=Stability, col=PC_Normalize)) +
  geom_line(size=1.1) +
  scale_color_d3() +
  theme_ipsum() +
  theme(legend.position='none')

xplot.enrich <- lapply(enrichments, score_fun) %>%
  reshape2::melt() %>%
  mutate(K=kvals[L1]) %>%
  mutate(PC_Normalize=TRUE)

p3 <- xplot.enrich %>% 
  ggplot(aes(x=K, y=value, col=PC_Normalize)) +
  geom_line(size=1.1) +
  scale_color_d3() +
  theme_ipsum() +
  theme(legend.position='none')

gridExtra::grid.arrange(p1, p2, p3, ncol=1)
