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

count_categories <- function(x) {
  # Count number of unique values in vector of categories
  x <- na.omit(x)
  return(length(unique(x)))
}

# Initialize analysis parameters
intensity.normalize <- TRUE
n.core <- 16
min.cat <- 5
n.subsample <- 100

# Initialize color palettes
heat.pal <- viridis::viridis(10)#c('#FFFFFF', pal_material("indigo")(10))
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


xreference <- filter(x, Compound_Usage == 'reference_cpd') %>%
  select(Compound_Usage, Compound_Category) %>%
  filter(Compound_Category != '')

################################################################################
# Category analysis
################################################################################
#' Figures below summarize compound counts by category.
#+ count_summaries, fig.height=8, fig.width=12
xcat.tab <- filter(x, !is.na(Compound_Category)) %>%
  filter(Compound_Category != '') %>%
  dplyr::select(Compound_ID, Compound_Category) %>%
  distinct() %>%
  group_by(Compound_Category) %>%
  summarize(Count=n())

filter(xcat.tab, Count > 1) %>%
  ggplot(aes(x=Count)) +
  geom_histogram() +
  theme_bw() +
  ggtitle('Distribution of compound counts per category') +
  geom_vline(xintercept=min.cat, col='red') +
  xlab('# compounds') +
  theme(tex=element_text(size=24))

xcat.tab.top <- filter(xcat.tab, Count >= min.cat)
cat.keep <- xcat.tab.top$Compound_Category

xcat.tab.top %>%
  ggplot(aes(x=reorder(Compound_Category, Count), y=Count)) +
  geom_bar(stat='identity', position=position_dodge(preserve="single")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Compounds per category counts', str_c('top ', nrow(xcat.tab.top), ' categories')) +
  scale_y_log10()

#' # Overview
#' This notebook considers the problems of bioactivity detection, compound 
#' classification, and clustering based on a panel of 6 cell lines: OVCAR4, 
#' A549, DU145, ALS-WT, HEPG2, and786-0.
#' 
#' # Bioactivity scoring
#' We compute bioactivity scores as follows:
#'
#' 1. Subsample DMSO wells — 2/3 of DMSO wells, approximating the proportion 
#' of unique wells from a bootstrap sample.
#' 2. Compute the center of the DMSO point cloud among subsampled wells (i.e. 
#' average feature values across each well in the subsample) 
#' 3. Generate a null distribution corresponding to the l2 distance between each 
#' DMSO well in the subsample and the DMSO point cloud center. 
#' 4. Define bioactivity within a subsample as any well that is further than the 
#' 90th percentile of the DMSO null distribution from (3)
#' 
#' Bioactivity p-values (for a given well) are defined as the proportion of 
#' subsamples for which a given well is defined as bioactive in (4).
#' 
#' **Note:** Bioactivity scores are computed within plate and subsequently 
#' aggregated when looking at treatments (i.e. dose/compound pairs).
#' 
#' **Note:** Features are $l2$ normalized prior to computing distances to 
#' prevent scaling differences among features from skewing results.
#' 
#' **Note:** We assess bioactivity relative to the highest dose level of a given 
#' compound. Some compounds appear under multiple usage categories (e.g. 
#' reference_cpd & positive_ctrl_cpd). In such cases, we compute bioactivity 
#' at the maximum dose over all use categories.
#+ bioactivity, fig.height=8, fig.width=14, echo=FALSE

################################################################################
# Compute bioactivity scores
################################################################################
# Normalize data matrix
l2norm <- function(z) z / sqrt(sum(z ^ 2))
qt90 <- function(x) quantile(x, 0.9)
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

#' ### Fig 1a. Bioactivity by cell line. 
#' Bioactivity calls by compound, cell line. Compounds with p-value = 0 are 
#' defined as bioactive, while compounds with p-value > 0 are defined as 
#' inactive.
#' 
#' **Note:** treatment (compound/dose pairs) replicates are grouped and p-values 
#' averaged before calling bioactivity. Thus, bioactive treatments correspond to 
#' those with p-value = 0 across all replicates.
#+ bioactivity_p, fig.height=12, fig.width=18
# Format cell x compound bioactivity for visualization
xplot <- matrix(as.numeric(xgroup$pval == 0), nrow=n.cell.line)
rownames(xplot) <- unique(xgroup$Cell_Line)
colnames(xplot) <- unique(str_c(xgroup$Compound_ID))
category <- matrix(xgroup$Compound_Category, nrow=n.cell.line)[1,]

# Initialize column ordering and bioactivity proportions
prop.bioactive <- rowMeans(xplot)
bioactive.order <- order(prop.bioactive, decreasing=TRUE)

col.order <- order(
  xplot[bioactive.order[1],],
  xplot[bioactive.order[2],],
  xplot[bioactive.order[3],],
  xplot[bioactive.order[4],],
  xplot[bioactive.order[5],],
  xplot[bioactive.order[6],]
)

# Plot bioactivity by treatment, cell line
superheat(
  xplot[bioactive.order, col.order],
  heat.pal=heat.pal,
  yr=prop.bioactive[bioactive.order],
  yr.plot.type='bar',
  yr.axis.name='Proportion bioactive',
  yr.axis.name.size=18,
  heat.pal.values=seq(0, 1, by=0.1),
  title='Bioactivity by cell line, compound\nall compounds'
)

#' ### Fig. 1b. Bioactivity generalist
#' Optimal cell line sets for a given budget of cell lines.
#+ bioactivity_generalist
opt_bioactive_generalist <- function(x, k=1) {
  # Evaluate cell line sets with optimal bioactivity across all compounds
  cell.combs <- combn(rownames(x), k, simplify=FALSE)
  bioactive <- lapply(cell.combs, function(z) colSums(x[z,,drop=FALSE]) > 0)
  bioactive.score <- sapply(bioactive, mean)
  names(bioactive.score) <- sapply(cell.combs, str_c, collapse=', ')
  return(list(opt=bioactive.score[which.max(bioactive.score)], score=bioactive.score))
}

x.generalist <- lapply(1:n.cell.line, opt_bioactive_generalist, x=xplot)

xopt <- sapply(x.generalist, function(z) z$opt)
xopt <- data.frame(CellSet=names(xopt), Score=xopt, K=1:n.cell.line) %>%
  mutate(Score=round(Score, 3)) %>%
  mutate(K=as.factor(K))

xdist <- unlist(lapply(x.generalist, function(z) z$score))
xdist <- data.frame(CellSet=names(xdist), Score=xdist) %>%
  mutate(K=(str_count(CellSet, ',') + 1)) %>%
  mutate(K=as.factor(K))

# Plot optimal bioactivity scores
ggplot(xdist, aes(x=K, y=Score)) +
  geom_boxplot() +
  geom_point(data=xopt) +
  geom_text(data=xopt, aes(label=CellSet), nudge_y=0.025, size=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  ylab('Proportion bioactive') +
  ggtitle('Proportion bioactive compounds in select categories by optimal cell set')

#' ### Fig 2a. Bioactivity by cell line, compound category.
#' Bioactivity calls by cell line, compound category. Treatments (compound dose 
#' pairs) with p-value = 0 are defined as bioactive, while treatments with 
#' p-value < 0 are defined as inactive. Top: bioactivity calls for each 
#' treatment, bottom: proportion of bioactivity calls within each compound 
#' category. Compound categories with fewer than `R min.cat` treatments are 
#' dropped for visualization.
#' 
#' **Note:** treatment replicates are grouped and p-values averaged before 
#' calling bioactivity. Thus, bioactive treatments correspond to those with 
#' p-value = 0 across all replicates.
#+ bioactivity_category, fig.height=12, fig.width=18, warning=FALSE
# Plot bioactivity by cell line, compound category
cat.keep.plot <- filter(xcat.tab, Count >= 10)$Compound_Category
xplot.cat <- sapply(cat.keep.plot, function(ctg) {
  rowMeans(xplot[,which(category == ctg)])
})

# reorder for visualization
xplot.cat <- xplot.cat[,order(colMeans(xplot.cat))]

avg.bioactive <- rowMeans(xplot.cat)
row.order <- order(avg.bioactive)

col.order <- order(xplot.cat[tail(row.order, 1),])
n.compounds <- xcat.tab.top$Count

superheat(
  xplot.cat[row.order, col.order],
  yt=n.compounds[col.order],
  yt.plot.type='bar',
  yt.axis.name='# compounds',
  yt.axis.name.size=16,
  yr=avg.bioactive[row.order],
  yr.plot.type='bar',
  yr.axis.name='Average proportion biactivity',
  yr.axis.name.size=16,
  bottom.label.text.angle=90,
  bottom.label.size=0.75,
  heat.pal=heat.pal, 
  heat.pal.values=seq(0, 1, by=0.1),
  title='Bioactivity by cell line, compound\nprevalent compound categories'
)

#' ### Fig 2b. Bioactivity specialist.
#' Optimal cell line sets for detecting bioactivity in select compound 
#' categories.
#+ bioactivity_specialist
opt_bioactive_specialist <- function(x, categories, categories.select, k=1) {
  # Evaluate cell line sets with optimal bioactivity across all compounds
  cell.combs <- combn(rownames(x), k, simplify=FALSE)
  
  # Compute proportion of bioactive calls by category
  bioactive <- lapply(cell.combs, function(z) colSums(x[z,,drop=FALSE]) > 0)
  
  bioactive.score.cat <- lapply(bioactive, function(z) {
    return(c(by(z, categories, mean))[categories.select])
  })
  
  bioactive.score <- sapply(bioactive.score.cat, mean)
  
  names(bioactive.score) <- sapply(cell.combs, str_c, collapse=', ')
  return(list(opt=bioactive.score[which.max(bioactive.score)], score=bioactive.score))
}

# Compute objective functions scores over cell line combinations
categories.select <- unique(xreference$Compound_Category)
x.specialist <- lapply(1:n.cell.line, function(k) {
  opt_bioactive_specialist(xplot, category, categories.select, k)
})

xopt <- sapply(x.specialist, function(z) z$opt)
xopt <- data.frame(CellSet=names(xopt), Score=xopt, K=1:n.cell.line) %>%
  mutate(Score=round(Score, 3)) %>%
  mutate(K=as.factor(K))

xdist <- unlist(lapply(x.specialist, function(z) z$score))
xdist <- data.frame(CellSet=names(xdist), Score=xdist) %>%
  mutate(K=(str_count(CellSet, ',') + 1)) %>%
  mutate(K=as.factor(K))

# Plot optimal bioactivity scores
ggplot(xdist, aes(x=K, y=Score)) +
  geom_boxplot() +
  geom_point(data=xopt) +
  geom_text(data=xopt, aes(label=CellSet), nudge_y=0.025, size=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  ylab('Proportion bioactive') +
  ggtitle('Proportion bioactive compounds in select categories by optimal cell set')


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
cell.lines <- unique(x$Cell_Line)

# Initialize normalization and scoring functions for clustering analysis
normalization <- function(z) {
  return(rank(z) / length(z))
  #z[z > quantile(z, 0.95)] <- quantile(z, 0.95)
  #z[z < quantile(z, 0.05)] <- quantile(z, 0.05)
  #return(z / sqrt(sum(z ^ 2)))
}

score_fun <- function(z) {
 # Compute enrichment scores from cluster enrichment matrix
 min.p <- apply(z, MAR=1, min)
 min.k <- apply(z, MAR=1, which.min)
 
 out <- data.frame(p=min.p, k=min.k) %>%
   group_by(k) %>%
   summarize(p=mean(p))
 
 return(mean(out$p))
}

cluster_stability <- function(x, k) {
  # Compute kmeans cluster stability
  c1 <- kmeans(x, centers=k, nstart=3, iter.max=25)$cluster
  c1 <- sapply(1:k, function(kk) as.numeric(c1 == kk))
  
  c2 <- kmeans(x, centers=k, nstart=3, iter.max=25)$cluster
  c2 <- sapply(1:k, function(kk) as.numeric(c2 == kk))
  
  c.cor <- cor(c1, c2)
  row.max <- apply(c.cor, MAR=2, max)
  col.max <- apply(c.cor, MAR=2, max)
  return((mean(row.max) + mean(col.max)) / 2)
}

#' ## Case study for A549
#' 
#' Analysis summary:
#' - Features normalized (rank v. thresholded and l2 normalized)
#' - PCA features computed (unnormalized, normalized)
#' - Feature set = PCs that explain up to `r pct.var` variance
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

#+ pca_heatmap, fig.height=8, fig.width=12
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
  bottom.label.size=0.75,
  heat.pal.values=seq(0, 1, by=0.1),
  title='PCA features'
)  

# Normalize PCA features
xc.feat.pca.norm <- apply(xc.pca$x[,id.select], MAR=2, normalization)

superheat(
  xc.feat.pca.norm,
  row.dendrogram=TRUE,
  pretty.order.rows=TRUE,
  yt=var.explained[id.select],
  yt.plot.type='bar',
  yt.axis.name='Variance',
  heat.pal=heat.pal,
  bottom.label.text.angle=90,
  bottom.label.text.size=3,
  bottom.label.size=0.75,
  heat.pal.values=seq(0, 1, by=0.1),
  title='PCA features, normalized'
) 

#' ### PC qualitative assessment
#' Below we consider whether PC components that capture low variance in the 
#' data are enriched for compound categories. In other words, whether these 
#' components are biologicially relevant or noise.
#+ pca_Assessment, fig.height=8, fig.width=12
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

#' ### PCA clustering unnormalized
#' Below we evaluate the strength of clusters (wrt silhouette score, enrichment)
#' relative to unnormalized PCs
#+ clustering_unnormalized, fig.height=8, fig.width=12
# Cluster compounds for select cell line
# Cluster compounds for select cell line
cluster.analysis <- lapply(kvals, function(k) {
  cluster <- kmeans(xc.feat.pca, centers=k, nstart=10, iter.max=25)$cluster
  
  # Cluster stability analysis
  stability <- replicate(n.stab, cluster_stability(xc.feat.pca, k))
  return(list(cluster=cluster, stability=stability))
})

clusters <- lapply(cluster.analysis, function(z) z$cluster)
stability <- lapply(cluster.analysis, function(z) z$stability)

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

xplot <- data.frame(K=kvals) %>%
  mutate(Silhouette=colMeans(silhouettes)) %>%
  mutate(Enrich=-log(sapply(enrichments, score_fun), base=10)) %>%
  mutate(Cell_Line='A549')

p1 <- ggplot(xplot, aes(x=K, y=Silhouette)) +
  geom_line(size=1.1, col=pal_jco()(3)[1]) +
  theme_ipsum() +
  theme(legend.position='none')

p2 <- ggplot(xplot, aes(x=K, y=Enrich)) +
  geom_line(size=1.1, col=pal_jco()(3)[2]) +
  theme_ipsum() +
  theme(legend.position='none')

p3 <- reshape2::melt(stability) %>%
  rename(K=L1, Stability=value) %>%
  mutate(K=kvals[K]) %>%
  group_by(K) %>%
  summarize(Stability=mean(Stability)) %>%
  ggplot(aes(x=K, y=Stability)) +
  geom_line(size=1.1, col=pal_jco()(3)[3]) +
  theme_ipsum() +
  theme(legend.position='none')

gridExtra::grid.arrange(p1, p2, p3, ncol=1)

#' ### PCA clustering normalized
#' Below we evaluate the strength of clusters (wrt silhouette score, enrichment)
#' relative to normalized PCs
#+ clustering_normalized, fig.height=8, fig.width=12


# Cluster compounds for select cell line
cluster.analysis <- lapply(kvals, function(k) {
  cluster <- kmeans(xc.feat.pca.norm, centers=k, nstart=10, iter.max=25)$cluster
  
  # Cluster stability analysis
  stability <- replicate(25, cluster_stability(xc.feat.pca.norm, k))
  return(list(cluster=cluster, stability=stability))
})

clusters <- lapply(cluster.analysis, function(z) z$cluster)
stability <- lapply(cluster.analysis, function(z) z$stability)

# Compute clusters + silhouettes for range of k
xc.dist <- dist(xc.feat.pca.norm)

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

xplot <- data.frame(K=kvals) %>%
  mutate(Silhouette=colMeans(silhouettes)) %>%
  mutate(Enrich=-log(sapply(enrichments, score_fun), base=10)) %>%
  mutate(Cell_Line='A549')

# Initialize visualization parameters
p1 <- ggplot(xplot, aes(x=K, y=Silhouette)) +
  geom_line(size=1.1, col=pal_jco()(3)[1]) +
  theme_ipsum() +
  theme(legend.position='none')

p2 <- ggplot(xplot, aes(x=K, y=Enrich)) +
  geom_line(size=1.1, col=pal_jco()(3)[2]) +
  theme_ipsum() +
  theme(legend.position='none')

p3 <- reshape2::melt(stability) %>%
  rename(K=L1, Stability=value) %>%
  mutate(K=kvals[K]) %>%
  group_by(K) %>%
  summarize(Stability=mean(Stability)) %>%
  ggplot(aes(x=K, y=Stability)) +
  geom_line(size=1.1, col=pal_jco()(3)[3]) +
  theme_ipsum() +
  theme(legend.position='none')

gridExtra::grid.arrange(p1, p2, p3, ncol=1)



#' ## Full cell line analysis
#+ cell_line_analysis
# Initialize tracking for scoring/enrichment/silhouettes
scores <- vector('list', length(cell.lines))
names(scores) <- cell.lines

clusters <- vector('list', length(cell.lines))
names(clusters) <- cell.lines

sil.plots <- vector('list', length(cell.lines))
names(sil.plots) <- cell.lines

enrichments <- vector('list', length(cell.lines))
names(enrichments) <- cell.lines

for (cl in cell.lines) {
  
  # Initialize feature set for select cell line
  xc <- filter(x, Cell_Line == cl)
  xc.feat <- select(xc, matches('^non')) %>% mutate_all(normalization)
  colnames(xc.feat) <- str_remove_all(colnames(xc.feat), '^nonborder\\.\\.\\.')
  
  # Compute principle component features for select cell line
  xc.pca <- prcomp(xc.feat)
  pct.var <- cumsum(xc.pca$sdev ^ 2) / sum(xc.pca$sdev ^ 2)
  id.select <- pct.var < pct.var.thresh
  xc.feat.pca <- xc.pca$x[,id.select]
  
  # Cluster compounds for select cell line
  clusters[[cl]] <- lapply(kvals, function(k) {
    kmeans(xc.feat.pca, centers=k, nstart=10, iter.max=25)$cluster
  })
  
  superheat(
    xc.feat.pca,
    pretty.order.rows=TRUE,
    membership.rows=clusters[[cl]][[4]],
    heat.pal=heat.pal,
    bottom.label.text.angle=90,
    bottom.label.text.size=3,
    bottom.label.size=0.75,
    heat.pal.values=seq(0, 1, by=0.1)
  )  
  
  # Compute clusters + silhouettes for range of k
  xc.dist <- dist(xc.feat.pca)
  
  silhouettes <- sapply(clusters[[cl]], function(cluster) {
    return(cluster::silhouette(cluster, xc.dist)[,'sil_width'])
  })
 
  sil.plots[[cl]] <- data.frame(K=kvals) %>%
    mutate(Silhouette=colMeans(silhouettes)) %>%
    mutate(Cell_Line=cl)
  
  # Compute category enrichment at select resolution
  enrichments[[cl]] <- lapply(clusters[[cl]], function(cluster) {
    
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
  
  scores[[cl]] <- sapply(enrichments[[cl]], score_fun)
}

################################################################################
# Clustering summary plots
################################################################################
#+ cluster_summary_plots, fig.height=8, fig.width=12
# Plot silhouette score v. enrichment
p1 <- rbindlist(sil.plots) %>%
  ggplot(aes(x=K, y=Silhouette, color=Cell_Line)) +
  geom_line(size=1.25) +
  geom_point() +
  scale_color_nejm() +
  theme_bw() +
  facet_grid(~Cell_Line) +
  theme(legend.position='none') +
  ggtitle("Silhouette")

p2 <- reshape2::melt(scores) %>%
  group_by(L1) %>%
  mutate(K=kvals[1:n()]) %>%
  ungroup() %>%
  ggplot(aes(x=K, y=value, color=L1)) +
  geom_line(size=1.25) +
  geom_point() +
  scale_color_nejm() +
  theme_bw() +
  facet_grid(~L1) +
  theme(legend.position='none') +
  ggtitle('Enrichment') +
  ylab('Mean enrichment')

gridExtra::grid.arrange(p1, p2, ncol=1)

################################################################################
# Plot cluster enrichment for select cluster size
################################################################################
#+ cluster_heatmaps, fig.height=12, fig.width=20
k <- 7
cl <- 'A549'
plot.thresh <- 1e-5

# Initialize compound x cell lline enrichment matrix
enrichments.k <- sapply(cell.lines, function(cl) {
  enrichment.cl <- enrichments[[cl]][[k - 1]]
  return(apply(enrichment.cl, MAR=1, min))
})

enrichments.k[enrichments.k < plot.thresh] <- plot.thresh
xplot <- -log(enrichments.k, base=10)

# Initialize row counts
compound.counts <- compound.table$Count
names(compound.counts) <- compound.table$Compound_Category
compound.counts <- compound.counts[unique(names(compound.counts))]

# Filter out compounds with no enrichment
id.drop <- apply(xplot, MAR=1, function(z) max(z) < 3)

# Initialize row ordering
col.means <- colMeans(xplot[!id.drop,])
max.col <- which.max(colMeans(xplot[!id.drop,]))
col.order <- order(col.means)
row.order <- order(xplot[!id.drop, max.col])

# Initialize 
superheat(
  xplot[!id.drop,][row.order, col.order],
  left.label.text.size=3,
  left.label.size=0.5,
  yr=compound.counts[rownames(enrichments.k)][!id.drop][row.order],
  yr.plot.type='bar',
  yr.axis.name='# compounds', 
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
)

# Plot cluster enrichment for select cell line
clusters.cl <- clusters[[cl]][[k - 1]]
xplot <- enrichments[[cl]][[k - 1]]
xplot[xplot < plot.thresh] <- plot.thresh
xplot <- -log(xplot, base=10)
xplot[xplot < 2] <- 0

# Filter out categories without enrichment
id.drop <- apply(xplot, MAR=1, function(z) max(z) < 3)

# Initialize row ordering
row.order <- order(
  xplot[!id.drop, 1],
  xplot[!id.drop, 2],
  xplot[!id.drop, 3],  
  xplot[!id.drop, 4],
  xplot[!id.drop, 5],
  xplot[!id.drop, 6],
  xplot[!id.drop, 7]
)

# Initialize cluster counts
superheat(
  xplot[!id.drop,][row.order,],
  left.label.text.size=3,
  yt=c(table(clusters.cl)),
  yt.plot.type='bar',
  yt.axis.name='Cluster size',
  yr=compound.counts[rownames(xplot)][!id.drop][row.order],
  yr.plot.type='bar',
  yr.axis.name='# compounds',
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  heat.lim=c(0, -log(plot.thresh, base=10)),
  left.label.size=0.6
)
