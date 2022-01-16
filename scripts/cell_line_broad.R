#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)

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
heat.pal <- c('#FFFFFF', pal_material("indigo")(10))
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
opt.scores <- sapply(x.generalist, function(z) z$opt)
data.frame(CellSet=names(opt.scores), PropBioactive=opt.scores) %>%
  mutate(CellSet=factor(CellSet, levels=CellSet)) %>%
  mutate(Label=round(PropBioactive, 3)) %>%
  ggplot(aes(x=CellSet, y=PropBioactive)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=Label), nudge_y=0.025) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  ylim(c(0, 1)) +
  ggtitle('Proportion bioactive compounds by optimal cell set')

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

categories.select <- c("VEGFR inhibitor", "EGFR inhibitor" )
x.specialist <- lapply(1:n.cell.line, function(k) {
  opt_bioactive_specialist(xplot, category, categories.select, k)
})

opt.scores <- sapply(x.specialist, function(z) z$opt)
data.frame(CellSet=names(opt.scores), PropBioactive=opt.scores) %>%
  mutate(CellSet=factor(CellSet, levels=CellSet)) %>%
  mutate(Label=round(PropBioactive, 3)) %>%
  ggplot(aes(x=CellSet, y=PropBioactive)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=Label), nudge_y=0.025) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
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
kvals <- 2:25
min.p <- 2e-16
cell.lines <- unique(x$Cell_Line)

normalization <- function(z) rank(z)

score_fun <- function(z) {
  z[z < min.p] <- min.p
  logz <- -log(z, base=10)
 return(mean(apply(logz, MAR=1, max)))
}

scores <- vector('list', length(cell.lines))
names(scores) <- cell.lines

clusters <- vector('list', length(cell.lines))
names(clusters) <- cell.lines

sil.plots <- vector('list', length(cell.lines))
names(sil.plots) <- cell.lines

enrichments <- vector('list', length(cell.lines))
names(enrichments) <- cell.lines

for (cl in cell.lines) {
  
  # Compute pairwise distances for select cell line
  xc <- filter(x, Cell_Line == cl)
  xc.feat <- select(xc, matches('^non')) %>% mutate_all(normalization)
  xc.dist <- dist(xc.feat)
  
  # Cluster compounds for select cell line
  xc.cluster <- hclust(xc.dist)
  clusters[[cl]] <- xc.feat # kmeans(xc.feat, centers=k.select, nstart=10)$cluster#xc.cluster
  
  # Compute clusters + silhouettes for range of k
  silhouettes <- sapply(kvals, function(k) {
    clusters <- cutree(xc.cluster, k=k)
    return(cluster::silhouette(clusters, xc.dist)[,'sil_width'])
  })
 
  sil.plots[[cl]] <- data.frame(K=kvals) %>%
    mutate(Silhouette=colMeans(silhouettes)) %>%
    mutate(Cell_Line=cl)
  
  # Compute category enrichment at select resolution
  enrichment.k <- lapply(kvals, function(k.select) {
    
    # Initialize clusters for select k  
    clusters <- kmeans(xc.feat, centers=k.select, nstart=10)$cluster#cutree(xc.cluster, k.select)
    
    # Compute enrichment of compound categories by cluster
    enrichment <- sapply(1:k.select, function(k) {
      hyper.test <- sapply(unique(xc$Compound_Category), function(cat) {
        x <- sum(xc$Compound_Category == cat & clusters == k)
        m <- sum(xc$Compound_Category == cat)
        n <- sum(xc$Compound_Category != cat)
        return(1 - phyper(x, m, n, sum(clusters == k)))
      })
      
      return(hyper.test)
    })
    
    # Normalize enrichment within clusters
    enrichment <- apply(enrichment, MAR=2, function(z) z / max(z))
    return(data.frame(enrichment))
  })
  
  enrichments[[cl]] <- enrichment.k
  scores[[cl]] <- sapply(enrichment.k, score_fun)
}

################################################################################
# Clustering summary plots
################################################################################
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

# Plot cluster enrichment for select cluster size
k <- 15
plot.thresh <- 1e-8

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

# Initialize 
superheat(
  xplot[,order(yt)],
  pretty.order.rows=TRUE,
  left.label.text.size=3,
  yr=compound.counts[rownames(enrichments.k)],
  yr.plot.type='bar',
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
)

# Plot cluster enrichment for select cell line
cl <- 'A549'
clusters.cl <- kmeans(clusters[[cl]], k, nstart=10)$cluster
xplot <- enrichments[[cl]][[k - 1]]
xplot[xplot < plot.thresh] <- plot.thresh

superheat(
  -log(xplot, base=10),
  pretty.order.rows=TRUE,
  left.label.text.size=3,
  yr=compound.counts[rownames(xplot)],
  yr.plot.type='bar',
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1)
)



################################################################################
#' # Compound predictions
#' To assess the similarity between different compound classes, we compute 
#' distances between compound phenotypic profiles and assess the similarity of 
#' compounds within the same category relative to a k-nearest neighbor classifier. 
#+ compound_prediction fig.height=12, fig.width=18, warning=FALSE
compound_prediction <- function(xfeat, categories, k=5) {
  # Compute within/between class similarity by compound category
  
  xdist.m <- as.matrix(dist(xfeat))
  diag(xdist.m) <- Inf
  
  cat.pred <- apply(xdist.m, MAR=1, function(z) {
    cat.table <- table(categories[order(z)[1:k]])
    return(names(cat.table[cat.table == max(cat.table)]))
  })
  
  out <- data.table(CategoryPred=cat.pred) %>%
   mutate(Count=sapply(CategoryPred, length)) %>%
   mutate(CategoryTrue=categories)
  
  return(out)
}

prediction_analysis <- function(x, min.cat=5, k=5) {
  
  # Filter categories for class balance
  xcat <- select(x, Compound_Category, Compound_ID) %>% 
    distinct() %>%
    group_by(Compound_Category) %>%
    sample_n(min.cat) %>%
    ungroup()
  
  # Generate feature matrix and response vector
  xc <- filter(x, Compound_ID %in% xcat$Compound_ID) 
  xfeat <- select(xc, matches('^nonborder'))
  categories <- xc$Compound_Category
  
  predict.cell <- lapply(unique(xc$Cell_Line), function(cl) {
    id.cl <- xc$Cell_Line == cl
    out <- compound_prediction(xfeat[id.cl,], categories[id.cl], k)
    out$Cell_Line <- cl
    out$Compound_ID <- xc$Compound_ID[id.cl]
    return(out)
  })
  
  return(rbindlist(predict.cell))
}

aggregate_prediction <- function(ypred, cell.line) {
  # Aggregate predictions from select cell lines and compute accuracy
  ypred.cl <- filter(ypred, Cell_Line %in% cell.line) %>%
    group_by(Compound_ID) %>%
    summarize(CategoryPred=list(unlist(CategoryPred)), CategoryTrue=unique(CategoryTrue)) %>%
    mutate(CategoryPred=lapply(CategoryPred, select_max)) %>%
    mutate(Count=sapply(CategoryPred, function(z) length(unique(z)))) %>%
    mutate(Acc=mapply(function(yp, yt) yt %in% yp, CategoryPred, CategoryTrue))
}

select_max <- function(x) {
  xtab <- table(x)
  return(names(xtab[xtab == max(xtab)]))
}

cell_set_accuracy <- function(ypred, cell.lines, categories.select=NULL, confidence.thresh=2) {
  ypred.agg <- aggregate_prediction(ypred, cell.lines)
  if (!is.null(categories.select)) 
    ypred.agg <- filter(ypred.agg, Compound_Category %in% categories.select)
    
  ypred.agg.highconf <- filter(ypred.agg, Count < confidence.thresh )
  return(c(mean(ypred.agg$Acc), mean(ypred.agg.highconf$Acc)))
}

################################################################################
# Compound  prediction analysis
################################################################################
set.seed(47)

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


ypred <- replicate(n.subsample, prediction_analysis(x), simplify=FALSE)

#' ### Fig 3a Compound prediction generalist
#+ prediction_generalist
opt_classify_generalist <- function(ypred, cell.lines, k=1) {
  # Evaluate cell line sets with optimal bioactivity across all compounds
  cell.combs <- combn(cell.lines, k, simplify=FALSE)
  accuracy <- lapply(cell.combs, cell_set_accuracy, ypred=ypred)
  accuracy <- data.frame(do.call(rbind, accuracy))
  accuracy$Cell_Set <- sapply(cell.combs, str_c, collapse=', ')
  accuracy$k <- k
  return(accuracy)
}

cell.lines <- unique(x$Cell_Line)
ypred.opt <- lapply(1:n.cell.line, function(k) {
  return(opt_classify_generalist(ypred[[1]], cell.lines, k))
})

ypred.opt <- rbindlist(ypred.opt) %>%
  group_by(k) %>%
  top_n(1, X1)

################################################################################
# Visualize full compound set clustering
################################################################################
# Generalist accuracy
reshape2::melt(accuracy) %>%
  dplyr::rename(Accuracy=value, CellLine=L1) %>%
  ggplot(aes(x=reorder(CellLine, Accuracy, median), y=Accuracy)) +
  geom_boxplot(aes(fill=CellLine, col=CellLine), alpha=0.5) +
  scale_color_nejm() +
  scale_fill_nejm() +
  theme_bw() +
  theme(legend.position='none') +
  theme(text=element_text(size=24)) +
  xlab(NULL)

reshape2::melt(accuracy.high.conf) %>%
  dplyr::rename(Accuracy=value, CellLine=L1) %>%
  ggplot(aes(x=reorder(CellLine, Accuracy, median), y=Accuracy)) +
  geom_boxplot(aes(fill=CellLine, col=CellLine), alpha=0.5) +
  scale_color_nejm() +
  scale_fill_nejm() +
  theme_bw() +
  theme(legend.position='none') +
  theme(text=element_text(size=24)) +
  xlab(NULL)

# Specialist accuracy
accuracy.category <- reshape2::melt(accuracy.category) %>%
  dplyr::rename(Accuracy=value, CellLine=L1) %>%
  group_by(CellLine, Category) %>%
  summarize(Accuracy=mean(Accuracy))

categories <- unique(accuracy.category$Category)
cell.line <- unique(accuracy.category$CellLine)

xplot <- matrix(accuracy.category$Accuracy, nrow=length(categories))
rownames(xplot) <- categories
colnames(xplot) <- cell.line

# Filter out compound categories that cannot be classified
id.drop <- apply(xplot, MAR=1, max) < 0.25
xplot <- xplot[!id.drop,]

# Initialize column ordering
yt <- colMeans(xplot)
category.table <- table(compound.table$Compound_Category)
  
superheat(
  xplot[,order(yt)],
  yr=c(category.table)[rownames(xplot)],
  yr.plot.type='bar',
  yr.axis.name='# compounds',
  yr.axis.name.size=16,
  yt=sort(yt),
  yt.plot.type='bar',
  yt.axis.name='Average accuracy\nselect categories',
  yt.axis.name.size=16,
  pretty.order.rows=TRUE,
  heat.pal=c('#FFFFFF', heat.pal),
  heat.pal.values=seq(0, 1, by=0.1),
  left.label.text.size=3,
  left.label.size=0.25
)
