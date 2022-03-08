#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)

intensity.normalize <- TRUE
n.core <- 16
min.cat <- 5

reference.compounds <- c(
  'HSP90',
  'MT',
  'mTOR',
  'DNA',
  'Proteasome',
  'HDAC'
)

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

silhouette_med <- function(x, labels) {
  # Compute silhouette from class medioids
  x <- as.matrix(x)
  s <- numeric(length(labels))
  
  for (i in 1:length(labels)) {
    idx <- setdiff(which(labels == labels[i]), i)
    ai <- median(x[i, idx])
    
    bi <- sapply(setdiff(labels, labels[i]), function(l) {
      return(median(x[i, labels == l]))
    })
    
    bi <- min(bi)
    
    if (ai == bi) {
      s[i] <- 0
      next
    }
    
    s[i] <- ifelse(ai > bi, bi / ai - 1, 1 - ai / bi)
  }
  
  return(s)
}

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

# Summarize compound category counts
xcat.tab <- group_by(x, Compound_ID) %>%
  summarize(Compound_Category=list(Compound_Category)) %>%
  mutate(Count=sapply(Compound_Category, count_categories)) %>%
  mutate(Compound_Category=sapply(Compound_Category, select_category))

mutate(xcat.tab, Count=as.factor(Count)) %>%
  ggplot(aes(x=Count)) +
  geom_bar() +
  theme_bw() +
  ggtitle('# compound category labels')

filter(xcat.tab, Compound_Category %in% reference.compounds) %>%
  mutate(Count=as.factor(Count)) %>%
  ggplot(aes(x=Count)) +
  geom_bar() +
  theme_bw() +
  ggtitle('# compound category labels', 'reference compounds')

filter(xcat.tab, !is.na(Compound_Category)) %>%
  group_by(Compound_Category) %>%
  summarize(Count=n()) %>%
  mutate(Reference=!Compound_Category %in% reference.compounds) %>%
  ggplot(aes(x=reorder(Compound_Category, Count), y=Count, fill=Reference)) +
  geom_bar(stat='identity') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90), legend.position='none') +
  ggtitle('Compound category counts') +
  scale_fill_hue(l=50)

filter(xcat.tab, Compound_Category %in% reference.compounds) %>%
  group_by(Compound_Category) %>%
  summarize(Count=n()) %>%
  ggplot(aes(x=reorder(Compound_Category, Count), y=Count)) +
  geom_bar(stat='identity') +
  theme_bw() +
  ggtitle('Compound category counts', 'reference compounds')

category.table <- table(xcat.tab$Compound_Category)
cat.keep <- names(category.table[category.table >= min.cat])
cat.keep <- setdiff(cat.keep, 'Others')

# Clean compound category/target/pathways
x <- group_by(x, Compound_ID) %>%
  mutate(Compound_Category=select_category(Compound_Category)) %>%
  mutate(Target=select_category(Target)) %>%
  mutate(Pathway=select_category(Pathway)) %>%
  ungroup()

xref <- filter(x, Compound_Category %in% reference.compounds) %>%
  select(Compound_ID, Compound_Category) %>%
  distinct()

table(xref$Compound_Category)

# Initialize color palettes
heat.pal <- c('#FFFFFF', pal_material("indigo")(10))
heat.pal.signed <- RColorBrewer::brewer.pal(9, 'RdBu')

#' # Overview
#' This notebook considers the problems of identifying bioactive compounds 
#' based on a panel of 6 cell lines: OVCAR4, A549, DU145, ALS-WT, HEPG2, and 
#' 786-0. In addition to individual cell lines, we consider whether bioactivity 
#' calls can be improved by aggregating information across cell lines.
#' 
#' # Bioactivity scoring
#' We compute bioactivity scores as follows:
#'
#' 1. Subsample DMSO wells — 2/3 of DMSO wells, approximating the proportion 
#' of unique wells from a bootstrap sample.
#' 2. Compute the center of the DMSO point cloud among subsampled wells (i.e. 
#' average feature values across each well in the subsample) 
#' 3. Compute the l2 distance between each DMSO well in the subsample and the 
#' DMSO point cloud center 
#' 4. Define the DMSO maximal distance as the maximum distance (from 3) over all 
#' subsampled wells.
#' 
#' We then ask whether a given well compound/dose is further from the DMSO
#' point cloud center than the maximal DMSO distance. Repeating this process
#' across many subsamples allows us to generate (i) bioactivity p-value and
#' (ii) bioactivity relative (to maximum DMSO well) distance. 
#' 
#' **Note:** Bioactivity scores are computed within plate and subsequently 
#' aggregated when looking at treatments (i.e. dose/compound pairs).
#' 
#' **Note:** Features are $l2$ normalized prior to computing distances to 
#' prevent scaling differences among features from skewing results.
#+ bioactivity, fig.height=8, fig.width=14, echo=FALSE

################################################################################
# Compute bioactivity scores
################################################################################
setwd(analysis.dir)

# Normalize data matrix
l2norm <- function(z) z / sqrt(sum(z ^ 2))
x <- mutate_if(x, is.numeric, l2norm)

# Compute bioactivity scores for each well
xdist <- mclapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(bioactivity(out))
}, mc.cores=n.core)

xtab <- rbindlist(xdist) %>%
  mutate(ID=str_c(PlateID, '_', WellID)) %>%
  group_by(Compound_ID) %>%
  group_by(Cell_Line, Compound_ID, Dose, Dose_Category, Compound_Category, Target, Pathway) %>% #, SMILES) %>% 
  summarize(Bioactive=(mean(pval) == 0), Dist=mean(DistNorm), ID=list(ID), .groups='drop') %>%
  mutate(NID=sapply(ID, length)) %>%
  filter(NID < 10) %>%
  select(-NID) %>%
  mutate(ID=sapply(ID, str_c, collapse=', ')) %>%
  arrange(desc(Dist))

#' Here we assess bioactivity calls by compound/cell line. We take the minimum
#' p-value across all dose levels of a given treatment, effectively asking
#' whether a cell line detects bioactivity at some dose level.
# Group bioactivity scores by treatment, cell line
n.cell.line <- length(unique(x$Cell_Line))
xgroup <- rbindlist(xdist) %>%
  group_by(Cell_Line, Compound_ID, Dose_Category, Compound_Category) %>% 
  summarize(DistNorm=mean(DistNorm), pval=mean(pval), .groups='drop') %>%
  group_by(Cell_Line, Compound_ID, Compound_Category) %>%
  summarize(DistNorm=max(DistNorm), pval=min(pval), .groups='drop') %>%
  group_by(Compound_ID) %>%
  mutate(Count=n()) %>%
  filter(Count == n.cell.line) %>%
  ungroup() %>%
  arrange(Compound_ID, Cell_Line)

#' ### Fig 1. Bioactivity by cell line. 
#' Bioactivity calls by compound), cell line. Compounds with p-value = 0 (for > 
#' 0 doses) are defined as bioactive, while compounds with p-value = 0 (for all 
#' doses) are defined as inactive.
#' 
#' **Note:** treatment (compound/dose pairs) replicates are grouped and p-values 
#' averaged before calling bioactivity. Thus, bioactive treatments correspond to 
#' those with p-value = 0 across all replicates.
#+ bioactivity_p, fig.height=12, fig.width=18
# Format cell x compound bioactivity for visualization
xplot <- matrix(as.numeric(xgroup$pval == 0), nrow=n.cell.line)
rownames(xplot) <- unique(xgroup$Cell_Line)
colnames(xplot) <- unique(str_c(xgroup$Compound_ID))#, ', ', xgroup$Dose_Category))
category <- matrix(xgroup$Compound_Category, nrow=n.cell.line)[1,]

# Filter to compounds with bioactive calls in > 0 cell lines
id.drop <- colMeans(xplot == 0) == 1
xplot <- xplot[,!id.drop]
category <- category[!id.drop]

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

#' ### Fig 2. Bioactivity by cell line, compound category.
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
id.keep <- (category %in% cat.keep)

xplot.cat <- sapply(unique(category[id.keep]), function(ctg) {
  rowMeans(xplot[,which(category == ctg)])
})

# reorder for visualization
xplot.cat <- xplot.cat[,order(colMeans(xplot.cat))]
xplot.cat[xplot.cat < 0.25] <- 0.25

avg.bioactive <- rowMeans(xplot.cat)
row.order <- order(avg.bioactive)
col.order <- order(xplot.cat[tail(row.order, 1),])

superheat(
  xplot.cat[row.order, col.order],
  yt=c(category.table[colnames(xplot.cat)][col.order]),
  yt.plot.type='bar',
  yt.axis.name='# of compounds',
  yt.axis.name.size=12,
  yr=avg.bioactive[row.order],
  yr.plot.type='bar',
  yr.axis.name='Average proportion biactivity',
  yr.axis.name.size=12,
  bottom.label.text.angle=90,
  bottom.label.size=0.75,
  heat.pal=heat.pal, 
  heat.pal.values=seq(0, 1, by=0.1),
  title='Bioactivity by cell line, compound\nprevalent compound categories'
)

#' # Compound clustering
#' To assess the similarity between different compound classes, we compute 
#' distances between compounds based on their phenotypic profiles and assess
#' the similarity of compounds within the same category versus between 
#' categories. Prior to assessing compound similarity, we filter to treatments
#' (compound/dose pairs) identified as bioactive.
#' 
#' ### Fig 3 Compound clustering
#+ compound_clustering, fig.height=12, fig.width=18, warning=FALSE
compound_similarity <- function(xdist, categories) {
  # Compute within/between class similarity by compound category
  
  categories.num <- as.numeric(as.factor(categories))
  score <- cluster::silhouette(categories.num, xdist)
  
  if (FALSE) {
    xdist <- as.matrix(xdist)
    score <- sapply(unique(categories), function(ctg) {
      dist.within <- mean(xdist[categories == ctg, categories == ctg])
      
      dist.between <- sapply(setdiff(unique(categories), ctg), function(btw) {
        mean(xdist[categories == ctg, categories != btw])
      })
      
      return(dist.within / min(dist.between))
    })
  }
  
  out <- score[,'sil_width']
  names(out) <- categories
  return(out)
}

################################################################################
# Compound similarity analysis
################################################################################
# Initialize bioactive compound set
compounds.select <- xtab %>%
  mutate(Treat=str_c(Compound_ID, ', ', Dose_Category)) %>%
  group_by(Treat) %>%
  mutate(Count=n()) %>%
  filter(Count == n.cell.line) %>%
  filter(Bioactive)

# Initialize cluster score variables
scores <- vector('list', length(unique(x$Cell_Line)))
scores.ref <- vector('list', length(unique(x$Cell_Line)))

names(scores) <- unique(x$Cell_Line)
names(scores.ref) <- unique(x$Cell_Line)

for (cl in unique(x$Cell_Line)) {
  
  # Filter to bioactive compounds, select cell line
  xc <- filter(x, Cell_Line == cl) %>%
    mutate(Treat=str_c(Compound_ID, ', ', Dose_Category)) %>%
    filter(Treat %in% compounds.select$Treat) %>%
    filter(Compound_Category %in% cat.keep) %>%
    select(-matches('^PC')) %>%
    group_by(Compound_ID, Dose_Category, Compound_Category) %>%
    summarize_if(is.numeric, mean) %>%
    ungroup()
  
  ##############################################################################
  # Full compound set similarity
  ##############################################################################
  # Compute within/between class similarity
  xfeat <- select(xc, matches('^nonborder'))
  colnames(xfeat) <- NULL
  categories <- xc$Compound_Category
  
  # TODO: why some categories below selection threshold?
  scores[[cl]] <- compound_similarity(dist(xfeat), categories)
  
  ncat <- length(unique(categories))
  cat.pal <- scales::hue_pal()(ncat)
  row.cols <- cat.pal[as.factor(categories)]
  
  # Cluster compounds by class
  heatmap(
    as.matrix(xfeat),
    col=heat.pal.signed,
    RowSideColors=row.cols,
    main=cl
  )
  
  legend('topright', fill=cat.pal, legend=levels(as.factor(categories)))
  
  ##############################################################################
  # Reference compound set similarity
  ##############################################################################
  xfeat.ref <- xfeat[categories %in% reference.compounds,]
  categories.ref <- categories[categories %in% reference.compounds]
  
  scores.ref[[cl]] <- compound_similarity(dist(xfeat.ref), categories.ref)
  
  ncat <- length(unique(categories.ref))
  cat.pal <- pal_simpsons()(ncat)
  row.cols <- cat.pal[as.factor(categories.ref)]
  
  # Cluster compounds by class
  heatmap(
    as.matrix(xfeat.ref),
    col=heat.pal.signed,
    RowSideColors=row.cols,
    main=cl
  )
  
  legend('topright', fill=cat.pal, legend=levels(as.factor(categories.ref)))
  
  # PCA plot for reference compounds
  xfeat.pca <- prcomp(xfeat.ref)
  
  p <- data.frame(PctVar=c(0, cumsum(xfeat.pca$sdev ^ 2))) %>%
    mutate(PctVar=PctVar / sum(xfeat.pca$sdev ^ 2)) %>%
    mutate(NPC=(1:n() - 1)) %>%
    ggplot(aes(x=NPC, y=PctVar)) +
    geom_point() +
    geom_line() +
    theme_bw()
  plot(p)
  
  p <- data.frame(xfeat.pca$x, Category=categories.ref) %>%
    ggplot(aes(x=PC1, y=PC2, col=Category)) +
    geom_point(size=3) +
    theme_bw() +
    scale_color_simpsons() +
    theme(text=element_text(size=12)) +
    ggtitle(cl)
  plot(p)
}

################################################################################
# Visualize full compound set clustering
################################################################################
# initialize data matrix for heatmap
scores <- lapply(scores, as.matrix)
scores <- reshape2::melt(scores) %>% 
  group_by(Var1, L1) %>%
  summarize(value=mean(value))

category <- unique(scores$Var1)
cell.line <- unique(scores$L1)

xplot <- matrix(scores$value, nrow=length(category))
rownames(xplot) <- category
colnames(xplot) <- cell.line

# Plot bioactivity by treatment, cell line
pal.range <- quantile(abs(xplot), 0.95)
xplot[xplot < -pal.range] <- -pal.range
xplot[xplot > pal.range] <- pal.range
pal.med <- mean(xplot <= 0)

superheat(
  xplot,
  yr=c(category.table)[rownames(xplot)],
  yr.plot.type='bar',
  yr.axis.name='# compounds',
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  heat.pal=heat.pal.signed,
  heat.pal.values=c(0, 0.5, 1),
  heat.lim=c(-pal.range, pal.range)
)


################################################################################
# Visualize reference compound set clustering
################################################################################
# initialize data matrix for heatmap
scores.ref <- lapply(scores.ref, as.matrix)
scores.ref <- reshape2::melt(scores.ref) %>% 
  group_by(Var1, L1) %>%
  summarize(value=mean(value))

category <- unique(scores.ref$Var1)
cell.line <- unique(scores.ref$L1)

xplot <- matrix(scores.ref$value, nrow=length(category))
rownames(xplot) <- category
colnames(xplot) <- cell.line

# Plot bioactivity by treatment, cell line
pal.range <- quantile(abs(xplot), 0.95)
xplot[xplot < -pal.range] <- -pal.range
xplot[xplot > pal.range] <- pal.range
pal.med <- mean(xplot <= 0)

superheat(
  xplot,
  yr=c(category.table)[rownames(xplot)],
  yr.plot.type='bar',
  yr.axis.name='# compounds',
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  heat.pal=heat.pal.signed,
  heat.pal.values=c(0, 0.5, 1),
  heat.lim=c(-pal.range, pal.range)
)

