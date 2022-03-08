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
  filter(n == 6) %>%
  mutate(ID=str_c(Compound_ID, Dose_Category, Compound_Usage))

x <- mutate(x, ID=str_c(Compound_ID, Dose_Category, Compound_Usage)) %>%
  filter(ID %in% x.treat$ID) %>%
  filter(!is.na(Compound_Usage))

xref <- filter(x, Compound_Usage == 'reference_cpd') %>%
  select(Compound_ID, Compound_Category) %>%
  distinct()

################################################################################
# Category analysis
################################################################################
#' Figures below summarize compound counts by category. We highlight reference 
#' compounds in particular as they provide a baseline in our clustering analysis 
#' below.
#+ count_summaries, fig.height=8, fig.width=12
# Summarize compound category counts
xcat.tab <- group_by(x, Compound_ID, Compound_Usage) %>%
  summarize(Compound_Category=list(Compound_Category)) %>%
  mutate(Count=sapply(Compound_Category, count_categories)) %>%
  mutate(Compound_Category=sapply(Compound_Category, select_category))

mutate(xcat.tab, Count=factor(Count, levels=0:2)) %>%
  ggplot(aes(x=Count)) +
  geom_bar() +
  theme_bw() +
  ggtitle('# compound category labels')

filter(xcat.tab, !is.na(Compound_Category)) %>%
  group_by(Compound_Category, Compound_Usage) %>%
  summarize(Count=n()) %>%
  mutate(Reference=Compound_Usage == 'reference_cpd') %>%
  ggplot(aes(x=reorder(Compound_Category, Count), y=Count, fill=Reference)) +
  geom_bar(stat='identity') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90), legend.position='none') +
  ggtitle('Compound category counts') +
  scale_fill_hue(l=50) +
  scale_y_log10()

filter(xcat.tab, Compound_Usage == 'reference_cpd') %>%
  group_by(Compound_Category, Compound_Usage) %>%
  summarize(Count=n()) %>%
  ggplot(aes(x=reorder(Compound_Category, Count), y=Count)) +
  geom_bar(stat='identity') +
  theme_bw() +
  ggtitle('Compound category counts', 'reference compounds')

xcat.tab <- select(xcat.tab, Compound_ID, Compound_Category) %>% distinct()
category.table <- table(xcat.tab$Compound_Category)
cat.keep <- names(category.table[category.table >= min.cat])
cat.keep <- setdiff(cat.keep, 'Others')

#' # Overview
#' This notebook considers the problems of identifying bioactive compounds 
#' based on a panel of 6 cell lines: OVCAR4, A549, DU145, ALS-WT, HEPG2, and 
#' 786-0. In addition to individual cell lines, we consider whether bioactivity 
#' calls can be improved by aggregating information across cell lines.
#' 
#' # Analysis updates
#' 
#' 1. Assessing bioactivity at the highest dose level rather than all dose levels.
#' 2. Clustering analysis is performed for highest dose level rather than all dose levels.
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
x <- mutate_if(x, is.numeric, l2norm)

# Compute bioactivity scores for each well
xdist <- mclapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(bioactivity(out))
}, mc.cores=n.core)

# Group bioactivity scores by treatment, cell line
xgroup <- rbindlist(xdist) %>%
  group_by(Compound_ID) %>% 
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup() %>%
  group_by(Cell_Line, Compound_ID, Compound_Category) %>% 
  summarize(DistNorm=mean(DistNorm), pval=mean(pval), .groups='drop') %>%
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
  
  out <- score[,'sil_width']
  names(out) <- categories
  return(out)
}

################################################################################
# Compound similarity analysis
################################################################################
# Initialize bioactive compound set
compounds.select <- filter(xgroup, pval == 0)

# Initialize cluster score variables
scores <- vector('list', length(unique(x$Cell_Line)))
scores.ref <- vector('list', length(unique(x$Cell_Line)))

names(scores) <- unique(x$Cell_Line)
names(scores.ref) <- unique(x$Cell_Line)

for (cl in unique(x$Cell_Line)) {
  
  # Filter to bioactive compounds, select cell line
  xc <- filter(x, Cell_Line == cl) %>%
    filter(Compound_ID %in% compounds.select$Compound_ID) %>%
    filter(Compound_Category %in% cat.keep) %>%
    select(-matches('^PC')) %>%
    group_by(Compound_ID) %>% 
    filter(Dose == max(as.numeric(Dose))) %>%
    group_by(Compound_ID, Compound_Category, Compound_Usage) %>%
    summarize_if(is.numeric, mean) %>%
    ungroup()
  
  ##############################################################################
  # Full compound set similarity
  ##############################################################################
  # Compute within/between class similarity
  xfeat <- select(xc, matches('^nonborder'))
  colnames(xfeat) <- NULL
  categories <- xc$Compound_Category
  
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
  xfeat.ref <- xfeat[xc$Compound_ID %in% xref$Compound_ID,]
  categories.ref <- categories[xc$Compound_ID %in% xref$Compound_ID]
  
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
category.table <- table(compounds.select$Compound_Category)

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
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  heat.pal=heat.pal.signed,
  heat.pal.values=c(0, 0.5, 1),
  heat.lim=c(-pal.range, pal.range)
)

