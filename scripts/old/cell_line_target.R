#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)

################################################################################
# Helper functions for analysis
################################################################################
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

target_bioactivity <- function(target, compounds, xtarget, xbioactive) {
  # Compute proportion of compounds within target category called bioactive
  
  # Initialize compound set for select target
  cpd.target <- unlist(select(xtarget, one_of(target)))
  names(cpd.target) <- compounds
  cpd.target <- cpd.target[cpd.target != 0] %>% names
  
  xbioactive.target <- filter(xbioactive, Compound_ID %in% cpd.target) %>%
    group_by(Cell_Line) %>%
    summarize(PropBioactive=mean(pval == 0))
  
  out <- xbioactive.target$PropBioactive
  names(out) <- xbioactive.target$Cell_Line
  return(out)
}

binarize <- function(z) as.numeric(z != 0)

zero <- function(z) ifelse(is.na(z), 0, z)

group_cpd <- function(z) as.numeric(any(z > 0))

l2norm <- function(z) z / sqrt(sum(z ^ 2))

################################################################################
# Setup
################################################################################
intensity.normalize <- TRUE
n.core <- 6
min.cpd.target <- 10

# Initialize input/output directories and load KS profile data
analysis.dir <- '~/github/cancer_translator/'
setwd(analysis.dir)
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

output.dir <- 'results/cell_line/'
dir.create(output.dir, showWarnings=FALSE)
output.file <- str_c(output.dir, 'bioactivity.Rdata') 

# Load phenotypic profiling data
intensity.normalize <- TRUE
load(str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata'))

x.treat <- select(x, Cell_Line, Compound_ID, Dose_Category, Compound_Usage) %>%
  distinct() %>%
  group_by(Compound_ID, Dose_Category, Compound_Usage) %>% 
  count() %>%
  filter(n == 6) %>%
  mutate(ID=str_c(Compound_ID, Dose_Category, Compound_Usage))

x <- mutate(x, ID=str_c(Compound_ID, Dose_Category, Compound_Usage)) %>%
  filter(ID %in% x.treat$ID) %>%
  filter(!is.na(Compound_Usage))

n.cell.line <- length(unique(x$Cell_Line))

# Load compound target table
target.file <- 'data/Target_and_Activity/' %>%
  str_c('Screened_Library_Activity_Combined_on_Compound.csv')

x.target <- fread(target.file) %>%
  mutate_if(is.numeric, zero) %>%
  #mutate_if(is.numeric, binarize) %>%
  group_by(Compound) %>%
  summarize_all(mean)
  #summarize_all(group_cpd)

# Initialize color palettes
heat.pal <- c('#FFFFFF', pal_material("indigo")(10))
heat.pal.signed <- RColorBrewer::brewer.pal(9, 'RdBu')

#' # Analysis updates
#' 
#' 1. Assessing bioactivity at the highest dose level rather than all dose levels.
#' 2. Assessing bioactivity within protein targets rather than compound categories.
#' 
#' # Notes/questions/concerns
#' 
#' 1. For compounds w/ single dose, why are dose categories not all "1"?
#' 2. Why do some compounds not appear in target table — remapped names?
#' 3. Should replicate compounds in target table be treated as same compound? We
#' are currently grouping based on bioactivity for any replicate. 
#'
#' # Compound targets
#' Figures below summarize compound activity by target. Compounds with no
#' reported target are dropped for visualization.
#+ target_summaries, fig.height=8, fig.width=12

# Filter out compounds with no target data
xplot.target <- select(x.target, -Compound) %>% as.data.frame
rownames(xplot.target) <- x.target$Compound

id.drop <- rowMeans(xplot.target == 0) == 1
xplot.target <- xplot.target[!id.drop,]
compounds <- x.target$Compound[!id.drop]
print(dim(xplot.target))

data.frame(Activity=unlist(xplot.target)) %>%
  filter(Activity != 0) %>%
  ggplot(aes(x=Activity)) +
  geom_histogram() +
  theme_bw() +
  xlab('Activity (non-zero)')

# Compound x target
superheat(
  xplot.target,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  dist.method='binary',
  title='Compound x target'
)

#+ target_summaries_targeted, fig.height=8, fig.width=15
# Compound x target, highly targeted proteins
target.counts <- colSums(xplot.target != 0)
targets.select <- sort(target.counts) %>% tail(50) %>% names()
xplot.target.f <- select(xplot.target, one_of(targets.select))

superheat(
  xplot.target.f[rowSums(xplot.target.f) != 0,],
  pretty.order.cols=TRUE,
  pretty.order.rows=TRUE,
  yt=target.counts[targets.select],
  yt.plot.type='bar',
  yt.axis.name='# compounds targeting',
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label.text.angle=90,
  bottom.label.text.size=3,
  dist.method='binary',
  title='Compound x target, highly targeted proteins'
)

#+ target_summaries_targeting, fig.height=16, fig.width=16
# Compound x target, highly targeting proteins
compound.counts <- rowSums(xplot.target != 0)
compound.select <- sort(compound.counts) %>% tail(50) %>% names()
xplot.target.f <- xplot.target[compound.select,]

superheat(
  xplot.target.f[,colSums(xplot.target.f) != 0],
  pretty.order.cols=TRUE,
  pretty.order.rows=TRUE,
  yr=compound.counts[compound.select],
  yr.plot.type='bar',
  yr.axis.name='# targets',
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1),
  left.label.text.size=3,
  dist.method='binary',
  title='Compound x target, highly targeting compounds'
)

#' Figures below report summary statistics of compound/target counts.
#+ target_summaries_stats, fig.height=8, fig.width=12
# Summarize # of compounds by target
data.frame(N_Compound=colSums(xplot.target != 0)) %>%
  ggplot(aes(x=N_Compound)) +
  geom_histogram() +
  theme_bw() +
  ggtitle('Distribution # compounds targeting protein')

data.frame(N_Target=rowSums(xplot.target != 0)) %>%
  ggplot(aes(x=N_Target)) +
  geom_histogram() +
  theme_bw() +
  ggtitle('Distribution # protein targets per compound')

# Summarize # compounds by identical target
target.group <- apply(xplot.target, MAR=2, str_c, collapse='')
data.frame(N_Compound=c(table(target.group))) %>%
  ggplot(aes(x=N_Compound)) +
  geom_histogram() +
  theme_bw() +
  ggtitle('Distribution # compounds by identical target set')


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
#' 
#' **Note:** We assess bioactivity relative to the highest dose level of a given 
#' compound. Some compounds appear under multiple usage categories (e.g. 
#' reference_cpd & positive_ctrl_cpd). In such cases, we compute bioactivity 
#+ bioactivity, fig.height=8, fig.width=14, echo=FALSE

################################################################################
# Compute bioactivity scores
################################################################################
# Normalize data matrix
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
  group_by(Cell_Line, Compound_ID) %>% 
  summarize(DistNorm=mean(DistNorm), pval=mean(pval), .groups='drop') %>%
  arrange(Compound_ID, Cell_Line) 

# Filter out compounds without target info
print(mean(xgroup$Compound_ID %in% x.target$Compound))
xgroup <- filter(xgroup, Compound_ID %in% x.target$Compound)

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

#' ### Fig 2. Bioactivity by cell line, target
#' Bioactivity calls by cell line, target. Treatments (compound/dose 
#' pairs) with p-value = 0 are defined as bioactive. Top: bioactivity calls for 
#' each treatment, bottom: proportion of bioactivity calls within each compound 
#' category.
#' 
#' **Note:** treatment replicates are grouped and p-values averaged before 
#' calling bioactivity. Thus, bioactive treatments correspond to those with 
#' p-value = 0 across all replicates.
#+ bioactivity_category, fig.height=12, fig.width=18, warning=FALSE
# Plot bioactivity by cell line, compound category
targets <- colnames(xplot.target)
xplot.bioactive.target <- sapply(targets, function(tgt) {
  target_bioactivity(tgt, compounds, xplot.target, xgroup)
})

# Drop targets with no bioactive compounds
id.drop <- sapply(xplot.bioactive.target, length) != n.cell.line
xplot.bioactive.target <- do.call(cbind, xplot.bioactive.target[!id.drop])

# filter out targets with < min.cpd.target compounds
targets <- colnames(xplot.bioactive.target)
ncpd <- colSums(xplot.target[,targets])

# reorder for visualization
avg.bioactive <- rowMeans(xplot.bioactive.target)
row.order <- order(avg.bioactive)
col.drop <- colSums(xplot.bioactive.target) == 0

superheat(
  xplot.bioactive.target[row.order, !col.drop],
  pretty.order.cols=TRUE,
  yr=avg.bioactive[row.order],
  yr.plot.type='bar',
  yr.axis.name='Average proportion biactivity',
  yr.axis.name.size=12,
  bottom.label.text.angle=90,
  bottom.label.size=0.75,
  heat.pal=heat.pal,
  heat.pal.values=c(0, seq(0.5, 1, by=0.1)),
  title='Bioactivity by cell line, target'
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

################################################################################
# Compound similarity analysis
################################################################################
# Initialize bioactive compound set
compounds.select <- filter(xgroup, pval == 0)
  
# Initialize cluster score variables
scores <- numeric(length(unique(x$Cell_Line)))
scores.high <- numeric(length(unique(x$Cell_Line)))
names(scores) <- unique(x$Cell_Line)
names(scores.high) <- unique(x$Cell_Line)

distances <- vector('list', length(unique(x$Cell_Line)))
names(distances) <- unique(x$Cell_Line)

# Initialize target distance matrix
target.dist <- dist(apply(xplot.target, MAR=2, l2norm))#, method='binary')
target.dist <- as.matrix(target.dist)


for (cl in unique(x$Cell_Line)) {
  
  xgroup.cl <- filter(xgroup, Cell_Line == cl)
  
  # Filter to bioactive compounds, select cell line
  xc <- filter(x, Cell_Line == cl) %>%
    filter(Compound_ID %in% compounds.select$Compound_ID) %>%
    filter(Compound_ID %in% rownames(target.dist)) %>%
    mutate(ID=str_c(Cell_Line, ', ', Compound_ID)) %>%
    select(-matches('^PC')) %>%
    group_by(Compound_ID) %>% 
    filter(Dose == max(as.numeric(Dose))) %>%
    group_by(Compound_ID, Compound_Usage) %>%
    summarize_if(is.numeric, mean) %>%
    ungroup() %>%
    left_join(xgroup.cl, by='Compound_ID')
  
  # Initialize phenotypic profiles for current cell line
  xfeat <- select(xc, matches('^nonborder')) %>% mutate_all(l2norm)
  colnames(xfeat) <- NULL
  
  # Compute distances in PP and target space
  phenotypic.dist <- dist(xfeat)
  target.dist.f <- as.dist(target.dist[xc$Compound_ID, xc$Compound_ID])
  
  # Compute average bioactivity for all compound pairs
  bioactivity.cl <- xc$DistNorm
  names(bioactivity.cl) <- xc$Compound_ID
  
  compound.pairs <- expand.grid(xc$Compound_ID, xc$Compound_ID)
  
  bioactivity.pairs <- pmin(
    bioactivity.cl[compound.pairs$Var1],
    bioactivity.cl[compound.pairs$Var2]
  )
  
  bioactivity.pairs <- matrix(bioactivity.pairs, nrow=nrow(xc))
  bioactivity.pairs <- bioactivity.pairs[lower.tri(bioactivity.pairs, diag=FALSE)]
  
  p <- data.frame(Target=c(target.dist.f), Phenotype=c(phenotypic.dist)) %>%
    mutate(Bioactivity=bioactivity.pairs) %>%
    mutate(Bioactivity=ifelse(Bioactivity > 2, 2, Bioactivity)) %>%
    filter(Target != 0) %>%
    sample_frac(0.01) %>%
    ggplot(aes(x=Target, y=Phenotype, color=Bioactivity)) +
    geom_point(alpha=0.5) +
    theme_bw() +
    scale_color_material('light-green') +
    ggtitle(cl)
  plot(p)
  
  p <- data.frame(Target=c(target.dist.f), Phenotype=c(phenotypic.dist)) %>%
    mutate(Bioactivity=bioactivity.pairs) %>%
    filter(Bioactivity > 1.5) %>%
    sample_frac(0.01) %>%
    ggplot(aes(x=Target, y=Phenotype, color=Bioactivity)) +
    geom_point() +
    theme_bw() +
    scale_color_material('light-green') +
    ggtitle(cl, 'high bioactivity')
  plot(p)
  
  p <- data.frame(Target=c(target.dist.f), Phenotype=c(phenotypic.dist)) %>%
    filter(Target != 0) %>%
    mutate(Target=rank(Target) / n()) %>%
    mutate(TargetT=cut(Target, seq(0, 1, by=0.25), include.lowest=TRUE)) %>%
    ggplot(aes(x=TargetT, y=Phenotype)) +
    geom_boxplot(alpha=0.6, fill=heat.pal[4]) +
    theme_bw() + 
    ggtitle(cl)
  plot(p)
  
  
  id.drop <- bioactivity.pairs < 1.5
  scores[[cl]] <- cor(c(phenotypic.dist), c(target.dist.f))
  scores.high[[cl]] <- cor(c(phenotypic.dist)[!id.drop], c(target.dist.f)[!id.drop])
  
  print(cl)
  print(scores[[cl]])
  print(scores.high[[cl]])
  
  
  phenotypic.dist.mat <- as.matrix(phenotypic.dist)
  rownames(phenotypic.dist.mat) <- xc$Compound_ID
  distances[[cl]] <- phenotypic.dist.mat
}

# Initialize aggregate distance matrix
phenotypic.dist <- Reduce('+', distances) / length(distances)
target.dist.f <- target.dist[rownames(phenotypic.dist), rownames(phenotypic.dist)]
phenotypic.dist <- as.dist(phenotypic.dist)
target.dist.f <- as.dist(target.dist.f)

data.frame(Target=1 - c(target.dist.f), Phenotype=1 - c(phenotypic.dist)) %>%
  filter(Target != 0) %>%
  sample_frac(0.05) %>%
  ggplot(aes(x=Target, y=Phenotype)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  ggtitle('Aggregate')

data.frame(Target=1 - c(target.dist.f), Phenotype=1 - c(phenotypic.dist)) %>%
  filter(Target != 0) %>%
  mutate(Target=rank(Target) / n()) %>%
  mutate(TargetT=cut(Target, seq(0, 1, by=0.25), include.lowest=TRUE)) %>%
  group_by(TargetT) %>%
  ungroup() %>%
  ggplot(aes(x=TargetT, y=Phenotype)) +
  geom_boxplot(alpha=0.6, fill=heat.pal[4]) +
  theme_bw() + 
  ggtitle('Aggregate')

id.drop <- target.dist.f == 1
aggregate <- cor(c(phenotypic.dist)[!id.drop], c(target.dist.f)[!id.drop])
scores <- c(scores, aggregate)
names(scores)[7] <- 'Aggregate'

data.frame(Cor=scores, CellLine=names(scores)) %>%
  ggplot(aes(x=reorder(CellLine, Cor), y=Cor)) +
  geom_bar(stat='identity') +
  theme_bw()
