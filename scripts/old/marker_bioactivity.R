#+ setup, echo=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(parallel)
library(superheat)
library(ggsci)

analysis.dir <- '~/github/cancer_translator/'
setwd(analysis.dir)
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

################################################################################
# Initialize analysis parameters
################################################################################
intensity.normalize <- TRUE
n.core <- 16
cell.line <- 'A549'
score.threshold <- 4

fin <- str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata')
load(fin)

# Initialize output files
output.dir <- 'results/marker/'
dir.create(output.dir, showWarnings=FALSE)
output.file <- str_c(output.dir, 'bioactivity_', cell.line, '.Rdata') 

# Initialize color palettes
heat.pal <- c('#FFFFFF', pal_material('light-green')(10))
group.pal <- pal_jco()(10)[c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10)]

# Filter to selected cell line
x <- filter(x, Cell_Line == cell.line)

# Clean compound category names
select_category <- function(x) na.omit(x)[1]

xcat.key <- select(x, Compound_ID, Compound_Category) %>%
  group_by(Compound_ID) %>%
  summarize(Compound_Category=select_category(Compound_Category))

################################################################################
# Initialize marker sets evaluated in analysis
################################################################################
# Initialize feature sets for each marker
colnames(x) <- str_replace_all(colnames(x), 'MitoTracker', 'Mitotracker')

markers <- c('HOECHST', 'Mitotracker', 'SYTO14', 'Alexa')
markers.re <- str_c('(', str_c(markers, collapse='|'), ')')
features.marker <- lapply(markers, function(m) str_subset(colnames(x), m))

id.morphology <- !str_detect(colnames(x), markers.re)
features.morphology <- colnames(x)[id.morphology] %>% str_subset('nonborder')

features <- c(features.marker, list(features.morphology))
names(features) <- c(markers, 'morphology')

# Initialize marker set combinations
n.marker <- length(markers)
marker.comb <- lapply(1:n.marker, function(k) combn(markers, k, simplify=FALSE))
marker.comb <- unlist(marker.comb, recursive=FALSE)
n.markers <- sapply(marker.comb, length)

marker.comb.m <- lapply(marker.comb, function(m) c(m, 'morphology'))
marker.comb <- c(marker.comb, marker.comb.m, 'morphology')

# Initialize marker set groups
groups <- str_c(n.markers, ' markers')
groups <- c(groups, str_c(n.markers, ' markers, morphology'), 'morphology')
names(groups) <- sapply(marker.comb, str_c, collapse='_')

#' # Overview
#' This notebook considers the problems of identifying bioactive compounds 
#' based on the cell painting marker set. In addition to the full panel of 
#' markers, we consider whether bioactivity calls require all markers or some 
#' subset.
#'
#' # Key takeaways
#' 
#' 1. Morphology features alone are sufficient to detect bioactivity over the 
#' full compound set. In the reference compound set, we observe a small increase
#' in the proportion of bioactive calls by including marker-specific features
#' in addition to morphology features.
#' 
#' 2. We do not see strong evidence that particular markers are better at
#' identifying bioactivity among particular compound classes in the reference 
#' set. Rather, particular compound classes are harder/easier to detect 
#' bioactivity in. I.e. difference between compound categories dominate 
#' differences between markers.
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
#+ bioactivity, fig.height=8, fig.width=14

################################################################################
# Compute bioactivity scores
################################################################################
setwd(analysis.dir)

# Normalize data matrix
l2norm <- function(z) z / sqrt(sum(z ^ 2))
x <- mutate_if(x, is.numeric, l2norm)

# Compute bioactivity scores for each well by plate
grid <- expand.grid(unique(x$PlateID), 1:length(marker.comb))
xdist <- mclapply(1:nrow(grid), function(i) {
  
  # Initialize features for current marker set
  m <- marker.comb[[grid$Var2[i]]]
  features.m <- unlist(features[m])
  
  feat.drop <- str_subset(colnames(x), '^nonborder') 
  feat.drop <- setdiff(feat.drop, features.m)
  
  # Subset data to select features/plate
  xf <- filter(x, PlateID == grid$Var1[i]) %>%
    dplyr::select(-one_of(feat.drop)) %>%
    dplyr::select(-matches('^PC'))

  # Compute bioactivity
  out <- bioactivity(xf) %>% mutate(Markerset=str_c(m, collapse='_'))
  
  return(out)
}, mc.cores=n.core)

xdist <- rbindlist(xdist)
save(file=output.file, xdist)

#' To assess differential activity across cell lines, we report bioactivity, 
#' measured as the relative distance from DMSO, for each compound/dose 
#' combination. Figures below show this distance for (i) all compounds (ii) 
#' compound categories with more than 20 instances (not including “Other”), 
#' called bioactive wrt > 0 markersets.
#' 
#' **Note:** a bioactivity score (distance relative to DMSO) of 1 below implies 
#' that a treatment (compound/dose pair) is not distinguishable from DMSO in 
#' phenotypic space. We threshold bioactivity scores at `r score.threshold` for 
#' visualization.
#+ bioactivity_grouped, fig.height=8, fig.width=14

# Group bioactivity scores by markerset
n.markerset <- length(unique(xdist$Markerset))
xgroup <- group_by(xdist, Markerset, Compound_ID, Dose_Category) %>% 
  summarize(DistNorm=mean(DistNorm), pval=mean(pval), .groups='drop') %>%
  group_by(Compound_ID, Dose_Category) %>%
  mutate(Count=n()) %>%
  filter(Count == n.markerset) %>%
  ungroup() %>%
  arrange(Compound_ID, Markerset) %>%
  left_join(xcat.key, by='Compound_ID')

################################################################################
# Bioactivity over all compounds
################################################################################
# Format data for heatmap
xplot <- matrix(xgroup$DistNorm, nrow=n.markerset)
rownames(xplot) <- unique(xgroup$Markerset)
colnames(xplot) <- unique(str_c(xgroup$Compound_ID, ', ', xgroup$Dose_Category))
category <- matrix(xgroup$Compound_Category, nrow=n.markerset)[1,]

# Threshold for visualization
xplot.t <- xplot
xplot.t[xplot.t < 1] <- 1
xplot.t[xplot.t > score.threshold] <- score.threshold

# Filter to treatments bioactive in > 0 markersets
id.drop <- colMeans(xplot.t == 1) == 1
xplot.t <- xplot.t[,!id.drop]
category <- category[!id.drop]
xplot <- xplot[,!id.drop]

# Initialize row scores
row.avg <- rowMeans(xplot)
row.order <- order(row.avg)

groups.x <- groups[rownames(xplot)]
row.pal <- group.pal[as.factor(groups.x)]

superheat(
  xplot.t[row.order,],
  pretty.order.cols=TRUE,
  membership.rows=groups.x[row.order],
  left.label='variable',
  left.label.col=row.pal[row.order],
  left.label.text.size=3,
  left.label.size=0.35,
  yr=sort(row.avg),
  yr.obs.col=row.pal[row.order],
  yr.axis.name='Mean relative distance',
  yr.plot.type='bar',
  heat.pal=heat.pal, 
  heat.pal.values=seq(0, 1, by=0.1),
  title='Relative distance from DMSO\nall compounds'
)

################################################################################
# Bioactivity over top compound categories
################################################################################
# Plot bioactivity by cell line, compound category
cat.table <- table(category)
cat.keep <- setdiff(names(cat.table[cat.table >= 20]), 'Others')
id.keep <- category %in% cat.keep
xplot <- xplot[,id.keep]
category <- category[id.keep]

# Initialize column ordering
id.drop <- colMeans(xplot <= 1) == 1
col.order <- order(colMeans(xplot))
col.order <- col.order[!id.drop]

# Initialize row scores
row.avg <- rowMeans(xplot)
row.order <- order(row.avg)

groups.x <- groups[rownames(xplot)]
row.pal <- group.pal[as.factor(groups.x)]

# Threshold for visualization
xplot.t <- xplot
xplot.t[xplot.t < 1] <- 1
xplot.t[xplot.t > score.threshold] <- score.threshold

# Filter to treatments bioactive in > 0 markersets
id.drop <- colMeans(xplot.t == 1) == 1
xplot.t <- xplot.t[,!id.drop]
category <- category[!id.drop]
xplot <- xplot[,!id.drop]

superheat(
  xplot.t[row.order, col.order],
  membership.cols=category[col.order],
  bottom.label.text.angle=90,
  bottom.label.size=0.5,
  bottom.label.text.size=3,
  membership.rows=groups.x[row.order],
  left.label='variable',
  left.label.col=row.pal[row.order],
  left.label.text.size=3,
  left.label.size=0.35,
  yr=sort(row.avg),
  yr.obs.col=row.pal[row.order],
  yr.axis.name='Mean relative distance',
  yr.plot.type='bar',
  heat.pal=heat.pal, 
  heat.pal.values=seq(0, 1, by=0.1),
  title='Relative distance from DMSO\nprevalent categories, bioactive compounds'
)

#' Finally, we report the proportion of bioactive treatment (compound/dose 
#' pairs) calls for reference compounds by markerset. For each markerset, we 
#' consider 2 definitions of bioactivity (i)p-value = 0 (ii) relative distance 
#' from DMSO center > 2. 
#' 
#' **Note:** Distance and p-values for replicates (treatment, cell line) are
#' averaged before calling bioactivity.
#' 
#' The top figures reports the proportion of bioactive calls over all compounds.
#' The bottom figure reports proportion of bioactivity calls within reference 
#' compound categories.
#+ bioactivity_treat, fig.height=10, fig.width=18, warning=FALSE

################################################################################
# Proportion bioactivity
################################################################################
# Compute bioactive relative to each markerset
bioactive.cat <- xdist %>%
  group_by(Markerset, Compound_ID, Compound_Category, Dose_Category, Compound_Usage) %>%
  summarize(pval=mean(pval), DistNorm=mean(DistNorm), count=n(), .groups='drop') %>%
  mutate(BioactiveP=(pval == 0), BioactiveDist=(DistNorm > 2))
  
# P-value = 0 bioactivity
group_by(bioactive.cat, Markerset) %>%
  summarize(PropBioactive=mean(BioactiveP)) %>%
  mutate(Group=groups[Markerset]) %>%
  mutate(PropText=round(PropBioactive, 2)) %>%
  ggplot(aes(x=reorder(Markerset, PropBioactive), y=PropBioactive)) +
  geom_bar(stat='identity', aes(fill=Group)) +
  geom_text(aes(label=PropText), nudge_y=0.01, size=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  ggtitle('Proportion of bioactive calls by markerset, all compounds', 'p-value = 0')

filter(bioactive.cat, Compound_Usage == 'reference_cpd') %>%
  group_by(Markerset) %>%
  summarize(PropBioactive=mean(BioactiveP)) %>%
  mutate(Group=groups[Markerset]) %>%
  mutate(PropText=round(PropBioactive, 2)) %>%
  ggplot(aes(x=reorder(Markerset, PropBioactive), y=PropBioactive)) +
  geom_bar(stat='identity', aes(fill=Group)) +
  geom_text(aes(label=PropText), nudge_y=0.01, size=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  ggtitle('Proportion of bioactive calls by markerset, reference compounds', 'p-value = 0')


# Distance > 2 bioactivity
group_by(bioactive.cat, Markerset) %>%
  summarize(PropBioactive=mean(BioactiveDist)) %>%
  mutate(Group=groups[Markerset]) %>%
  mutate(PropText=round(PropBioactive, 2)) %>%
  ggplot(aes(x=reorder(Markerset, PropBioactive), y=PropBioactive)) +
  geom_bar(stat='identity', aes(fill=Group)) +
  geom_text(aes(label=PropText), nudge_y=0.005, size=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  ggtitle('Proportion of bioactive calls by markerset, all compounds', 'Distance > 2')

filter(bioactive.cat, Compound_Usage == 'reference_cpd') %>%
  group_by(Markerset) %>%
  summarize(PropBioactive=mean(BioactiveDist)) %>%
  mutate(Group=groups[Markerset]) %>%
  mutate(PropText=round(PropBioactive, 2)) %>%
  ggplot(aes(x=reorder(Markerset, PropBioactive), y=PropBioactive)) +
  geom_bar(stat='identity', aes(fill=Group)) +
  geom_text(aes(label=PropText), nudge_y=0.01, size=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=group.pal) +
  ggtitle('Proportion of bioactive calls by markerset, reference compounds', 'Distance > 2')

#' Figures below report bioactivity calls by compound category for reference
#' set compounds. Top figure corresponds to bioactivity defined as p-value = 0,
#' bottom figure as distance > 2.
#+ bioactivity_treat_heat, fig.height=12, fig.width=24, warning=FALSE
# Initialize bioactivity table
bioactive.cat.t <- group_by(bioactive.cat, Compound_Category, Markerset, Compound_Usage) %>%
  summarize(BioactiveP=mean(BioactiveP), BioactiveDist=mean(BioactiveDist), .groups='drop') %>%
  arrange(Markerset, Compound_Category) %>%
  filter(Compound_Usage == 'reference_cpd')

################################################################################
# P-value bioactivity
################################################################################
# Format data for heatmap
n.markerset <- length(unique(bioactive.cat.t$Markerset))
xplot <- matrix(bioactive.cat.t$BioactiveP, nrow=n.markerset, byrow=TRUE)

rownames(xplot) <- unique(bioactive.cat.t$Markerset)
colnames(xplot) <- unique(bioactive.cat.t$Compound_Category)

# Initialize row attributes
row.avg <- rowMeans(xplot)
row.order <- order(row.avg)

groups.x <- groups[rownames(xplot)]
row.pal <- group.pal[as.factor(groups.x)]

# Plot bioactivity heatmap
superheat(
  xplot[row.order,],
  pretty.order.cols=TRUE,
  membership.rows=groups.x[row.order],
  left.label='variable', 
  left.label.text.size=3,
  left.label.size=0.3,
  left.label.col=row.pal[row.order],
  yr=row.avg[row.order],
  yr.plot.type='bar',
  yr.obs.col=row.pal[row.order],
  yr.axis.name='Average proportion bioactive',
  heat.pal=heat.pal
)

################################################################################
# Distance bioactivity
################################################################################
# Format data for heatmap
n.markerset <- length(unique(bioactive.cat.t$Markerset))
xplot <- matrix(bioactive.cat.t$BioactiveDist, nrow=n.markerset, byrow=TRUE)

rownames(xplot) <- unique(bioactive.cat.t$Markerset)
colnames(xplot) <- unique(bioactive.cat.t$Compound_Category)

# Initialize row attributes
row.avg <- rowMeans(xplot)
row.order <- order(row.avg)

groups.x <- groups[rownames(xplot)]
row.pal <- group.pal[as.factor(groups.x)]

# Plot bioactivity heatmap
superheat(
  xplot[row.order,],
  pretty.order.cols=TRUE,
  membership.rows=groups.x[row.order],
  left.label='variable', 
  left.label.text.size=3,
  left.label.size=0.3,
  left.label.col=row.pal[row.order],
  yr=row.avg[row.order],
  yr.plot.type='bar',
  yr.obs.col=row.pal[row.order],
  yr.axis.name='Average proportion bioactive',
  heat.pal=heat.pal
)

