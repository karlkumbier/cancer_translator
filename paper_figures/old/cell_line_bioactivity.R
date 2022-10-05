#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)

intensity.normalize <- TRUE
n.core <- 16
score.threshold <- 4

analysis.dir <- '~/github/cancer_translator/'
setwd(analysis.dir)
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

output.dir <- 'results/cell_line/'
dir.create(output.dir, showWarnings=FALSE)
output.file <- str_c(output.dir, 'bioactivity.Rdata') 

intensity.normalize <- TRUE
fin <- str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata')
load(fin)

# Clean compound category names
select_category <- function(x) na.omit(x)[1]

xcat.key <- select(x, Compound_ID, Compound_Category) %>%
  group_by(Compound_ID) %>%
  summarize(Compound_Category=select_category(Compound_Category))

options(knitr.table.format = function() {
  ifelse(knitr::is_latex_output(), 'latex', 'pipe')
})

# Initialize color palettes
heat.pal <- c('#FFFFFF', pal_material("light-green")(10))
ncell.pal <- pal_material("purple")(10)
compound.pal <- pal_jco()(10)
cell.pal <- pal_nejm()(8)

#' # Overview
#' This notebook considers the problems of identifying bioactive compounds 
#' based on a panel of 6 cell lines: OVCAR4, A549, DU145, ALS-WT, HEPG2, and 
#' 786-0. In addition to individual cell lines, we consider whether bioactivity 
#' calls can be improved by aggregating information across cell lines.
#' 
#' # Key takeaways
#'
#' 1. With the exception of DU145, cell lines show similar patterns of 
#' bioactivity across treatments (compound/dose pairs) for the most bioactive 
#' treatments.
#' 
#' 2. Individual cell lines have a wide range in the proportion of bioactive 
#' treatments identified for both the full compound set and reference set. 
#' OVCAR4 performs best for the reference set (and well for the full compound 
#' set), ALS-WT performs best for the full compound set (but poorly for the 
#' reference set).
#' 
#' 3. For compounds with moderate bioactivity, we observe variation wrt which
#' cell lines detect bioactivity. Compounds that vary across cell lines tend to 
#' fall outside the reference compound & most prevelant categories.
#' 
#' 4. Increasing the number of cell lines leads to a steady increase in the 
#' proportion of bioactive calls among the full compound set (0.34, 1 cell line
#' to 0.57, 6 cell lines). Gains are less pronounced in the reference compound 
#' set (0.6, 1 cell line to 0.67, 6 cell lines).
#' 
#' 5. Trends are robust to how we define bioactivity (distance-based v. 
#' stablity-based).
#' 
#' 6. We do not see strong evidence that particular cell lines are better at
#' identifying bioactivity among particular compound classes in the reference 
#' set. Rather, particular compound classes are harder/easier to detect 
#' bioactivity in. I.e. difference between compound categories dominate 
#' differences between cell lines.
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

# Compute bioactivity scores for each well
xdist <- mclapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(bioactivity(out))
}, mc.cores=n.core)

xdist <- rbindlist(xdist)
save(file=output.file, xdist)

#' The figure below reports normalized distance from DMSO against p-value for
#' each treatment (i.e. compound/dose pairs) evaluated across all cell lines. 
#' Each point corresponds to a single treatment condition.
#+ dist_v_p, fig.height=8, fig.width=12
# Group bioactivity scores by cell line/compound
n.cell.line <- length(unique(x$Cell_Line))
xgroup <- group_by(xdist, Cell_Line, Compound_ID, Dose_Category) %>% 
  summarize(DistNorm=mean(DistNorm), pval=mean(pval), .groups='drop') %>%
  group_by(Compound_ID, Dose_Category) %>%
  mutate(Count=n()) %>%
  filter(Count == n.cell.line) %>%
  ungroup() %>%
  arrange(Compound_ID, Cell_Line) %>%
  left_join(xcat.key, by='Compound_ID')

# Plot pvalue x dist
ggplot(xgroup, aes(x=DistNorm, y=pval, col=Cell_Line)) +
  geom_jitter(alpha=0.5, width=0, height=0.025) +
  facet_wrap(~Cell_Line) +
  geom_vline(xintercept=1, color='grey', linetype=2) +
  theme_bw() +
  theme(legend.position='none') +
  scale_color_manual(values=cell.pal) +
  xlab('Distance relative to DMSO') +
  ylab('p-value')

#' To assess differential activity across cell lines, we report bioactivity, 
#' measured as the relative distance from DMSO, for each compound/dose 
#' combination. Figures below show this distance for (i) all compounds (ii) 
#' compound categories with more than 20 instances (not including "Other"). 
#' 
#' **Note:** a bioactivity score (distance relative to DMSO) of 1 below implies 
#' that a treatment (compound/dose pair) is not distinguishable from DMSO in 
#' phenotypic space. We threshold bioactivity scores at `r score.threshold` for 
#' visualization.
#+ bioactivity_dist, fig.height=12, fig.width=18
# Format cell x compound bioactivity for visualization
xplot <- matrix(xgroup$DistNorm, nrow=n.cell.line)
rownames(xplot) <- unique(xgroup$Cell_Line)
colnames(xplot) <- unique(str_c(xgroup$Compound_ID, ', ', xgroup$Dose_Category))
category <- matrix(xgroup$Compound_Category, nrow=n.cell.line)[1,]

# Threshold for visualization
xplot.t <- xplot
xplot.t[xplot.t < 1] <- 1
xplot.t[xplot.t > score.threshold] <- score.threshold

# Filter to compounds bioactive in > 0 cell lines
id.drop <- colMeans(xplot.t == 1) == 1
xplot.t <- xplot.t[,!id.drop]
category <- category[!id.drop]
xplot <- xplot[,!id.drop]

superheat(
  xplot.t,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  heat.pal=heat.pal, 
  heat.pal.values=seq(0, 1, by=0.1),
  title='Bioactivity by cell line, compound\nall compounds'
)

# Plot bioactivity by cell line, compound category
cat.table <- table(category)
cat.keep <- setdiff(names(cat.table[cat.table >= 20]), 'Others')
id.keep <- category %in% cat.keep

# Initialize column ordering
col.order <- order(colMeans(xplot[,id.keep]))

superheat(
  xplot.t[,id.keep][,col.order],
  pretty.order.rows=TRUE,
  membership.cols=category[id.keep][col.order],
  bottom.label.text.angle=90,
  bottom.label.size=0.75,
  heat.pal=heat.pal, 
  heat.pal.values=seq(0, 1, by=0.1),
  title='Bioactivity by cell line, compound\nprevalent compound categories'
)

#' Finally, we report the proportion of "bioactive" treatments (compound/dose 
#' pairs)  by cell line and across subsets of the full cell line set. For an 
#' individual cell line, we consider 2 definitions of bioactivity (i) 
#' p-value = 0 (ii) relative distance from DMSO center > 2. For cell line sets, 
#' we define any treatment with that is bioactive in at least one cell line as 
#' bioactive (i.e. minimum p-value = 0, maximum relative distance > 2).
#' 
#' **Note:** Distance and p-values for replicates (treatment, cell line) are
#' averaged before calling bioactivity.
#' 
#' **Note:** By definition, the proportion of bioactive treatments/compounds for 
#' a set of cell lines will be greater than or equal to the proportion of 
#' bioactive treatments for any subset of cell lines. 
#+ bioactivity_treat, fig.height=10, fig.width=18, warning=FALSE
# Initialize cell line sets
cell.lines <- unique(xdist$Cell_Line)     

cell.sets <- lapply(1:n.cell.line, function(k) {
  combn(cell.lines, k, simplify=FALSE)
})

cell.sets <- unlist(cell.sets, recursive=FALSE)

# Compute bioactive by treatment relative to each cell line set
bioactive.cell <- lapply(cell.sets, function(s) {
  filter(xdist, Cell_Line %in% s) %>%
    group_by(Compound_ID, Dose_Category, Compound_Usage, Cell_Line) %>%
    summarize(pval=mean(pval), DistNorm=mean(DistNorm), .groups='drop') %>%    
    group_by(Compound_ID, Dose_Category, Compound_Usage) %>%
    summarize(pval=min(pval), DistNorm=max(DistNorm), .groups='drop') %>%
    mutate(BioactiveP=(pval == 0), BioactiveDist=(DistNorm > 2)) %>%
    mutate(Cell_Line=str_c(s, collapse=', ')) %>%
    mutate(Ncells=length(s))
})

# P-value = 0 bioactivity
rbindlist(bioactive.cell) %>%
  group_by(Cell_Line) %>%
  summarize(PropBioactive=mean(BioactiveP), Ncells=mean(Ncells)) %>%
  mutate(PropText=round(PropBioactive, 2)) %>%
  ggplot(aes(x=reorder(Cell_Line, PropBioactive), y=PropBioactive)) +
  geom_bar(stat='identity', aes(fill=Ncells)) +
  geom_text(aes(label=PropText), nudge_y=0.01, size=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  ylab('Proportion bioactive treatments') +
  ggtitle('Bioactivity proportion by cell line set, all compounds', 'p-value = 0')

rbindlist(bioactive.cell) %>%
  filter(Compound_Usage == 'reference_cpd') %>%
  group_by(Cell_Line) %>%
  summarize(PropBioactive=mean(BioactiveP), Ncells=mean(Ncells)) %>%
  mutate(PropText=round(PropBioactive, 2)) %>%
  ggplot(aes(x=reorder(Cell_Line, PropBioactive), y=PropBioactive)) +
  geom_bar(stat='identity', aes(fill=Ncells)) +
  geom_text(aes(label=PropText), nudge_y=0.01, size=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  ylab('Proportion bioactive treatments') +
  ggtitle('Bioactivity proportion by cell line set, reference compounds', 'p-value = 0')

# Distance > 2 bioactivity
rbindlist(bioactive.cell) %>%
  group_by(Cell_Line) %>%
  summarize(PropBioactive=mean(BioactiveDist), Ncells=mean(Ncells)) %>%
  mutate(PropText=round(PropBioactive, 2)) %>%
  ggplot(aes(x=reorder(Cell_Line, PropBioactive), y=PropBioactive)) +
  geom_bar(stat='identity', aes(fill=Ncells)) +
  geom_text(aes(label=PropText), nudge_y=0.005, size=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  ylab('Proportion bioactive treatments') +
  ggtitle('Bioactivity proportion by cell line set, all compounds', 'Distance > 2')

rbindlist(bioactive.cell) %>%
  filter(Compound_Usage == 'reference_cpd') %>%
  group_by(Cell_Line) %>%
  summarize(PropBioactive=mean(BioactiveDist), Ncells=mean(Ncells)) %>%
  mutate(PropText=round(PropBioactive, 2)) %>%
  ggplot(aes(x=reorder(Cell_Line, PropBioactive), y=PropBioactive)) +
  geom_bar(stat='identity', aes(fill=Ncells)) +
  geom_text(aes(label=PropText), nudge_y=0.01, size=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  ylab('Proportion bioactive treatments') +
  ggtitle('Bioactivity proportion by cell line set, reference compounds', 'Distance > 2')


#' Figures below report bioactivity by compound category, filtered to reference 
#' compounds. Bioactivity proportions include all doses of a given compound.
#+ bioactivity_treat_heat, fig.height=12, fig.width=24, warning=FALSE, message=FALSE
# Compute bioactive relative to each cell line set
bioactive.cat <- lapply(cell.sets, function(s) {
  filter(xdist, Cell_Line %in% s) %>%
    group_by(Compound_ID, Dose_Category, Compound_Category, Compound_Usage, Cell_Line) %>%
    summarize(pval=mean(pval), DistNorm=mean(DistNorm), .groups='drop') %>%    
    group_by(Compound_ID, Dose_Category, Compound_Category, Compound_Usage) %>%
    summarize(pval=min(pval), DistNorm=max(DistNorm), .groups='drop') %>%
    mutate(BioactiveP=(pval == 0), BioactiveDist=(DistNorm >= 2)) %>%
    mutate(Cell_Line=str_c(s, collapse=', ')) %>%
    mutate(Ncells=length(s))
})

# Initialize bioactivity table
bioactive.cat.t <- rbindlist(bioactive.cat) %>% 
  filter(Compound_Usage == 'reference_cpd') %>%
  group_by(Compound_Category, Compound_Usage, Cell_Line) %>%
  summarize(BioactiveP=mean(BioactiveP), BioactiveDist=mean(BioactiveDist), Count=n(), Ncells=mean(Ncells)) %>%
  arrange(Cell_Line, Compound_Category) %>%
  ungroup()

################################################################################
# Bioactivity p-value = 0
################################################################################
# Format data for heatmap
n.cell.line <- length(unique(bioactive.cat.t$Cell_Line))
xplot <- matrix(bioactive.cat.t$BioactiveP, nrow=n.cell.line, byrow=TRUE)

rownames(xplot) <- unique(bioactive.cat.t$Cell_Line)
colnames(xplot) <- unique(bioactive.cat.t$Compound_Category)

# Initialize row attributes
ncell <- select(bioactive.cat.t, Cell_Line, Ncells) %>% 
  distinct() %>%
  arrange(match(Cell_Line, rownames(xplot)))

ncell.col <- ncell.pal[-1][floor(ncell$Ncells * 1.5)]

avg.bioactive <- rowMeans(xplot)
row.order <- order(avg.bioactive)

superheat(
  xplot[row.order,],
  pretty.order.cols=TRUE,
  membership.rows=ncell$Ncells[row.order],
  left.label='variable', 
  left.label.text.size=3,
  left.label.col=ncell.col[row.order],
  yr=avg.bioactive[row.order],
  yr.plot.type='bar',
  yr.axis.name='Average proportion bioactive',
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1)
)

################################################################################
# Bioactivity distance > 2
################################################################################
# Format data for heatmap
n.cell.line <- length(unique(bioactive.cat.t$Cell_Line))
xplot <- matrix(bioactive.cat.t$BioactiveDist, nrow=n.cell.line, byrow=TRUE)

rownames(xplot) <- unique(bioactive.cat.t$Cell_Line)
colnames(xplot) <- unique(bioactive.cat.t$Compound_Category)

# Initialize row attributes
ncell <- select(bioactive.cat.t, Cell_Line, Ncells) %>% 
  distinct() %>%
  arrange(match(Cell_Line, rownames(xplot)))

ncell.col <- ncell.pal[-1][floor(ncell$Ncells * 1.5)]

avg.bioactive <- rowMeans(xplot)
row.order <- order(avg.bioactive)

superheat(
  xplot[row.order,],
  pretty.order.cols=TRUE,
  membership.rows=ncell$Ncells[row.order],
  left.label='variable', 
  left.label.text.size=3,
  left.label.col=ncell.col[row.order],
  yr=avg.bioactive[row.order],
  yr.plot.type='bar',
  yr.axis.name='Average proportion bioactive',
  heat.pal=heat.pal,
  heat.pal.values=seq(0, 1, by=0.1)
)

# Bioactivity table by cell line/category
setwd(analysis.dir)
write.csv(file='results/bioactivity_cell_category.csv', bioactive.cat.t, quote=FALSE)
