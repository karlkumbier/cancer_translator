#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)
library(patchwork)
library(hrbrthemes)
library(umap)
library(readr)
library(see)
library(twosamples)


theme_set(
  theme_ipsum(
    axis_title_size=18, 
    strip_text_size=18, 
    axis_text_size=14,
    base_size=18, 
    base_family='sans'
  )
)

################################################################################
# Analysis parameters
################################################################################
# Initialize analysis parameters
intensity.normalize <- TRUE
n.core <- 16
min.cat <- 5
n.subsample <- 100
save.fig <- TRUE

# Initialize color palettes
heat.pal <- viridis::viridis(10)
heat.pal.blue <- c('#FFFFFF', pal_material("blue")(10))
heat.pal.signed <- RColorBrewer::brewer.pal(9, 'RdBu')
col.pal <- pal_nejm()(8)

################################################################################
# Load dataset
################################################################################
setwd('~/github/cancer_translator/')
source('scripts/utilities.R')

# Load KS profile data
data.dir <- 'data/screens/LH_CDC_1/'
load(str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata'))
x <- dplyr::select(x, -matches('^PC'))

# Initialize # of cell lines
n.cell.line <- length(unique(x$Cell_Line))

# Initialize output directories
output.dir <- 'results/cell_line/'
fig.dir <- 'results/figures/fig1/'

dir.create(output.dir, showWarnings=FALSE)
output.file <- str_c(output.dir, 'bioactivity.Rdata') 

# Clean compound categories, taking maximum vote across replicates
x <- group_by(x, Compound_ID) %>%
  mutate(Dose=as.numeric(Dose)) %>%
  mutate(Compound_Category=select_category(Compound_Category)) %>%
  ungroup()

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

# Load compound target table
x.target <- 'data/Bioactivity_broad_label/' %>%
  str_c('bioactivity_broad_institute_label_20220108.csv') %>%
  fread() %>%
  dplyr::select(Compound_ID, MOA_broad_institute, Target_broad_institute) %>%
  distinct() %>%
  dplyr::rename(Compound_Category=MOA_broad_institute)

x <- left_join(x, x.target, by='Compound_ID')

# Initialize plate x cell line key
plate.key <- dplyr::select(x, PlateID, Cell_Line) %>% distinct()

#' # Category summary
#' Figures below summarize compound counts by category across the full compound 
#' library as well as categories with the largest number of compounds 
#' represented.
#+ count_summaries, fig.height=8, fig.width=12
xcat.tab <- filter(x, !is.na(Compound_Category)) %>%
  dplyr::select(Compound_ID, Compound_Category) %>%
  distinct() %>%
  group_by(Compound_Category) %>%
  summarize(Count=n())

filter(xcat.tab, Count > 1) %>%
  filter(Compound_Category != '') %>%
  ggplot(aes(x=Count)) +
  geom_histogram() +
  ggtitle('Distribution of compound counts per category') +
  xlab('# compounds') +
  theme(text=element_text(size=24))

# Summarize prevalent compound categories
xcat.tab.top <- filter(xcat.tab, Count >= min.cat, Compound_Category != '')
cat.keep <- xcat.tab.top$Compound_Category

xcat.tab.top %>%
  ggplot(aes(x=reorder(Compound_Category, Count), y=Count)) +
  geom_bar(stat='identity', position=position_dodge(preserve="single")) +
  geom_text(aes(label=Count), nudge_y=0.025) +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Compounds per category counts', str_c('top categories')) +
  scale_y_log10()

# Write table for supplement
category.table <- dplyr::select(x, Compound_ID, Compound_Category, Dose, Cell_Line) %>%
  dplyr::rename(MOA=Compound_Category) %>%
  group_by(Compound_ID, MOA, Dose) %>%
  summarize(N_Cell_Line=length(unique(Cell_Line))) %>%
  mutate(MOA=str_replace_all(MOA, ',', ';')) %>%
  arrange(Compound_ID)

fout <- str_c(output.dir, 'compound_library.csv')
write_excel_csv(file=fout, x=category.table)

#' # Overview
#' This notebook considers the problems of bioactivity detection based on 
#' a panel of 6 cell lines: OVCAR4, A549, DU145, ALS-WT, HEPG2, and 786-0.
#' 
#' # Bioactivity scoring
#' We compute bioactivity scores within eact plate as follows:
#'
#' 1. Compute the center of the DMSO point cloud  (i.e. median feature values 
#' across each DMSO well).
#' 2. Generate a null distribution corresponding to the weighted l2 distance 
#' between each DMSO well and the DMSO point cloud center. 
#' 3. Define bioactivity as any well that is further than 1 IQR from the 
#' DMSO centroid
#' 
#' 
#' **Note:** Bioactivity scores are computed within plate and subsequently 
#' aggregated when looking at treatments (i.e. dose/compound pairs).
#' 
#' **Note:** We assess bioactivity relative to the highest dose level of a given 
#' compound. Some compounds appear under multiple usage categories (e.g. 
#' reference_cpd & positive_ctrl_cpd). In such cases, we compute bioactivity 
#' at the maximum dose over all use categories.
#+ bioactivity, fig.height=8, fig.width=14, echo=FALSE

# Mean-center features by cell line for subsetquent
x <- correlation_weight(x, normalize=TRUE)
  
################################################################################
# Initialize parameters for bioactivity analysis
################################################################################
# Initialize DMSO distance summary function
null_summary <- function(x) {
 iqr.thresh <- IQR(x)
 return(max(x[x < median(x) + iqr.thresh]))
}

################################################################################
# Summarize null distribution
################################################################################
xnull <- lapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(null_dist(out))
})

xthresh <- reshape2::melt(xnull) %>%
  group_by(L1) %>%
  summarize(value=null_summary(value)) %>%
  mutate(PlateID=unique(x$PlateID)[L1]) %>%
  left_join(plate.key, by='PlateID')

reshape2::melt(xnull) %>%
  mutate(PlateID=unique(x$PlateID)[L1]) %>%
  left_join(plate.key, by='PlateID') %>%
  ggplot(aes(x=reorder(PlateID, value), y=value, fill=Cell_Line, col=Cell_Line)) +
  geom_boxplot(alpha=0.6) +
  geom_point(data=xthresh, shape=8) +
  theme(axis.text.x=element_blank()) +
  scale_fill_nejm() +
  scale_color_nejm() +
  facet_wrap(~Cell_Line, scales='free_x') +
  theme(legend.position='none') +
  ylab('Distance from centroid') +
  xlab(NULL) +
  ggtitle('DMSO wells')

################################################################################
# Compute bioactivity scores
################################################################################
xdist <- lapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(bioactivity(out, null_summary))
})

names(xdist) <- unique(x$PlateID)

rbindlist(xdist) %>%
  filter(!str_detect(Compound_Usage, 'query')) %>%
  filter(!str_detect(Compound_Usage, 'reference')) %>%
  ggplot(aes(x=reorder(PlateID, DistNorm), y=DistNorm)) +
  geom_boxplot(alpha=0.6, aes(fill=Compound_Usage, col=Compound_Usage)) +
  theme(axis.text.x=element_blank()) +
  scale_fill_d3() +
  scale_color_d3() +
  facet_wrap(~Cell_Line, scales='free_x') +
  ylab('Distance from DMSO centroid')

rbindlist(xdist) %>%
  mutate(Type=ifelse(Compound_ID == 'DMSO', 'DMSO', 'Other')) %>%
  ggplot(aes(x=reorder(Cell_Line, DistNorm), y=DistNorm)) +
  geom_violinhalf(scale='width', aes(fill=Type, color=Type)) +
  xlab(NULL) +
  scale_fill_d3(name='') +
  scale_color_d3(name='') +
  ylab('Distance from DMSO centroid') +
  theme(legend.position=c(0.1, 0.9))


#' ### Fig 1a UMAP projection by cell line
#' UMAP projections of phenotypic profiles by cell line. Points colored dark 
#' blue represent DMSO-treated wells. Points colored orange and light blue 
#' represent bioactive and inactive calls respsectively.
#+ umap_projections, fig.height=8, fig.width=12

# Select one plate per cell line for visualization
xumap <- lapply(unique(x$Cell_Line), function(cl) {
  
  xdist.merge <- rbindlist(xdist)  %>% 
    filter(Cell_Line == cl) %>%
    mutate(ID=str_c(PlateID, '_', WellID))
  
  # Compute pairwise distances
  xp <- filter(x, Cell_Line == cl)
  xp.feat <- dplyr::select(xp, matches('^non')) 
  
  # Filter outliers
  umap.config <- umap.defaults
  xumap <- umap(xp.feat, config=umap.config)
  colnames(xumap$layout) <- c('UMAP1', 'UMAP2')
  
  out <- data.frame(xumap$layout) %>% 
    mutate(ID=str_c(xp$PlateID, '_', xp$WellID)) %>%
    left_join(xdist.merge, by='ID')
  
  return(out)
})

rbindlist(xumap) %>%
  filter(Compound_ID == 'DMSO') %>%
  group_by(Compound_ID) %>%
  mutate(ColID=as.factor(ColID)) %>%
  ggplot(aes(x=UMAP1, y=UMAP2, col=ColID)) +
  geom_jitter(height=0.25, width=0.25, alpha=0.5) +
  facet_wrap(~Cell_Line) +
  scale_color_nejm() +
  theme(legend.position='none')

p <- rbindlist(xumap) %>%
  mutate(Group=ifelse(Compound_ID == 'DMSO', 'DMSO', 'Non-bioactive')) %>%
  mutate(Group=ifelse(DistNorm > 1, 'Bioactive', Group)) %>%
  mutate(Group=factor(Group, levels=c('DMSO', 'Bioactive', 'Non-bioactive'))) %>%
  mutate(DistNorm=ifelse(Compound_ID == 'DMSO', NA, DistNorm)) %>%
  filter(Cell_Line %in% c('OVCAR4', 'HEPG2')) %>%
  mutate(Compound_Category=ifelse(Compound_ID == 'DMSO', 'DMSO', Compound_Category)) %>%
  filter(!is.na(Compound_Category)) %>%
  mutate(Shape=Compound_Category == 'glucocorticoid receptor agonist') %>%
  mutate(Size=ifelse(Shape, 1, 0.75)) %>%
  group_by(Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose)) | Compound_ID == 'DMSO') %>%
  ggplot(aes(x=UMAP1, y=UMAP2, col=Group, shape=Shape, size=Size, alpha=Size)) +
  geom_jitter(height=0.05, width=0.05) +
  facet_wrap(~Cell_Line, ncol=1) +
  scale_color_jama() +
  labs(col=NULL) +
  scale_shape_manual(breaks=c('TRUE', 'FALSE'), values=c(8, 20)) +
  scale_alpha(range=c(0.3, 1)) +
  scale_size(range=c(1, 1.5)) +
  guides(shape='none') +
  guides(size='none') +
  guides(alpha='none') +
  theme(legend.position=c(0.15, 0.95))

if (save.fig) pdf(str_c(fig.dir, 'umap.pdf'), height=12, width=8)
plot(p)  
if(save.fig) dev.off()

#' ### Fig 1b. Bioactivity by cell line. 
#' Bioactivity calls by compound, cell line. Compounds with p-value = 0 are 
#' defined as bioactive, while compounds with p-value > 0 are defined as 
#' inactive.
#' 
#' **Note:** treatment (compound/dose pairs) replicates are grouped and p-values 
#' averaged before calling bioactivity. Thus, bioactive treatments correspond to 
#' those with p-value = 0 across all replicates.
#+ bioactivity_p, fig.height=12, fig.width=18
# Group bioactivity scores by treatment, cell line, usage
xdmso <- rbindlist(xdist) %>%
  filter(Compound_ID == 'DMSO') %>%
  mutate(Compound_Category='DMSO') %>%
  dplyr::select(Cell_Line, Compound_ID, Compound_Category, DistNorm)

xgroup <- rbindlist(xdist) %>%
  filter(Compound_ID != 'DMSO') %>%
  group_by(Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup() %>%
  group_by(Cell_Line, Compound_ID, Compound_Category) %>% 
  summarize(DistNorm=mean(DistNorm), .groups='drop') %>%
  arrange(Compound_ID, Cell_Line)
  

rbind(xdmso, xgroup) %>%
  mutate(Type=ifelse(Compound_ID == 'DMSO', 'DMSO', 'Other')) %>%
  ggplot(aes(x=DistNorm)) +
  geom_density(aes(fill=Type, color=Type), alpha=0.5) +
  scale_fill_d3(name='') +
  scale_color_d3(name='') +
  xlab('Distance from DMSO centroid') +
  ylab(NULL) +
  facet_wrap(~Cell_Line) +
  theme(legend.position=c(0.9, 0.9))

#' # Stability analysis on distributions
#' We compute the difference in these ECDF using signed CVM. We resample 
#' compounds with probability proortional to given weight — for prioritizing MOA.
xmerge <- rbind(xdmso, xgroup)

bootstrap_bioactive_difference <- function(x, wt, nbs=5) {
  # Compute bioactive difference by cell line
  bdiff <- replicate(nbs, bioactive_difference(bootstrap_compounds(x, wt)))
  # Question: how do distinct/non-distinct values vary?
}

bioactive_difference <- function(x) {
  # Compute bioactive difference by cell line
  cell.lines <- unique(x$Cell_Line)
  distance <- sapply(cell.lines, function(cl) {
    xcl <- filter(x, Cell_Line == cl)
    id.dmso <- xcl$Compound_ID == 'DMSO'
    cvm_stat(xcl$DistNorm[id.dmso], xcl$DistNorm[!id.dmso], power=1)
  })
}

bootstrap_compounds <- function(x, wt) {
  # Weighted bootstrap sample of compounds, weights = MOA sampling probability
  xf <- filter(x, Compound_ID != 'DMSO')
  xdmso <- filter(x, Compound_ID == 'DMSO')
  
  n <- nrow(xf) / length(unique(xf$Cell_Line))
  wt <- wt / sum(wt)
  
  # Initialize ID by category table
  xcat.tab <- dplyr::select(xf, Compound_ID, Compound_Category) %>% distinct()
  
  id.bs.sample <- lapply(names(wt), function(ctg) {
    filter(xcat.tab, Compound_Category == ctg) %>% 
      sample_n(wt[ctg] * n, replace=TRUE) %>%
      dplyr::select(Compound_ID) %>%
      unlist()
  })
  
  id.bs.sample <- unlist(id.bs.sample) %>% unname
  
  # Initialize bootstrap sampled data table
  x.bs.sample <- lapply(unique(xf$Cell_Line), function(cl) {
    xfc <- filter(xf, Cell_Line == cl)
    return(xfc[fastmatch::fmatch(id.bs.sample, xfc$Compound_ID),])
  })
  
  return(rbind(rbindlist(x.bs.sample), xdmso))
}

# Explanation:
# 1. We compute reference distributions in phenotypic space
# 2. We score bioactivity / compound similarity relative to different reference 
#   distributions.
# 3. Uncertainty / stability analysis based on weighted bootstrap sampling, 
#   weights define the priority for different MOAs

# TODO (1) make this figure for Compound similarity
# TODO (2) run stability analysis

# TODO: cell line by category matrix based on distributional difference

# Format cell x compound bioactivity for visualization
xplot <- matrix(as.numeric(xgroup$DistNorm > 1), nrow=n.cell.line)
rownames(xplot) <- unique(xgroup$Cell_Line)
colnames(xplot) <- unique(xgroup$Compound_ID)
category <- matrix(xgroup$Compound_Category, nrow=n.cell.line)[1,]
category[is.na(category)] <- ''

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
  heat.pal=heat.pal.blue,
  yr=prop.bioactive[bioactive.order],
  yr.plot.type='bar',
  yr.axis.name='Proportion\nbioactive',
  yr.axis.name.size=14,
  heat.pal.values=seq(0, 1, by=0.1),
  title='Bioactivity by cell line, compound'
)

#' ### Fig 1c. Bioactivity by cell line, compound category.
#' Bioactivity calls by cell line, compound category. Treatments (compound dose 
#' pairs) with normalized distance from centroid > 1 are defined as bioactive, 
#' while remaining treatments are defined as non-bionactive. Top: bioactivity 
#' calls for each treatment, bottom: proportion of bioactivity calls within each 
#' compound category. Compound categories with fewer than `R min.cat` treatments 
#' are dropped for visualization.
#' 
#' **Note:** treatment replicates are grouped and distances averaged before 
#' calling bioactivity.
#+ bioactivity_category, fig.height=12, fig.width=18, warning=FALSE
# Plot bioactivity by cell line, compound category
cat.keep <- filter(xcat.tab, Count >= 10)$Compound_Category
cat.keep <- setdiff(cat.keep, '')

xplot.cat <- sapply(cat.keep, function(ctg) {
  rowMeans(xplot[,which(category == ctg)])
})

# reorder for visualization
xplot.cat <- xplot.cat[,order(colMeans(xplot.cat))]
col.order <- order(xplot.cat[head(bioactive.order, 1),])
n.compounds <- xcat.tab.top$Count
colnames(xplot.cat) <- str_remove_all(colnames(xplot.cat), ',.*$')

superheat(
  xplot.cat[bioactive.order, col.order],
  yt=n.compounds[col.order],
  yt.plot.type='bar',
  yt.axis.name='# compounds',
  yt.axis.name.size=16,
  bottom.label.text.angle=90,
  bottom.label.size=0.75,
  heat.pal=heat.pal, 
  heat.pal.values=seq(0, 1, by=0.1),
  title='Bioactivity by cell line, compound category'
)

#' ### Fig. 1d. Bioactivity generalist
#' Optimal cell line sets over all compounds — equivalent to setting category 
#' weights proportional to # of compounds in each category.
#+ bioactivity_generalist
# Initialize category weights
specialist.cat <- c('CDK inhibitor', 'MEK inhibitor', 'JAK inhibitor')
wt.gen <- table(category)
wt.spec <- setNames(rep(0, length(wt.gen)), names(wt.gen))
wt.spec[specialist.cat] <- 1

# Compute optimal cell line set by budget
x.generalist <- lapply(1:3, function(k) {
  opt_bioactive(xplot, category, wt.gen, k)
})

x.specialist <- lapply(1:3, function(k) {
  opt_bioactive(xplot, category, wt.spec, k)
})

# Initialize table for figures
xopt.gen <- sapply(x.generalist, function(z) z$opt)
xopt.spec <- sapply(x.specialist, function(z) z$opt)

xopt <- data.frame(CellSet=c(names(xopt.gen), names(xopt.spec))) %>%
  mutate(Score=c(xopt.gen, xopt.spec)) %>%
  mutate(K=c(1:3, 1:3)) %>%
  mutate(Type=c(
    rep('Generalist', length(xopt.gen)), 
    rep('Specialist', length(xopt.gen))
  )) %>%
  mutate(Score=round(Score, 3)) %>%
  mutate(K=as.factor(K))

xdist.gen <- unlist(lapply(x.generalist, function(z) z$score))
xdist.spec <- unlist(lapply(x.specialist, function(z) z$score))

xdist <- data.frame(CellSet=c(names(xdist.gen), names(xdist.spec))) %>%
  mutate(Score=c(xdist.gen, xdist.spec)) %>%
  mutate(K=str_count(CellSet, ',') + 1) %>%
  mutate(Type=c(
    rep('Generalist', length(xdist.gen)), 
    rep('Specialist', length(xdist.spec))
  )) %>%
  mutate(Score=round(Score, 3)) %>%
  mutate(K=as.factor(K))

# Plot optimal bioactivity scores
p <- ggplot(xdist, aes(x=K, y=Score, col=Type)) +
  geom_boxplot(aes(fill=Type), alpha=0.7) +
  geom_point() +
  geom_point(data=xopt, shape=8, size=3) +
  geom_text(data=xopt, aes(label=CellSet), nudge_y=0.025, size=5) +
  xlab('# cell lines') +
  theme(axis.text.x=element_text(angle=90), text=element_text(size=24)) +
  ylab('Proportion bioactive') +
  scale_color_nejm() +
  scale_fill_nejm() +
  theme(legend.position='none') +
  facet_wrap(~Type, ncol=1)

if (save.fig) pdf(str_c(fig.dir, 'selection.pdf'), height=12, width=8)
plot(p)  
if(save.fig) dev.off()
