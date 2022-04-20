#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)
library(RColorBrewer)
library(patchwork)
library(hrbrthemes)
library(umap)
library(rprofiler)
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
# Initialize constant parameters to be used in analysis
################################################################################
# Initialize analysis parameters
intensity.normalize <- TRUE
n.core <- 16
min.cat <- 5
ncell.thresh <- 50
n.bootstrap <- 25
umap.cell <- 'A549'
save.fig <- TRUE

# Visualization parameters
heat.pal <- viridis::viridis(10)
heat.pal.signed <- brewer.pal(9, 'RdBu')
col.pal <- pal_nejm()(8)
thresh.plot <- 0.3

# Function for normalizing distances within plate relative to DMSO + 
# thresholding bioactive/inactive compounds
null_summary <- function(x) {
  iqr.thresh <- IQR(x)
  return(max(x[x < median(x) + iqr.thresh]))
}

# Function for comparing ECDFs
summary_fun <- function(x, y) cvm_stat(x, y, power=1)

################################################################################
# Load dataset
################################################################################
analysis.dir <- '~/github/cancer_translator/'
setwd(analysis.dir)
source('scripts/utilities.R')

# Load KS profile data
data.dir <- 'data/screens/LH_CDC_1/'
load(str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata'))
x <- dplyr::select(x, -matches('^PC'))

# Initialize # of cell lines
n.cell.line <- length(unique(x$Cell_Line))

# initialize output directories
output.dir <- '~/github/cancer_translator/results/distance/'
fig.dir <- '~/github/cancer_translator/results/figures/fig2/'
dir.create(output.dir, showWarnings=FALSE)

# Clean compound categories
x <- group_by(x, Compound_ID) %>%
  mutate(Compound_Category=select_category(Compound_Category)) %>%
  ungroup()

# Filter to compound/dose combinations evaluated in all cell lines
x.treat <- dplyr::select(x, Cell_Line, Compound_ID, Dose_Category, Compound_Usage) %>%
  distinct() %>%
  group_by(Compound_ID, Dose_Category, Compound_Usage) %>% 
  count() %>%
  filter(n == n.cell.line) %>%
  mutate(ID=str_c(Compound_ID, Dose_Category, Compound_Usage))

x <- mutate(x, ID=str_c(Compound_ID, Dose_Category, Compound_Usage)) %>%
  filter(ID %in% x.treat$ID) %>%
  filter(!is.na(Compound_Usage)) %>%
  dplyr::rename(Category_Vendor=Compound_Category)

# Load cleaned MOA table
xmoa <- fread('data/Unique_Screened_Compound_MOA.csv') %>%
  dplyr::rename(Compound_ID=Compound_ID_1) %>%
  dplyr::select(-InchIKeys_All) %>%
  dplyr::rename(Category=MOA)

# Clean compound ID names 
id2 <- match(x$Compound_ID, xmoa$Compound_ID_2)
x$Compound_ID[!is.na(id2)] <- xmoa$Compound_ID[na.omit(id2)]

id3 <- match(x$Compound_ID, xmoa$Compound_ID_3)
x$Compound_ID[!is.na(id3)] <- xmoa$Compound_ID[na.omit(id3)]
xmoa <- dplyr::select(xmoa, -Compound_ID_2, Compound_ID_3)

# Merge data with MOA table
x <- left_join(x, xmoa, by='Compound_ID') %>%
  mutate(Category=ifelse(Compound_ID == 'DMSO', 'DMSO', Category))

# Filter to compound/dose combinations evaluated in all cell lines
x.treat <- select(x, Cell_Line, Compound_ID, Dose_Category, Compound_Usage) %>%
  distinct() %>%
  group_by(Compound_ID, Dose_Category, Compound_Usage) %>% 
  count() %>%
  filter(n == n.cell.line) %>%
  mutate(ID=str_c(Compound_ID, Dose_Category, Compound_Usage))

x <- mutate(x, ID=str_c(Compound_ID, Dose_Category, Compound_Usage)) %>%
  filter(ID %in% x.treat$ID) %>%
  filter(!is.na(Compound_Usage))

# Initialize plate x cell line key
plate.key <- dplyr::select(x, PlateID, Cell_Line) %>% distinct()

#' # Overview
#' This notebook considers the problem of compound similarity based on 
#' phenotypic profiles. Compounds are first filtered based on 
#' bioactivity — as specified in the `bioactivity` notebook. At a high level, 
#' bioactivity is based on the distance (evaluated over phenotypic profiles) 
#' between a target well and a collection of DMSO control wells. 
#+ bioactivity, fig.height=8, fig.width=12, echo=FALSE
# Mean-center + correlation weight features by cell line for subsequent analysis
x <- correlation_weight(x, TRUE)

# Compute distance to DMSO centroid
xdist <- lapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(dmso_distance(out, null_summary))
})

#' ### Fig. S3 Cell count v. distance from DMSO by cell line
#' Number of cells versus normalized distance from DMSO centroid by cell line.
#' Points called bioactive (centroid distance > 1 IQR median DMSO sample 
#' distances) are colored in orange. Each point corresponds to a well replicate.
#+ bioactivity_count, fig.height=8, fig.width=12, echo=FALSE
# Average bioactivity scores within treatment/cell line groups
xgroup <- rbindlist(xdist) %>%
  group_by(Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup() %>%
  group_by(Cell_Line, Compound_ID, Category) %>% 
  summarize(DistNorm=mean(DistNorm), NCells=mean(NCells), .groups='drop') %>%
  arrange(Compound_ID, Cell_Line) 

mutate(xgroup, Bioactive=DistNorm > 1) %>%
  mutate(Bioactive=ifelse(Bioactive, 'Bioactive', 'Non-bioactive')) %>%
  mutate(Bioactive=factor(Bioactive, levels=c('Non-bioactive', 'Bioactive'))) %>%
  ggplot(aes(x=NCells, y=DistNorm)) +
  geom_point(aes(col=Bioactive), alpha=0.4) +
  facet_wrap(~Cell_Line) +
  scale_x_log10() +
  scale_color_jama(name=element_blank()) +
  geom_vline(xintercept=ncell.thresh, lty=2, col='grey') +
  ylab('Distance from DMSO')

#' # Similarity analysis
#' The following case study consider the problem of assessing phenotypic 
#' similarity  relative to compound MOA. We filter to compounds as follows:
#' 
#' 1. MOA categories with at least `r min.cat` compounds 
#' 2. Compounds that are called as bioactive (distance > 1)
#' 3. Compounds that do not result in cell death phenotype (Ncells >= `r ncell.thresh`)
#+ clustering_analysis
################################################################################
# Initialize cluster analysis compounds
################################################################################

# Get compound counts for bioactive set and filter based on minimum category size
compound.table <- select(xgroup, Compound_ID, Category) %>%
  distinct() %>%
  group_by(Category) %>%
  mutate(Count=n()) %>%
  filter(Count >= min.cat | Compound_ID == 'DMSO') %>%
  filter(!is.na(Category)) %>%
  filter(Category != 'Others') %>%
  filter(Category != '')

x <- filter(x, Compound_ID %in% compound.table$Compound_ID) %>%
  group_by(Compound_ID, Cell_Line) %>% 
  mutate(Dose=as.numeric(Dose)) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup()

#' ## Similarity by MOA
#' To assess the similarity of compounds with the same MOA in phenotypic space, 
#' we consider both within and between category similarity. Within category 
#' similarity compares the distributions of: (i) pairwise distances between 
#' samples with the same MOA and (ii) pairwise distances between DMSO samples. 
#' Intuitively, this asks whether points with the same MOA are closer in 
#' phenotypic space than DMSO controls. Between category similarity compares the 
#' distributions of: (i) pairwise distances between samples with the same MOA
#' (ii) pairwise distances between samples with the given MOA and their nearest 
#' neighbors. Intuitively, this asks how close points with the same MOA are to
#' one another relative to neighboring compounds.
#' 
#' We compare distributions using a signed CVM statistic. For within category
#' similarity, values < 0 imply that the DMSO samples are closer than than 
#' samples with a given MOA. For between category similarity, points are by 
#' definition closer to their nearest neighbors. Thus the signed CVM statistic
#' will always be $\le 0$ with equality when all neighbors are from the same MOA
#' category. We therefore take 1 + the CVM statistic as our measure of between 
#' category similarity.
#+ analysis, fig.height=8, fig.width=16
# Initialize feature set for select cell line

# Compute within/between MOA similarity for each cell line
moa.dist <- lapply(unique(x$Cell_Line), function(cl) {
  
  # Filter to select cell line
  xc <- filter(x, Cell_Line == cl)
  xc.feat <- select(xc, matches('^non'))
  
  # Compute within/between MOA distance
  out <- category_distance(xc.feat, xc$Category, summary_fun) %>%
    mutate(DistBetween= 1 + DistBetween) %>%
    mutate(Count=table(xc$Category)[Category]) %>%
    mutate(Cell_Line=cl)
  
  return(out)
})

moa.dist <- rbindlist(moa.dist)

# Bootstrap within/between MOA similarity for each cell line
moa.dist.bs <- mclapply(1:n.bootstrap, function(i) {
  set.seed(i)
  
  out <- lapply(unique(x$Cell_Line), function(cl) {
  
    # Filter to select cell line
    xc <- filter(x, Cell_Line == cl)
    xc.feat <- select(xc, matches('^non'))
    xc.dist <- as.matrix(dist(xc.feat))
    
    # Compute bootstrap within/between MOA distance
    out <- category_distance(xc.feat, xc$Category, summary_fun, xc.dist, TRUE) %>%
      mutate(DistBetween= 1 + DistBetween) %>%
      mutate(Count=table(xc$Category)[Category]) %>%
      mutate(Cell_Line=cl)
  })
  
  return(rbindlist(out) %>% mutate(Rep=i))
}, mc.cores=n.core)


#' ### Fig 2a. UMAP projection of bioactive + DMSO compounds.
#' Phenotypic profiles projected into umap space. Compounds with high within and
#' between MOA similarity are shown in different colors. Note that UMAP 
#' projection groups strongly clustering compounds (high between MOA similarity)
#' while high within MOA similarity compounds are near the DMSO point cloud but
#' more tightly clustered within.
#+ umap, fig.height=12, fig.width=12
################################################################################
# UMAP visualization
################################################################################
# Subset to selected cell line for UMAP visualization
xc <- filter(x, Cell_Line == umap.cell)
xc.feat <- select(xc, matches('^non'))
dist.category <- filter(moa.dist, Cell_Line == umap.cell) %>% 
  arrange(desc(DistBetween)) %>%
  filter(Category != 'DMSO')

# Initialize umap embedding for visualization
configs <- umap.defaults
configs$n_neighbors <- min.cat
xumap <- umap(xc.feat)

# Initialize data for visualization
xplot <- data.frame(Category=xc$Category) %>%
  mutate(X1=xumap$layout[,1]) %>%
  mutate(X2=xumap$layout[,2])

# Select compounds with high within/between compound similarity
between.select <- dist.category$Category[1:3]
category.select <- c(between.select, 'DMSO', 'HDAC inhibitor', 'DNA inhibitor')

xplot.select <- filter(xplot, Category %in% category.select)
xplot.other <- filter(xplot, !Category %in% category.select) %>%
  mutate(Category='Other')

# Plot select compound category in UMAP space, top neighborhood
col.pal.umap <- pal_jama()(7)
col.pal.umap <- col.pal.umap[c(1, 3, 2, 4:7)]

p <- rbind(xplot.select, xplot.other) %>%
  mutate(Size=(Category != 'Other') + 1) %>%
  mutate(Size=Size + (!Category %in% c('DMSO', 'Other'))) %>%
  ggplot(aes(x=X1, y=X2, col=Category, size=Size, alpha=Size)) +
  geom_jitter(height=0.05, width=0.05) +
  xlab('UMAP1') +
  ylab('UMAP2') +
  scale_color_manual(values=col.pal.umap) +
  scale_size(guide='none', range=c(1.5, 2.5)) +
  scale_alpha(guide='none', range=c(0.25, 0.9)) +
  labs(col=NULL) +
  theme(legend.position=c(0.15, 0.8))

if (save.fig) pdf(str_c(fig.dir, 'umap_', umap.cell, '.pdf'), height=12, width=12)
plot(p)
if(save.fig) dev.off()

#' ### Fig 2b. Query v reference densities
#' Kernel densities of pairwise distances for query (color) and reference (grey)
#' groups. Reference population defined as nearest neighbors of the query 
#' population.
#+ densities, fig.height=12, fig.width=12
  
# Compute within/between MOA distance distributions
category.select <- setdiff(category.select, 'DMSO')
moa.distn <- category_distance_distn(xc.feat, xc$Category)
moa.distn <- moa.distn[category.select]

# Format data table for figure
col.pal.umap <- pal_jama()(7)
col.pal.umap <- col.pal.umap[c(3,2,4,5,7,1)]
col.pal.umap[6] <- '#CCCCCC'

p <- reshape2::melt(moa.distn) %>%
  mutate(L2=ifelse(L2 == 'Query', L1, L2)) %>%
  ggplot(aes(x=value, fill=L2, col=L2)) +
  geom_density(alpha=0.6) +
  scale_fill_manual(values=col.pal.umap) +
  scale_color_manual(values=col.pal.umap) +
  facet_wrap(~L1, ncol=1) +
  xlab('Pairwise distance') +
  theme(legend.position='none')

if (save.fig) pdf(str_c(fig.dir, 'density_', umap.cell, '.pdf'), height=12, width=6)
plot(p)
if(save.fig) dev.off()

#' ### Fig 2c. Between MOA similarity
#' To assess the ability of different cell lines to phenotypically group 
#' compounds by MOA, we consider the distance between compounds (in phenotypic 
#' space) with the same MOA relative to nearby compounds with different MOA. The 
#' figure below reports this relative distance by MOA and cell line.
#+ moa_sim_between, fig.height=12, fig.width=16
################################################################################
# Visualize heatmap of clustering strength
################################################################################
moa.dist <- arrange(moa.dist, Category, Cell_Line)
categories <- unique(moa.dist$Category)
cell.lines <- unique(moa.dist$Cell_Line)

# Initialize data to be plotted
xplot <- matrix(moa.dist$DistBetween, nrow=length(cell.lines))
colnames(xplot) <- categories
rownames(xplot) <- cell.lines

# Rank compounds by average MOA similarity
avg.sim <- rowMeans(xplot)
row.order <- order(avg.sim, decreasing=TRUE)

# Filter to compounds in top quartile
thresh <- quantile(xplot, 0.75)
xplot <- xplot[,colSums(xplot > thresh) > 1]
xplot <- xplot[,colnames(xplot) != 'DMSO']
xplot <- xplot[,colnames(xplot) != 'Others']
between.moa <- colnames(xplot)

# Filter to categories enriched in at least one cell line
col.order <- order(
  xplot[row.order[1],],
  xplot[row.order[2],],
  xplot[row.order[3],],
  xplot[row.order[4],],
  xplot[row.order[5],],
  xplot[row.order[6],]
)

# Clean column names for figure
colnames(xplot) <- str_remove_all(colnames(xplot), ',.*$')

fout <- str_c(fig.dir, 'between_sim.png')
if (save.fig) png(fout, height=12.5, width=26, units='in', res=300)
superheat(
  xplot[row.order, col.order], 
  bottom.label.text.angle=90,
  bottom.label.size=0.5,
  bottom.label.text.size=5,
  heat.pal=heat.pal
)
if (save.fig) dev.off()


#' ### Fig 2d.  Optimal cell line selection
#' Below we compute the optimal cell lines for detecting within and between 
#' compound similarity. For between similarity, we use MOA weights to consider 
#' only "strongly clustering compounds" — at least one cell line > 75th 
#' percentile. For Within compound similarity, we consider all MOAs weighted 
#' equally.
#+ moa_opt, fig.height=12, fig.width=16
# Initialize cell line sets
cell.pairs <- combn(cell.lines, 2, simplify=FALSE)
cell.sets <- c(cell.pairs, cell.lines)
cell.sets.str <- sapply(cell.sets, str_c, collapse=', ')
kcells <- sapply(cell.sets, length)

# Initialize weights for within/between MOA scoring
wt.between <- setNames(rep(1, length(between.moa)), between.moa)

# Score cell line sets by between MOA similarity
moa.score.between <- sapply(moa.dist.bs, function(z) {
  
  # Aggregate between similarity scores for cell line sets
  xagg <- sapply(cell.sets, function(cs) {
    out <- aggregate_similarity(z, cs) %>% arrange(desc(Category))
    return(setNames(out$DistBetween, out$Category))
  })
  
  colnames(xagg) <- cell.sets.str
  return(score_bioactive(t(xagg), wt.between))
})

xopt.between <- reshape2::melt(moa.score.between) %>%
  mutate(CellSet=cell.sets.str[Var1]) %>%
  mutate(KCells=as.factor(kcells[Var1])) %>%
  mutate(Score=value) %>%
  mutate(Type='Between MOA')


# Plot optimal bioactivity scores
p <- xopt.between %>%
  mutate(Alpha=as.numeric(KCells)) %>%
  ggplot(aes(x=reorder(CellSet, Score), y=Score)) +
  geom_boxplot(aes(alpha=Alpha), fill=col.pal[2], color=col.pal[2]) +
  ylab('MOA similarity') +
  scale_color_nejm(drop=FALSE) +
  scale_fill_nejm(drop=FALSE) +
  scale_alpha(range=c(0.25, 0.75)) +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  xlab(NULL)

if (save.fig) pdf(str_c(fig.dir, 'selection.pdf'), height=8, width=8)
plot(p)
if(save.fig) dev.off()

#' ### Fig 2d.  Bioactivity v. cell similarity
#' Joint distribution of bioactivity/MOA similarity scores
#+ bioactive_similarity
setwd(analysis.dir)
load('results/cell_line/category_bioactivity.Rdata')

bioactive.dist <- reshape2::melt(xplot) %>%
  mutate(ID=str_c(Var1, Var2, sep='_')) %>%
  mutate(Bioactivity=value)

moa.dist <- mutate(moa.dist, ID=str_c(Cell_Line, Category, sep='_')) %>%
  left_join(bioactive.dist, by='ID')

xplot <- filter(moa.dist, Cell_Line == umap.cell)
xplot.text <- filter(xplot, Category %in% category.select)

p <- ggplot(xplot, aes(x=Bioactivity, y=DistBetween)) +
  geom_point() +
  geom_text(data=xplot.text, aes(label=Category), nudge_y=0.0075, nudge_x=-0.01) +
  xlab('Bioactivity score') +
  ylab('MOA similarity')

if (save.fig) pdf(str_c(fig.dir, 'bioactive_v_moa.pdf'), height=12, width=12)
plot(p)
if(save.fig) dev.off()
