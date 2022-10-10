#' This script runs bioactivity analyses described in figure 3. Compound/MOA 
#' table and normalized KS profiles are loaded by a call to `load_ks_data.R`,
#' which requires (i) that normalized KS profiles have been generated using 
#' `make_ks_profiles.R` and `normalize_ks_profiles.R` and (ii) that a clean 
#' compound/MOA table has been generated using `make_compound_table.R`.
#' 
#' Last updated 10/5/2022, Karl Kumbier
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
    plot_title_size=28,
    axis_title_size=24, 
    strip_text_size=28, 
    axis_text_size=12,
    base_size=12, 
    base_family='sans'
  )
)

################################################################################
# Initialize constant parameters to be used in analysis
################################################################################
# Initialize analysis parameters
n.core <- 16
min.cat <- 5
n.bootstrap <- 50
umap.cells <- c('A549', 'FB')
save.fig <- TRUE

# Visualization parameters
heat.pal <- viridis::viridis(10)
heat.pal.signed <- brewer.pal(9, 'RdBu')
col.pal <- pal_jama()(7)
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
# Initialize analysis directories
ccl.dir <- Sys.getenv('CCL_DIR')
fig.dir <- str_c(ccl.dir, 'paper_figures/Figure3/fig3/')
phenosimilarity.file <- str_c(ccl.dir, 'data/Table4_phenosimilarity.csv')
dir.create(fig.dir, showWarnings=FALSE)

source(str_c(ccl.dir, 'paper_figures/utilities.R'))
source(str_c(ccl.dir, 'paper_figures/load_ks_data.R'))

cell.lines <- unique(x$Cell_Line)

#' # Overview
#' This notebook considers the problem of compound similarity based on 
#' phenotypic profiles. Similarities are computed over correlation re-weighted
#' phenotypic profiles.
#+ bioactivity, fig.height=8, fig.width=12, echo=FALSE
# Mean-center + correlation weight features by cell line for subsequent analysis
x <- correlation_weight(x, normalize=TRUE)

# Compute distance to DMSO centroid
xdist <- lapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(dmso_distance(out, null_summary))
})

#' # Similarity analysis
#' The following case study consider the problem of assessing phenotypic 
#' similarity relative to compound category. Prior to computing similarities, 
#' compounds are filtered as follows:
#' 
#' 1. MOA categories with at least `r min.cat` compounds 
#' 2. Filtering out compounds with no category / category == Other
#' 3. Filtering out "cell death" compounds - # cells below 100
#+ clustering_analysis
################################################################################
# Initialize cluster analysis compounds
################################################################################

# Filter to categories exceeding minimum number bioactive compounds
compound.table <- select(x, Compound_ID, Category, Cell_Line) %>%
  distinct() %>%
  group_by(Category, Cell_Line) %>%
  mutate(Count=n()) %>%
  filter(Count >= min.cat | Compound_ID == 'DMSO') %>%
  filter(!is.na(Category)) %>%
  filter(Category != 'Others')

x <- filter(x, Compound_ID %in% compound.table$Compound_ID)

# Filter out cell death compounds
x <- filter(x, NCells >= 100)

#' ## Similarity by MOA
#' To assess the similarity of compounds with the same MOA in phenotypic space, 
#' we consider the distance between samples with the same MOA relative to their 
#' nearest neighbors.
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
moa.dist <- lapply(cell.lines, function(cl) {
  
  # Filter to select cell line
  xc <- filter(x, Cell_Line == cl)
  xc.feat <- select(xc, matches('^non'))
  
  # Compute within/between MOA distance
  xdist <- phenosimilarity(xc.feat, xc$Category, summary_fun)
  out <- data.frame(Category=xdist$Category, Dist=xdist$Dist) %>%
    mutate(Prop=xdist$Prop) %>%
    mutate(Dist=1 + (Dist / 0.5)) %>%
    mutate(Count=table(xc$Category)[Category]) %>%
    mutate(Cell_Line=cl)
  
  return(out)
})

moa.dist <- rbindlist(moa.dist)

# Compute within/between MOA similarity for each cell line
cpd.dist <- lapply(unique(x$Cell_Line), function(cl) {
  
  # Filter to select cell line
  xc <- filter(x, Cell_Line == cl, Compound_ID != 'DMSO')
  xc.feat <- select(xc, matches('^non'))
  
  # Compute within/between MOA distance
  out <- phenosimilarity_cpd(xc.feat, xc$Category, xc$Compound_ID) %>%
    mutate(Cell_Line=cl)
    
  return(out)
})


# Write table of distance from DMSO by compound
xtable <- rbindlist(cpd.dist) %>%
  dplyr::select(Cell_Line, Compound_ID, Category, Dist) %>%
  mutate(Key=str_c(Cell_Line, Compound_ID)) %>%
  left_join(xkey, by='Key') %>%
  dplyr::select(-Key) %>%
  arrange(Cell_Line, Category, desc(Dist)) %>%
  mutate(Compound_ID=str_replace_all(Compound_ID, ',', ';')) %>%
  mutate(Category=str_replace_all(Category, ',', ';')) %>%
  dplyr::rename(Dist2Neighbors=Dist) %>%
  dplyr::rename(ExperimentalID=ID) %>%
  group_by(Cell_Line) %>%
  arrange(desc(Dist2Neighbors))
  
write.csv(file=phenosimilarity.file, xtable, row.names=FALSE)


# Bootstrap within/between MOA similarity for each cell line
moa.dist.bs <- lapply(unique(x$Cell_Line), function(cl) {
  
  # Filter to select cell line
  xc <- filter(x, Cell_Line == cl)
  xc.feat <- select(xc, matches('^non'))
  xc.dist <- as.matrix(dist(xc.feat))
  
  out <- mclapply(1:n.bootstrap, function(i) {

    set.seed(i)
    
    # Compute bootstrap within/between MOA distance
    xdist <- phenosimilarity(xc.feat, xc$Category, summary_fun, xc.dist, TRUE)
    out <- data.frame(Category=xdist$Category, Dist=xdist$Dist) %>%
      mutate(Prop=xdist$Prop) %>%
      mutate(Dist=Dist / 0.5) %>%
      mutate(Dist= 1 + Dist) %>%
      mutate(Count=table(xc$Category)[Category]) %>%
      mutate(Cell_Line=cl) %>%
      mutate(Rep=i)
    
    return(out)
  }, mc.cores=n.core)
  
  return(rbindlist(out))
})


#' ### Fig 3a. UMAP projection of bioactive + DMSO compounds.
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
xumap <- lapply(umap.cells, function(cl) {
 
  set.seed(47)
  
  # Filter to data from select cell line
  xc <- filter(x, Cell_Line == cl)
  xc.feat <- select(xc, matches('^non'))
  dist.category <- filter(moa.dist, Cell_Line == cl) %>% 
    arrange(desc(Dist)) %>%
    filter(Category != 'DMSO')

  # Initialize umap embedding for visualization
  configs <- umap.defaults
  configs$n_neighbors <- min.cat
  xumap <- umap(xc.feat, config=configs)

  # Initialize data for visualization
  out <- data.frame(Category=xc$Category) %>%
    mutate(Compound_ID=xc$Compound_ID) %>%
    mutate(Cell_Line=cl) %>%
    mutate(X1=xumap$layout[,1]) %>%
    mutate(X2=xumap$layout[,2]) 
  
  return(out)
})
  

# Select compounds with high within/between compound similarity
category.select <- c(
  'HSP inhibitor',
  'glycogen synthase kinase inhibitor',
  'SYK inhibitor',
  'mTOR inhibitor',
  'DMSO'
)

xplot <- rbindlist(xumap)
xplot.select <- filter(xplot, Category %in% category.select)
xplot.other <- xplot %>%
  filter(!Category %in% category.select) %>%
  mutate(Category='Other')

# Plot select compound category in UMAP space, top neighborhood
col.pal.umap <- pal_jama()(7)

p <- rbind(xplot.select, xplot.other) %>%
  mutate(Size=(Category != 'Other') + 1) %>%
  mutate(Size=Size + (!Category %in% c('DMSO', 'Other'))) %>%
  mutate(Category=str_replace_all(Category, 'glycogen synthase kinase', 'GSK')) %>%
  ggplot(aes(x=X1, y=X2, col=Category, size=Size, alpha=Size)) +
  geom_jitter(height=0.25, width=0.25) +
  xlab('UMAP1') +
  ylab('UMAP2') +
  scale_color_manual(values=col.pal.umap) +
  scale_size(guide='none', range=c(1.5, 2.5)) +
  scale_alpha(guide='none', range=c(0.25, 0.9)) +
  labs(col=NULL) +
  facet_wrap(~Cell_Line, ncol=1) +
  theme(legend.position='none')

if (save.fig) pdf(str_c(fig.dir, 'fig3a.pdf'), height=18, width=12)
plot(p)
if(save.fig) dev.off()

################################################################################
# mTOR comparison
################################################################################
xf <- filter(x, Category == 'mTOR inhibitor', Cell_Line == 'A549')
xfeat <- dplyr::select(xf, matches('^nonborder')) 

rownames(xfeat) <- str_remove_all(xf$Compound_ID, '\\(.*')
pal <- RColorBrewer::brewer.pal(8, 'RdBu')
pal <- c(pal[1:4], '#FFFFFF', pal[5:8])

fout <- str_c(fig.dir, 'fig3e.png')
if (save.fig) png(fout, height=20, width=12, units='in', res=300)
superheat(
  xfeat,
  row.dendrogram=TRUE,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  bottom.label='none',
  left.label='variable',
  left.label.text.size=9,
  left.label.size=0.7,
  heat.pal=pal,
  heat.lim=c(-1, 1)
)
if (save.fig) dev.off()

#' ### Fig 3b. Query v reference densities
#' Kernel densities of pairwise distances for query (color) and reference (grey)
#' groups. Reference population defined as nearest neighbors of the query 
#' population.
#+ densities, fig.height=12, fig.width=12
  
# Compute within/between MOA distance distributions
moa.distn <- lapply(umap.cells, function(cl) {
  xc <- filter(x, Cell_Line == cl)
  xc.feat <- select(xc, matches('^non'))

  category.select <- setdiff(category.select, 'DMSO')
  out <- phenosimilarity_distn(xc.feat, xc$Category)
  return(out[setdiff(category.select, 'DMSO')])
})

# Format data table for figure
col.pal.umap <- pal_jama()(7)
col.pal.umap <- c(col.pal.umap[c(2, 3, 4, 6)], '#CCCCCC')

# Initialize text for scores
category.select <- setdiff(category.select, 'DMSO')

xplot.text <- filter(moa.dist, Category %in% category.select) %>%
  filter(Cell_Line %in% umap.cells) %>%
  mutate(Dist=str_c('Phenosimilarity = ', round(Dist, 3))) %>%
  mutate(Category=str_replace_all(Category, 'glycogen synthase kinase', 'GSK')) %>%
  mutate(G=Category)
  
xplot <- reshape2::melt(moa.distn) %>%
  mutate(Cell_Line=umap.cells[L1]) %>%
  mutate(L2=str_replace_all(L2, 'glycogen synthase kinase', 'GSK')) %>%
  mutate(L3=str_replace_all(L3, 'glycogen synthase kinase', 'GSK')) %>%
  mutate(G=ifelse(L3 == 'Query', L2, 'ZZZ')) %>%
  dplyr::rename(Category=L2)
  
p <- ggplot(xplot, aes(x=value, fill=G, col=G)) +
  geom_density(alpha=0.6) +
  geom_text(data=xplot.text, x=2.75, y=2, size=8, aes(label=Dist)) +
  facet_grid(Cell_Line~Category) +
  scale_fill_manual(values=col.pal.umap) +
  scale_color_manual(values=col.pal.umap) +
  xlab('Pairwise distance') +
  ylab('Density') +
  theme(legend.position='none')

if (save.fig) pdf(str_c(fig.dir, 'fig3b.pdf'), height=14, width=24)
plot(p)
if(save.fig) dev.off()

#' ### Fig 3c. Between MOA similarity
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
xplot <- matrix(moa.dist$Dist, nrow=length(cell.lines))
colnames(xplot) <- categories
rownames(xplot) <- cell.lines

rownames(xplot)
tt <- table(apply(xplot, MAR=2, which.max))
sum(tt) - max(tt)

# Rank compounds by average MOA similarity
avg.sim <- rowMeans(xplot)
row.order <- order(avg.sim, decreasing=TRUE)

# Filter to compounds in top quartile
xplot <- xplot[,colnames(xplot) != 'DMSO']
xplot <- xplot[,colnames(xplot) != 'Others']
categories.moa <- colnames(xplot)

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
xplot[xplot < 0] <- 0

fout <- str_c(fig.dir, 'fig3c.png')
if (save.fig) png(fout, height=20, width=8, units='in', res=300)
superheat(
  t(xplot[row.order, col.order]), 
  bottom.label.text.angle=90,
  bottom.label.text.size=9,
  heat.pal=heat.pal,
  heat.pal.values=c(seq(0, 0.45, by=0.05), seq(0.5, 1, by=0.1)),
  heat.lim=c(0, max(xplot))
)
if (save.fig) dev.off()


#' ### Fig 3d.  Optimal cell line selection
#' Below we compute the optimal cell lines for detecting within and between 
#' compound similarity. For between similarity, we use MOA weights to consider 
#' only "strongly clustering compounds" — at least one cell line > 75th 
#' percentile. For Within compound similarity, we consider all MOAs weighted 
#' equally.
#+ moa_opt, fig.height=12, fig.width=16
# Initialize cell line sets
cell.pairs <- combn(cell.lines, 2, simplify=FALSE)
cell.sets <- c(list(cell.lines), cell.pairs, cell.lines)
cell.sets.str <- sapply(cell.sets, str_c, collapse=', ')
cell.sets.str[1] <- str_replace_all(cell.sets.str[1], '145,', '145,\n')

kcells <- sapply(cell.sets, length)

# Initialize weights for within/between MOA scoring
weights <- setNames(rep(1, length(categories.moa)), categories.moa)

# Score cell line sets by between MOA similarity
moa.dist.bs <- rbindlist(moa.dist.bs)

moa.score <- sapply(1:n.bootstrap, function(i) {
  
  # Aggregate between similarity scores for cell line sets
  z <- filter(moa.dist.bs, Rep == i)
  xagg <- sapply(cell.sets, function(cs) {
    out <- aggregate_similarity(z, cs) %>% arrange(desc(Category))
    return(setNames(out$Dist, out$Category))
  })
  
  colnames(xagg) <- cell.sets.str
  scores <- apply(xagg, MAR=2, score_bioactive, weights=weights)
  return(scores)
})

xopt <- reshape2::melt(moa.score) %>%
  mutate(CellSet=cell.sets.str[Var1]) %>%
  mutate(KCells=as.factor(kcells[Var1])) %>%
  mutate(Score=value) %>%
  mutate(Type='MOA')


# Plot optimal bioactivity scores
p <- xopt %>%
  mutate(Alpha=as.numeric(KCells)) %>%
  ggplot(aes(x=reorder(CellSet, Score), y=Score)) +
  geom_boxplot(aes(alpha=Alpha), fill=col.pal[1], color=col.pal[1]) +
  ylab('Phenosimilarity') +
  scale_color_nejm(drop=FALSE) +
  scale_fill_nejm(drop=FALSE) +
  scale_alpha(range=c(0.25, 0.75)) +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  xlab(NULL)

if (save.fig) pdf(str_c(fig.dir, 'fig3d.pdf'), height=8, width=6)
plot(p)
if(save.fig) dev.off()
