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
n.bootstrap <- 25
save.fig <- FALSE #TRUE

# Function for normalizing distances within plate relative to DMSO
null_summary <- function(x) {
  iqr.thresh <- IQR(x)
  return(max(x[x < median(x) + iqr.thresh]))
}

# Initialize color palettes
heat.pal <- viridis::viridis(10)
col.pal <- pal_nejm()(8)

################################################################################
# Load dataset
################################################################################
setwd('~/github/cancer_translator/')
source('scripts/utilities.R')

# Load KS profile data
data.dir <- 'data/screens/LH_CDC_1/'
#load(str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata'))
load(str_c(data.dir, 'profiles_normalized.Rdata'))
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
  dplyr::select(-Compound_Category) %>%
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

#' # Category summary
#' Figures below summarize compound counts by category across the full compound 
#' library as well as categories with the largest number of compounds 
#' represented.
#+ count_summaries, fig.height=8, fig.width=12
xcat.tab <- filter(x, !is.na(Category)) %>%
  dplyr::select(Compound_ID, Category) %>%
  distinct() %>%
  group_by(Category) %>%
  summarize(Count=n())

filter(xcat.tab, Count > 1) %>%
  filter(Category != '') %>%
  ggplot(aes(x=Count)) +
  geom_histogram() +
  ggtitle('Distribution of compound counts per category') +
  xlab('# compounds') +
  theme(text=element_text(size=24))

# Summarize prevalent compound categories
xcat.tab.top <- filter(xcat.tab, Count >= min.cat, Category != '')
cat.keep <- xcat.tab.top$Category

xcat.tab.top %>%
  ggplot(aes(x=reorder(Category, Count), y=Count)) +
  geom_bar(stat='identity', position=position_dodge(preserve="single")) +
  geom_text(aes(label=Count), nudge_y=0.025) +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Compounds per category counts', str_c('top categories')) +
  scale_y_log10()

# Write table for supplement
category.table <- dplyr::select(x, Compound_ID, Category, Dose, Cell_Line) %>%
  dplyr::rename(MOA=Category) %>%
  group_by(Compound_ID, MOA, Dose) %>%
  summarize(N_Cell_Line=length(unique(Cell_Line)), .groups='drop') %>%
  mutate(MOA=str_replace_all(MOA, ',', ';')) %>%
  arrange(Compound_ID)

fout <- str_c(output.dir, 'compound_library.csv')
write_excel_csv(file=fout, x=category.table)

#' # Overview
#' This notebook considers the problems of bioactivity detection based on 
#' a panel of 6 cell lines: OVCAR4, A549, DU145, ALS-WT, HEPG2, and 786-0.
#' 
#' # Bioactivity scoring
#' We compute bioactivity scores as follows:
#'
#' 1. Compute the center of the DMSO point cloud  (i.e. median feature values 
#' across each DMSO well) by plate.
#' 2. Generate a null distribution corresponding to the weighted l2 distance 
#' between each DMSO well and the DMSO point cloud center by plate.
#' 3. Normalize distances relative to the DMSO point cloud by plate. 
#' 4. Define bioactivity for a collection of interest (e.g. MOA) based on the 
#' deviation between MOA-centroid and DMSO-centroid ECDFs
#' 
#' 
#' **Note:** Distances are computed within plate, normalized relative to the 
#' DMSO null distribution, and aggregated by treatment (i.e. dose/compound pairs).
#' 
#' **Note:** We assess bioactivity relative to the highest dose level of a given 
#' compound.
#+ bioactivity, fig.height=8, fig.width=14, echo=FALSE

# Mean-center and correlation weight features by cell line
x <- correlation_weight(x, normalize=TRUE)


#' ### Fig. S1 positive control bioactivity by plate
#' Distribution of distances from DMSO centroids for positive + negative 
#' controls by plate, cell line.
#+ bioactivity_control, fig.height=8, fig.width=14, echo=FALSE
# Compute distances between compounds and DMSO centroids
xdist <- lapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(dmso_distance(out, null_summary))
})

names(xdist) <- unique(x$PlateID)

# Plot distance distributions for positive/negative controls
levels <- c('DMSO', 'Bortezomib', 'Gemcitabine')

p <- rbindlist(xdist) %>%
  filter(!str_detect(Compound_Usage, 'query')) %>%
  filter(!str_detect(Compound_Usage, 'reference')) %>%
  mutate(Compound_ID=factor(Compound_ID, level=levels)) %>%
  ggplot(aes(x=reorder(PlateID, DistNorm), y=DistNorm, fill=Compound_ID)) +
  geom_boxplot(alpha=0.6, aes(col=Compound_ID), outlier.size=0.5) +
  theme(axis.text.x=element_blank()) +
  scale_fill_jama() +
  scale_color_jama() +
  theme(legend.position=c(0.9, 0.9), title=element_blank()) +
  facet_wrap(~Cell_Line, scales='free_x') +
  ylab('Distance from DMSO centroid') +
  xlab(NULL)

if (save.fig) pdf(str_c(fig.dir, 'controls.pdf'), height=8, width=12)
plot(p)  
if(save.fig) dev.off()

#' ### Fig 1a UMAP projection by cell line
#' UMAP projections of phenotypic profiles by cell line. Points colored dark 
#' blue represent DMSO-treated wells. Points colored orange and light blue 
#' represent bioactive and inactive calls respsectively.
#+ umap_projections, fig.height=8, fig.width=12

# Select one plate per cell line for visualization
xumap <- lapply(c('OVCAR4', 'HEPG2'), function(cl) {
  
  # Filter DMSO distance table to select cell line
  xdist.cl <- rbindlist(xdist)  %>% 
    filter(Cell_Line == cl) %>%
    mutate(ID=str_c(PlateID, '_', WellID))
  
  # Filter features to select cell line
  xc <- filter(x, Cell_Line == cl)
  xc.feat <- dplyr::select(xc, matches('^non')) 
  
  # Compute UMAP projection
  xumap <- umap(xc.feat)
  colnames(xumap$layout) <- c('UMAP1', 'UMAP2')
  
  out <- data.frame(xumap$layout) %>% 
    mutate(ID=str_c(xc$PlateID, '_', xc$WellID)) %>%
    left_join(xdist.cl, by='ID')
  
  return(out)
})

# Plot samples projected into UMAP space
rbindlist(xumap) %>%
  filter(Compound_ID == 'DMSO') %>%
  group_by(Compound_ID) %>%
  mutate(ColID=as.factor(ColID)) %>%
  ggplot(aes(x=UMAP1, y=UMAP2, col=ColID)) +
  geom_jitter(height=0.25, width=0.25, alpha=0.5) +
  facet_wrap(~Cell_Line) +
  scale_color_nejm() +
  theme(legend.position='none') +
  ggtitle('DMSO samples by column ID')

rbindlist(xumap) %>%
  filter(Compound_ID == 'DMSO') %>%
  group_by(Compound_ID) %>%
  mutate(RowID=as.factor(RowID)) %>%
  ggplot(aes(x=UMAP1, y=UMAP2, col=RowID)) +
  geom_jitter(height=0.25, width=0.25, alpha=0.5) +
  facet_wrap(~Cell_Line) +
  theme(legend.position='none') +
  ggtitle('DMSO samples by column ID')

rbindlist(xumap) %>%
  filter(Compound_ID == 'DMSO') %>%
  ggplot(aes(x=UMAP1, y=UMAP2, col=PlateID)) +
  geom_jitter(height=0.25, width=0.25, alpha=0.5) +
  facet_wrap(~Cell_Line) +
  scale_color_hue(l=60) +
  theme(legend.position='none') +
  ggtitle('DMSO samples by plate ID')


# Initialize parameters for visualization
cat.umap <- c('glucocorticoid receptor agonist', 'HDAC inhibitor')
col.pal.umap <- pal_jama()(5)[c(1, 2, 4)]

p <- rbindlist(xumap) %>%
  filter(Cell_Line %in% c('OVCAR4', 'HEPG2')) %>%
  mutate(Color=ifelse(Compound_ID == 'DMSO', 'DMSO', 'Other')) %>%
  mutate(Color=ifelse(Category == cat.umap[1], cat.umap[1], Color)) %>%
  mutate(Color=ifelse(Category == cat.umap[2], cat.umap[2], Color)) %>%
  mutate(Alpha=as.numeric(Color != 'DMSO')) %>%
  filter(!is.na(Category)) %>%
  group_by(Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose)) | Compound_ID == 'DMSO') %>%
  filter(Color != 'Other') %>%
  ggplot(aes(x=UMAP1, y=UMAP2, col=Color, alpha=Alpha)) +
  geom_jitter(height=0.05, width=0.05) +
  facet_wrap(~Cell_Line, ncol=1) +
  scale_color_manual(values=col.pal.umap) +
  labs(col=NULL) +
  scale_alpha(range=c(0.3, 1)) +
  guides(shape='none') +
  guides(size='none') +
  guides(alpha='none') +
  theme(legend.position=c(0.25, 0.1))

if (save.fig) pdf(str_c(fig.dir, 'umap.pdf'), height=12, width=12)
plot(p)  
if(save.fig) dev.off()

#' # Bioactivity comparisons
#' To assess the degree to which different cell lines detect bioactivity, we 
#' compare the distributions of DMSO well to DMSO centroid with compound to
#' DMSO centroid. For non-DMSO compounds, we take the maximum dose and average
#' DMSO centroid distances across replicates.
#+ bioactivity_merge
# Group bioactivity scores by treatment, cell line, usage
xdmso <- rbindlist(xdist) %>%
  filter(Compound_ID == 'DMSO') %>%
  dplyr::select(Cell_Line, Compound_ID, Category, DistNorm)

xgroup <- rbindlist(xdist) %>%
  filter(Compound_ID != 'DMSO') %>%
  group_by(Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup() %>%
  group_by(Cell_Line, Compound_ID, Category) %>% 
  summarize(DistNorm=mean(DistNorm), .groups='drop') %>%
  arrange(Compound_ID, Cell_Line)

#' ### Fig 1b. DMSO centroid distance distributions
#' Density and ECDFs of distance to DMSO centroid by cell line, DMSO/compound. 
#' Larger distances correspond to compounds that induce phenotypes that are
#' more dissimilar to DMSO, representing stronger bioactivity. Deviations 
#' between the DMSO/compound distributions represent cell line sensitivity.
#+ bioactivity__distn, fig.height=8, fig.width=12
xmerge <- rbind(xdmso, xgroup)
  
p <- xmerge %>%
  mutate(Type=ifelse(Compound_ID == 'DMSO', 'DMSO', 'Other')) %>%
  ggplot(aes(x=DistNorm)) +
  geom_density(aes(fill=Type, color=Type), alpha=0.5) +
  scale_fill_d3(name='') +
  scale_color_d3(name='') +
  xlab('Distance from DMSO centroid') +
  facet_wrap(~Cell_Line) +
  theme(legend.position=c(0.9, 0.9))

if (save.fig) pdf(str_c(fig.dir, 'density.pdf'), height=6, width=14)
plot(p)  
if(save.fig) dev.off()


xplot1 <- xmerge %>%
  filter(Category %in% c(cat.umap[1], 'DMSO')) %>%
  filter(Cell_Line %in% c('OVCAR4', 'HEPG2')) %>%
  mutate(Type=cat.umap[1])

xplot2 <- xmerge %>%
  filter(Category %in% c(cat.umap[2], 'DMSO')) %>%
  filter(Cell_Line %in% c('OVCAR4', 'HEPG2')) %>%
  mutate(Type=cat.umap[2])

p <- rbind(xplot1, xplot2) %>%
  ggplot(aes(x=DistNorm)) +
  geom_density(aes(fill=Category, color=Category), alpha=0.5, bw='nrd') +
  scale_fill_manual(values=col.pal.umap, name='') +
  scale_color_manual(values=col.pal.umap, name='') +
  xlab('Distance from DMSO centroid') +
  facet_grid(Type~Cell_Line) +
  theme(legend.position='none')

if (save.fig) pdf(str_c(fig.dir, 'density_select.pdf'), height=12, width=12)
plot(p)  
if(save.fig) dev.off()

p <- xmerge %>%
  mutate(Type=ifelse(Compound_ID == 'DMSO', 'DMSO', 'Other')) %>%
  group_by(Type, Cell_Line) %>%
  arrange(DistNorm) %>%
  mutate(ECDF=(0:(n() - 1)) / n()) %>%
  ggplot(aes(x=DistNorm, y=ECDF)) +
  geom_line(aes(color=Type), lwd=1) +
  scale_color_d3(name='') +
  xlab('Distance from DMSO centroid') +
  ylab(NULL) +
  facet_wrap(~Cell_Line) +
  theme(legend.position=c(0.9, 0.9))

if (save.fig) pdf(str_c(fig.dir, 'ecdf.pdf'), height=8, width=12)
plot(p)  
if(save.fig) dev.off()

p <- xmerge %>%
  filter(Category %in% c(cat.umap, 'DMSO')) %>%
  filter(Cell_Line %in% c('OVCAR4', 'HEPG2')) %>%
  group_by(Category, Cell_Line) %>%
  arrange(DistNorm) %>%
  mutate(ECDF=(0:(n() - 1)) / n()) %>%
  ggplot(aes(x=DistNorm, y=ECDF)) +
  geom_line(aes(color=Category), lwd=1) +
  scale_color_manual(values=col.pal.umap, name='') +
  xlab('Distance from DMSO centroid') +
  ylab(NULL) +
  facet_wrap(~Cell_Line) +
  theme(legend.position=c(0.9, 0.15))

if (save.fig) pdf(str_c(fig.dir, 'ecdf_select.pdf'), height=8, width=12)
plot(p)  
if(save.fig) dev.off()

#' ### Fig 1c. Bioactivity by cell line, compound category 
#' 
#' Bioactivity scores by cell line compound category. Scores for a given MOA
#' category are computed by comparing DMSO centroid distances between DMSO wells
#' and wells with a select MOA.
#+ bioactivity_category, fig.height=12, fig.width=18
# Compute bioactivity score by compound category
cell.lines <- unique(xmerge$Cell_Line)
xplot <- bioactive_difference_ctg(xmerge, cell.lines, cat.keep)
save(file=str_c(output.dir, 'category_bioactivity.Rdata'), xplot)

# Initialize column ordering and bioactivity proportions
bioactive.score <- rowMeans(xplot)
bioactive.order <- order(bioactive.score, decreasing=TRUE)

# Filter to compounds with high bioactivity in > 1 cell line
id.keep <- colSums(xplot > quantile(xplot, 0.75)) > 0
xplot <- xplot[,id.keep]
xplot[xplot < 0] <- 0

col.order <- order(
  xplot[bioactive.order[1],],
  xplot[bioactive.order[2],],
  xplot[bioactive.order[3],],
  xplot[bioactive.order[4],],
  xplot[bioactive.order[5],],
  xplot[bioactive.order[6],]
)

# Plot bioactivity by treatment, cell line
colnames(xplot) <- str_remove_all(colnames(xplot), ',.*$')

fout <- str_c(fig.dir, 'bioactive_ctg.png')
if (save.fig) png(fout, height=16, width=20, units='in', res=300)
superheat(
  xplot[bioactive.order, col.order],
  heat.pal=heat.pal,
  yr=bioactive.score[bioactive.order],
  yr.plot.type='bar',
  yr.axis.name='Average bioactivity score',
  yr.axis.name.size=14,
  heat.pal.values=seq(0, 1, by=0.1),
  bottom.label.text.angle=90,
  bottom.label.size=0.4,
  bottom.label.text.size=5
)
if (save.fig) dev.off()


#' ### Fig. 1d. Bioactivity optimization
#' Optimal cell line sets over all compounds — equivalent to setting category 
#' weights proportional to # of compounds in each category.
#+ bioactivity_generalist
# Initialize cell line sets
cell.pairs <- combn(cell.lines, 2, simplify=FALSE)
cell.sets <- c(cell.pairs, cell.lines)
cell.sets.str <- sapply(cell.sets, str_c, collapse=', ')
kcells <- sapply(cell.sets, length)

# Compute bioactivity for cell sets
xcat.bootstrap <- mclapply(1:n.bootstrap, function(i) {
  set.seed(i)
  xplot <- bioactive_difference_ctg(xmerge, cell.sets, cat.keep, TRUE)
}, mc.cores=n.core)

# Initialize weights for generalist/specialist
specialist.cat <- c('CDK inhibitor', 'MEK inhibitor', 'JAK inhibitor')
category <- xmerge$Category[xmerge$Category %in% cat.keep]

wt.gen <- table(category) %>% c
wt.spec <- setNames(rep(1, length(specialist.cat)), specialist.cat)

# Compute score for generalist/specialist
xopt.gen <- sapply(xcat.bootstrap, score_bioactive, weights=wt.gen)
xopt.spec <- sapply(xcat.bootstrap, score_bioactive, weights=wt.spec)

# Initialize table for figures
xopt.gen <- sapply(xcat.bootstrap, score_bioactive, weights=wt.gen)
xopt.gen <- reshape2::melt(xopt.gen) %>%
  mutate(CellSet=cell.sets.str[Var1]) %>%
  mutate(KCells=as.factor(kcells[Var1])) %>%
  mutate(Score=value) %>%
  mutate(Type='Generalist')
  
xopt.spec <- sapply(xcat.bootstrap, score_bioactive, weights=wt.spec)
xopt.spec <- reshape2::melt(xopt.spec) %>%
  mutate(CellSet=cell.sets.str[Var1]) %>%
  mutate(KCells=as.factor(kcells[Var1])) %>%
  mutate(Score=value) %>%
  mutate(Type='Specialist')

# Plot optimal bioactivity scores
p1 <- xopt.gen %>%
  mutate(Alpha=as.numeric(KCells)) %>%
  ggplot(aes(x=reorder(CellSet, Score), y=Score)) +
  geom_boxplot(aes(alpha=Alpha), fill=col.pal[2], color=col.pal[2]) +
  ylab('Bioactivity score') +
  scale_color_nejm(drop=FALSE) +
  scale_fill_nejm(drop=FALSE) +
  scale_alpha(range=c(0.25, 0.75)) +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Generalist') +
  xlab(NULL)

p2 <- xopt.spec %>%
  mutate(Alpha=as.numeric(KCells)) %>%
  ggplot(aes(x=reorder(CellSet, Score), y=Score)) +
  geom_boxplot(aes(alpha=Alpha), fill=col.pal[2], color=col.pal[2]) +
  ylab('Bioactivity score') +
  scale_color_nejm(drop=FALSE) +
  scale_fill_nejm(drop=FALSE) +
  scale_alpha(range=c(0.25, 0.75)) +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Specialist', str_c(specialist.cat, collapse=', ')) +
  xlab(NULL)

if (save.fig) pdf(str_c(fig.dir, 'selection.pdf'), height=12, width=12)
gridExtra::grid.arrange(p1, p2, ncol=1)
if(save.fig) dev.off()
