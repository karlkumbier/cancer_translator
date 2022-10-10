#' This script runs bioactivity analyses described in figure 2. Compound/MOA 
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
library(patchwork)
library(hrbrthemes)
library(umap)
library(readr)
library(twosamples)
library(see)
library(ggridges)

theme_set(
  theme_ipsum(
    plot_title_size=28,
    axis_title_size=24, 
    strip_text_size=28, 
    axis_text_size=18,
    base_size=18, 
    base_family='sans'
  )
)

set.seed(47)

################################################################################
# Analysis parameters
################################################################################
# Initialize analysis parameters
n.core <- 16
min.cat <- 5 # Minimum category size threshold
n.bootstrap <- 50 # number of bootstrap samples for optimization
save.fig <- TRUE

# Within plate normalization function 
null_summary <- function(x) {
  iqr.thresh <- IQR(x)
  return(max(x[x < median(x) + iqr.thresh]))
}

# Initialize color palettes
col.pal <- pal_jama()(7)
heat.pal <- viridis::viridis(10)

################################################################################
# Load dataset
################################################################################
# Initialize analysis directories
ccl.dir <- Sys.getenv('CCL_DIR')
fig.dir <- str_c(ccl.dir, 'paper_figures/Figure2/fig2/')
bioactivity.file <- str_c(fig.dir, 'bioactivity.Rdata') 
cat.bioactivity.file <- str_c(fig.dir, 'category_bioactivity.Rdata') 
feature.file <- str_c(ccl.dir, 'data/Table2_features.csv')
phenoactivity.file <- str_c(ccl.dir, 'data/Table3_dist2dmso.csv')
analysis.file <-  str_c(fig.dir, 'analysis.Rdata')
dir.create(fig.dir, showWarnings=FALSE)

source(str_c(ccl.dir, 'paper_figures/utilities.R'))
source(str_c(ccl.dir, 'paper_figures/load_ks_data.R'))

# Initialize prevalent compound categories for figures
cat.keep <- group_by(x, Category, Cell_Line) %>%
  count() %>%
  filter(n >= min.cat)

cat.keep <- unique(cat.keep$Category)

# Save feature table
features <- colnames(x) %>% 
  str_subset('^nonborder') %>%
  str_remove_all('nonborder\\.\\.\\.')

write.csv(file=feature.file, features, row.names=FALSE)

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
#' DMSO null distribution.
#' 
#' **Note:** We assess bioactivity relative to the highest dose level of a given 
#' compound (processing performed in `load_ks_data.R`).
#' 
#' To account for redundancy among features, each feature is re-weighted 
#' relative to it's correlationwith other features in the data set. As a 
#' result, features that are highly correlated with others are down-weighted
#+ bioactivity, fig.height=8, fig.width=14, echo=FALSE

# Mean-center and correlation weight features by cell line
x <- correlation_weight(x, normalize=TRUE)

#' ### Fig. S1 positive control bioactivity by plate
#' Distribution of distances from DMSO centroids for positive + negative 
#' controls by plate, cell line.
#+ bioactivity_control, fig.height=8, fig.width=14, echo=FALSE
# Compute distances between compounds and DMSO centroids
xdist <- lapply(unique(x$Cell_Line), function(cl) {
  xc <- filter(x, Cell_Line == cl)
  return(dmso_distance(xc, null_summary))
})

# Write table of distance from DMSO by compound
xtable <- rbindlist(xdist) %>%
  filter(Compound_ID != 'DMSO') %>%
  dplyr::select(Cell_Line, Compound_ID, Category, DistNorm) %>%
  mutate(Key=str_c(Cell_Line, Compound_ID)) %>%
  left_join(xkey, by='Key') %>%
  dplyr::select(-Key) %>%
  dplyr::rename(Dist=DistNorm) %>%
  arrange(Cell_Line, Category, desc(Dist)) %>%
  dplyr::rename(Distance2DMSO=Dist) %>%
  dplyr::rename(ExperimentalID=ID) %>%
  group_by(Cell_Line) %>%
  arrange(desc(Distance2DMSO)) %>%
  ungroup()

write.csv(file=phenoactivity.file, xtable, row.names=FALSE)

# Plot distance distributions for positive/negative controls
levels <- c('DMSO', 'Bortezomib', 'Gemcitabine')

xplot <- rbindlist(xdist)

p <- xplot %>%
  mutate(XCpd=ifelse(Compound_ID == 'DMSO', 'DMSO', 'Other')) %>%
  ggplot(aes(x=DistNorm)) +
  geom_density(alpha=0.5, aes(fill=XCpd, col=XCpd)) +
  geom_segment(data=filter(xplot, Compound_ID == levels[2]), 
               aes(x=DistNorm, xend=DistNorm, y=0, yend=0.1), 
               col=pal_jama()(4)[3], size=2) +
  geom_segment(data=filter(xplot, Compound_ID == levels[3]), 
               aes(x=DistNorm, xend=DistNorm, y=0, yend=0.1), 
               col=pal_jama()(4)[4], size=2) +
  scale_fill_jama() +
  scale_color_jama() +
  theme(legend.position='none', title=element_blank()) +
  theme(axis.text.x=element_text(angle=90)) +
  facet_wrap(~Cell_Line, ncol=3) +
  xlab('Distance from DMSO centroid') +
  ylab('Density')
  

if (save.fig) pdf(str_c(fig.dir, 'supplemental_fig_ctrls.pdf'), height=8, width=12)
plot(p)  
if(save.fig) dev.off()


#' ### Fig 2a UMAP projection by cell line
#' UMAP projections of phenotypic profiles by cell line. Points colored dark 
#' blue represent DMSO-treated wells. Points colored orange and light blue 
#' represent bioactive and inactive calls respectively.
#+ umap_projections, fig.height=8, fig.width=12
# Select one plate per cell line for visualization
xumap <- lapply(c('OVCAR4', 'HEPG2'), function(cl) {
  
  # Filter DMSO distance table to select cell line
  xdist.cl <- rbindlist(xdist)  %>% 
    filter(Cell_Line == cl)
  
  # Filter features to select cell line
  xc <- filter(x, Cell_Line == cl)
  xc.feat <- dplyr::select(xc, matches('^non')) 
  
  # Compute UMAP projection
  xumap <- umap(xc.feat)
  colnames(xumap$layout) <- c('UMAP1', 'UMAP2')
  return(data.frame(xumap$layout, xdist.cl) )
})


# Initialize parameters for visualization
cat.umap <- c('glucocorticoid receptor agonist', 'HDAC inhibitor')
lines.umap <- c('OVCAR4', 'HEPG2')
col.pal.umap <- pal_jama()(5)[c(1, 2, 4)]

xplot <- rbindlist(xumap) %>%
  filter(Cell_Line %in% c('OVCAR4', 'HEPG2')) %>%
  mutate(Color=ifelse(Compound_ID == 'DMSO', 'DMSO', 'Other')) %>%
  mutate(Color=ifelse(Category == cat.umap[1], cat.umap[1], Color)) %>%
  mutate(Color=ifelse(Category == cat.umap[2], cat.umap[2], Color)) %>%
  mutate(Alpha=as.numeric(Color != 'DMSO')) %>%
  filter(!is.na(Category)) %>%
  group_by(Compound_ID) %>%
  filter(Color != 'Other') %>%
  ungroup()

p <- ggplot(xplot, aes(x=UMAP1, y=UMAP2, col=Color, alpha=Alpha)) +
  geom_jitter(height=0.05, width=0.05) +
  facet_wrap(~Cell_Line, ncol=1) +
  scale_color_manual(values=col.pal.umap) +
  scale_alpha(range=c(0.3, 1)) +
  theme(legend.position='none')

if (save.fig) pdf(str_c(fig.dir, 'fig2a.pdf'), height=12, width=12)
plot(p)  
if(save.fig) dev.off()

# Summary statistics for select MOAs
xplot <- rbindlist(xdist) %>%
  filter(Category %in% cat.umap, Cell_Line %in% lines.umap) %>%
  group_by(Cell_Line, Category) %>%
  summarize(PropActive=mean(DistNorm > 1), Count=sum(DistNorm > 1), N=n())

p <- ggplot(xplot, aes(x=Cell_Line, y=PropActive)) +
  facet_grid(. ~ Category) +
  geom_bar(stat='identity') +
  scale_fill_jama() +
  ylab('Proportion compounds > 1 IQR from DMSO centroid') +
  xlab(NULL)

if (save.fig) pdf(str_c(fig.dir, 'supplemental_fig_summary.pdf'), height=12, width=15)
plot(p)  
if(save.fig) dev.off()

#' # Pheno-activity comparisons
#' ### Fig 2b. DMSO centroid distance distributions
#' Density and ECDFs of distance to DMSO centroid by cell line, DMSO/compound. 
#' Larger distances correspond to compounds that induce phenotypes that are
#' more dissimilar to DMSO, representing stronger bioactivity. Deviations 
#' between the DMSO/compound distributions represent cell line sensitivity.
#+ bioactivity_distn, fig.height=8, fig.width=12
# Compute pheno-activity scores
xmerge <- rbindlist(xdist) %>% arrange(Compound_ID, Cell_Line)

# Compute bioactivity score by compound category
cell.lines <- unique(xmerge$Cell_Line)
xpa <- sapply(cell.lines, phenoactivity, x=xmerge, categories=cat.keep)
xpa <- t(xpa) / 0.5

save(file=cat.bioactivity.file, xpa)
  
# Plot densities of selected compound/cell line pairs
lines.select <- c('OVCAR4', 'HEPG2')

xplot1 <- xmerge %>%
  filter(Category %in% c(cat.umap[1], 'DMSO')) %>%
  filter(Cell_Line %in% lines.select) %>%
  mutate(Type=cat.umap[1])

xplot2 <- xmerge %>%
  filter(Category %in% c(cat.umap[2], 'DMSO')) %>%
  filter(Cell_Line %in% lines.select) %>%
  mutate(Type=cat.umap[2])

xtext <- data.frame(Cell_Line=rep(lines.select, 2)) %>%
  mutate(Type=rep(cat.umap, each=2)) %>%
  mutate(Type=factor(Type, levels=cat.umap[2:1])) %>%
  mutate(Category=Type) %>%
  mutate(Score=c(xpa[lines.select, cat.umap[1]], xpa[lines.select, cat.umap[2]])) %>%
  mutate(Score=str_c('Phenoactivity = ', round(Score, 3)))

p <- rbind(xplot1, xplot2) %>%
  mutate(Type=factor(Type, levels=cat.umap[2:1])) %>%
  ggplot(aes(fill=Category, color=Category)) +
  geom_density(aes(x=DistNorm), alpha=0.5, bw='nrd') +
  geom_text(data=xtext, x=3, y=1, size=8, aes(label=Score)) +
  scale_fill_manual(values=col.pal.umap, name='') +
  scale_color_manual(values=col.pal.umap, name='') +
  xlab('Distance from DMSO centroid') +
  ylab('Density') +
  facet_grid(Cell_Line~Type) +
  theme(legend.position='none')

if (save.fig) pdf(str_c(fig.dir, 'fig2b.pdf'), height=15, width=15)
plot(p)  
if(save.fig) dev.off()

#' ### Fig 2c. Bioactivity by cell line, compound category 
#' 
#' Bioactivity scores by cell line compound category. Scores for a given MOA
#' category are computed by comparing DMSO centroid distances between DMSO wells
#' and wells with a select MOA.
#+ bioactivity_category, fig.height=18, fig.width=24
ctg.table <- table(rownames(xpa)[apply(xpa, MAR=2, which.max)])

print(ctg.table)
print(sum(ctg.table))

# Initialize column ordering and bioactivity proportions
bioactive.score <- rowMeans(xpa)
bioactive.order <- order(bioactive.score, decreasing=TRUE)

# Filter to compounds with high bioactivity in > 1 cell line
#id.keep <- colSums(xpa > quantile(xpa, 0.75)) > 1
xpa.plot <- xpa#[,]#id.keep]
xpa.plot[xpa.plot < 0] <- 0

col.order <- order(
  xpa.plot[bioactive.order[1],],
  xpa.plot[bioactive.order[2],],
  xpa.plot[bioactive.order[3],],
  xpa.plot[bioactive.order[4],],
  xpa.plot[bioactive.order[5],],
  xpa.plot[bioactive.order[6],]
)

# Plot bioactivity by treatment, cell line
colnames(xpa) <- str_remove_all(colnames(xpa), ',.*$')

fout <- str_c(fig.dir, 'fig2c.png')
if (save.fig) png(fout, height=20, width=8, units='in', res=300)

xpa.plot[xpa.plot < 0.25] <- 0.25

superheat(
  t(xpa.plot[bioactive.order, col.order]),
  heat.pal=heat.pal,
  yt=bioactive.score[bioactive.order],
  yt.plot.type='bar',
  yt.axis.name='Mean\nphenoactivity',
  yt.axis.name.size=24,
  yt.axis.size=12,
  yt.plot.size=0.1,
  heat.pal.values=seq(0, 1, by=0.1),
  left.label.size=0.5,
  left.label.text.size=5,
  bottom.label.text.angle=90,
  bottom.label.text.size=10
)

if (save.fig) dev.off()


#' ### Fig. 1d. Bioactivity optimization
#' Optimal cell line sets over all compounds — equivalent to setting category 
#' weights proportional to # of compounds in each category.
#+ bioactivity_generalist, fig.height=12, fig.width=12
# Initialize cell line sets
cell.pairs <- combn(cell.lines, 2, simplify=FALSE)
cell.sets <- c(list(cell.lines), cell.pairs, cell.lines)
cell.sets.str <- sapply(cell.sets, str_c, collapse=', ')
cell.sets.str[1] <- str_replace_all(cell.sets.str[1], '145,', '145,\n')
kcells <- sapply(cell.sets, length)

# Compute bioactivity for cell sets
xcat.bootstrap <- mclapply(1:n.bootstrap, function(i) {
  set.seed(i)
  out <- sapply(cell.lines, phenoactivity, x=xmerge, categories=cat.keep, bootstrap=TRUE)
  return(t(out) / 0.5)
}, mc.cores=n.core)

# Initialize weights for generalist/specialist
specialist.cat <- c('JAK inhibitor', 'VEGFR inhibitor', 'HDAC inhibitor', 'PI3K inhibitor')
category <- xmerge$Category[xmerge$Category %in% cat.keep]

wt.gen <- table(category) %>% c
wt.spec <- setNames(rep(1, length(specialist.cat)), specialist.cat)

# Initialize table for figures
xopt.gen <- sapply(xcat.bootstrap, function(z) {
  sapply(cell.sets, function(s) {
    zmax <- apply(z[s,,drop=FALSE], MAR=2, max)
    return(score_bioactive(zmax, weights=wt.gen))
  })
})
  
xopt.gen <- reshape2::melt(xopt.gen) %>%
  mutate(CellSet=cell.sets.str[Var1]) %>%
  mutate(KCells=as.factor(kcells[Var1])) %>%
  mutate(Score=value) %>%
  mutate(Type='Generalist')
  
xopt.spec <- sapply(xcat.bootstrap, function(z) {
  sapply(cell.sets, function(s) {
    zmax <- apply(z[s,,drop=FALSE], MAR=2, max)
    return(score_bioactive(zmax, weights=wt.spec))
  })
})

xopt.spec <- reshape2::melt(xopt.spec) %>%
  mutate(CellSet=cell.sets.str[Var1]) %>%
  mutate(KCells=as.factor(kcells[Var1])) %>%
  mutate(Score=value) %>%
  mutate(Type='Specialist')

# Plot optimal bioactivity scores
group_by(xopt.gen, CellSet) %>%
  summarize(Score=mean(Score)) %>%
  arrange(desc(Score)) %>%
  as.data.frame

p <- xopt.gen %>%
  mutate(Alpha=as.numeric(KCells)) %>%
  ggplot(aes(x=reorder(CellSet, Score), y=Score)) +
  geom_boxplot(aes(alpha=Alpha), fill=col.pal[1], color=col.pal[1]) +
  ylab('Phenoactivity') +
  scale_alpha(range=c(0.25, 0.75)) +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Generalist') +
  xlab(NULL)


if (save.fig) pdf(str_c(fig.dir, 'fig2d.pdf'), height=12, width=12)
plot(p)
if(save.fig) dev.off()

save.image(file=analysis.file)
