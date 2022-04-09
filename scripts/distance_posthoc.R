#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)
library(RColorBrewer)
library(hrbrthemes)

theme_set(
  theme_ipsum(
    axis_title_size=18, 
    strip_text_size=18, 
    axis_text_size=14,
    base_size=18, 
    base_family='sans'
  )
)

setwd('~/github/cancer_translator/')
source('scripts/utilities.R')

################################################################################
# Initialize constant parameters to be used in analysis
################################################################################
# Initialize color palettes
heat.pal <- viridis::viridis(10)
heat.pal <- brewer.pal(9, 'RdBu')
thresh.plot <- 0.3

# Initialize input/output directories and load KS profile data
analysis.dir <- '~/github/cancer_translator/'
result.dir <- str_c(analysis.dir, 'results/distance/')
fig.dir <- str_c(analysis.dir, 'results/figures/fig2/')

################################################################################
# Load results
################################################################################
files <- list.files(result.dir)

x.similarity <- lapply(files, function(f) {
  load(str_c(result.dir, f))
  return(dist.category)
})

# Sort by category/cell line
x.similarity <- rbindlist(x.similarity) %>% arrange(Category, Cell_Line)

#' ### Fig 2b. Between MOA similarity
#' To assess the ability of different cell lines to phenotypically group 
#' compounds by MOA, we consider the distance between compounds (in phenotypic 
#' space) with the same MOA relative to nearby compounds with different MOA. The 
#' figure below reports this relative distance by MOA and cell line.
#+ moa_sim_between, fig.height=12, fig.width=16
################################################################################
# Visualize heatmap of clustering strength
################################################################################
categories <- unique(x.similarity$Category)
cell.lines <- unique(x.similarity$Cell_Line)

# Initialize data to be plotted
xplot <- matrix(x.similarity$DistBetween, nrow=length(cell.lines))
colnames(xplot) <- categories
rownames(xplot) <- cell.lines
colnames(xplot) <- str_remove_all(colnames(xplot), ',.*$')

# Rank compounds by average MOA similarity
avg.sim <- rowMeans(xplot)
row.order <- order(avg.sim, decreasing=TRUE)

# Filter to compounds in top quartile
thresh <- quantile(xplot, 0.75)
xplot <- xplot[,colSums(xplot > thresh) > 1]
xplot <- xplot[,colnames(xplot) != 'DMSO']
xplot <- xplot[,colnames(xplot) != 'Others']

# Filter to categories enriched in at least one cell line
col.order <- order(
  xplot[row.order[1],],
  xplot[row.order[2],],
  xplot[row.order[3],],
  xplot[row.order[4],],
  xplot[row.order[5],],
  xplot[row.order[6],]
)

superheat(
  xplot[row.order, col.order], 
  bottom.label.text.angle=90,
  bottom.label.size=0.75
)

#' ### Fig 2c. Within MOA similarity
#' To assess the ability of different cell lines to phenotypically group 
#' compounds by MOA, we consider the distance between compounds (in phenotypic 
#' space) with the same MOA relative to the distance between DMSO compounds. The 
#' figure below reports this relative distance by MOA and cell line.
#+ moa_sim_within, fig.height=12, fig.width=16

# Initialize data to be plotted
xplot <- matrix(x.similarity$DistWithin, nrow=length(cell.lines))
colnames(xplot) <- categories
rownames(xplot) <- cell.lines

# Filter to compounds in top quartile
thresh <- quantile(xplot, 0.75)
xplot <- xplot[,colSums(xplot > thresh) > 1]
xplot <- xplot[,colnames(xplot) != 'DMSO']
xplot[xplot < -thresh.plot] <- -thresh.plot
xplot[xplot > thresh.plot] <- thresh.plot

# Filter to categories enriched in at least one cell line
col.order <- order(
  xplot[row.order[1],],
  xplot[row.order[2],],
  xplot[row.order[3],],
  xplot[row.order[4],],
  xplot[row.order[5],],
  xplot[row.order[6],]
)

superheat(
  xplot[row.order, col.order], 
  bottom.label.text.angle=90,
  bottom.label.size=0.75,
  heat.pal=heat.pal,
  heat.lim=c(-thresh.plot, thresh.plot)
)

#' ### Fig 2d.  Optimal cell line selection
#' Below we compute the optimal cell lines for detecting within and between 
#' compound similarity. For between similarity, we use MOA weights to consider 
#' only "strongly clustering compounds" — at least one cell line > 75th 
#' percentile. For Within compound similarity, we consider all MOAs weighted 
#' equally.
#+ moa_opt, fig.height=12, fig.width=16
xplot.between <- matrix(x.similarity$DistBetween, nrow=length(cell.lines))
xplot.within <- matrix(x.similarity$DistWithin, nrow=length(cell.lines))
xplot <- (xplot.between + xplot.within) / 2
colnames(xplot) <- categories
rownames(xplot) <- cell.lines

wt <- rep(1, ncol(xplot))
names(wt) <- colnames(xplot)

top.category <- setdiff(categories[colSums(xplot.between > 0.75) > 0], 'DMSO')
wt.specialist <- wt
wt.specialist[!names(wt) %in% top.category] <- 0

# Compute optimal cell line set by budget
x.generalist <- lapply(1:3, function(k) {
  opt_similarity(xplot, names(wt), wt, k)
})

x.specialist <- lapply(1:3, function(k) {
  opt_similarity(xplot, names(wt.specialist), wt.specialist, k)
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
  ylab('Similarity score') +
  scale_color_nejm() +
  scale_fill_nejm() +
  theme(legend.position='none') +
  facet_wrap(~Type, ncol=1)

if (save.fig) pdf(str_c(fig.dir, 'selection.pdf'), height=8, width=10)
plot(p)  
if(save.fig) dev.off()
