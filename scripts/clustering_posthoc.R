#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)
library(hrbrthemes)

theme_set(theme_ipsum(axis_title_size=14, plot_title_size=14, strip_text_size=14))

################################################################################
# Initialize constant parameters to be used in analysis
################################################################################

# Initialize color palettes
heat.pal <- viridis::viridis(10)

################################################################################
# Load results
################################################################################
# Initialize input/output directories and load KS profile data
analysis.dir <- '~/github/cancer_translator/'
result.dir <- str_c(analysis.dir, 'results/clustering/')
files <- list.files(result.dir)

x.cluster <- lapply(files, function(f) {
  load(str_c(result.dir, f))
  return(left_join(x.cluster.tab, loss, by='Cluster'))
})

x.cluster.full <- rbindlist(x.cluster) %>%
  mutate(Loss=(Silhouette * Prec * Rec) ^ (1/3))

x.cluster <- group_by(x.cluster.full, CellLine, Compound_Category) %>%
  summarize(Silhouette=max(Silhouette), Prec=max(Prec), Rec=max(Rec), Loss=max(Loss)) %>%
  ungroup()

#' # Loss by compound category, cell line
#+ cluster_loss, fig.height=12, fig.width=16
################################################################################
# Visualize heatmap of clustering strength
################################################################################
compounds <- unique(x.cluster$Compound_Category)
cell.lines <- unique(x.cluster$CellLine)

# Loss plot - geometric mean of 3 categories
xplot <- t(matrix(x.cluster$Loss, nrow=length(compounds)))
colnames(xplot) <- compounds
rownames(xplot) <- cell.lines

avg.loss <- rowMeans(xplot)
max.loss <- apply(xplot, MAR=2, max)
xplot <- xplot[,max.loss > (0.05) ^ (1 / 3)]

row.order <- sort(avg.loss, decreasing=TRUE) %>% names

col.order <- order(
  xplot[row.order[1],],
  xplot[row.order[2],],
  xplot[row.order[3],],
  xplot[row.order[4],],
  xplot[row.order[5],],
  xplot[row.order[6],]
)

superheat(
  xplot[rev(row.order), col.order], 
  yr=sort(avg.loss),
  yr.plot.type='bar',
  yr.axis.name='Average loss',
  yr.axis.name.size=18,
  bottom.label.text.angle=90,
  bottom.label.size=0.75,
  heat.pal.values=seq(0, 1, by=0.1)
)

#' # Loss and metrics by cluster, cell line
#+ loss_metrics, fig.height=12, fig.width=16
x.cluster.full <- x.cluster.full %>% 
  group_by(CellLine, Cluster) %>%
  summarize(
    Silhouette=max(Silhouette), 
    Prec=max(Prec), 
    Rec=max(Rec), 
    Loss=max(Loss), 
    Category=Compound_Category[which.max(Loss)]
  ) %>%
  ungroup() %>%
  mutate(Cluster=str_c(Cluster, ', ', Category))

ggplot(x.cluster.full, aes(x=Cluster, y=Loss, fill=Category)) + 
  geom_bar(stat='identity')  +
  facet_wrap(~CellLine, scales='free_y') +
  theme(legend.position='none') +
  coord_flip() +
  xlab(NULL)

reshape2::melt(x.cluster.full, id.vars=c('CellLine', 'Cluster', 'Category')) %>%
  ggplot(aes(x=Cluster, y=value, fill=variable)) + 
  geom_bar(stat='identity', position='dodge')  +
  facet_wrap(~CellLine, scales='free_y') +
  scale_fill_nejm() +
  coord_flip() +
  xlab(NULL)
  