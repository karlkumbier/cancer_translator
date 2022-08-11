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
enrich.thresh <- 0.3

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
  mutate(Loss=((Silhouette > 0) * Prec * Rec) ^ (1/2))

# Compute maximal enrichment by cell line, category
x.cluster.cl <- x.cluster %>%
  mutate(Enrichment=Enrichment * (Silhouette > 0)) %>%
  group_by(CellLine, Compound_Category) %>%
  summarize(
    Precision=max(Precision), 
    Recall=max(Recall), 
    Enrichment=max(Enrichment),
    .groups='drop'
  )

#' # Enrichment by compound category, cell line
#' To assess the ability of different cell lines to detect different compound 
#' categories, we consider our previously defined enrichment score — the 
#' geometric mean of precision, recall over clusters with silhouette > 0 — by 
#' compound category and cell line. We enforce the silhouette > 0 constraint by 
#' including an indcator term in the geometric mean.
#' 
#' The figure below reports the maximum enrichment over all clusters within each 
#' cell line and compound category pair (i.e. is the category detected in any 
#' "strong" cluster). We filter to compound categories with enrichment > 0.5
#' in at least 1 cell line.
#+ cluster_enrich, fig.height=12, fig.width=16
################################################################################
# Visualize heatmap of clustering strength
################################################################################
compounds <- unique(x.cluster.cl$Compound_Category)
cell.lines <- unique(x.cluster.cl$CellLine)

# Initialize data to be plotted
xplot <- t(matrix(x.cluster.cl$Enrichment, nrow=length(compounds)))
colnames(xplot) <- compounds
rownames(xplot) <- cell.lines

avg.loss <- rowMeans(xplot)
max.loss <- apply(xplot, MAR=2, max)
xplot <- xplot[,max.loss > (0.25) ^ (1 / 2)]

# Filter to categories enriched in at least one cell line
max.enrich <- apply(xplot, MAR=2, max)
xplot <- xplot[,max.enrich > enrich.thresh]

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
  yr=sort(avg.enrich),
  yr.plot.type='bar',
  yr.axis.name='Average enrichment',
  yr.axis.name.size=18,
  bottom.label.text.angle=90,
  bottom.label.size=0.75,
  heat.pal.values=c(seq(0, 0.1, by=0.1), seq(0.15, 1, by=0.05)),
  title='Enrichment by cell line, compound category'
)

#' # Loss and metrics by cluster, cell line
#' The top figure below reports loss by cell line and cluster. For each cluster,
#' we report the category associated with each cluster — i.e. the maximal 
#' cluster. The bottom figure breaks down the cluster/cell line loss by each
#' metric (precision, recall, silhouette) included in the geometric mean.
#+ loss_metrics, fig.height=12, fig.width=16
# Initizlie cluster table
x.cluster <- lapply(files, function(f) {
  load(str_c(result.dir, f))
  return(mutate(loss, CellLine=str_remove_all(f, '.Rdata')))
})

x.cluster <- rbindlist(x.cluster) %>%
  mutate(Cluster=str_c(Cluster, ', ', Category)) %>%
  mutate(Enrichment=Enrichment * (Silhouette > 0))

filter(x.cluster, Silhouette > 0) %>%
  ggplot(aes(x=reorder_within(Cluster, Enrichment, CellLine), y=Enrichment)) + 
  geom_bar(stat='identity', aes(fill=Category))  +
  geom_hline(yintercept=enrich.thresh, lty=2, col='grey') +
  facet_wrap(~CellLine, scales='free_y') +
  theme(legend.position='none') +
  coord_flip() +
  xlab(NULL)

# Plot loss by each sub component
id.vars <- c('CellLine', 'Cluster', 'Category')
xplot <- mutate(x.cluster, Cluster=str_c(Cluster, '__', CellLine)) %>%
  group_by(CellLine) %>%
  mutate(Cluster=factor(Cluster, levels=Cluster[order(Enrichment)])) %>%
  ungroup() %>%
  dplyr::select(-PrecisionFull, -RecallFull) %>%
  filter(Enrichment > 0) %>%
  reshape2::melt(id.vars=id.vars)

metrics <- c('Enrichment', 'Precision', 'Recall', 'Silhouette')
mutate(xplot, variable=factor(variable, levels=metrics)) %>%
  ggplot(aes(x=Cluster, y=value, fill=variable)) + 
  geom_bar(stat='identity', position='dodge')  +
  facet_grid(CellLine~variable, scales='free_y') +
  scale_fill_nejm() +
  coord_flip() +
  xlab(NULL) +
  theme(legend.position='none')
  
#' # Aggregating
#' Finally, we ask to whether combining information across cell lines allows
#' one to identify a larger set of phenotypically similar compounds associated
#' with specific compound categories. We filter to cell line/cluster
#' combinations with enrichment > `0.5`r enrich.thresh`` and select compounds 
#' belonging to these clusters.
#+ aggregate_cluster, fig.height=8, fig.width=12
enriched <- filter(x.cluster, Enrichment > enrich.thresh) %>%
  mutate(Cluster=as.numeric(str_remove_all(Cluster, ', .*$')))

x.compound <- lapply(files, function(f) {
  
  load(str_c(result.dir, f))
  
  # Filter to enriched clusters for cell line
  enriched.f <- filter(enriched, CellLine == x.cluster.tab$CellLine[1])
  
  # Filter to compounds in enriched clusters
  out <- filter(x.cluster.tab, Cluster %in% enriched.f$Cluster)
  return(out)
})

# Aggregate selected compounds across cell lines
x.compound <- rbindlist(x.compound)
cell.order <- sort(avg.enrich, decreasing=TRUE) %>% names

compound.set <- vector('list', length(cell.order))
names(compound.set) <- cell.order
for (cl in cell.order) {
  x.compound.f <- filter(x.compound, CellLine == cl)
  cpds <- setdiff(x.compound.f$Compound_ID, unlist(compound.set))
  compound.set[[cl]] <- cpds
}

# Aggregate compound count for visualization
xplot <- reshape2::melt(compound.set) %>%
  rename(CellLine=L1) %>%
  group_by(CellLine) %>%
  count() %>%
  arrange(desc(n)) %>%
  ungroup()

mutate(xplot, CellLine=factor(CellLine, levels=CellLine)) %>%
  mutate(NCpd=cumsum(xplot$n)) %>%
  ggplot(aes(x=CellLine, y=NCpd)) +
  geom_bar(stat='identity')
