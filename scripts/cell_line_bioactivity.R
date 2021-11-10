#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)

intensity.normalize <- TRUE
n.core <- 6

setwd('~/github/cancer_translator/')
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

output.dir <- 'results/cell_line/'
dir.create(output.dir, showWarnings=TRUE)
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
#' # Bioactivity filtering
#' We compute bioactivity scores as follows:
#'
#' 1. Subsample DMSO wells 
#' 2. Compute the center of the DMSO point cloud (i.e. average feature values 
#' across each well in the subsample) 
#' 3. Compute the l2 distance between each DMSO well in the subsample and the 
#' DMSO point cloud center 
#' 4. Define the DMSO maximal distance as the maximum distance (from 3) over all 
#' subsampled wells.
#' 
#' We then ask whether a given well compound/dose is further from the DMSO
#' point cloud center than the maximal DMSO distance. Repeating this process
#' across many subsamples allows us to generate a bioactivity p-value. 
#' 
#' #' **Note:** Bioactivity scores are computed within plate.
#+ bioactivity, fig.height=8, fig.width=14

################################################################################
# Compute bioactivity scores
################################################################################
setwd('~/github/cancer_translator/')

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
#' compound/dose pairs evaluated across all cell lines. Each point 
#' corresponds to a single compound/dose.
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
  scale_color_manual(values=cell.pal)

#' To assess differential activity across cell lines, we report bioactivity, 
#' measured as the relative distance from DMSO, for each compound/dose 
#' combination. Figures below show this distance for (i) all compounds (ii) 
#' compound categories with more than 20 instances (not including "Other").
#+ bioactivity_dist, fig.height=12, fig.width=18
# Format cell x compound bioactivity for visualization
xplot <- matrix(xgroup$DistNorm, nrow=n.cell.line)
rownames(xplot) <- unique(xgroup$Cell_Line)
colnames(xplot) <- unique(str_c(xgroup$Compound_ID, ', ', xgroup$Dose_Category))
category <- matrix(xgroup$Compound_Category, nrow=n.cell.line)[1,]

# Threshold for visualization
xplot.t <- xplot
xplot.t[xplot.t < 1] <- 1
xplot.t[xplot.t > 3] <- 3

superheat(
  xplot.t,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  heat.pal=heat.pal, 
  heat.pal.values=seq(0, 1, by=0.1),
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
  title='Bioactivity by cell line, compound (select compounds)\nl2 distance from DMSO center relative to DMSO wells'
)

#' Finally, we report the proportion of bioactive treatments (compound/dose 
#' pairs) for reference compounds by cell line and across the full cell line 
#' set. For an individual cell line, we define treatments with p-value = 0 as 
#' bioactive. For the "full" cell line set, we define any treatment with that is 
#' bioactive in at least one cell line as bioactive (i.e. minimum p-value = 0).
#+ bioactivity_treat, fig.height=10, fig.width=18, warning=FALSE
# Initialize cell line sets
cell.lines <- unique(xdist$Cell_Line)

cell.sets <- lapply(1:n.cell.line, function(k) {
  combn(cell.lines, k, simplify=FALSE)
})

cell.sets <- unlist(cell.sets, recursive=FALSE)

# Compute bioactive relative to each cell line set
bioactive.cat <- lapply(cell.sets, function(s) {
  out <- filter(xdist, Cell_Line %in% s) %>%
    group_by(Compound_ID, Compound_Category, Dose_Category, Compound_Usage) %>%
    summarize(pval=min(pval), .groups='drop') %>%
    group_by(Compound_Category, Compound_Usage) %>%
    summarize(PropBioactive=mean(pval == 0), Count=n(), .groups='drop') %>%
    mutate(Cell_Line=str_c(s, collapse=', ')) %>%
    select(Compound_Category, Cell_Line, PropBioactive, Compound_Usage, Count) %>%
    mutate(Ncells=length(s))
})
   
cpd.count <- filter(rbindlist(bioactive.cat), Cell_Line =='A549')$Count
cpd.count <- sum(cpd.count)
print(str_c('# Treatments: ', cpd.count))
  
rbindlist(bioactive.cat) %>%
  group_by(Cell_Line) %>%
  summarize(PropBioactive=mean(PropBioactive), Ncells=mean(Ncells)) %>%
  mutate(PropText=round(PropBioactive, 2)) %>%
  ggplot(aes(x=reorder(Cell_Line, PropBioactive), y=PropBioactive)) +
  geom_bar(stat='identity', aes(fill=Ncells)) +
  geom_text(aes(label=PropText), nudge_y=0.01, size=3) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  scale_fill_gradientn(colors=ncell.pal[-1]) +
  ggtitle('Proportion of bioactive calls by cell line set')

#+ bioactivity_treat_heat, fig.height=12, fig.width=24, warning=FALSE
################################################################################
# Bioactivity heatmap by cell line/category
################################################################################
# Initialize bioactivity table
bioactive.cat.t <- rbindlist(bioactive.cat) %>% 
  arrange(Cell_Line, Compound_Category) %>%
  filter(Compound_Usage == 'reference_cpd')

# Format data for heatmap
n.cell.line <- length(unique(bioactive.cat.t$Cell_Line))
xplot <- matrix(bioactive.cat.t$PropBioactive, nrow=n.cell.line, byrow=TRUE)

rownames(xplot) <- unique(bioactive.cat.t$Cell_Line)
colnames(xplot) <- unique(bioactive.cat.t$Compound_Category)

# Initialize row attributes
ncell <- select(bioactive.cat.t, Cell_Line, Ncells) %>% distinct()
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
  heat.pal=heat.pal
)

# Bioactivity table by cell line/category
setwd('~/github/cancer_translator/')
write.csv(file='results/bioactivity_cell_category.csv', bioactive.cat, quote=FALSE)

# Generate table of bioactivity by plate
bioactive.plate <- group_by(xdist, PlateID) %>% 
  summarize(PropBioactive=mean(pval == 0))

bioactive.plate %>%
  ggplot(aes(x=reorder(PlateID, PropBioactive), y=PropBioactive)) +
  geom_bar(stat='identity') +
  theme_bw() +
  theme(legend.position='none') +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Bioactivity by plate')

knitr::kable(bioactive.plate)
write.csv(file='results/bioactivity_plate.csv', bioactive.plate, quote=FALSE)
save(file=output.file, xdist)