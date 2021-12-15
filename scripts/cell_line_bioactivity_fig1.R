#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)
library(umap)

intensity.normalize <- TRUE
n.core <- 16
min.cat <- 5

analysis.dir <- '~/github/cancer_translator/'
setwd(analysis.dir)
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'

output.dir <- 'results/cell_line/'
dir.create(output.dir, showWarnings=FALSE)
output.file <- str_c(output.dir, 'bioactivity.Rdata') 

intensity.normalize <- TRUE
load(str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata'))
load(str_c(analysis.dir, 'data/smiles_table.Rdata'))

# Clean compound category names
select_category <- function(x) na.omit(x)[1]

# Clean compound category/target/pathways
x <- group_by(x, Compound_ID) %>%
  mutate(Compound_Category=select_category(Compound_Category)) %>%
  mutate(Target=select_category(Target)) %>%
  mutate(Pathway=select_category(Pathway)) %>%
  ungroup()

x <- left_join(x, xsmiles, by='Compound_ID')

# Initialize color palettes
heat.pal <- c('#FFFFFF', pal_material("indigo")(10))
heat.pal.signed <- RColorBrewer::brewer.pal(9, 'RdBu')

#' # Overview
#' This notebook considers the problems of identifying bioactive compounds 
#' based on a panel of 6 cell lines: OVCAR4, A549, DU145, ALS-WT, HEPG2, and 
#' 786-0. In addition to individual cell lines, we consider whether bioactivity 
#' calls can be improved by aggregating information across cell lines.
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

# Generate compound table
# Note: for compounds with multiple category/pathway/target annotations, 
# we take the first non-na value.
xtab <- rbindlist(xdist) %>%
  group_by(Compound_ID) %>%
  group_by(Cell_Line, Compound_ID, Dose, Dose_Category, Compound_Category, Target, Pathway, SMILES) %>% 
  summarize(Bioactive=(mean(pval) == 0), Dist=mean(DistNorm), ID=list(ID), .groups='drop') %>%
  mutate(NID=sapply(ID, length)) %>%
  filter(NID < 10) %>%
  select(-NID) %>%
  mutate(ID=sapply(ID, str_c, collapse=', ')) %>%
  arrange(desc(Dist))

write.csv(file='~/Desktop/bioactivity.csv', xtab, quote=FALSE, row.names=FALSE)

# Summarize bioactivity by compound category
xtab.cat <- filter(xtab, !is.na(Compound_Category)) %>%
  group_by(Compound_Category) %>%
  summarize(PropBioactive=mean(Bioactive), NumBioactive=sum(Bioactive), N=n()) %>%
  arrange(desc(PropBioactive))

write.csv(file='~/Desktop/bioactivity_category.csv', xtab.cat, quote=FALSE, row.names=FALSE)

# Summarize bioactivity by pathway
xtab.pw <- filter(xtab, !is.na(Pathway)) %>%
  group_by(Pathway) %>%
  summarize(PropBioactive=mean(Bioactive), NumBioactive=sum(Bioactive), N=n()) %>%
  arrange(desc(PropBioactive))

write.csv(file='~/Desktop/bioactivity_pathway.csv', xtab.pw, quote=FALSE, row.names=FALSE)

# Summarize bioactivity by target
xtab.target <- filter(xtab, !is.na(Target)) %>%
  group_by(Target) %>%
  summarize(PropBioactive=mean(Bioactive), NumBioactive=sum(Bioactive), N=n()) %>%
  arrange(desc(PropBioactive))

write.csv(file='~/Desktop/bioactivity_target.csv', xtab.target, quote=FALSE, row.names=FALSE)


#' Here we assess bioactivity calls by compound/cell line. We take the minimum
#' p-value across all dose levels of a given treatment, effectively asking
#' whether a cell line detects bioactivity at some dose level.
# Group bioactivity scores by treatment, cell line
n.cell.line <- length(unique(x$Cell_Line))
xgroup <- rbindlist(xdist) %>%
  group_by(Cell_Line, Compound_ID, Dose_Category, Compound_Category) %>% 
  summarize(DistNorm=mean(DistNorm), pval=mean(pval), .groups='drop') %>%
  group_by(Cell_Line, Compound_ID, Compound_Category) %>%
  summarize(DistNorm=max(DistNorm), pval=min(pval), .groups='drop') %>%
  group_by(Compound_ID) %>% #, Dose_Category) %>%
  mutate(Count=n()) %>%
  filter(Count == n.cell.line) %>%
  ungroup() %>%
  arrange(Compound_ID, Cell_Line)

#' ### Fig 1. Bioactivity by cell line. 
#' Bioactivity calls by compound), cell line. Compounds with p-value = 0 (for > 
#' 0 doses) are defined as bioactive, while compounds with p-value = 0 (for all 
#' doses) are defined as inactive.
#' 
#' **Note:** treatment (compound/dose pairs) replicates are grouped and p-values 
#' averaged before calling bioactivity. Thus, bioactive treatments correspond to 
#' those with p-value = 0 across all replicates.
#+ bioactivity_p, fig.height=12, fig.width=18
# Format cell x compound bioactivity for visualization
xplot <- matrix(as.numeric(xgroup$pval == 0), nrow=n.cell.line)
rownames(xplot) <- unique(xgroup$Cell_Line)
colnames(xplot) <- unique(str_c(xgroup$Compound_ID))#, ', ', xgroup$Dose_Category))
category <- matrix(xgroup$Compound_Category, nrow=n.cell.line)[1,]

# Filter to compounds with bioactive calls in > 0 cell lines
id.drop <- colMeans(xplot == 0) == 1
xplot <- xplot[,!id.drop]
category <- category[!id.drop]

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
  heat.pal=heat.pal,
  yr=prop.bioactive[bioactive.order],
  yr.plot.type='bar',
  yr.axis.name='Proportion bioactive',
  yr.axis.name.size=18,
  heat.pal.values=seq(0, 1, by=0.1),
  title='Bioactivity by cell line, compound\nall compounds'
)

#' ### Fig 2. Bioactivity by cell line, compound category.
#' Bioactivity calls by cell line, compound category. Treatments (compound dose 
#' pairs) with p-value = 0 are defined as bioactive, while treatments with 
#' p-value < 0 are defined as inactive. Top: bioactivity calls for each 
#' treatment, bottom: proportion of bioactivity calls within each compound 
#' category. Compound categories with fewer than `$r min.cat` treatments are 
#' dropped for visualization.
#' 
#' **Note:** treatment replicates are grouped and p-values averaged before 
#' calling bioactivity. Thus, bioactive treatments correspond to those with 
#' p-value = 0 across all replicates.
#+ bioactivity_category, fig.height=12, fig.width=18, warning=FALSE
# Plot bioactivity by cell line, compound category
cat.table <- table(category)
cat.keep <- setdiff(names(cat.table[cat.table >= min.cat]), 'Others')
id.keep <- (category %in% cat.keep)

xplot.cat <- sapply(unique(category[id.keep]), function(ctg) {
  rowMeans(xplot[,which(category == ctg)])
})

# reorder for visualization
xplot.cat <- xplot.cat[,order(colMeans(xplot.cat))]
xplot.cat[xplot.cat < 0.25] <- 0.25

avg.bioactive <- rowMeans(xplot.cat)
row.order <- order(avg.bioactive)

superheat(
  xplot.cat[row.order,],
  yt=c(cat.table[colnames(xplot.cat)]),
  yt.plot.type='bar',
  yt.axis.name='# of compounds',
  yt.axis.name.size=12,
  yr=avg.bioactive[row.order],
  yr.plot.type='bar',
  yr.axis.name='Average proportion biactivity',
  yr.axis.name.size=12,
  bottom.label.text.angle=90,
  bottom.label.size=0.75,
  heat.pal=heat.pal, 
  heat.pal.values=seq(0, 1, by=0.1),
  title='Bioactivity by cell line, compound\nprevalent compound categories'
)

#' # Compound clustering
#' To assess the similarity between different compound classes, we compute 
#' distances between compounds based on their phenotypic profiles and assess
#' the similarity of compounds within the same category versus between 
#' categories. Prior to assessing compound similarity, we filter to treatments
#' (compound/dose pairs) identified as bioactive.
### Fig 3 Compound clustering
#+ compound_clustering, fig.height=12, fig.width=18, warning=FALSE
compound_similarity <- function(xdist, categories) {
  # Compute within/between class similarity by compound category
  xdist <- as.matrix(xdist)
  score <- sapply(unique(categories), function(ctg) {
    dist.within <- mean(xdist[categories == ctg, categories == ctg])
   
    dist.between <- sapply(setdiff(unique(categories), ctg), function(btw) {
      mean(xdist[categories == ctg, categories != btw])
    })
    
    return(dist.within / min(dist.between))
  })
  
  return(score)
}

# Initialize bioactive compound set
compounds.select <- filter(xtab, Bioactive) %>%
  mutate(Compound=str_c(Compound_ID, ', ', Dose_Category))

scores <- vector('list', length(unique(x$Cell_Line)))
names(scores) <- unique(x$Cell_Line)

for (cl in unique(x$Cell_Line)) {

  # Filter to bioactive compounds, select cell line
  xc <- filter(x, Cell_Line == cl) %>%
    mutate(Compound = str_c(Compound_ID, ', ', Dose_Category)) %>%
    filter(Compound %in% compounds.select$Compound) %>%
    filter(Compound_Category %in% cat.keep) %>%
    select(-matches('^PC')) %>%
    group_by(Compound_ID, Dose_Category, Compound_Category, Pathway, Target) %>%
    summarize_if(is.numeric, mean) %>%
    ungroup()

  # Compute within/between class similarity
  xfeat <- select(xc, matches('^nonborder'))
  colnames(xfeat) <- NULL
  scores[[cl]] <- compound_similarity(dist(xfeat), xc$Compound_Category)
  
  ncat <- length(unique(xc$Compound_Category))
  cat.pal <- scales::hue_pal(l=60)(ncat)
  row.cols <- cat.pal[as.factor(xc$Compound_Category)]

  # Cluster compounds by class
  heatmap(
    as.matrix(xfeat),
    col=heat.pal.signed,
    RowSideColors=row.cols,
    main=cl
  )
  
  cat <- levels(as.factor(xc$Compound_Category))
  legend('topright', fill=cat.pal, legend=cat)
}


# initialize data matrix for heatmap
scores <- lapply(scores, as.matrix)
scores <- reshape2::melt(scores)

category <- unique(scores$Var1)
cell.line <- unique(scores$L1)

xplot <- matrix(scores$value, nrow=length(category))
rownames(xplot) <- category
colnames(xplot) <- cell.line
xplot[xplot > 1] <- 1

# Plot bioactivity by treatment, cell line
superheat(
  xplot,
  pretty.order.rows=TRUE,
  pretty.order.cols=TRUE,
  heat.pal=heat.pal,
  heat.pal.values=c(0, 0.75, seq(0.75, 1, by=0.05))
)
