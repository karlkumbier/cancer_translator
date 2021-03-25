library(data.table)
library(tidyverse)
library(xpose4)
library(rcellminer)
library(rcellminerData)
library(glmnet)
library(PRROC)
library(caret)
library(viridis)
library(parallel)
library(iRF)
library(superheat)
setwd('~/github/cancer_translator/')
source('scripts/utilities.R')

train.prop <- 0.8

#' The NCI60 data, describing bioactivity for 60 CCLs are available through 
#' rcellminer. Below we consider a subset of drug classes from the ORACL paper,
#' derived from matching MoA annotations in rcellminer against class annotations
#' in the ORACL data.  
#' 
#' *Note:* MOA data are available for 1040 compounds. The data we consider 
#' represent a subset of 261 of these compounds.
#' 
#' For each of these compounds, we load bioactivity data from the NCI60 cell 
#' line set as `x.activity`. We impute NA values as the median within cell lines.
#+ load_data

# Set compound key for select classes
class.key <- c(
  MT='Tu',
  DNA='(TOP2|Ds|Df|DHFR|AM|GARTF)',
  VEGFR='VEGFR',
  A7='A7',
  YK='YK',
  PIK3='PIK3',
  HDAC='HDAC',
  mTOR='MTOR',
  AKT='AKT'
)

# Read in MOA data from NCI60
nsc.set <- unique(c(getMoaToCompounds(), recursive=TRUE))
moa <- getMoaStr(nsc.set)

moa.split <- unlist(str_split(moa, ','))
moa.table <- data.frame(MoA=moa.split) %>%
  group_by(MoA) %>%
  summarize(Count=n()) %>%
  arrange(desc(Count))

write.csv(file='~/Desktop/moa_table.csv', moa.table, row.names=F, quote=F)

# Map MOA to compound class
class <- sapply(class.key, function(k) str_detect(moa, k))

# Drop multi- and 0- target compounds 
id.drop <- rowSums(class) != 1
class.id <- names(class.key)[apply(class[!id.drop,], MAR=1, which)]
class.id <- as.factor(class.id)
nsc.set <- nsc.set[!id.drop]

# Read in bioactivity data for compound set
x.activity <- getDrugActivityData(nsc.set)

# Impute NAs as median within cell line
x.activity <- apply(x.activity, MAR=2, median_impute)
print(dim(x.activity))

# Upsample for class imbalance
# If upsampling, we need to split train/test to not including training samples 
# in test set. Discuss up/down sampling - should models be giving equal weight 
# across drug classes?
xbalanced <- downSample(x.activity, class.id)
class.id <- xbalanced$Class
x.activity <- select(xbalanced, -Class)
print(dim(x.activity))

#' We visualize bioactivity patterns within our compound set below. Heatmap rows 
#' correspond to different compounds and columns cell lines. Compound classes 
#' cluster relative to PC1 and PC2 but show high variability.
#+ bioactivity_heatmap, fig.height=10, fig.width=15
# Bioactivity heatmap visualization
xplot <- apply(x.activity, MAR=2, rank)
type <- as.factor(str_remove_all(colnames(xplot), ':.*'))
type.pal <- scales::hue_pal()(length(unique(type)))

superheat(xplot, 
          membership.rows=class.id,
          membership.cols=type,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=inferno(100),
          bottom.label='variable',
          bottom.label.text.col=type.pal[type],
          bottom.label.size=0.5,
          bottom.label.text.angle=90)

# Bioactivity PCA plot visualization
xpca <- prcomp(apply(x.activity, MAR=2, rank))
data.frame(PctVar=cumsum(xpca$sdev ^ 2) / sum(xpca$sdev ^ 2)) %>%
  mutate(NPC=1:n()) %>%
  ggplot(aes(x=NPC, y=PctVar)) +
  geom_point() +
  geom_line() +
  theme_bw()

data.frame(xpca$x, Class=class.id) %>%
  ggplot(aes(x=PC1, y=PC2, col=Class)) +
  geom_point() +
  theme_bw()

#' ## Predictive modeling
#' Below, we consider the problem of predicting which compound class a given 
#' compound belongs to, and how many cell lines are required for prediction.
#+ supervised_modeling, warning=FALSE
y <- as.numeric(class.id) - 1
x <- x.activity
ncell <- 2:10
n.replicate <- 50

fits <- mclapply(1:n.replicate, function(i) {
  set.seed(i)
  
  # Generate sample split
  id.train <- createDataPartition(y, p=train.prop)$Resample1
  id.test <- setdiff(1:nrow(x), id.train)
  
  # Fit RF to predict moa from bioactivity
  out <- irf_fit(x, y, id.train, id.test, ncell)
  return(list(accuracy=out$accuracy, selected=out$selected))
}, mc.cores=10)

accuracy <- lapply(fits, function(f) as.matrix(f$accuracy))
selected <- lapply(fits, function(f) f$selected)

#+ supervised_modeling_accuracy, warning=FALSE, fig.height=10, fig.width=10
reshape2::melt(accuracy) %>%
  mutate(Ncell=as.factor(ncell[Var1])) %>%
  ggplot(aes(x=Ncell, y=value)) +
  geom_boxplot(fill='#0388D1') +
  theme_bw() +
  ylab("Top-1 classification accuracy") +
  ggtitle('Model accuracy by cell lines used')

#+ supervised_modeling_beta, warning=FALSE, fig.height=18, fig.width=12
reshape2::melt(selected) %>%
  rename(CellLine=Var1, Ncell=Var2, Rep=L1) %>%
  mutate(Ncell=ncell[Ncell]) %>%
  group_by(CellLine, Ncell) %>%
  summarize(Proportion=mean(value)) %>%
  filter(Proportion != 0) %>%
  ggplot(aes(x=reorder(CellLine, Proportion), y=Proportion)) +
  facet_wrap(~Ncell) +
  geom_bar(stat='identity', aes(fill=CellLine)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position='none') +
  ggtitle('Cell line selection proportion')

#+ supervised_modeling_heatmap, warning=FALSE, fig.height=12, fig.width=12
# Plot distribution of compounds from select cell lines
ncell.plot <- 5
cell.select <- reshape2::melt(selected) %>%
  rename(CellLine=Var1, Ncell=Var2) %>%
  mutate(Ncell=ncell[Ncell])
  
cell.select.k <- cell.select %>%
  filter(Ncell == ncell.plot) %>%
  group_by(CellLine, Ncell) %>%
  summarize(Proportion=mean(value)) %>%
  ungroup() %>%
  top_n(ncell.plot, Proportion)

selected <- c(as.character(cell.select.k$CellLine), 'LC:A549/ATCC')
xplot <- apply(x[,selected], MAR=2, rank)
id.highlight <- colnames(xplot) == 'LC:A549/ATCC'

superheat(xplot, 
          membership.rows=class.id,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=inferno(100),
          bottom.label.size=1,
          bottom.label.text.angle=90,
          bottom.label.text.col=ifelse(id.highlight, 'red', 'black'))

#' ## Null predictive modeling - permuted responses
#' To assess performance against random chance, we repeat the experiment above 
#' using randomly permuted response labels.
#+ random_permute, fig.height=10, fig.width=15
fits.null <- mclapply(1:n.replicate, function(i) {
  set.seed(i)
  
  # Generate sample split
  id.train <- createDataPartition(y, p=train.prop)$Resample1
  id.test <- setdiff(1:nrow(x), id.train)
  
  # Fit lasso to predict moa from bioactivity
  yy <- sample(y)
  out <- irf_fit(x, yy, id.train, id.test, ncell)
  return(out$accuracy)
}, mc.cores=10)

accuracy <- lapply(fits.null, as.matrix)
reshape2::melt(accuracy) %>%
  mutate(Ncell=as.factor(ncell[Var1])) %>%
  ggplot(aes(x=Ncell, y=value)) +
  geom_boxplot(fill='#0388D1') +
  theme_bw() +
  ylab("Top-1 classification accuracy") +
  ggtitle('Model accuracy by cell lines used - permuted responses')

#' ## Null predictive modeling - random cell lines
#' To assess the performance against randomly selected cell lines, we repeat 
#' the experiment above using a random sample of $k$ cell lines for $k = 1, ..., 
#' 20$.
#+ random cell line
# Drop single cell due to modeling constraints
fits.null <- mclapply(1:n.replicate, function(i) {
  set.seed(i)
  
  # Generate sample split
  id.train <- createDataPartition(y, p=train.prop)$Resample1
  id.test <- setdiff(1:nrow(x), id.train)
  
  # Fit RF to predict moa from bioactivity
  out.null <- irf_fit(x, y, id.train, id.test, ncell, null.cells=TRUE)
  out.opt <-  irf_fit(x, y, id.train, id.test, ncell)
  
  return(rbind(out.null$accuracy, out.opt$accuracy))
}, mc.cores=10)


reshape2::melt(fits.null) %>%
  rename(Selection=Var1, Ncell=Var2, Rep=L1) %>%
  mutate(Selection=ifelse(Selection == 1, 'Random', 'Optimal')) %>%
  mutate(Ncell=Ncell + 1) %>%
  mutate(Ncell=as.factor(Ncell)) %>%
  group_by(Ncell, Rep) %>%
  summarize(Diff=value[Selection == 'Optimal'] - value[Selection == 'Random']) %>%
  ggplot(aes(x=Ncell, y=Diff)) +
  geom_boxplot(fill='#0388D1') +
  theme_bw()
