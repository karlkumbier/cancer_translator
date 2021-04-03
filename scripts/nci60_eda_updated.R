#+ setup, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(readxl)
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
model <- irf_fit

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

# Read in MOA data from NCI60
nsc.set <- unique(c(getMoaToCompounds(), recursive=TRUE))
moa <- getMoaStr(nsc.set)

# Generate MoA table for Louise annotations
moa.split <- unlist(str_split(moa, ','))
moa.table <- data.frame(MoA=moa.split) %>%
  group_by(MoA) %>%
  summarize(Count=n()) %>%
  arrange(desc(Count))

# Read in Louise compound annotation table and clean for class ID
setwd('~/github/cancer_translator/')
moa.louise <- read_excel('data/nci60/moa_table_LH.xlsx') %>%
  filter(!(is.na(comment) & is.na(reference_category_setof6))) %>%
  filter(!(comment == 'cant ID' & is.na(reference_category_setof6))) %>%
  mutate(CleanMoA=str_remove_all(MoA, 'PK:')) %>%
  group_by(CleanMoA) %>%
  summarize(Class=list(reference_category_setof6), MoA=list(MoA), Count=sum(Count))

# Group by compound class, using clean MoA when no class indicated
moa.louise <- mutate(moa.louise, Class=sapply(Class, clean_class))

moa.class <- filter(moa.louise, !is.na(Class)) %>%
  group_by(Class) %>%
  summarize(CleanMoA=list(CleanMoA), MoA=list(unlist(MoA)), Count=sum(Count))

moa.noclass <- filter(moa.louise, is.na(Class)) %>%
  mutate(Class=CleanMoA) %>%
  group_by(Class) %>%
  summarize(CleanMoA=list(CleanMoA), MoA=list(unlist(MoA)), Count=sum(Count))

moa.table <- rbind(moa.class, moa.noclass) %>% arrange(desc(Count)) %>%
  filter(Class != 'Not Available')

# Clean MoA table for output
moa.table.clean <- moa.table %>%
  mutate(MoA=sapply(MoA, function(m) str_c(unique(m), collapse='|'))) %>%
  mutate(CleanMoA=sapply(CleanMoA, function(m) str_c(unique(m), collapse='|')))

# Map MOA to compound class
class.key <- str_c('(', moa.table.clean$MoA, ')')
class <- sapply(class.key, function(k) str_detect(moa, k))

# Drop multi- and 0- target compounds 
id.drop <- rowSums(class) != 1
class.id <- moa.table.clean$Class[apply(class[!id.drop,], MAR=1, which)]
class.id <- as.factor(class.id)
nsc.set <- nsc.set[!id.drop]

# Read in bioactivity data for compound set
x.activity <- getDrugActivityData(nsc.set)

# Impute NAs as median within cell line
x.activity <- apply(x.activity, MAR=2, median_impute)
print(dim(x.activity))

# Filter out compounds with fewer than 10 instances
class.table <- data.frame(Class=class.id) %>%
  group_by(Class) %>%
  summarize(Count=n()) %>%
  arrange(desc(Count))

class.f <- filter(class.table, Count >= 10)
x.activity.f <- x.activity[class.id %in% class.f$Class,]
class.id.f <- class.id[class.id %in% class.f$Class]

#' We visualize bioactivity patterns within our compound set below. Heatmap rows 
#' correspond to different compounds and columns cell lines. Compound classes 
#' cluster relative to PC1 and PC2 but show high variability.
#+ bioactivity_heatmap, fig.height=10, fig.width=15
# Bioactivity heatmap visualization
xplot <- apply(x.activity.f, MAR=2, rank)
type <- as.factor(str_remove_all(colnames(xplot), ':.*'))
type.pal <- scales::hue_pal()(length(unique(type)))

superheat(xplot, 
          membership.rows=class.id.f,
          membership.cols=type,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=inferno(100),
          bottom.label='variable',
          bottom.label.text.col=type.pal[type],
          bottom.label.size=0.5,
          bottom.label.text.angle=90)

# Bioactivity PCA plot visualization
xpca <- prcomp(apply(x.activity.f, MAR=2, rank))
data.frame(PctVar=cumsum(xpca$sdev ^ 2) / sum(xpca$sdev ^ 2)) %>%
  mutate(NPC=1:n()) %>%
  ggplot(aes(x=NPC, y=PctVar)) +
  geom_point() +
  geom_line() +
  theme_bw()

data.frame(xpca$x, Class=class.id.f) %>%
  ggplot(aes(x=PC1, y=PC2, col=Class)) +
  geom_point(alpha=0.5) +
  theme_bw()

#' ## Predictive modeling
#' Below, we consider the problem of predicting which compound class a given 
#' compound belongs to, and how many cell lines are required for prediction.
#+ supervised_modeling, warning=FALSE
# Convery y indexing to range 0:nclass
y <- as.numeric(as.factor(as.numeric(class.id.f))) - 1
x <- x.activity.f
ncell <- 2:10
n.replicate <- 50

fits <- mclapply(1:n.replicate, function(i) {
  set.seed(i)
  
  # Generate sample split
  id.train <- createDataPartition(y, p=train.prop)$Resample1
  id.test <- setdiff(1:nrow(x), id.train)
  
  # Rebalance classes within training/test sets
  xtrain <- upSample(x[id.train,], as.factor(y[id.train]))
  xtest <- upSample(x[id.test,], as.factor(y[id.test]))
  
  # Combine balanced training/test sets
  x <- rbind(xtrain, xtest)
  y <- x$Class
  x <- select(x, -Class)
  
  id.train <- 1:nrow(xtrain)
  id.test <- setdiff(1:nrow(x), id.train)
  
  # Fit RF to predict moa from bioactivity
  out <- model(x, y, id.train, id.test, ncell)
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
ncell.plot <- 10
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
          membership.rows=class.id.f,
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
  
  # Rebalance classes within training/test sets
  xtrain <- upSample(x[id.train,], as.factor(y[id.train]))
  xtest <- upSample(x[id.test,], as.factor(y[id.test]))
  
  # Combine balanced training/test sets
  x <- rbind(xtrain, xtest)
  y <- x$Class
  x <- select(x, -Class)
  
  id.train <- 1:nrow(xtrain)
  id.test <- setdiff(1:nrow(x), id.train)
  
  # Fit lasso to predict moa from bioactivity
  yy <- sample(y)
  out <- model(x, yy, id.train, id.test, ncell)
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
  
  # Rebalance classes within training/test sets
  xtrain <- upSample(x[id.train,], as.factor(y[id.train]))
  xtest <- upSample(x[id.test,], as.factor(y[id.test]))
  
  # Combine balanced training/test sets
  x <- rbind(xtrain, xtest)
  y <- x$Class
  x <- select(x, -Class)
  
  id.train <- 1:nrow(xtrain)
  id.test <- setdiff(1:nrow(x), id.train)
  
  # Fit RF to predict moa from bioactivity
  out.null <- model(x, y, id.train, id.test, ncell, null.cells=TRUE)
  out.opt <-  model(x, y, id.train, id.test, ncell)
  
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
