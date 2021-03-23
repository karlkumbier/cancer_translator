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
setwd('~/github/cancer_translator/')
source('scripts/utilities.R')

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
  HDAC='HDAC',
  mTOR='MTOR',
  AKT='AKT'
)

# Read in MOA data from NCI60
nsc.set <- unique(c(getMoaToCompounds(), recursive=TRUE))
moa <- getMoaStr(nsc.set)

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

#' We visualize bioactivity patterns within our compound set below. Heatmap rows 
#' correspond to different compounds and columns cell lines. Compound classes 
#' cluster relative to PC1 and PC2 but show high variability.
#+ bioactivity_heatmap, fig.height=10, fig.width=15
# Bioactivity heatmap visualization
col.pal <- RColorBrewer::brewer.pal(length(class.key), 'Set1')
heatmap(x.activity, 
        col=inferno(100),
        RowSideColors=col.pal[as.numeric(class.id)])
legend('topleft', legend=names(class.key), fill=col.pal)

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
#+ supervised_modeling, warning=FALSE, fig.height=10, fig.width=15
y <- as.numeric(class.id) - 1
x <- x.activity
ncell <- 1:20

fits <- mclapply(1:50, function(i) {
  set.seed(i)
  
  # Generate sample split
  id.train <- createDataPartition(y, p=0.8)$Resample1
  id.test <- setdiff(1:nrow(x), id.train)
  
  # Fit lasso to predict moa from bioactivity
  fit <- cv.glmnet(x=x[id.train,], y=y[id.train], family='multinomial')
  out <- lapply(ncell, model_accuracy, fit=fit, x=x[id.test,], y=y[id.test])

  accuracy <- sapply(out, function(z) z$accuracy)
  betas <- sapply(out, function(z) z$betas)
  ypred <- data.frame(ID=as.factor(id.test), sapply(out, function(z) z$ypred))
  
  return(list(accuracy=accuracy, betas=betas, ypred=ypred))
}, mc.cores=10)


accuracy <- lapply(fits, function(f) as.matrix(f$accuracy))
betas <- lapply(fits, function(f) f$betas)
ypred <- rbindlist(lapply(fits, function(f) f$ypred)) %>%
  group_by(ID) %>%
  summarize_if(is.numeric, mean) %>%
  arrange(match(ID, 1:nrow(x)))

reshape2::melt(betas) %>%
  rename(CellLine=Var1, Ncell=Var2, Rep=L1) %>%
  mutate(Ncell=ncell[Ncell]) %>%
  group_by(CellLine) %>%
  mutate(Drop=all(value == 0)) %>%
  filter(!Drop) %>%
  ggplot(aes(x=reorder(CellLine, value, median), y=value)) +
  facet_wrap(~Ncell) +
  geom_boxplot(aes(fill=CellLine)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position='none') +
  ggtitle('Distribution of model coefficients')

reshape2::melt(accuracy) %>%
  mutate(Ncell=as.factor(ncell[Var1])) %>%
  ggplot(aes(x=Ncell, y=value)) +
  geom_boxplot(fill='#0388D1') +
  theme_bw() +
  ylab("Top-1 classification accuracy") +
  ggtitle('Model accuracy by cell lines used')

# Plot PCA distribution of compounds from select cell lines
ncell.plot <- 5
cell.select <- reshape2::melt(betas) %>%
  filter(L1 == ncell.plot) %>%
  group_by(Var1) %>%
  summarize(value=mean(value)) %>%
  top_n(ncell.plot, value)

xpca <- prcomp(apply(x[,cell.select$Var1], MAR=2, rank))
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

xplot <- apply(x[,cell.select$Var1], MAR=2, rank)
yplot <- select(ypred, one_of(str_c('X', ncell.plot)))[[1]]
superheat(xplot, 
          membership.rows=class.id,
          yr=yplot,
          yr.plot.size=1,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=inferno(100),
          bottom.label.text.angle=90)

# Compute random forest performance using select cell lines
xcell <- x[,cell.select$Var1]
rf.accuracy <- mclapply(1:1000, function(i) {
  set.seed(i)
  
  # Generate sample split
  id.train <- createDataPartition(y, p=0.8)$Resample1
  id.test <- setdiff(1:nrow(x), id.train)
  
  # Fit lasso to predict moa from bioactivity
  fit <- iRF(x=xcell[id.train,], y=as.factor(y[id.train]), type='ranger', n.iter=1)
  ypred <- predict(fit$rf.list, xcell[id.test,])
  accuracy <- mean(ypred$predictions == y[id.test])
  
  
  return(accuracy)
}, mc.cores=10)

data.frame(Accuracy=unlist(rf.accuracy)) %>%
  ggplot(aes(x=Accuracy)) +
  geom_histogram(bins=15) +
  theme_bw() +
  ggtitle('RF prediction accuracy', str_c(ncell.plot, ' cell lines'))

#' ## Null predictive modeling
#' To assess performance against random chance, we repeat the experiment above 
#' using randomly permuted response labels.
#+ random_permute, fig.height=10, fig.width=15
fits.null <- mclapply(1:50, function(i) {
  set.seed(i)
  
  # Generate sample split
  id.train <- createDataPartition(y, p=0.8)$Resample1
  id.test <- setdiff(1:nrow(x), id.train)
  
  # Fit lasso to predict moa from bioactivity
  yy <- sample(y)
  fit <- cv.glmnet(x=x[id.train,], y=yy[id.train], family='multinomial')
  out <- lapply(ncell, model_accuracy, fit=fit, x=x[id.test,], y=yy[id.test])
  
  accuracy <- sapply(out, function(z) z$accuracy)
  betas <- sapply(out, function(z) z$betas)
  
  return(list(accuracy=accuracy, betas=betas))
}, mc.cores=10)

accuracy <- lapply(fits.null, function(f) as.matrix(f$accuracy))

reshape2::melt(accuracy) %>%
  mutate(Ncell=as.factor(ncell[Var1])) %>%
  ggplot(aes(x=Ncell, y=value)) +
  geom_boxplot(fill='#0388D1') +
  theme_bw() +
  ylab("Top-1 classification accuracy") +
  ggtitle('Model accuracy by cell lines used - permuted responses')
