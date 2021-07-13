#+ setup, echo=FALSE
library(readxl)
library(tidyverse)
library(tidytext)
library(superheat)
library(RColorBrewer)
library(glmnet)
library(caret)
library(Matrix)
library(data.table)

mean_impute <- function(x) {
  x[is.na(x)] <- mean(x, na.rm=TRUE)
  return(x)
}

clean_target <- function(x) {
  # Note: issues with e.g "aurora kinases A" "B"
  x <- str_remove_all(x, '^.* of ') %>% str_split('(, | and )')
  x <- unlist(x) %>% str_remove_all('and ') %>% str_remove_all('^the ')
  return(x)
}

get_beta <- function(betas, n.cell) {
  # Wrapper to get beta across bootstrap replicates
  return(sapply(betas, get_beta_, n.cell=n.cell))
}

get_beta_ <- function(betas, n.cell) {
  # Get betas for selected # of cell lines
  s <- colSums(betas != 0) - 1
  id.select <- max(which(s <= n.cell))
  beta <- betas[,id.select][-1] != 0
  return(beta)
}

setwd('~/github/cancer_translator/')
xcell <- read_excel('data/41589_2016_BFnchembio1986_MOESM242_ESM.xlsx')
xcpd <- read_excel('data/41589_2016_BFnchembio1986_MOESM243_ESM.xlsx')
xbio <- read_excel('data/41589_2016_BFnchembio1986_MOESM244_ESM.xlsx')
bioactive.thresh <- 10 # AUC threshold for calling bioactivity

xtarget.clean <- fread('compound_key_LH.csv', skip=1) %>% 
  rename(ID=V1, Target=V2) %>%
  filter(Target != '') %>%
  arrange(ID)

# Filter to cell lines contained in both datasets
ccl <- intersect(xcell$master_ccl_id, xbio$master_ccl_id)
xcell <- filter(xcell, master_ccl_id %in% ccl)
xbio <- filter(xbio, master_ccl_id %in% ccl)

# Generate data matrix for <cell line> x <compound>
n <- length(unique(xbio$master_ccl_id))
p <- length(unique(xbio$master_cpd_id))
x <- matrix(NA, nrow=n, ncol=p)

rownames(x) <- unique(xbio$master_ccl_id)
colnames(x) <- unique(xbio$master_cpd_id)

for (i in 1:nrow(xbio)) {
  ccl <- as.character(xbio$master_ccl_id[i])
  cpd <- as.character(xbio$master_cpd_id[i])
  x[ccl, cpd] <- xbio$area_under_curve[i]
}

# Reorder data matrix to match cell metadata
x <- x[match(xcell$master_ccl_id, rownames(x)),]

#' # Exploratory data analysis
#' Below, we conduct a preliminary assessment of cancer cell line bioactivity 
#' data available through the cancer therapeutics response portal (CTRP). 
#' Specifically, we consider the availabe cell lines (with metadata), compounds 
#' (with metadata), and bioactivity data available for optimal cell line 
#' selection.
#' 
#' ## Cell lines
#' The CTRP dataset contain information on 842 cancer cell lines. Each cell line
#' is associated with features indicating (1) name (2) identifier (3) 
#' site/histology (4) growth conditions.
#' 
#' **Note:** We filter out all suspension cell lines per Louise's suggestion.
#+ cell_data, fig.width=12, fig.height=7
dim(xcell)
head(xcell)

group_by(xcell, ccle_primary_site) %>%
  summarize(Count=n()) %>%
  ggplot(aes(x=reorder(ccle_primary_site, Count), y=Count)) +
  geom_bar(stat='identity') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Cell line counts by site')

group_by(xcell, ccle_primary_hist) %>%
  summarize(Count=n()) %>%
  ggplot(aes(x=reorder(ccle_primary_hist, Count), y=Count)) +
  geom_bar(stat='identity') +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Cell line counts by histology')

# Filter out suspension cells
id.drop <- xcell$growth_mode == 'suspension'
xcell <- xcell[!id.drop,]
x <- x[!id.drop,]

#' ## Compounds
#' The CTRP dataset measures bioactivity for 481 compounds. In addition to 
#' identifiers for the compound/source, the dataset contains information on 
#' protein target (gene symbol), target, and smiles representation. Many of the 
#' compounds in the dataset target multiple proteins (see figure below). 
#' In total, the CTRP compounds target 322 unique compounds (+1 NA).
#+ compound_data, warning=FALSE, fig.height=7, fig.width=12
dim(xcpd)
head(xcpd)

# Generate binary matrix of cpd x target
compound.targets <- str_split(xcpd$gene_symbol_of_protein_target, ';')
unique.targets <- unique(unlist(compound.targets))

xtarget <- matrix(0, nrow=length(compound.targets), ncol=length(unique.targets))
rownames(xtarget) <- xcpd$cpd_name
colnames(xtarget) <- unique.targets
for (i in 1:length(compound.targets)) {
  j <- which(unique.targets %in% compound.targets[[i]])
  xtarget[i, j] <- 1
}

# List out targets
targets <- str_split(xcpd$target_or_activity_of_compound, '; ')
data.frame(Target=unlist(targets)) %>%
  group_by(Target) %>%
  summarize(Count=n()) %>% 
  filter(Count > 1) %>%
  ggplot(aes(x=reorder(Target, Count), y=Count)) +
  geom_bar(stat='identity') +
  theme_bw() +
  coord_flip() +
  ggtitle('Target counts - full target')

# Clean targets
drop <- c('diversity oriented synthesis', 'screening hit', 'natural product')
targets.cl <- lapply(targets, clean_target)

# Merge target IDs with Louise's labels
cpd.table <- data.frame(Targets=sapply(targets, str_c, collapse='; ')) %>%
  mutate(TargetsClean=sapply(targets.cl, str_c, collapse='; ')) %>%
  mutate(Compound=xcpd$cpd_name) %>%
  mutate(ID=1:n()) %>%
  left_join(xtarget.clean, by='ID') %>%
  mutate(TargetsClean=ifelse(is.na(Target), TargetsClean, Target)) %>%
  select(-Target)

group_by(cpd.table, TargetsClean) %>%
  summarize(Count=n(), .groups='drop') %>% 
  filter(Count > 1) %>%
  filter(!TargetsClean %in% drop) %>%
  ggplot(aes(x=reorder(TargetsClean, Count), y=Count)) +
  geom_bar(stat='identity') +
  theme_bw() +
  coord_flip() +
  ggtitle('Target counts - consolidated target')

# Drop unique target and unknown target compounds
id.unknown <- cpd.table$TargetsClean == 'unknown'
target.table <- table(cpd.table$TargetsClean)
id.single <- cpd.table$TargetsClean %in% names(target.table[target.table == 1])

x <- x[,!(id.unknown | id.single)]
cpd.table <- cpd.table[!(id.unknown | id.single),]

#' ## Bioactivity
#' The bioactivity data measures area under the curve (AUC) — evaluated at 
#' multiple concentrations — for cell line/compound pairs. We convert these 
#' data into a <n cell lines> x <n compounds> matrix.
#' 
#' On average, 20% of all cell line/compound pairs are missing AUC data. However,
#' missingness varies considerably by compound. For instance, a subset 
#' of ~100 compounds are missing data for > 50% of cell lines. We handle 
#' missingness by imputing AUC as mean within each compound.
#+ bioactivity_data, fig.width=12, fig.height=7
dim(x)
mean(is.na(x))

# Plot missingness proportion by compound
data.frame(PropNA=colMeans(is.na(x))) %>%
  ggplot(aes(x=PropNA)) +
  geom_histogram() +
  theme_bw() +
  ggtitle('Distribution of missingness by compound')

# Drop compounds with high missingness
id.drop <- colMeans(is.na(x)) > 0.3
x <- x[,!id.drop]
cpd.table <- cpd.table[!id.drop,]

# Plot missingness proportion by cell line
data.frame(PropNA=rowMeans(is.na(x))) %>%
  ggplot(aes(x=PropNA)) +
  geom_histogram() +
  theme_bw() +
  ggtitle('Distribution of missingness by cell line')

# Drop cell lines with high missingness
id.drop <- rowMeans(is.na(x)) > 0.3
x <- x[!id.drop,]
xcell <- xcell[!id.drop,]

# Impute missing AUC as mean within compound
x <- apply(x, MAR=2, mean_impute)
reshape2::melt(x) %>%
  ggplot(aes(x=value)) +
  geom_histogram() +
  theme_bw() +
  xlab('Bioactivity (AUC)') +
  ggtitle('Distribution of compound bioactivity')

#' We first consider the maximum bioactivity (minimum area under the dose 
#' response curve) for each compound in the dataset.
#+ compound_bioactivity, fig.height=7, fig.width=7
data.frame(AUC=apply(x, MAR=2, min)) %>%
  ggplot(aes(x=AUC)) +
  geom_histogram() +
  theme_bw() +
  ggtitle('Distribution of minimum AUC across compounds')

data.frame(AUC=colSums(x <= bioactive.thresh)) %>%
  group_by(AUC) %>%
  summarize(Count=n()) %>%
  ggplot(aes(x=AUC, y=Count)) +
  geom_point() +
  theme_bw() +
  ylab('# Compounds') +
  xlab('# Cell lines with AUC < 10 by compound')

# Filter out compounds w/ AUC < 10 across all cell lines
id.drop <- colSums(x <= bioactive.thresh) == 0
x <- x[,!id.drop]
cpd.table <- cpd.table[!id.drop,]

group_by(cpd.table, TargetsClean) %>%
  summarize(Count=n(), .groups='drop') %>% 
  filter(Count > 1) %>%
  filter(!TargetsClean %in% drop) %>%
  ggplot(aes(x=reorder(TargetsClean, Count), y=Count)) +
  geom_bar(stat='identity') +
  theme_bw() +
  coord_flip() +
  ggtitle('Target counts - consolidated target, filtered')

#' Next, we plot raw bioactivity data — area under the dose response curve 
#' (AUC) — for each cell line/compound pair.
#' 
#' Below, we plot bioactivity for each cell line/compound pair. The majority of 
#' compounds do 
#' not exhibit differential bioactivity by cell — i.e. nearly all cell lines
#' exhibit bioactivity at or near the maximum bioactivity across cell lines.
#' To focus in on compounds with differential bioactivity, we subset based
#' on the standard deviation of compound bioactivity across cell lines.

#+ bioactivity_plot_1, fig.height=12, fig.width=18
# Threshold bioactivity and scale to 0-1
xplot <- x
xplot[xplot > 15] <- 15

# Compute average bioactivity within each compound
xmean <- colMeans(xplot)

cols <- c(rev(brewer.pal(9, 'Blues')), brewer.pal(9, 'Reds'))
superheat(xplot[,order(xmean)],
          pretty.order.rows=TRUE,
          membership.cols=cpd.table$TargetsClean[order(xmean)],
          yt=xmean[order(xmean)],
          yt.plot.size=0.5,
          yt.axis.name='Average AUC',
          heat.pal=cols,
          heat.pal.values=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
          bottom.label.text.angle=90,
          bottom.label.size=1,
          bottom.label.text.size=3, 
          title='Compound bioactivity by cell line')

superheat(xplot[,head(order(xmean), 50)],
          pretty.order.rows=TRUE,
          membership.cols=cpd.table$TargetsClean[head(order(xmean), 50)],
          yt=xmean[head(order(xmean), 50)],
          yt.plot.size=0.5,
          yt.axis.name='Average AUC',
          heat.pal=cols,
          heat.pal.values=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
          bottom.label.text.angle=90,
          bottom.label.size=1,
          title='Compound bioactivity by cell line')

#' # Modeling
#' To account for the long tail of target counts, we restrict our initial 
#' analysis to predicting the most prevalent target classes in the CTRP dataset.
#' Toward that end, we filter out all compounds with fewer than 5 targets, 
#' leaving us with a total of only 6 compounds.
#+ modeling, warning=FALSE
################################################################################
# Fit model
################################################################################
# Set grid for # of cell lines to be used
n.cells <- c(5, 10, 15, 20, 25, 50)

# Filter to prevalent compounds
cpd.counts <- group_by(cpd.table, TargetsClean) %>%
  summarize(Count=n(), .groups='drop') %>%
  filter(Count >= 5)

id.select <- cpd.table$TargetsClean %in% cpd.counts$TargetsClean
y <- as.factor(cpd.table$TargetsClean[id.select])
target.levels <- levels(y)
x <- t(x)[id.select,]

# Generate key for selected compounds/id
cpd.key <- rename(cpd.table, master_cpd_id=ID)

# Create cell key for plot
cell.key <- select(xcell, ccl_name, master_ccl_id)

# Generate bootstrap samples for stability analysis
bootstrap <- createResample(y, times=200)

# Fit models across bootstrap replicates
fits <- lapply(bootstrap, function(i) {
  if (any(table(y[i]) < 4)) return(NULL)
  fit <- glmnet(x=x[i,], y=as.numeric(y)[i], family='multinomial')
  
  betas <- lapply(coef(fit), abs)
  betas <- Reduce('+', betas)
  s <- colSums(betas != 0) - 1
  
  id.test <- setdiff(1:nrow(x), i)
  ypred <- predict(fit, x[id.test,], type='response')
  ypred.n <- lapply(n.cells, function(n) {
    out <- ypred[,,max(which(s <= n))]
    rownames(out) <- id.test
    return(out)
  })
  
  return(list(fit=fit, ypred=ypred.n))
})

# Filter replicates with small class size in bootstrap sample
id.drop <- sapply(fits, is.null)
ypreds <- lapply(fits[!id.drop], function(f) f$ypred)
fits <- lapply(fits[!id.drop], function(f) f$fit)

# Reformat predictions as data table and aggregate across replciates
ypreds <- reshape2::melt(ypreds) %>%
  rename(ID=Var1, PredictedClass=Var2) %>%
  mutate(NCells=n.cells[L2]) %>%
  mutate(PredictedClass=target.levels[PredictedClass]) %>%
  group_by(ID, PredictedClass, NCells) %>%
  summarize(Probability=mean(value)) %>%
  mutate(TrueClass=y[ID])

################################################################################
# Posthoc analysis on selected cell lines
################################################################################
# Get selected cell lines across replicates
betas <- lapply(fits, coef)
betas <- lapply(betas, function(b) Reduce('+', lapply(b, abs)))

# Compute proportion of times each cell line is selected across sample splits
betas.select <- lapply(n.cells, get_beta, betas=betas)
xplot <- reshape2::melt(betas.select) %>%
  mutate(NCells=n.cells[L1]) %>%
  rename(master_ccl_id=Var1, Selected=value) %>%
  left_join(cell.key, by='master_ccl_id') %>%
  group_by(NCells, master_ccl_id, ccl_name) %>%
  summarize(Proportion=mean(Selected), .groups='drop') %>%
  group_by(NCells) %>%
  top_n(10, Proportion)

# Raw selection proportion
ggplot(xplot, aes(x=reorder_within(ccl_name, Proportion, NCells), y=Proportion)) +
  geom_bar(stat='identity', aes(fill=ccl_name)) +
  facet_wrap(~NCells, scales='free_y') +
  coord_flip() +
  theme_bw() +
  theme(legend.position='none') +
  ylab('Proportion times selected across resampling replicates')

# Selection proportion relative to expected
mutate(xplot, Expected=choose(ncol(x), NCells - 1) / choose(ncol(x), NCells)) %>%
  mutate(Enrichment=log(Proportion, base=2) - log(Expected, base=2)) %>%
  ggplot(aes(x=reorder_within(ccl_name, Enrichment, NCells), y=Enrichment)) +
  geom_bar(stat='identity', aes(fill=ccl_name)) +
  facet_wrap(~NCells, scales='free_y') +
  coord_flip() +
  theme_bw() +
  theme(legend.position='none') +
  ylab('Selection enrichment across resampling replicates')

# Plot heatmap of top selected cell lines
x.select <- filter(xplot, NCells == 5) %>% top_n(5, Proportion)
id.select <- colnames(x) %in% x.select$master_ccl_id

xplot <- x[,id.select]
xplot[xplot > 15] <- 15
id.match <- match(colnames(xplot), x.select$master_ccl_id)
colnames(xplot) <- x.select$ccl_name[id.match]

superheat(xplot,
          membership.rows=y,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=cols,
          heat.pal.values=c(0, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
          bottom.label.text.angle=90,
          bottom.label.size=1,
          bottom.label.text.size=4, 
          title='Compound bioactivity by selected cell lines')
################################################################################
# Posthoc analysis on model accuracy
################################################################################
# Summarize model predictions using different #s of cell lines
accuracy <- c()
ncells <- c()
for (i in n.cells) {

  ypred.i <- filter(ypreds, NCells == i)
  xplot <- matrix(ypred.i$Probability, ncol=max(ypred.i$ID))
  ypred.top <- apply(xplot, MAR=2, which.max)
  accuracy <- c(accuracy, mean(ypred.top == as.numeric(y)))
  ncells <- c(ncells, i)
  
  xplot <- apply(xplot, MAR=2, function(z) z / max(z))
  superheat(t(xplot), 
            membership.rows=y,
            membership.cols=target.levels,
            pretty.order.rows=TRUE,
            pretty.order.cols=TRUE,
            title=str_c('Predictions: # cell lines = ', i),
            heat.pal=cols,
            heat.pal.values=seq(0, 1, by=0.1))
}

data.frame(Ncells=as.factor(ncells), Accuracy=accuracy) %>%
  ggplot(aes(x=Ncells, y=Accuracy)) +
  geom_bar(stat='identity') +
  theme_bw() +
  ggtitle('Top 1 prediction accuracy by # of cell lines used') +
  theme(text=element_text(size=14))

