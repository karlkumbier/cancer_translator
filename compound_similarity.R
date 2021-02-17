#+ setup, echo=FALSE
library(readxl)
library(data.table)
library(tidyverse)
library(superheat)
library(RColorBrewer)
library(ChemmineR)

clean_target <- function(x) {
  # Note: issues with e.g "aurora kinases A" "B"
  x <- str_remove_all(x, '^.* of ') %>% str_split('(, | and )')
  x <- unlist(x) %>% str_remove_all('and ') %>% str_remove_all('^the ')
  return(x)
}

gg_color_hue <- function(n) {
  # Generate ggplot color palette
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

setwd('~/github/cancer_translator/')
xcell <- read_excel('data/41589_2016_BFnchembio1986_MOESM242_ESM.xlsx')
xcpd <- read_excel('data/41589_2016_BFnchembio1986_MOESM243_ESM.xlsx')
xbio <- read_excel('data/41589_2016_BFnchembio1986_MOESM244_ESM.xlsx')

xtarget.clean <- fread('compound_key_LH.csv', skip=1) %>% 
  rename(ID=V1, Target=V2) %>%
  filter(Target != '') %>%
  arrange(ID)

# Filter to cell lines contained in both datasets
ccl <- intersect(xcell$master_ccl_id, xbio$master_ccl_id)
xcell <- filter(xcell, master_ccl_id %in% ccl)
xbio <- filter(xbio, master_ccl_id %in% ccl)

# List out and clean targets
drop <- c('diversity oriented synthesis', 'screening hit', 'natural product')
targets <- str_split(xcpd$target_or_activity_of_compound, '; ')
targets.cl <- lapply(targets, clean_target)

# Merge target IDs with Louise's labels
cpd.table <- data.frame(Targets=sapply(targets, str_c, collapse='; ')) %>%
  mutate(TargetsClean=sapply(targets.cl, str_c, collapse='; ')) %>%
  mutate(Compound=xcpd$cpd_name) %>%
  mutate(ID=1:n()) %>%
  left_join(xtarget.clean, by='ID') %>%
  mutate(TargetsClean=ifelse(is.na(Target), TargetsClean, Target)) %>%
  select(-Target)

################################################################################
# Compute compound similarity
################################################################################
# Convert cpd smiles to atom pair
smiles <- xcpd$cpd_smiles
sdfset <- smiles2sdf(smiles)
cid(sdfset) <- xcpd$cpd_name
apset <- sdf2ap(sdfset)

# Compute fingerprint
fpset <- desc2fp(apset, descnames=128, type="FPset")

dmat.tanimoto <- sapply(1:length(fpset), function(i) {
  fpSim(fpset[i], fpset, method="Tanimoto", sorted=FALSE)
})

dmat.euclidean <- sapply(1:length(fpset), function(i) {
  fpSim(fpset[i], fpset, method="Euclidean", sorted=FALSE)
})

dmat.dice <- sapply(1:length(fpset), function(i) {
  fpSim(fpset[i], fpset, method="Dice", sorted=FALSE)
})

################################################################################
# Validate similarity against compounds with high prevalence
################################################################################
# Filter to target classes with > 5 compounds
cpd.counts <- group_by(cpd.table, TargetsClean) %>% 
  summarize(Count=n()) %>%
  filter(Count > 5, TargetsClean != 'unknown')

# Set color palette for selected target classes
id <- cpd.table$TargetsClean %in% cpd.counts$TargetsClean
cpd <- as.numeric(as.factor(cpd.table$TargetsClean[id]))
cols <- gg_color_hue(nrow(cpd.counts))

colnames(dmat.tanimoto) <- rownames(dmat.tanimoto)
superheat(dmat.tanimoto,#[id, id], 
          pretty.order.rows=TRUE, 
          pretty.order.cols=TRUE,
          bottom.label.text.angle=90,
          bottom.label.text.size=4,
          bottom.label.text.col=cols,#[cpd],
          left.label.text.size=4,
          left.label.text.col=cols)#[cpd])

colnames(dmat.euclidean) <- rownames(dmat.euclidean)
superheat(dmat.euclidean,#[id, id], 
          pretty.order.rows=TRUE, 
          pretty.order.cols=TRUE,
          bottom.label.text.angle=90,
          bottom.label.text.size=4,
          bottom.label.text.col=cols,#[cpd],
          left.label.text.size=4,
          left.label.text.col=cols)#[cpd])

colnames(dmat.dice) <- rownames(dmat.dice)
superheat(dmat.dice,#[id, id], 
          pretty.order.rows=TRUE, 
          pretty.order.cols=TRUE,
          bottom.label.text.angle=90,
          bottom.label.text.size=4,
          bottom.label.text.col=cols,#[cpd],
          left.label.text.size=4,
          left.label.text.col=cols)#[cpd])

################################################################################
# Clustering stability analysis
################################################################################
stability_score <- function(id1, id2) {
  # Compute stability between two clusterings as defined in Wu et al. 2016
  clusters1 <- unique(id1)
  clusters2 <- unique(id2) 
  stopifnot(length(clusters1) == length(clusters2))
  
  x1 <- sapply(clusters1, function(i) as.numeric(id1 == i))
  x2 <- sapply(clusters2, function(i) as.numeric(id2 == i))
  
  x.intersect <- t(x1) %*% x2
  x.union <- nrow(x1) - t(1 - x1) %*% (1 - x2)
  dmat <- x.intersect / x.union

  k <- length(clusters1)
  maxr <- sum(apply(dmat, MAR=1, max))
  maxc <- sum(apply(dmat, MAR=2, max))
  return(1 / (2 * k) * (2 * k - maxr - maxc))
}

stab_k <- function(x) {
  # Wrapper to compute stability scores for select k
  n.rep <- length(x)
  rep.pairs <- expand.grid(1:n.rep, 1:n.rep)
  
  out <- sapply(1:nrow(rep.pairs), function(i) {
    stability_score(x[[rep.pairs$Var1[i]]], x[[rep.pairs$Var2[i]]])
  })
  
  return(out)
}


# Run replicates of kmeans across grid of k
n.cluster <- 2:10
n.rep <- 25

run_kmeans <- function(k) kmeans(as.matrix(fpset), centers=k)$cluster
clusters <- lapply(n.cluster, function(k) {
  replicate(n.rep, run_kmeans(k), simplify=FALSE)
})

stab.scores <- sapply(clusters, stab_k)
reshape2::melt(stab.scores) %>%
  rename(K=Var2, Instability=value) %>%
  mutate(K=n.cluster[K]) %>%
  group_by(K) %>%
  summarize(SD=sd(Instability), Instability=mean(Instability)) %>%
  ggplot(aes(x=K, y=Instability)) +
  geom_line() +
  geom_point() +
  theme_bw()

# Generate clusters from optimal k
clusters <- kmeans(as.matrix(fpset), centers=6, nstart=10)$cluster

# Rescale distance matrix for visualization
colnames(dmat.euclidean) <- rownames(dmat.euclidean)
dmat.euclidean <- (dmat.euclidean - min(dmat.euclidean))
dmat.euclidean <- dmat.euclidean / max(dmat.euclidean)

superheat(dmat.euclidean, 
          pretty.order.rows=TRUE, 
          pretty.order.cols=TRUE,
          membership.rows=clusters,
          membership.cols=clusters)

# TODO: check how compound classes map to stable clusters
