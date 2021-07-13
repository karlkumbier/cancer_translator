#+ setup, echo=FALSE
library(readxl)
library(data.table)
library(tidyverse)
library(superheat)
library(RColorBrewer)
library(ChemmineR)
library(tidytext)

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

hypergeometric_test <- function(cluster.table, full.table) {
  # Test for target enrichment across all compounds in a given cluster
  targets.test <- names(cluster.table)
  pval <- sapply(targets.test, hypergeometric_test_target,
                 cluster.table=cluster.table,
                 full.table=full.table)
  
  pval <- as.matrix(pval)
  rownames(pval) <- targets.test
  return(pval)
}

hypergeometric_test_target <- function(target, cluster.table, full.table) {
  # Test enrichment of specified target using hypergeometric test
  x <- cluster.table[target]
  m <- full.table[target]
  n <- sum(full.table) - m
  k <- sum(cluster.table)
  return(phyper(x, m, n, k))
}

filter_target <- function(x, target.keep)  {
  x <- x[names(x) %in% target.keep]
  return(x)
}

################################################################################
# Load data
################################################################################
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
fpset <- desc2fp(apset, descnames=1024, type="FPset")

dmat.euclidean <- sapply(1:length(fpset), function(i) {
  fpSim(fpset[i], fpset, method="Euclidean", sorted=FALSE)
})

################################################################################
# Validate similarity against compounds with high prevalence
################################################################################
# Filter to target classes with > 10 compounds
cpd.counts <- group_by(cpd.table, TargetsClean) %>% 
  summarize(Count=n()) %>%
  filter(Count > 10, TargetsClean != 'unknown')

# Set color palette for selected target classes
id <- cpd.table$TargetsClean %in% cpd.counts$TargetsClean
cpd <- as.numeric(as.factor(cpd.table$TargetsClean[id]))
cols <- gg_color_hue(nrow(cpd.counts))

colnames(dmat.euclidean) <- rownames(dmat.euclidean)
superheat(dmat.euclidean[id, id], 
          pretty.order.rows=TRUE, 
          pretty.order.cols=TRUE,
          bottom.label.text.angle=90,
          bottom.label.text.size=4,
          bottom.label.text.col=cols[cpd],
          left.label.text.size=4,
          left.label.text.col=cols[cpd])

################################################################################
# Clustering stability analysis
################################################################################
# Run replicates of kmeans across grid of k
n.cluster <- 2:20
n.rep <- 50

run_kmeans <- function(k) kmeans(as.matrix(fpset), centers=k)$cluster
clusters <- mclapply(n.cluster, function(k) {
  replicate(n.rep, run_kmeans(k), simplify=FALSE)
}, mc.cores=min(length(n.cluster), 12))

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
k <- n.cluster[which.min(colMeans(stab.scores))]
clusters <- kmeans(as.matrix(fpset), centers=k, nstart=10)$cluster

# Rescale distance matrix for visualization
colnames(dmat.euclidean) <- rownames(dmat.euclidean)
superheat(dmat.euclidean, 
          pretty.order.rows=TRUE, 
          pretty.order.cols=TRUE,
          membership.rows=clusters,
          membership.cols=clusters)

################################################################################
# Cluster target enrichment analysis
################################################################################
# Compute counts of targets across full compound set
cpds <- str_split(cpd.table$TargetsClean, '; ') 
cpd.table.full <- table(unlist(cpds))

# Compute counts of targets withing each cluster
cpd.table.cluster <- lapply(1:k, function(kk) table(unlist(cpds[clusters == kk])))

# Drop select targets
targets.keep <- names(cpd.table.full[cpd.table.full > 5])
targets.keep <- setdiff(targets.keep, 'unknown')

cpd.table.full <- filter_target(cpd.table.full, targets.keep)
cpd.table.cluster <- lapply(cpd.table.cluster, filter_target, target.keep=targets.keep)

# Compute hypergeomteric enrichment scores for each cluster
enrichment <- lapply(cpd.table.cluster, hypergeometric_test, full.table=cpd.table.full)

# Plot top targets for each cluster
reshape2::melt(enrichment) %>%
  rename(Target=Var1, Cluster=L1, Pval=value) %>%
  mutate(Pval=unlist(Pval)) %>%
  mutate(Logp=-log(Pval, base=10)) %>%
  mutate(Cluster=as.factor(Cluster)) %>%
  group_by(Cluster) %>%
  top_n(5, Logp) %>%
  ggplot(aes(x=reorder_within(Target, Logp, Cluster), y=Logp, fill=Target)) +
  geom_bar(stat='identity') +
  facet_wrap(~Cluster, scales='free_y') +
  theme_bw() +
  theme(legend.position='none') +
  coord_flip()
