library(twosamples)
library(MASS)
library(tidyverse)
library(hrbrthemes)
library(ggsci)
library(data.table)

theme_set(theme_ipsum(
  plot_title_size=18, 
  subtitle_size=14,
  strip_text_size=16,
  axis_title_size=18 
))

################################################################################
# Initialize functions for phenoactivity/phenosimilarity
################################################################################
phenoactivity <- function(x, h, id.null=1, features='^X') {
  # Computes distance for each well to DMSO centroid
  
  # Select data
  xdata <- select(x, matches(features))
  
  # Initialize control point cloud centroid
  xctl <- apply(filter(xdata, x$Label == id.null), MAR=2, median)
  
  # Compute distance to centroid
  xdist <- sqrt(colMeans((t(xdata) - xctl) ^ 2))
  
  # Generate null distribution of maximum distance DMSO to DMSO center
  scores <- sapply(unique(x$Label), function(l) {
    h(xdist[x$Label == id.null], xdist[x$Label == l])
  })
  
  out <- mutate(x, Dist=xdist)
  return(list(scores=scores, x=out))
}

phenosimilarity <- function(x, h, features='^X') {
  
  # Initialize category set
  categories <- x$Label
  x <- select(x, matches(features))
  
  # Compute pairwise distances
  dist.mat <- dist(x) %>% as.matrix
  
  dist.summary <- lapply(unique(categories), function(ctg) {
    
    # Compute within category distance distribution
    id.ctg <- which(categories == ctg)
    ctg.dist <- dist.mat[id.ctg, id.ctg]
    ctg.dist <- ctg.dist[upper.tri(ctg.dist)]
    
    # Compute category neighborhood distance distribution
    nbhd.dist <- apply(dist.mat[id.ctg, ], MAR=1, function(z) sort(z))
    nbhd.dist <- nbhd.dist[2:length(unique(id.ctg)),]
    

    score <- 1 + h(ctg.dist, nbhd.dist)
    
    out <- rbind(
      data.frame(Dist=c(ctg.dist), Group='Query', Label=ctg),
      data.frame(Dist=c(nbhd.dist), Group='Reference', Label=ctg)
    )
    
    return(list(score=score, dist=out))
  })
  
  scores <- sapply(dist.summary, function(z) z$score)
  out <- rbindlist(lapply(dist.summary, function(z) z$dist))

  return(list(scores=scores, dist=out))
}

#' # Case studies
#' The following case studies examine the performance of our phenoactivity and
#' phenosimilarity under a variety of settings. We initialize the summary
#' function for comparing ECDFs below.
#+ summary_function

# Initialize summary function
#h <- function(x, y) cvm_stat(x, y, 1)
h <- function(x, y) wass_stat(x, y)
#h <- function(x, y) ks_stat(x, y)


#' ### Increasing effect size
#' Our first case study examines the impact of increasing effect size between
#' reference and query populations. We take reference and control populations 
#' to be multivariate Gaussian with increasing difference in means and identity
#' covariance.

# Initialize grid of mean vectors
mus <- list(
  c(0, 0),
  c(0.5, 0.5),
  c(1, 1), 
  c(2, 2),
  c(5, 5),
  c(10, 10)
)

# Initialize sample/reference sizes
n <- c(100, rep(25, 5))

# Initialize covariance
sigmas <- lapply(mus, function(z) diag(length(z)))

# Compute effect size
dist.mu <- as.matrix(dist(do.call(rbind, mus)))[,1]

# Initialize multivariate normal data
x <- mapply(function(n, m, s) mvrnorm(n, m, s), n, mus, sigmas, SIMPLIFY=FALSE)

x <- do.call(rbind, x)
x <- data.frame(x, Label=as.factor(rep(1:length(n), n)))

# Compute phenoactivity and phenosimilarity
xpa <- phenoactivity(x, h)
xps <- phenosimilarity(x, h)

# Summary plots
ggplot(xpa$x, aes(x=X1, y=X2, col=Label)) +
  geom_point() +
  scale_color_jama(name='') +
  theme(legend.position=c(0.9, 0.15))

group_by(xpa$x, Label) %>%
  mutate(Score=round(xpa$scores[as.numeric(Label)], 2)) %>%
  mutate(Label=ifelse(Label == 1, 'Reference', str_c('Query, ', Label))) %>%
  mutate(Label=str_c(Label, ', ', Score)) %>%
  arrange(Dist) %>%
  mutate(ECDF=(0:(n() - 1)) / n()) %>%
  ggplot(aes(x=Dist, y=ECDF)) +
  geom_line(aes(color=Label), lwd=1) +
  scale_color_jama(name='') +
  xlab('Distance from DMSO centroid') +
  ylab(NULL) +
  theme(legend.position=c(0.9, 0.15))

group_by(xpa$x, Label) %>%
  mutate(Score=round(xpa$scores[as.numeric(Label)], 2)) %>%
  mutate(Label=ifelse(Label == 1, 'Reference', str_c('Query, ', Label))) %>%
  mutate(Label=str_c(Label, ', ', Score)) %>%
  ggplot(aes(x=Dist)) +
  geom_density(aes(color=Label, fill=Label), alpha=0.5) +
  scale_color_jama(name='') +
  scale_fill_jama(name='') +
  xlab('Distance from DMSO centroid') +
  ylab(NULL) +
  theme(legend.position=c(0.9, 0.15))

group_by(xps$dist, Label, Group) %>%
  mutate(Label=str_c(Label, ', ', round(xps$scores[as.numeric(Label)], 2))) %>%
  arrange(Dist) %>%
  mutate(ECDF=(0:(n() - 1)) / n()) %>%
  ggplot(aes(x=Dist, y=ECDF)) +
  geom_line(aes(color=Group), lwd=1) +
  scale_color_jama(name='') +
  xlab('Pairwise distance') +
  ylab(NULL) +
  facet_wrap(~Label) +
  theme(legend.position=c(0.9, 0.1))

data.frame(MuDist=dist.mu, Score=xpa$scores) %>%
  ggplot(aes(x=MuDist, y=Score)) +
  geom_line(lwd=1) +
  xlab('Effect size') +
  ylab('Phenoactivity')

#' ### Heterogeneity
#' Our second case study examines the impact label heterogeneity. We set the 
#' reference population as a mulivariate gaussian and the query population as
#' a mixture of multivariate gaussians. The mixture includes one component near 
#' control and one far from control. Figures below compare results as the 
#' proportion of samples near control varies.

# Initialize grid of mean vectors
mus <- list(
  c(0, 0),
  c(0, 0),
  c(-3, -3)
)

for (i in c(5, 10, 20, 30, 40, 45)) {
  
  # Initialize sample/reference sizes
  n <- c(100, i, 50-i)
  
  # Initialize covariance
  sigmas <- lapply(mus, function(z) diag(length(z)))
  
  # Compute effect size
  dist.mu <- as.matrix(dist(do.call(rbind, mus)))[,1]
  
  # Initialize multivariate normal data
  x <- mapply(function(n, m, s) mvrnorm(n, m, s), n, mus, sigmas, SIMPLIFY=FALSE)
  
  x <- do.call(rbind, x)
  x <- data.frame(x, Label=as.factor(rep(c(1, 2, 2), n)))
  
  # Compute phenoactivity and phenosimilarity
  xpa <- phenoactivity(x, h)
  xps <- phenosimilarity(x, h)
  
  # Summary plots
  p1 <- ggplot(xpa$x, aes(x=X1, y=X2, col=Label)) +
    geom_point() +
    scale_color_jama(name='') +
    theme(legend.position=c(0.9, 0.15)) +
    ggtitle('Heterogeneity', str_c('p = ', i / 50))
  
  p2 <- group_by(xpa$x, Label) %>%
    mutate(Score=round(xpa$scores[as.numeric(Label)], 2)) %>%
    mutate(Label=ifelse(Label == 1, 'Reference', str_c('Query, ', Label))) %>%
    mutate(Label=str_c(Label, ', ', Score)) %>%
    arrange(Dist) %>%
    mutate(ECDF=(0:(n() - 1)) / n()) %>%
    ggplot(aes(x=Dist, y=ECDF)) +
    geom_line(aes(color=Label), lwd=1) +
    scale_color_jama(name='') +
    xlab('Distance from DMSO centroid') +
    ylab(NULL) +
    ggtitle('Phenoactivity') +
    theme(legend.position=c(0.9, 0.2))
  
  p3 <- group_by(xps$dist, Label, Group) %>%
    filter(Label != 1) %>%
    mutate(Label=str_c(Label, ', ', round(xps$scores[as.numeric(Label)], 2))) %>%
    arrange(Dist) %>%
    mutate(ECDF=(0:(n() - 1)) / n()) %>%
    ggplot(aes(x=Dist, y=ECDF)) +
    geom_line(aes(color=Group), lwd=1) +
    scale_color_jama(name='') +
    xlab('Pairwise distance') +
    ylab(NULL) +
    ggtitle('Phenosimilarity') +
    facet_wrap(~Label) +
    theme(legend.position=c(0.85, 0.15))
  
  gridExtra::grid.arrange(p1, p2, p3, ncol=3)
}
