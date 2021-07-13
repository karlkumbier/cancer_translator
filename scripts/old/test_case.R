library(readxl)
setwd('~/github/cancer_translator/')
n.core <- 16

mean_impute <- function(x) {
  x[is.na(x)] <- mean(x, na.rm=TRUE)
  return(x)
}

greedy_search <- function(x, 
                          objective=objective_1,
                          depth=5,
                          nsample=sqrt(nrow(x)), 
                          psample=sqrt(ncol(x))) {

  s <- numeric(depth)
  for (i in 1:depth) {
    s[i] <- greedy_select(x, s, objective, nsample, psample)
  }
  
  return(data.table(s=list(s), val=objective(x[s,])))
}
 
greedy_select <- function(x, s, objective,
                          nsample=sqrt(nrow(x)), 
                          psample=sqrt(ncol(x))) {
  
  n.sample <- sample(nrow(x), nsample)
  p.sample <- sample(ncol(x), psample)
  out <- sapply(n.sample, function(i) objective(x[c(i, s), p.sample, drop=FALSE]))
  
  return(sample(n.sample[out == max(out)], 1))
}

objective_1 <- function(x) {
  active <- colMeans(x == 0) != 1
  return(sum(active))
}

xcell <- read_excel('data/41589_2016_BFnchembio1986_MOESM242_ESM.xlsx')
xcpd <- read_excel('data/41589_2016_BFnchembio1986_MOESM243_ESM.xlsx')
xbio <- read_excel('data/41589_2016_BFnchembio1986_MOESM244_ESM.xlsx')

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

# Check missingness in bioactivity data
n.na <- rowMeans(is.na(x))
hist(n.na, main='Proporiton missing by cell line')

p.na <- colMeans(is.na(x))
hist(p.na, main='Proportion missing by compound')

# Drop compounds missing in > 20% of samples
x <- x[,p.na < 0.2]
x <- apply(x, MAR=2, mean_impute)

x[x < 12] <- 0
x[x >= 12] <- 1
heatmap(x)

depth <- 3:20
out <- mclapply(depth, function(d) {
  out <- replicate(100, greedy_search(x, depth=d), simplify=FALSE)
  out <- rbindlist(out) %>% mutate(Depth=d)
  return(out)
}, mc.cores=n.core)

out <- rbindlist(out)
ggplot(out, aes(x=as.factor(Depth), y=val)) +
  geom_boxplot()
