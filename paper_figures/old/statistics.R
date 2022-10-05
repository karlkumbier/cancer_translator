library(tidyverse)
library(hrbrthemes)
library(ggsci)

n <- 1000
col.pal <- pal_jama()(7)[c(1, 4, 5)]

# Initialize distributions
d1 <- rnorm(n)
d2 <- c(rexp(n / 2) - 1, rnorm(n/2, mean=2))
d <- c(d1, d2)

e <- c(rep(1 / n, n), rep(0, n))
f <- c(rep(0, n), rep(1 / n, n))
dorder <- order(d)

# Initialize data for visualization
xplot <- data.frame(X=d[dorder]) %>%
  mutate(Y1=cumsum(e[dorder]), Y2=cumsum(f[dorder])) %>%
  mutate(Delta=Y1 - Y2)

id.max <- which.max(xplot$Delta)
xks <- xplot$X[id.max]
ymin.ks <- xplot$Y2[id.max]
ymax.ks <- xplot$Y1[id.max]
delta.ks <- round(xplot$Delta[id.max], 2)

ggplot(xplot, aes(x=X))  +
  geom_line(size=1, col=col.pal[1], aes(y=Y1)) +
  geom_line(size=1, col=col.pal[3], aes(y=Y2)) +
  geom_segment(x=xks, xend=xks, y=ymin.ks, yend=ymax.ks, col=col.pal[2]) +
  geom_text(x=xks + 0.2, y=ymin.ks + delta.ks / 2, label=delta.ks, col=col.pal[2], size=5) +
  ylab('P(X < x)') +
  xlab('x')
  

delta.cvm <- round(mean(xplot$Delta) * 2, 2)
ggplot(xplot, aes(x=X))  +
  geom_line(size=1, col=col.pal[1], aes(y=Y1)) +
  geom_line(size=1, col=col.pal[3], aes(y=Y2)) +
  geom_ribbon(aes(ymin=pmin(Y1, Y2), ymax=pmax(Y1, Y2)), alpha=0.25) +
  geom_text(x=xks + 0.2, y=ymin.ks + delta.ks / 2, label=delta.cvm, col=col.pal[2], size=5) +
  ylab('P(X < x)') +
  xlab('x')
  