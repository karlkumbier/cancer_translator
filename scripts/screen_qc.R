#+ setup, echo=FALSE
library(data.table)
library(tidyverse)
library(rprofiler)
library(superheat)
library(tidytext)
library(parallel)

################################################################################
# Setup
################################################################################
col.pal <- RColorBrewer::brewer.pal(11, 'RdYlBu')
col.pal[6] <- '#FFFFFF'
intensity.normalize <- TRUE
n.core <- 16

setwd('~/github/cancer_translator/')
source('scripts/utilities.R')
data.dir <- 'data/screens/LH_CDC_1/'
load(str_c(data.dir, 'ks_profiles.Rdata'))

xks <- mutate(xks, ID=str_c(PlateID, '_', WellID)) %>%
  dplyr::select(-PlateID, -WellID)

xmeta <- mutate(xmeta, ID=str_c(PlateID, '_', WellID))
x <- left_join(xks, xmeta, by='ID') %>%
  mutate(Control=Compound_ID == 'DMSO') %>%
  distinct() %>%
  filter(!is.na(PlateID)) %>%
  mutate(Col=as.numeric(sapply(str_split(WellID, '-'), tail, 1)))

# TODO: why do we need distinct?
# TODO: what are NA plates?

# Normalize intensity within each plate
if (intensity.normalize) {
  plates <- unique(x$PlateID)
  x <- lapply(plates, function(p) {
    intensity_normalize(filter(x, PlateID == p))
  })
  
  x <- rbindlist(x)
}

#+ eda, fig.height=8, fig.width=8
################################################################################
# Simple EDA - counts and intensity
################################################################################
# Cell count plots
filter(x, Control) %>%
  ggplot(aes(x=NCells)) +
  geom_histogram() +
  theme_bw() +
  ggtitle('Cell counts, DMSO wells')

filter(x, Control) %>%
  ggplot(aes(x=NCells, fill=Cell_Line)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~PlateID) +
  ggtitle('Cell counts, DMSO wells')


filter(x, Control) %>%
  ggplot(aes(x=NCells)) +
  geom_histogram() +
  theme_bw() +
  facet_wrap(~Cell_Line) +
  ggtitle('Cell counts, DMSO wells')

# Drop plate 2021018028 - strange cell count distribution
x <- filter(x, PlateID != '2021018028')

# Feature distribution plots
xcontrol <- filter(x, Control)
plate.id <- xcontrol$PlateID
cell.line <- xcontrol$Cell_Line
xcontrol <- dplyr::select(xcontrol, matches('^nonborder'))
colnames(xcontrol) <- str_remove_all(colnames(xcontrol), 'nonborder\\.\\.\\.')

superheat(xcontrol,
          membership.rows=cell.line,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=col.pal,
          heat.pal.values=seq(0, 1, by=0.1),
          bottom.label.text.angle=90,
          bottom.label.text.size=3,
          bottom.label.size=0.75)

xcontrol[xcontrol < -0.3] <- -0.3
xcontrol[xcontrol > 0.3] <- 0.3
superheat(xcontrol,
          membership.rows=cell.line,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=col.pal,
          heat.pal.values=seq(0, 1, by=0.1),
          bottom.label.text.angle=90,
          bottom.label.text.size=3,
          bottom.label.size=0.75)

# Intensity distribution plots
filter(x, Control) %>%
  dplyr::select(matches('(Intensity.*Mean|PlateID|WellID)')) %>%
  mutate(PlateID=as.factor(PlateID)) %>%
  mutate(Col=sapply(str_split(WellID, '-'), tail, 1)) %>%
  reshape2::melt() %>%
  mutate(variable=str_remove_all(variable, 'nonborder\\.\\.\\.')) %>%
  ggplot(aes(x=PlateID, y=value, fill=Col)) +
  geom_boxplot() +
  facet_wrap(~variable) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  ylab('KS')

filter(x, Control) %>%
  dplyr::select(matches('(Intensity.*Mean|PlateID|WellID|Cell_Line)')) %>%
  mutate(PlateID=as.factor(PlateID)) %>%
  mutate(Col=sapply(str_split(WellID, '-'), tail, 1)) %>%
  reshape2::melt() %>%
  mutate(variable=str_remove_all(variable, 'nonborder\\.\\.\\.')) %>%
  ggplot(aes(x=Cell_Line, y=value, fill=Col)) +
  geom_boxplot() +
  facet_wrap(~variable) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  ylab('KS')

#+ 
################################################################################
# Full dataset PCA
################################################################################
xfeat <- dplyr::select(x, matches('nonborder'))
xpca <- prcomp(xfeat)
x <- cbind(xpca$x, x)

data.frame(PctVar=cumsum(xpca$sdev ^ 2) / sum(xpca$sdev ^ 2)) %>%
  mutate(NPC=1:n()) %>%
  ggplot(aes(x=NPC, y=PctVar)) +
  geom_point() +
  geom_line() +
  theme_bw()

pc.list <- list(c('PC1', 'PC2'), c('PC3', 'PC4'), c('PC5', 'PC6'))
for (pcs in pc.list) {
  
  p <- ggplot(x, aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=log(NCells), shape=Control), alpha=0.5) +
    facet_wrap(~PlateID) +
    theme_bw() +
    scale_color_gradientn(colors=viridis::inferno(10)[1:9]) +
    ggtitle('Cell counts by plate')
  plot(p)
  
  p <- filter(x, !str_detect(Compound_Usage, '(query|reference)')) %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=Compound_Usage), alpha=0.5) +
    theme_bw() +
    ggtitle('Positive v. negative controls') +
    facet_wrap(~Cell_Line)
  plot(p)
  
  p <- filter(x, Compound_Usage == 'reference_cpd') %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=Compound_Category), alpha=0.5) +
    theme_bw() +
    ggtitle('Positive v. negative controls') +
    facet_wrap(~Cell_Line)
  plot(p)
  
  p <- filter(x, Compound_Usage == 'reference_cpd') %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=Pathway), alpha=0.5) +
    theme_bw() +
    ggtitle('Positive v. negative controls') +
    facet_wrap(~Cell_Line)
  plot(p)
  
  p <- filter(x, Compound_Usage == 'reference_cpd') %>%
    ggplot(aes_string(x=pcs[1], y=pcs[2])) +
    geom_point(aes(col=Target), alpha=0.5) +
    theme_bw() +
    ggtitle('Positive v. negative controls') +
    facet_wrap(~Cell_Line)
  plot(p)
}

################################################################################
# Bioactivity
################################################################################
# Normalize data matrix
l2norm <- function(z) z / sqrt(sum(z ^ 2))
x <- mutate_if(x, is.numeric, l2norm)

# Compute bioactivity scores for each well
xdist <- mclapply(unique(x$PlateID), function(p) {
    out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
    return(bioactivity(out))
}, mc.cores=n.core)

xdist <- rbindlist(xdist)

# Filter to compounds that have the same dose schedule for all cell lines
xdist.summary <- select(xdist, Cell_Line, Compound_ID, Dose_Category) %>%
  distinct() %>%
  group_by(Compound_ID, Cell_Line) %>%
  summarize(Dose=str_c(sort(Dose_Category), collapse='_'), .groups='drop') %>%
  group_by(Compound_ID) %>%
  mutate(Count=n()) %>%
  filter(Count == 3) %>%
  filter(Dose=='1_2_3_4_5')
  #filter(length(unique(Dose)) == 1)
  
xdist <- filter(xdist, Compound_ID %in% xdist.summary$Compound_ID)
#' # Bioactivity by compound category/target
#' To assess differential compound bioactivity by cell lines, we first compute 
#' the average distance between each `Compound_Category`/`Target` and DMSO center 
#' (separately within each cell line). **Note:** this average is taken across 
#' all dose categories. We restrict to compounds measured in all cell lines.
#+ distance_category, fig.heigh=6, fig.width=12
# Average distance by category
reference.category <- filter(xdist, Compound_Usage == 'reference_cpd')$Compound_Category
x.category <- group_by(xdist, Compound_Category, Cell_Line) %>%
  summarize(DistNorm=mean(DistNorm), .groups='drop') %>%
  group_by(Compound_Category) %>%
  mutate(Count=n()) %>%
  filter(Count == 3)

xplot <- matrix(x.category$DistNorm, nrow=3)
colnames(xplot) <- unique(x.category$Compound_Category)
rownames(xplot) <- unique(x.category$Cell_Line)

superheat(xplot,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=viridis::inferno(10),
          bottom.label.text.angle=90,
          bottom.label.size=0.5,
          title='Bioactivity by cell line, compound category')


# Average distance by target
reference.target <- filter(xdist, Compound_Usage == 'reference_cpd')$Target
x.target <- group_by(xdist, Target, Cell_Line) %>%
  summarize(DistNorm=mean(DistNorm), .groups='drop') %>%
  group_by(Target) %>%
  mutate(Count=n()) %>%
  filter(Count == 3) %>%
  filter(!is.na(Target))

xplot <- matrix(x.target$DistNorm, nrow=3)
colnames(xplot) <- unique(x.target$Target)
rownames(xplot) <- unique(x.target$Cell_Line)

superheat(xplot,
          pretty.order.rows=TRUE,
          pretty.order.cols=TRUE,
          heat.pal=viridis::inferno(10),
          bottom.label.text.angle=90,
          bottom.label.size=0.5)

#' ## Differential phenotypes
#' Below we consider examples of compopunds with targets that show differential 
#' activity across cell lines (above):
#' 1. Camptothecin (Topoisomerase) — active in A549, OVCAR4
#' 2. Gemcitabine (Nucleic Acids) — active in OVCAR, 786-0
#' 3. Pemetrexed (DHFR. TS. GARFT) — active in 786-0
#+ differential_activity, fig.height=8, fig.width=12
plot_compound(x, 'Camptothecin')
plot_compound(x, 'Gemcitabine')
plot_compound(x, 'Pemetrexed')

# TODO: classification, treat different doses independently of one another

################################################################################
# Modeling
################################################################################
# (1) Run bioactivity - drop bioinactive
# (2) Train multiclass RF on category (w/ holdout compounds) for each cell line
# (3) Compare single cell accuracy (average, optimal) v. bagged
# (4) RF distance for pathway, target

