#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(tidytext)
library(parallel)
library(superheat)
library(ggsci)
library(patchwork)
library(hrbrthemes)
library(dbscan)
library(umap)
library(rprofiler)
library(twosamples)

theme_set(
  theme_ipsum(
    axis_title_size=18, 
    strip_text_size=18, 
    axis_text_size=14,
    base_size=18, 
    base_family='sans'
  )
)

################################################################################
# Initialize constant parameters to be used in analysis
################################################################################
# Initialize analysis parameters
intensity.normalize <- TRUE
n.core <- 16
min.cat <- 5
pval.thresh <- 0
ncell.thresh <- 50
save.fig <- TRUE

# Initialize clustering parameters
#cell.line <- 'OVCAR4'

################################################################################
# Load dataset
################################################################################
setwd('~/github/cancer_translator/')
source('scripts/utilities.R')

# Load KS profile data
data.dir <- 'data/screens/LH_CDC_1/'
load(str_c(data.dir, 'profiles_qc_norm=', intensity.normalize, '.Rdata'))
x <- dplyr::select(x, -matches('^PC'))

# Initialize # of cell lines
n.cell.line <- length(unique(x$Cell_Line))

# initialize output directories
output.dir <- '~/github/cancer_translator/results/distance/'
fig.dir <- '~/github/cancer_translator/results/figures/fig2/'
dir.create(output.dir, showWarnings=FALSE)

# Clean compound categories
x <- group_by(x, Compound_ID) %>%
  mutate(Compound_Category=select_category(Compound_Category)) %>%
  ungroup()

# Filter to compound/dose combinations evaluated in all cell lines
x.treat <- dplyr::select(x, Cell_Line, Compound_ID, Dose_Category, Compound_Usage) %>%
  distinct() %>%
  group_by(Compound_ID, Dose_Category, Compound_Usage) %>% 
  count() %>%
  filter(n == n.cell.line) %>%
  mutate(ID=str_c(Compound_ID, Dose_Category, Compound_Usage))

x <- mutate(x, ID=str_c(Compound_ID, Dose_Category, Compound_Usage)) %>%
  filter(ID %in% x.treat$ID) %>%
  filter(!is.na(Compound_Usage)) %>%
  dplyr::rename(Category_Vendor=Compound_Category)

# Load compound target table
x.target <- 'data/Bioactivity_broad_label/' %>%
  str_c('bioactivity_broad_institute_label_20220108.csv') %>%
  fread() %>%
  dplyr::select(Compound_ID, MOA_broad_institute, Target_broad_institute) %>%
  distinct() %>%
  dplyr::rename(Category=MOA_broad_institute) 

# Merge Broad MOA labels
x <- left_join(x, x.target, by='Compound_ID') %>%
  mutate(Category=ifelse(Category == '', Category_Vendor, Category)) %>%
  mutate(Category=ifelse(Compound_ID == 'DMSO', 'DMSO', Category))
  

# Replace missing categories with vendor labels
xcat.tab <- filter(x, !is.na(Category)) %>%
  filter(Category != '') %>%
  dplyr::select(Compound_ID, Category) %>%
  distinct() %>%
  group_by(Category) %>%
  summarize(Count=n())

# Initialize reference set compounds
x.reference <- filter(x, Compound_Usage == 'reference_cpd')
x.reference.group <- group_by(x.reference, Category) %>%
  summarize(Compound_ID=list(Compound_ID))

#' # Overview
#' This notebook considers the problem of compound similarity based on 
#' `r cell.line` phenotypic profiles. Compounds are first filtered based on 
#' bioactivity — as specified in the `bioactivity` notebook. At a high level, 
#' bioactivity is based on the distance (evaluated over phenotypic profiles) 
#' between a target well and a collection of DMSO control wells. 
#+ bioactivity, fig.height=8, fig.width=12, echo=FALSE
# Mean-center + correlation weight features by cell line for subsequent analysis
x <- correlation_weight(x, TRUE)

################################################################################
# Initialize parameters for bioactivity analysis
################################################################################
# Initialize DMSO distance summary function
null_summary <- function(x) {
  iqr.thresh <- IQR(x)
  return(max(x[x < median(x) + iqr.thresh]))
}

################################################################################
# Compute bioactivity scores
################################################################################
xdist <- lapply(unique(x$PlateID), function(p) {
  out <- filter(x, PlateID == p) %>% dplyr::select(-matches('^PC'))
  return(bioactivity(out, null_summary))
})

# Average bioactivity scores within treatment/cell line groups
xgroup <- rbindlist(xdist) %>%
  group_by(Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup() %>%
  group_by(Cell_Line, Compound_ID, Category) %>% 
  summarize(DistNorm=mean(DistNorm), NCells=mean(NCells), .groups='drop') %>%
  arrange(Compound_ID, Cell_Line) 

mutate(xgroup, Bioactive=DistNorm > 1) %>%
  ggplot(aes(x=NCells, y=DistNorm)) +
  geom_point(aes(col=Bioactive), alpha=0.6) +
  facet_wrap(~Cell_Line) +
  scale_x_log10() +
  scale_color_d3() +
  geom_vline(xintercept=ncell.thresh, lty=2, col='grey') +
  ylab('Normalized distance from DMSO')

#' # Similarity analysis
#' The following case study consider the problem of clustering compounds based 
#' on phenotypic profiles. Prior to clustering we filter compounds as follows:
#' 
#' 1. Categories with at least `r min.cat` compounds 
#' 2. Compounds that are called as bioactive (distance > 1)
#' 3. Compounds that do not result in cell death phenotype (Ncells >= `r ncell.thresh`)
#+ clustering_analysis
################################################################################
# Initialize cluster analysis compounds
################################################################################
# Initialize bioactive compound set
compounds.select <- filter(xgroup, DistNorm > 1 | Compound_ID == 'DMSO') %>%
  filter(NCells >= ncell.thresh)

# Get compound counts for bioactive set and filter based on minimum category size
compound.table <- select(compounds.select, Compound_ID, Category) %>%
  distinct() %>%
  group_by(Category) %>%
  mutate(Count=n()) %>%
  filter(Count >= min.cat | Compound_ID == 'DMSO') %>%
  filter(!is.na(Category))

x <- filter(x, Compound_ID %in% compound.table$Compound_ID) %>%
  group_by(Compound_ID, Cell_Line) %>% 
  mutate(Dose=as.numeric(Dose)) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup()

#' ## Case study for `r cell.line`
#' To assess the similarity of compounds with the same MOA in phenotypic space, 
#' we consider both within and between category similarity. Within category 
#' similarity compares the distributions of: (i) pairwise distances between 
#' samples with the same MOA and (ii) pairwise distances between DMSO samples. 
#' Intuitively, this asks whether points with the same MOA are closer in 
#' phenotypic space than DMSO controls. Between category similarity compares the 
#' distributions of: (i) pairwise distances between samples with the same MOA
#' (ii) pairwise distances between samples with the given MOA and their nearest 
#' neighbors. Intuitively, this asks whether points with the same MOA are closer
#' to one another than other compounds.
#' 
#' We compare distributions using a signed KS statistic. For within category
#' similarity, values < 0 imply that the DMSO samples are closer than than 
#' samples with a given MOA. For between category similarity, values > 0 imply 
#' that samples with the same MOA are closer to one another than to their 
#' nearest neighbors.
#+ analysis, fig.height=8, fig.width=16
# Initialize feature set for select cell line
xc <- filter(x, Cell_Line == cell.line)
xc.feat <- select(xc, matches('^non'))
colnames(xc.feat) <- str_remove_all(colnames(xc.feat), '^non\\.\\.\\.')
category.table <- c(table(xc$Category))

summary_fun <- function(x, y) cvm_stat(x, y, power=1)
dist.category <- category_distance(xc.feat, xc$Category, summary_fun) %>%
  mutate(Score=DistBetween + DistWithin) %>%
  arrange(desc(Score)) %>%
  mutate(Count=category.table[Category])

#' Below we visualize phenotypic profiles projeted into umap space. We report 
#' the top 5 most distinguishable categories based on profile similarity.
#+ umap
################################################################################
# UMAP visualization
################################################################################
# Initialize umap embedding for visualization
configs <- umap.defaults
configs$n_neighbors <- min.cat
xumap <- umap(xc.feat)

# Initialize data for visualization
xplot <- data.frame(X1=xumap$layout[,1], X2=xumap$layout[,2]) %>%
  mutate(Category=xc$Category)

# Plot select compound category in UMAP space, top neighborhood
within.select <- arrange(dist.category, desc(DistWithin))$Category[1:2]
between.select <- arrange(dist.category, desc(DistBetween))$Category[1:2]
category.select <- c(within.select, between.select, 'DMSO')

xplot.select <- filter(xplot, Category %in% category.select)
xplot.other <- filter(xplot, !Category %in% category.select) %>%
  mutate(Category='Other')

p <- rbind(xplot.select, xplot.other) %>%
  mutate(Size=(Category != 'Other') + 1) %>%
  mutate(Size=Size + (!Category %in% c('DMSO', 'Other'))) %>%
  ggplot(aes(x=X1, y=X2, col=Category, size=Size, alpha=Size)) +
  geom_jitter(height=0.025, width=0.025) +
  xlab('UMAP1') +
  ylab('UMAP2') +
  scale_color_jama() +
  scale_size(guide='none', range=c(1.5, 2.5)) +
  scale_alpha(guide='none', range=c(0.25, 1)) +
  labs(col=NULL) +
  theme(legend.position=c(0.8, 0.85))

if (save.fig) pdf(str_c(fig.dir, cell.line, '.pdf'), height=12, width=12)
plot(p)
if(save.fig) dev.off()

dist.category <- mutate(dist.category, Cell_Line=cell.line)
save(file=str_c(output.dir, cell.line, '.Rdata'), dist.category)


