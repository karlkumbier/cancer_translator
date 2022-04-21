#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(superheat)
library(ggsci)
library(readr)
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

if (is.null(min.cat)) min.cat <- 5

# Initialize color palettes
heat.pal <- viridis::viridis(10)
col.pal <- pal_nejm()(8)

################################################################################
# Load dataset
################################################################################
base.dir <- '~/github/cancer_translator/'
source(str_c(base.dir, 'scripts/utilities.R'))

# Load KS profile data, drop PCA features
data.dir <- str_c(base.dir, 'data/screens/LH_CDC_1/')
load(str_c(data.dir, 'profiles_normalized.Rdata'))
x <- dplyr::select(x, -matches('^PC'))

# Initialize # of cell lines
n.cell.line <- length(unique(x$Cell_Line))

# Clean compound categories, taking maximum vote across replicates
x <- group_by(x, Compound_ID) %>%
  mutate(Dose=as.numeric(Dose)) %>%
  mutate(Compound_Category=select_category(Compound_Category)) %>%
  ungroup()

# Load cleaned MOA table
xmoa <- fread(str_c(base.dir, 'data/Unique_Screened_Compound_MOA.csv')) %>%
  dplyr::rename(Compound_ID=Compound_ID_1) %>%
  dplyr::select(-InchIKeys_All) %>%
  dplyr::rename(Category=MOA)

# Clean compound ID names 
id2 <- match(x$Compound_ID, xmoa$Compound_ID_2)
x$Compound_ID[!is.na(id2)] <- xmoa$Compound_ID[na.omit(id2)]

id3 <- match(x$Compound_ID, xmoa$Compound_ID_3)
x$Compound_ID[!is.na(id3)] <- xmoa$Compound_ID[na.omit(id3)]
xmoa <- dplyr::select(xmoa, -Compound_ID_2, Compound_ID_3)

# Merge data with MOA table
x <- left_join(x, xmoa, by='Compound_ID') %>%
  dplyr::select(-Compound_Category) %>%
  mutate(Category=ifelse(Compound_ID == 'DMSO', 'DMSO', Category))

# Filter to compound/dose combinations evaluated in all cell lines
x.treat <- select(x, Cell_Line, Compound_ID, Dose_Category, Compound_Usage) %>%
  distinct() %>%
  group_by(Compound_ID, Dose_Category, Compound_Usage) %>% 
  count() %>%
  filter(n == n.cell.line) %>%
  mutate(ID=str_c(Compound_ID, Dose_Category, Compound_Usage))

x <- mutate(x, ID=str_c(Compound_ID, Dose_Category, Compound_Usage)) %>%
  filter(ID %in% x.treat$ID) %>%
  filter(!is.na(Compound_Usage))

# Initialize plate x cell line key
plate.key <- dplyr::select(x, PlateID, Cell_Line) %>% distinct()

#' # Category summary
#' Figures below summarize compound counts by category across the full compound 
#' library as well as categories with the largest number of compounds 
#' represented.
#+ count_summaries, fig.height=8, fig.width=12
xcat.tab <- filter(x, !is.na(Category)) %>%
  dplyr::select(Compound_ID, Category) %>%
  distinct() %>%
  group_by(Category) %>%
  summarize(Count=n())

filter(xcat.tab, Count > 1) %>%
  filter(Category != '') %>%
  ggplot(aes(x=Count)) +
  geom_histogram() +
  ggtitle('Distribution of compound counts per category') +
  xlab('# compounds') +
  theme(text=element_text(size=24))

# Summarize prevalent compound categories
xcat.tab.top <- filter(xcat.tab, Count >= min.cat, Category != '')
cat.keep <- xcat.tab.top$Category

xcat.tab.top %>%
  ggplot(aes(x=reorder(Category, Count), y=Count)) +
  geom_bar(stat='identity', position=position_dodge(preserve="single")) +
  geom_text(aes(label=Count), nudge_y=0.025) +
  theme(axis.text.x=element_text(angle=90)) +
  ggtitle('Compounds per category counts', str_c('top categories')) +
  scale_y_log10()

# Initialize category table
category.table <- dplyr::select(x, Compound_ID, Category, Dose, Cell_Line) %>%
  dplyr::rename(MOA=Category) %>%
  group_by(Compound_ID, MOA, Dose) %>%
  summarize(N_Cell_Line=length(unique(Cell_Line)), .groups='drop') %>%
  mutate(MOA=str_replace_all(MOA, ',', ';')) %>%
  arrange(Compound_ID)

################################################################################
# Filter to highest dose
################################################################################
xctrl <- filter(x, str_detect(Compound_Usage, 'ctrl')) 

xtreat <- filter(x, !str_detect(Compound_Usage, 'ctrl')) %>%
  group_by(Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup()

x <- rbind(xctrl, xtreat) %>%
  mutate(Cell_Line=ifelse(Cell_Line == 'ALS-WT', 'FB', Cell_Line))
