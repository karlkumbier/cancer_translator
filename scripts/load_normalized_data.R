#+ setup, echo=FALSE, warning=FALSE, message=FALSE
library(data.table)
library(tidyverse)
library(superheat)
library(ggsci)
library(readr)
library(twosamples)

if (!exists('min.cat')) min.cat <- 5

# Initialize color palettes


################################################################################
# Load dataset
################################################################################
base.dir <- '~/github/cancer_translator/'
source(str_c(base.dir, 'scripts/utilities.R'))

# Load KS profile data, drop PCA features
data.dir <- str_c(base.dir, 'data/screens/LH_CDC_1/')
load(str_c(data.dir, 'profiles_normalized.Rdata'))
x <- dplyr::select(x, -matches('^PC'))

# Messy category table
xx <- group_by(x, Compound_ID) %>% 
  mutate(Category=Compound_Category) %>%
  summarize(Category=list(na.omit(unique(Category)))) %>% 
  mutate(NCat=sapply(Category, length)) %>%
  arrange(desc(NCat)) %>%
  mutate(Category=sapply(Category, str_c, collapse=', ')) %>%
  select(-NCat)

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

xmoa.clean <- fread(str_c(base.dir, 'data/Screened_Compound_ID_MOA_added_20200420.csv')) %>%
  dplyr::select(-InchIKeys_All) %>%
  dplyr::rename(Category=MOA) %>%
  dplyr::select(-Compound_Category, -Pathway, -Target)

x <- left_join(x, xmoa.clean, by='Compound_ID') %>%
  mutate(Category=ifelse(Compound_ID == 'DMSO', 'DMSO', Category))

# Clean compound ID names 
id2 <- match(x$Compound_ID, xmoa$Compound_ID_2)
x$Compound_ID[!is.na(id2)] <- xmoa$Compound_ID[na.omit(id2)]

id3 <- match(x$Compound_ID, xmoa$Compound_ID_3)
x$Compound_ID[!is.na(id3)] <- xmoa$Compound_ID[na.omit(id3)]
xmoa <- dplyr::select(xmoa, -Compound_ID_2, Compound_ID_3)

# Merge data with MOA table
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
xctrl <- filter(x, str_detect(Compound_Usage, 'ctrl')) %>%
  group_by(Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup()

xtreat <- filter(x, !str_detect(Compound_Usage, 'ctrl')) %>%
  filter(!Compound_ID %in% xctrl$Compound_ID) %>%
  group_by(Compound_ID) %>%
  filter(Dose == max(as.numeric(Dose))) %>%
  ungroup()

x <- rbind(xctrl, xtreat) %>%
  mutate(Cell_Line=ifelse(Cell_Line == 'ALS-WT', 'FB', Cell_Line))

################################################################################
# Clean inconsistent categories
################################################################################
xcat.clean <- data.frame(
  Compound_ID=c(
    '(-)-Menthol', 
    'Ademetionine', 
    'Berbamine', 
    'Camphor', 
    'Cisapride hydrate', 
    'Epothilone B', 
    'Harmine hydrochloride', 
    'Scopine HCl', 
    'Tetrahydropalmatine', 
    'Tilmicosin'
  ),
  Category=c(
    "opioid receptor antagonist", 
    "methyltransferase stimulant, phosphodiesterase inhibitor", 
    "Bcr-Abl kinase inhibitor", 
    "TRPV antagonist", 
    "serotonin receptor agonist", 
    "microtubule inhibitor", 
    "DYRK inhibitor", 
    "adrenergic receptor antagonist",
    "dopamine receptor antagonist",
    "bacterial cell wall synthesis inhibitor"
    )
)

id.match <- match(x$Compound_ID, xcat.clean$Compound_ID)
x$Category[!is.na(id.match)] <- xcat.clean$Category[na.omit(id.match)]

################################################################################
# Merge replicates
################################################################################
# Generate key for plate/well replicates by compound
xkey <- filter(x, Compound_ID != 'DMSO') %>%
  mutate(ID=str_c(PlateID, '_', WellID)) %>%
  group_by(Cell_Line, Compound_ID) %>%
  summarize(ID=list(ID)) %>%
  mutate(ID=sapply(ID, function(z) str_c(z, collapse='; '))) %>%
  mutate(Key=str_c(Cell_Line, Compound_ID)) %>%
  ungroup() %>%
  dplyr::select(-Cell_Line, -Compound_ID)

xmeta <- dplyr::select(x, Cell_Line, Compound_ID, Category)
x <- dplyr::select(x, matches('^nonborder'))
x <- cbind(xmeta, x)

# Initialize DMSO set
id.dmso <- xmeta$Category == 'DMSO'
xdmso <- x[id.dmso,]

# Merge replicates in non-dmso set
xtreat <- group_by(x[!id.dmso,], Cell_Line, Compound_ID, Category) %>%
  summarize_if(is.numeric, median) %>%
  ungroup()

x <- rbind(xdmso, xtreat) %>%
  mutate(Compound_Usage=ifelse(Compound_ID == 'DMSO', 'negative_ctrl_cpd', 'query_cpd')) %>%
  mutate(Compound_Usage=ifelse(Compound_ID == 'Gemcitabine', 'positive_ctrl_cpd', Compound_Usage)) %>%
  mutate(Compound_Usage=ifelse(Compound_ID == 'Bortezomib', 'positive_ctrl_cpd', Compound_Usage))
  