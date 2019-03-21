library(phytools); library(ape); library(Biostrings); library(seqinr);
library(lubridate); library(magrittr); library(dplyr); library(stringr);

pattern='15-18'

setwd('/Users/belie/Documents/Phylogeny/World_12-18/world_15-18')

meta = read.csv(paste('data_world_',pattern,'_subsamp.csv', sep=''))

sts <- decimal_date(ymd(meta$Collection_Date))
names(sts) <- meta$Isolate_Id
sts <- sts[tre$tip.label]

region <- meta$Region
names(region) <- meta$Isolate_Id

AA_Prefs <- read.csv(file= 'AA_prefs_avg.csv', header = TRUE)

C <- rep(0,81)