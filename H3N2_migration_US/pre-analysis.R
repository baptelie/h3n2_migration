library(Biostrings); library(lubridate); library(data.table); library(dplyr)

#load the cleaned data
setwd('~/Documents/Phylogeny/USA')
meta = read.csv('data_world_12-19.csv')
(full_seq <- readDNAStringSet('nt_world_13-19_clean.fasta'))
meta <- meta[ !duplicated(meta$Isolate_Id), ]
meta <- meta[meta$Isolate_Id %in% names(full_seq),]

###keep only the USA
#identify the country
list_loc <- lapply(meta$Location, function(loc) unlist(strsplit(as.character(loc), split=" / ", fixed=T)))
names(list_loc) <- meta$Isolate_Id
list_country <- lapply(list_loc, function(loc) loc[2])

rowsUS <- grepl('States',list_country)

meta_US <- meta[rowsUS,]
seq_US <- full_seq[meta_US$Isolate_Id]

###Sort the informations on the dataframe
#Add decimal date
meta_US$Decimal_Date <- decimal_date(ymd(meta_US$Collection_Date))

#Add State
list_loc <- lapply(meta_US$Location, function(loc) unlist(strsplit(as.character(loc), split=" / ", fixed=T)))
names(list_loc) <- meta_US$Isolate_Id
meta_US$State <- lapply(list_loc, function(loc) loc[3])
list_states <- unique(meta_US$State)[c(-2,-8,-30,-53)]
meta_US <- meta_US[meta_US$State %in% list_states, ]
seq_US <- seq_US[meta_US$Isolate_Id]

#count the number of sample/state/month
nb_per_month_reg <- matrix(0,nrow=51, ncol=(2019.25-2013.5-(2/12))*12)
rownames(nb_per_month_reg)<- list_states
colnames(nb_per_month_reg)<- round(seq(from=2013.5+(2/12), to=2019.25-(1/12), by=1/12), digits=2)
for(r in list_states){
  sub <- meta_US[(meta_US$State == r),]
  col<-0
  for(m in seq(from=2013.5+(3/12), to=2019.25, by=1/12)){
    col<-col+1
    sub_month <- sub[sub$Decimal_Date <m & sub$Decimal_Date >= m-(1/12), ]
    nb_per_month_reg[r,col]<- nrow(sub_month)
  }
}
write.csv(nb_per_month_reg, file='nb_seq_per_month_state.csv')
size_states <- c(7, 40, 10, 29, 6, 4, 20, 10, 3, 13, 6, 6, 1, 3, 12, 3, 21, 5, 13, 2, 11, 1, 8, 6, 3, 10, 2, 4, 1, 2, 6, 4, 7, 1, 5, 1, 2, 7, 5, 1, 1, 9, 1, 1, 7, 4, 3, 3, 1)
names(size_states)<- list_states

subsamp_meta <- meta_US[0,]
for (r in list_states ){
  sub <- meta_US[(meta_US$State == r),]
  col<-0
  thresh <- max(1, size_states[r]/2)
  for(m in seq(from=2013.5+(4/12), to=2019.25, by=1/12)){
    col<- col+1
    sub_month <- sub[sub$Decimal_Date <m & sub$Decimal_Date >= m-(1/12), ]
    nbs <- nrow(sub_month)
    if (nbs>thresh){
      order = sample.int(nbs, thresh)
      subsamp_meta <- bind_rows( subsamp_meta, sub_month[order, ] )
    } else subsamp_meta <- bind_rows(subsamp_meta, sub_month)
  }
}

nrow(subsamp_meta)

nb_per_month_reg <- matrix(0,nrow=49, ncol=(2019.25-2013.5-(2/12))*12)
rownames(nb_per_month_reg)<- list_states
colnames(nb_per_month_reg)<- round(seq(from=2013.5+(2/12), to=2019.25-(1/12), by=1/12), digits=2)
for(r in list_states){
  sub <- subsamp_meta[(subsamp_meta$State == r),]
  col<-0
  for(m in seq(from=2013.5+(3/12), to=2019.25, by=1/12)){
    col<-col+1
    sub_month <- sub[sub$Decimal_Date <m & sub$Decimal_Date >= m-(1/12), ]
    nb_per_month_reg[r,col]<- nrow(sub_month)
  }
}
write.csv(nb_per_month_reg, file='nb_seq_per_month_state_subsamp.csv')

subsamp_seq <- seq_US[subsamp_meta$Isolate_Id]
writeXStringSet(subsamp_seq, filepath='subsamp_seq_us.fasta')
