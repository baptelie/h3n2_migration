library(ape) ;  library(magrittr); library(lubridate); library(seqinr)
library(Biostrings); library(dplyr); library(data.table); library(DECIPHER)

prefix = "_12-19"

setwd(paste("~/Documents/Phylogeny/", prefix, sep=""))
meta = read.csv(paste("meta_12-19",'.csv', sep=""))
reg = read.csv("states_per_regBedford.csv")
(full_seq <- readDNAStringSet (paste("nt_12-19", ".fasta", sep="")))

### Remove duplicates
meta <- meta[ !duplicated(meta$Isolate_Id), ]
(full_seq <- full_seq[meta$Isolate_Id])

### Idfentify the incomplete dates
sts <- decimal_date(mdy(meta$Collection_Date))
names(sts) <- meta$Isolate_Id
hist( sts , main = 'Time of sequence sampling', breaks=200) 
incompl_date <- c( names(which(is.na(sts))), names(which(sts>2019 & sts<2012)) )
meta$Decimal_Date <- sts

###Identify the incomplete locations (country)
list_loc <- lapply(meta$Location, function(loc) unlist(strsplit(as.character(loc), split=" / ", fixed=T)))
names(list_loc) <- meta$Isolate_Id
list_country <- lapply(list_loc, function(loc) loc[2])

List_reg <- list()
for (country in list_country){
  for(r in 1:(ncol(reg)+1)){
    if(r==ncol(reg)+1) {List_reg <- c(List_reg, NA)
    } else if (country%in%reg[[r]]) {
      List_reg <- c(List_reg, names(reg[r]))
      break()
    }
  }
}

meta$Region <- List_reg
names(List_reg) <- meta$Isolate_Id
out_Reg <- names(which(is.na(List_reg)))

###Remove passaged cells
pattern=regex('original|o[tr]igial|cs|primary|PI|clinic|^(s|SIAT)\\s*[0-9]$|^SX$|^(s|SIAT)\\s*[0-9]\\s*[[:punct:]]\\s*(s|SIAT)[0-9]|direct|passage details\\:\\s*(s|SIAT)\\s*[0-9]$|SIATX\\s*[[:punct:]]\\s*(SIAT|s)[0-9]|(S|SIAT)[0-9]\\s*[[:punct:]]\\s*(SIAT|S)[1-9]|OR|SX\\/S[1-9]$|^SIAT$|SIATx\\s*[[:punct:]]\\s*SIAT[1-9]|^siat\\s*[0-9]$|SIAT 2 \\+SIAT1  |SIAT0 \\+\\s*[1-9]|SIATX*\\s*[[:punct:]]\\s*SIAT[1-9]|passage details:\\s*S[1-9]|no passage|Not passaged|CS|org|SIAT0 \\+[1-9]$|^S[1-9]\\+S[1-9]$|passage details: MDCK-SIAT1|SIAT1\\/SIAT2 \\+SIAT1|^S[0-9] \\(20|P0\\/SIAT1|^SIAT1SIAT1$|SIAT \\,SIAT[1-9]$|SIAT1-ori|MDCK-SIAT1 passage|MDCK-SIAT1 1 +SIAT1|SIAT\\-2\\, SIAT2|SIAT1\\/S1|SIAT, SIAT[1-9]|SIAT 3\\, SIAT1|passage details: MDCK-Siat/1|^P[1-9]*\\s*[[:punct:]]\\s*SIAT[1-9]$|^SIAT\\s*[1-9]\\s*[[:punct:]]\\s*[1-9]$|^S[1-9]\\s*[[:punct:]]\\s*SIAT[1-9]$|^MDCKSIAT1$|MDCK[1-9]SIAT[1-9]|^P[1-9]\\s*[[:punct:]]\\s*SIAT[1-9]$|initial|SIAT1\\s*[[:punct:]]\\s*MDCK1|^S\\s*[[:punct:]]\\s*SIAT[1-9]|SIAT2-ori|P0|NA extract|SIAT 0\\s*\\+[1-9]|pi|S2\\+2|C[0-9]\\/SIAT[0-9]|SIAT1\\,MDCK1|MDCK\\-SIAT1 2 \\+SIAT1|^C[0-9][[:punct:]]*S[0-9]$|^MDCK\\s*[0-9]\\s*[[:punct:]]\\s*SIAT[0-9]$|passage details\\: MDCK1\\,SIAT1')
Passage_history <- meta$Passage_History
keep <- grepl(pattern, x=Passage_history)

unpass_siat <- meta$Isolate_Id[keep]

duplicated.seq <- duplicated(clean_seq0)
dupl <- c()
for(seq in names(unique(clean_seq0[duplicated.seq])) ){
  Ids <- which(clean_seq0==clean_seq0[seq])
  Ids_names <- names(clean_seq0[Ids])
  dates <- clean_meta0$Decimal_Date[Ids]
  R <- list_loc[Ids_names]
  R <- sapply(R, function(r) r[min(3,length(r))])
  dupl.dates <- duplicated(dates)
  for (i in 1:length(Ids)) {
    for(j in (1:length(Ids))[-i]) if(R[[i]]==R[[j]] & dates[i]==dates[j] & duplicated.seq[i]) dupl <- c(dupl, Ids[i])
  }
}
dupl <- names(clean_seq0[dupl])
### Make a first fasta file cleaned to align in ugene
clean_meta0 <- meta[ !(meta$Isolate_Id %in% c(incompl_date, out_Reg, dupl)), ]
clean_meta0 <- clean_meta0[(clean_meta0$Isolate_Id %in%unpass_siat), ]
(clean_seq0 <- full_seq[clean_meta0$Isolate_Id])
writeXStringSet(clean_seq0, filepath=paste("nt_",prefix,"_clean0.fasta", sep=""))
cat('first cleaning saved as ',paste("nt_",prefix,"_clean0.fasta", sep=""))

#align with MAFFT on Ugene & trim extremities outside the CDS
(clean_seq <- readDNAStringSet (paste("nt",prefix,"_clean0.fasta", sep="")))
clean_seq <- RemoveGaps(clean_seq)

### Remove ambiguous sequences
#remove the sequences containing gaps
good <-c()
for(s in 1:length(clean_seq)){
  if(length(clean_seq[[s]])== 1701) good <- c(good, s)
}
clean_seq <- clean_seq[good]

#remove the sequences containing ambiguous terms
unambiguousl <- sapply(clean_seq, function(seq) grepl("R|Y|S|W|K|M|B|D|H|V|N", seq))
ambiguous <- names(clean_seq[unambiguousl])
clean_meta <- clean_meta0[clean_meta0$Isolate_Id %in% names(clean_seq),]

#remove identical sequences sampled the same day
duplicated.seq <- duplicated(clean_seq)
dupl <- c()
for(seq in names(unique(clean_seq[duplicated.seq])) ){
  Ids <- which(clean_seq==clean_seq[seq])
  Ids_names <- names(clean_seq[Ids])
  dates <- clean_meta$Decimal_Date[Ids]
  R <- list_loc[Ids_names]
  R <- sapply(R, function(r) r[min(3,length(r))])
  dupl.dates <- duplicated(dates)
  for (i in 1:length(Ids)) {
    for(j in (1:length(Ids))[-i]) if(R[[i]]==R[[j]] & dates[i]==dates[j] & duplicated.seq[i]) dupl <- c(dupl, Ids[i])
  }
}
dupl <- names(clean_seq[dupl])

(clean_seq <- clean_seq[!names(clean_seq)%in% c(dupl, ambiguous)])

### identify outliers
outliers = function(sq, n = 50, q = 0.9, score = 4){
  library(ape); library(dplyr)
  ##  For really large datasets, take instead a sample of 50 samples + each of the samples. 
  md = data.frame()
  sq_pos = 1:length(sq)  ##  all indices of each sequence.
  for(i in sq_pos){
    cat('\rWork in progress ... ',round(i/length(sq)*100),'%')
    flush.console()
    ix <- c(i, sample(sq_pos[-i], n - 1 ) )  ##  Sample of indices
    tmp = as.DNAbin( sq[ix] ) %>% dist.dna(model = 'raw') %>% as.matrix() %>% apply(1, median)
    md = rbind(md, data.frame(tip = names(tmp), median = tmp))
  }
  ## We can guarantee that each node is sampled at least 1 per 
  md2 = group_by(md, tip) %>% summarize(min = min(median),  mean = mean(median) ,max = max(median))
  mu = mean(md2$mean)
  s = sd(md2$mean)
  md2$z = with(md2, (mean - mean(mean)) / s)
  ##  Return outliers. 
  return( as.character(filter(md2, z >= score)$tip) )
}

spill <- outliers(clean_seq, score=5)
clean_seq <- clean_seq[!names(clean_seq)%in%spill]
clean_meta <- filter(clean_meta0, clean_meta0$Isolate_Id %in% names(clean_seq))

writeXStringSet(clean_seq, filepath=paste("nt_",prefix,"_clean.fasta", sep=""))
fwrite(clean_meta, paste("data_",prefix,"_clean.csv", sep="") )

clean_meta <- read.csv('data_world_13-19_clean.csv')
clean_seq <- readDNAStringSet (paste("nt_",prefix,"_clean.fasta", sep=""))

### Random subsample
nb_per_month_reg <- matrix(0,nrow=ncol(reg), ncol=88)
rownames(nb_per_month_reg)<- colnames(reg)
colnames(nb_per_month_reg)<- round(seq(from=2012, to=2019.25, by=1/12), digits=2)
for(r in colnames(reg)){
  sub <- clean_meta[(clean_meta$Region == r),]
  col<-0
  for(m in seq(from=2012+(1/12), to=2019.25+(1/12), by=1/12)){
    col<-col+1
    sub_month <- sub[sub$Decimal_Date <m & sub$Decimal_Date >= m-(1/12), ]
    nb_per_month_reg[r,col]<- nrow(sub_month)
  }
}
nb_per_month_reg
write.csv(nb_per_month_reg, 'nb_seq_per_month_reg_init.csv')

identify_seasons<-matrix(0,nrow=ncol(reg), ncol=7)
rownames(identify_seasons)<- colnames(reg)
for(r in colnames(reg)){
  for(y in 0:6){
    identify_seasons[r,y+1] <- names( which( nb_per_month_reg[r,((y*12)+1):min(((y+1)*12),ncol(nb_per_month_reg))] == min(nb_per_month_reg[r,((y*12)+1):min(((y+1)*12),ncol(nb_per_month_reg))]) ) ) [[1]]
  }
}
identify_seasons <- matrix(as.numeric(identify_seasons), ncol=7)
rownames(identify_seasons)<- colnames(reg)
identify_seasons <- identify_seasons[,2:7]

thresh <- 300
subsamp_meta <- clean_meta[0,]
for( r in colnames(reg)){
  sub <- clean_meta[(clean_meta$Region == r),]
  begin <- 2012.25-(1/12)
  for(m in c(identify_seasons[r,],2019.25)){
    sub_year <- sub[sub$Decimal_Date <m & sub$Decimal_Date >= begin, ]
    nbs <- nrow(sub_year)
    Tr <- round(thresh*(m-begin))
    if(nbs>Tr){
      order = sample.int(nbs, Tr)
      subsamp_meta <- bind_rows( subsamp_meta, sub_year[order, ] )
    } else subsamp_meta <- bind_rows(subsamp_meta, sub_year)
    begin <- m
  }
}
nrow(subsamp_meta)

nb_per_month_reg <- matrix(0,nrow=ncol(reg), ncol=88)
rownames(nb_per_month_reg)<- colnames(reg)
colnames(nb_per_month_reg)<- round(seq(from=2012, to=2019.25, by=1/12), digits=2)

for(r in colnames(reg)){
  sub <- subsamp_meta[(subsamp_meta$Region == r),]
  col<-0
  for(m in seq(from=2012+(1/12), to=2019.25+(1/12), by=1/12)){
    col<-col+1
    sub_month <- sub[sub$Decimal_Date <m & sub$Decimal_Date >= m-(1/12), ]
    nb_per_month_reg[r,col]<- nrow(sub_month)
  }
}

write.csv(nb_per_month_reg, 'nb_seq_per_month_reg_subsamp.csv')

hist(subsamp_meta$Decimal_Date, breaks=200)
Tr
(subsamp_seq <- clean_seq[subsamp_meta$Isolate_Id])

writeXStringSet(subsamp_seq, filepath=paste("nt_",prefix,"_subsamp.fasta", sep=""))
fwrite(subsamp_meta, paste("data_",prefix,"_subsamp.csv", sep="") )
