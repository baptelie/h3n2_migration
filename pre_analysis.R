library(ape) ;  library(magrittr); library(lubridate); library(seqinr)
library(Biostrings); library(dplyr); library(data.table)

prefix = "world_13-19"

setwd(paste("~/Documents/Phylogeny/", prefix, sep=""))
meta = read.csv(paste("data_",prefix,".csv", sep=""))
reg = read.csv("/Users/belie/Documents/Phylogeny/states_per_regBedford.csv")
(full_aln <- readDNAStringSet (paste("nt_", prefix, ".fasta", sep="")))

### Remove duplicates
meta <- meta[ !duplicated(meta$Isolate_Id), ]
(full_aln <- full_aln[meta$Isolate_Id])

### Remove ambiguous sequences
unambiguousl <- sapply(full_aln, function(seq) grepl("R|Y|S|W|K|M|B|D|H|V|N", seq))
ambiguous <- names(full_aln[unambiguousl])

### Idfentify the incomplete dates
sts <- decimal_date(ymd(meta$Collection_Date))
names(sts) <- meta$Isolate_Id
hist( sts , main = 'Time of sequence sampling', breaks=100) 
incompl_date <- c( names(which(is.na(sts))), names(which(sts>2019.5)) )
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
keep = list('primary','clinic','SIAT','Original','Clinical','direct','P0','S[1,2,3,4,5,6,7,8,9]','cs')
Passage_history <- meta$Passage_History
names(Passage_history) <- meta$Isolate_Id
unpass_siat <- unique(unlist(lapply(keep, function(x) names(Passage_history[grepl(x,Passage_history, useBytes=TRUE)]) )))

### Make a first fasta file cleaned to align in ugene
clean_meta0 <- meta[ !(meta$Isolate_Id %in% c(incompl_date, out_Reg, ambiguous)), ]
clean_meta0 <- clean_meta0[(clean_meta0$Isolate_Id %in%unpass_siat), ]
(clean_seq0 <- full_aln[clean_meta0$Isolate_Id])

L = lapply((1:length(clean_seq0)), function(i) as.character(clean_seq0[[i]]) )

write.fasta(L, names = names(clean_seq0), file.out=paste("nt_",prefix,"_clean0.fasta", sep=""))
#align with MAFFT on Ugene & trim extremities outside the CDS
tmp.aln <- readDNAStringSet (paste("nt_",prefix,"_clean0.fasta", sep=""))

### identify outliers
outliers = function(sq, n = 50, q = 0.9, score = 4){
  library(ape); library(dplyr)
  ##  For really large datasets, take instead a sample of 50 samples + each of the samples. 
  md = data.frame()
  sq_pos = 1:length(sq)  ##  all indices of each sequence.
  for(i in sq_pos){
    cat('\rWork in progress ... ',round(i/length(sq)*100),'%', sep='')
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

spill <- outliers(tmp.aln, n=100, score=3.5)

clean_meta <- clean_meta0[ !(clean_meta0$Isolate_Id %in% c(spill)), ]
(clean_seq <- tmp.aln[clean_meta$Isolate_Id])

L = lapply((1:length(clean_seq)), function(i) as.character(clean_seq[[i]]) )

write.fasta(L, names = names(clean_seq), file.out=paste("nt_",prefix,"_clean.fasta", sep=""))
fwrite(clean_meta, paste("data_",prefix,"_clean.csv", sep="") )


### Random subsample max 20 seq/month/reg
nb_per_month_reg <- matrix(0,nrow=ncol(reg), ncol=(2019.25-2013.5-(3/12))*12)
rownames(nb_per_month_reg)<- colnames(reg)
colnames(nb_per_month_reg)<- round(seq(from=2013.5+(3/12), to=2019.25-(1/12), by=1/12), digits=2)
for(r in colnames(reg)){
  sub <- clean_meta0[(clean_meta0$Region == r),]
  col<-0
  for(m in seq(from=2013.5+(4/12), to=2019.25, by=1/12)){
    col<-col+1
    sub_month <- sub[sub$Decimal_Date <m & sub$Decimal_Date >= m-(1/12), ]
    nb_per_month_reg[r,col]<- nrow(sub_month)
  }
}
nb_per_month_reg
fwrite(nb_per_month_reg, 'nb_seq_per_month_reg.csv')

subsamp_meta <- meta[0,]
for (r in colnames(reg) ){
  sub <- clean_meta[(clean_meta$Region == r),]
  col<-0
  for(m in seq(from=2013.5+(4/12), to=2019.25, by=1/12)){
    col<- col+1
    thresh <- max(c(round(quantile(unlist(nb_per_month_reg[, max(c(0,col-5)):min(c(ncol(nb_per_month_reg),col+5))]), probs=seq(from=0, to=1, by=1/9))[7]), 30))
    sub_month <- sub[sub$Decimal_Date <m & sub$Decimal_Date >= m-(1/12), ]
    nbs <- nrow(sub_month)
    if (nbs>thresh){
      order = sample.int(nbs, thresh)
      subsamp_meta <- bind_rows( subsamp_meta, sub_month[order, ] )
    } else subsamp_meta <- bind_rows(subsamp_meta, sub_month)
  }
}
 
hist(subsamp_meta$Decimal_Date, breaks=60)

(subsamp_seq <- clean_seq[subsamp_meta$Isolate_Id])

L = lapply((1:length(subsamp_seq)), function(i) as.character(subsamp_seq[[i]]) )
write.fasta(L, names = names(subsamp_seq), file.out=paste("nt_",prefix,"_subsamp.fasta", sep=""))
fwrite(subsamp_meta, paste("data_",prefix,"_subsamp.csv", sep="") )