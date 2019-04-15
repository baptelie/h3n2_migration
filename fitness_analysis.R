library(phytools); library(ape); library(Biostrings); library(seqinr);
library(lubridate); library(magrittr); library(dplyr); library(stringr);

pattern='12-18'

setwd('/Users/belie/Documents/Phylogeny/World_12-18/ML_analysis')

#Run treetime and load the tree
system(paste('/anaconda3/bin/treetime --aln nt_world_', pattern, '_subsamp.fasta --tree ', "nt_world_12-18_subsamp.fasta.treefile", ' --dates dates.csv --clock-rate 0.005 --outdir treetime2', sep=''), wait=TRUE)
tre.tt <- ape::read.nexus("./treetime2/timetree.nexus")
# tre.tt <- drop.tip(tre.tt, c('EPI_ISL_195549','EPI_ISL_197247'))

#Load and process metadata
meta = read.csv(paste('data_world_',pattern,'_subsamp.csv', sep=''))

sts <- decimal_date(ymd(meta$Collection_Date))
names(sts) <- meta$Isolate_Id
sts <- sts[tre.tt$tip.label]

region <- meta$Region
names(region) <- meta$Isolate_Id

meta_nodes <- data.frame(Isolate_Id = tre.tt$node.label)
tip.order <- c(1:length(tre.tt$tip.label))
names(tip.order) <- meta$Isolate_Id
tip.order <- tip.order[tre.tt$tip.label]
meta_tree <- bind_rows(meta[tip.order, ], meta_nodes) #build a new metadata folder with all computed informations from all nodes, in the order of the absolute node numbers
meta_tree$Decimal_Date <- timeNodes(tre.tt)
row.names(meta_tree) <- meta_tree$Isolate_Id
meta_tree$Isolate_Id <- NULL

#compute the fitness on each node
Alignment <- readDNAStringSet (paste("./treetime2", "/ancestral_sequences.fasta", sep=""))
AA_Prefs <- read.csv(file= '/Users/belie/Documents/Phylogeny/New\ folder/AA_prefs_avg.csv', header = TRUE)
Alignment <- Alignment[c(tre.tt$tip.label, tre.tt$node.label)]
Alignment_AA <- Biostrings::translate(Alignment)
Alignment_AA <- Alignment_AA[c(tre.tt$tip.label, tre.tt$node.label)] #reorder in the order of the absolute node number
list_mut_AA = apply(tre.tt$edge, 1, function(edge) mismatches(Alignment_AA[ edge[2] ], Alignment_AA[ edge[1] ]))
list_mut_nt = apply(tre.tt$edge, 1, function(edge) mismatches(Alignment[ edge[2] ], Alignment[ edge[1] ]))
nb_mut_ns = sapply(list_mut_AA, function(e) nrow(e))
nb_mut_s = sapply(list_mut_nt, function(e) nrow(e)) - nb_mut_ns
ka_ks = apply(tre.tt$edge, 1, function(edge) kaks(DNAbin2Aln(Alignment[c(edge[2], edge[1]) ]) ))
kadks = sapply(ka_ks, function(k) k$ka/k$ks)
list_ME = lapply(list_mut_AA, function(ListMut) Sum_ME(ListMut))
Fitness <- FitnessNodes(tre.tt) #compute the fitness of each node/tip, with 0 for the root
meta_tree$Fitness <- Fitness


#identify the homogeneous antigenic clades
meta_tree <- define_HAC(meta_tree)

#Run the MCMC stochastic mapping
dst<- tre.tt
dst$edge.length[dst$edge.length==0]<-max(nodeHeights(tre.tt))*1e-6 # set zero-length branches to be 1/1000000 total tree length
nsim=10

start=Sys.time()
trees.sm <- make.simmap(dst, region, model='ARD', nsim=nsim, Q='mcmc', burnin=10000, samplefreq=200, pi='estimated')
end=Sys.time()
print(end-start)

write.simmap(trees.sm, 'trees.sm')

#Analyze each stochastic mapping
C <- rep(0,64) #
D <- rep(0,64)
nreg=length(levels(region))
slope <- array(0,dim=c(nreg, nreg, 3), dimnames=list(levels(region), levels(region), c('inf','equal','sup')))

setwd("./fitness_evol")
for (t in 1:nsim){
  # pdf(paste('tree_',t,'.pdf',sep=''))
  cnt <- rep(0,64)
  tre.sm <- trees.sm[[t]]
  #identify the migration events
  listClades <- list.clades(tre.sm)
  membersClades <- getMembers(listClades, tre.sm)
  dt.edges <- dates.edges()
  PM<-pairsMigr(listClades)
  print(sum(PM))
  
  #analysis of the migrant fitness probability
  LPM <- listPairsMigr(listClades)
  x <- sapply(levels(region), function(d) sapply(levels(region), function(r) proba.obs(d,r,LPM)))
  cnt[!is.na(x)] <- 1
  x[is.na(x)] <- 0
  sup <- (x>0.1)*rep(1, length(x))
  C <- C+sup
  D <- D+cnt
  
  #plot the fitness evolution after migration
  # for(don in levels(region)) for(rec in levels(region)){
  #   if(don==rec) next()
  #   FE <- fitness_evol(don, rec, LPM)
  #   if(!is.na(FE[1])) {
  #     for(i in 1:length(FE$slope)){
  #       if (FE$slope[i]<0 & FE$CI95[2*i]<0) {slope[don,rec,'inf'] <- slope[don,rec,'inf'] + 1
  #       } else if (FE$slope[i]>0 & FE$CI95[2*i-1]>0) {slope[don,rec,'sup'] <- slope[don,rec,'sup'] + 1
  #       } else slope[don,rec,'equal'] <- slope[don,rec,'equal'] + 1
  #     }
  #     slope[don,rec,] <- slope[don,rec,] / length(FE$slope)
  #   }
  # }
  # dev.off()
}
C/D
LM <- matrix(C, nrow=length(unique(region)), byrow=TRUE)
rownames(LM) <- levels(region)
colnames(LM) <- levels(region)
LM
for(don in levels(region)) for(rec in levels(region)){
  slope[don,rec,] <- slope[don,rec,]/sum(slope[don,rec,])
}
round(slope, digits=2)

#compare fitness distri after migration vs sister within a same region
# fitness_migr <- c()
# fitness_within <- c()
# GS <-c(getStates(tre.sm, type='tips'), getStates(tre.sm, type='nodes')) #region of each tip/node in the order : tips and then nodes
# pdf('fitness_migr_within.pdf')
# for (n in 1:nrow(listClades)){
#   me <- listClades[n,]
#   FP <- Fitness[getParent(tre.tt,me$Node)]
#   nodes_rec <- membersClades[[nodenumb.to.lab(me$Node)]]
#   nodes_don <- nodes.sameReg(me$Parent_Reg, getParent(tre.tt, me$Node), GS, tre.sm)
#   if(length(nodes_don)<5) next()
#   f_rec <- Fitness[nodes_rec]-FP
#   f_don <- Fitness[nodes_don]-FP
#   plot.ecdf(f_don, col='red', xlim=c(-15,10), main=paste(me$Node, me$Parent_Reg, me$Reciep_Reg), verticals=TRUE)
#   par(new=TRUE)
#   plot.ecdf(f_rec, col='green', xlim=c(-15,10), main='', verticals=TRUE)
# }
# dev.off()
# 
# plot.ecdf(fitness_migr, xlim=c(-15,10), col='green')
# par(new=TRUE)
# plot.ecdf(fitness_within, xlim=c(-15,10), col='red')
# 
# hist(fitness_migr, breaks=100, col='green')
# hist(fitness_within, breaks=100, col='red', add=TRUE)

pdf('fitness_migr_vs_within_dates_inf05.pdf')
for (t in 1:nsim){
  print(t)
  # cnt <- rep(0,64)
  tre.sm <- trees.sm[[t]]
  #identify the migration events
  listClades <- list.clades(tre.sm)
  membersClades <- getMembers(listClades, tre.sm)

  fitness_migr <- c()
  dates_migr <- c()
  fitness_within <- c()
  dates_within <- c()
  count <- 0
  GS <-c(getStates(tre.sm, type='tips'), getStates(tre.sm, type='nodes')) #region of each tip/node in the order : tips and then nodes
  for (n in listClades$Node){
    # print(n)
    FP <- Fitness[getParent(tre.tt,n)]
    reg=GS[n]
    nodes_clade <- membersClades[[nodenumb.to.lab(n)]]
    fam <- gd_children(n)
    nodes_migr <- strtoi(check_within(unlist(fam), reg, GS))
    dates_m <- meta_tree[nodes_migr,]$Decimal_Date - meta_tree[getParent(tre.tt,n),]$Decimal_Date
    nodes_migr <- nodes_migr[dates_m<0.5]
    dates_migr <- c(dates_migr, dates_m[dates_m<0.5])
    fitness_migr <- c(fitness_migr, Fitness[nodes_migr]-FP)

    dates <- meta_tree[nodes_clade,]$Decimal_Date
    if(length(nodes_clade)>25 & max(dates)-min(dates)>0.5){
      for(gc in fam$gd_ch){
        children <- tre.tt$edge[which(tre.tt$edge[,1]==gc),2]
        for(c in children){
          f <- unlist(gd_children(node=c))
          if(length(f)<4) next()
          if(length(check_within(f, reg, GS))==length(f)){
            fitness_within <- c(fitness_within, Fitness[strtoi(f)]-Fitness[gc]);  #if all nodes in the same region
            # if(any(meta_tree[strtoi(f),]$Decimal_Date - meta_tree[gc,]$Decimal_Date<0)) print(c(f,gc))
            dates_within <- c(dates_within, meta_tree[strtoi(f),]$Decimal_Date - meta_tree[gc,]$Decimal_Date)
          } 
        }
      }
    }
  }
  count
  plot.ecdf(fitness_migr, xlim=c(-10,8), col='green')
  par(new=TRUE)
  plot.ecdf(fitness_within, xlim=c(-10,8), col='red')
  plot.ecdf(dates_migr, col='green', xlim=c(0,1))
  par(new=TRUE)
  plot.ecdf(dates_within, col='red', xlim=c(0,1))
}

dev.off()
