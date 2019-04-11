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
C <- rep(0,64)
D <- rep(0,64)
setwd("./fitness_evol")
for (t in 1:nsim){
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
  for(don in levels(region)) for(rec in levels(region)){
    if(don==rec) next()
    FE <- fitness_evol(don, rec, LPM)
    if(!is.na(FE[1])) {
      pdf(paste('fitness_evol_',don,'-', rec,'.pdf', sep=''))
      reg=lm(FE$fitness~FE$time)
      plot(FE$time, FE$fitness, main=paste('Fitness evolution during the migration from',don,'to',rec), sub=summary(reg)$r.squared)
      abline(reg)
      dev.off()
      }
  }
}
C/D
LM <- matrix(C, nrow=length(unique(region)), byrow=TRUE)
rownames(LM) <- levels(region)
colnames(LM) <- levels(region)
LM


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

# 
# pdf('fitness_migr_vs_within_per_reg.pdf')
# for (t in 1:nsim){
#   print(t)
#   # cnt <- rep(0,64)
#   tre.sm <- trees.sm[[t]]
#   #identify the migration events
#   listClades <- list.clades(tre.sm)
#   membersClades <- getMembers(listClades, tre.sm)
#   PM<-pairsMigr(listClades)
#   print(sum(PM))
#   LPM <- listPairsMigr(listClades)
#   GS <-c(getStates(tre.sm, type='tips'), getStates(tre.sm, type='nodes')) #region of each tip/node in the order : tips and then nodes
#   
#   for(don in levels(region)) for(rec in levels(region)){
#     if(don==rec) next()
#     fitness_migr <- c()
#     fitness_within <- c()
#     
#     for (n in LPM[toString(don), toString(rec), ][!is.na(LPM[toString(don), toString(rec), ])]){
#       n <- nodelab.to.numb(n)
#       FP <- Fitness[getParent(tre.tt,n)]
#       reg=GS[n]
#       nodes_clade <- membersClades[[nodenumb.to.lab(n)]]
#       fam <- gd_children(n)
#       nodes_migr <- strtoi(check_within(unlist(fam), reg, GS))
#       fitness_migr <- c(fitness_migr, Fitness[nodes_migr]-FP)
#       
#       dates <- meta_tree[nodes_clade,]$Decimal_Date
#       if(length(nodes_clade)>25 & max(dates)-min(dates)>0.5){
#         for(gc in fam$gd_ch){
#           children <- tre.tt$edge[which(tre.tt$edge[,1]==gc),2]
#           for(c in children){
#             f <- unlist(gd_children(node=c))
#             if(length(f)<4) next()
#             if(length(check_within(f, reg, GS))==length(f)) fitness_within <- c(fitness_within, Fitness[strtoi(f)]-Fitness[c]) #if all nodes in the same region
#           }
#         }
#       }
#     }
#     if(length(LPM[toString(don), toString(rec), ][!is.na(LPM[toString(don), toString(rec), ])])>0) plot.ecdf(fitness_migr, xlim=c(-10,8), col='green', main=paste(t,don,rec))
#     if(!is.null(fitness_within)){par(new=TRUE); plot.ecdf(fitness_within, xlim=c(-10,8), col='red', main='')}
#     
#   }
# }
# 
# dev.off()
# 
# pdf('10fitness_migr_vs_within_grandchildren.pdf')
# for (t in 1:nsim){
#   tre.sm <- trees.sm[[t]]
#   
#   #identify the migration events
#   listClades <- list.clades(tre.sm)
#   membersClades <- getMembers(listClades, tre.sm)
#   
#   rates_migr_np2 <- c()
#   rates_within_np2 <- c()
#   rates_migr_np1 <- c()
#   rates_within_np1 <- c()
#   rates_migr_n <- c()
#   rates_within_n <- c()
#   # fitness_migr <- c()
#   dates_migr_np2 <- c()
#   # fitness_within <- c()
#   dates_within_np2 <- c()
#   dates_migr_n <- c()
#   dates_within_n <- c()
#   dates_migr_np1 <- c()
#   dates_within_np1 <- c()
#   count <- 0
#   GS <-c(getStates(tre.sm, type='tips'), getStates(tre.sm, type='nodes')) #region of each tip/node in the order : tips and then nodes
#   print(length(listClades$Node))
#   for (n in listClades$Node){
#     FP <- Fitness[getParent(tre.tt,n)]
#     reg=GS[n]
#     DP <- meta_tree[getParent(tre.tt,n),]$Decimal_Date
#     nodes_clade <- membersClades[[nodenumb.to.lab(n)]]
#     fam <- gd_children(n)
#     nodes_migr_np2 <- strtoi(check_within(fam$gd_ch, reg, GS))
#     nodes_migr_np1 <- strtoi(check_within(fam$ch, reg, GS))
#     nodes_migr_n <- strtoi(check_within(fam$node, reg, GS))
#     
#     # fitness_migr <- c(fitness_migr, Fitness[nodes_migr]-FP)
#     rates_mnp2 <- (Fitness[nodes_migr_np2]-FP) / tre.tt$edge.length[nodes_migr_np2]
#     rates_migr_np2 <- c(rates_migr_np2, rates_mnp2)
#     dates_migr_np2 <- c(dates_migr_np2, tre.tt$edge.length[nodes_migr_np2])
#     
#     rates_mnp1 <- (Fitness[nodes_migr_np1]-FP) / tre.tt$edge.length[nodes_migr_np1]
#     rates_migr_np1 <- c(rates_migr_np1, rates_mnp1)
#     dates_migr_np1 <- c(dates_migr_np1, tre.tt$edge.length[nodes_migr_np1])
#     
#     rates_mn <- (Fitness[nodes_migr_n]-FP) / tre.tt$edge.length[nodes_migr_n]
#     rates_migr_n <- c(rates_migr_n, rates_mn)
#     dates_migr_n <- c(dates_migr_n, tre.tt$edge.length[nodes_migr_n])
# 
#     if(length(nodes_clade)>20){
#       for(gc in fam$gd_ch){
#         children <- tre.tt$edge[which(tre.tt$edge[,1]==gc),2]
#         for(c in children){
#           f <- gd_children(node=c)
#           if(length(f$gd_ch)==0) next()
#           if(length(check_within(unlist(f), reg, GS))==length(f)){ #if all the n, n+1 & n+2 are still within the same clade
# 
#             rates_wnp2 <- (Fitness[strtoi(f$gd_ch)]-Fitness[gc]) / tre.tt$edge.length[strtoi(f$gd_ch)]
#             rates_within_np2 <- c(rates_within_np2, rates_wnp2)
#             dates_within_np2 <- c(dates_within_np2, tre.tt$edge.length[strtoi(f$gd_ch)])
#             
#             rates_wnp1 <- (Fitness[strtoi(f$ch)]-Fitness[gc]) / tre.tt$edge.length[strtoi(f$ch)]
#             rates_within_np1 <- c(rates_within_np1, rates_wnp1)
#             dates_within_np1 <- c(dates_within_np1, tre.tt$edge.length[strtoi(f$ch)])
#             
#             rates_wn <- (Fitness[strtoi(f$node)]-Fitness[gc]) / tre.tt$edge.length[strtoi(f$node)]
#             rates_within_n <- c(rates_within_n, rates_wn)
#             dates_within_n <- c(dates_within_n, tre.tt$edge.length[strtoi(f$node)])
#             
#           } 
#         }
#       }
#     }
#   }
#   plot.ecdf(rates_migr_np2[!rates_migr_np2==0], xlim=c(-150,100), ylim=c(0,1), col='green', main=paste('fitness evolution rate, within a region vs during a migration, n+2, tree',t ))
#   par(new=TRUE)
#   try(plot.ecdf(rates_within_np2[!rates_within_np2==0], xlim=c(-150,100),ylim=c(0,1), col='red', main=''))
#   plot.ecdf(dates_migr_np2[!rates_migr_np2==0], col='green', xlim=c(0,0.8), main='distribution of the branch lengths, n+2')
#   par(new=TRUE)
#   plot.ecdf(dates_within_np2[!rates_within_np2==0], col='red', xlim=c(0,0.8), main='')
#   
#   plot.ecdf(rates_migr_np1[!rates_migr_np1==0], xlim=c(-150,100), ylim=c(0,1), col='green', main=paste('fitness evolution rate, within a region vs during a migration, n+1, tree',t ))
#   par(new=TRUE)
#   try(plot.ecdf(rates_within_np1[!rates_within_np1==0], xlim=c(-150,100),ylim=c(0,1), col='red', main=''))
#   plot.ecdf(dates_migr_np1[!rates_migr_np1==0], col='green', xlim=c(0,0.8), main='distribution of the branch lengths, n+1')
#   try(par(new=TRUE); plot.ecdf(dates_within_np1[!rates_within_np1==0], col='red', xlim=c(0,0.8), main=''))
#   
#   plot.ecdf(rates_migr_n[!rates_migr_n==0], xlim=c(-150,100), ylim=c(0,1), col='green', main=paste('fitness evolution rate, within a region vs during a migration, n, tree',t ))
#   try(par(new=TRUE);plot.ecdf(rates_within_n[!rates_within_n==0], xlim=c(-150,100),ylim=c(0,1), col='red', main=''))
#   plot.ecdf(dates_migr_n[!rates_migr_n==0], col='green', xlim=c(0,0.8), main='distribution of the branch lengths, n')
#   try(par(new=TRUE); plot.ecdf(dates_within_n[!rates_within_n==0], col='red', xlim=c(0,0.8), main=''))
#   # plot.ecdf(rates_migr[!rates_migr==0], xlim=c(-120,70), ylim=c(0,1), col='green')
#   # par(new=TRUE)
#   # plot.ecdf(rates_within[!rates_within==0], xlim=c(-120,70),ylim=c(0,1), col='red')
# }
# 
# dev.off()
# 
# pdf('plotSimmap.pdf')
# cols<-setNames(palette()[1:length(unique(region))],sort(unique(region))) #choose the colors for each region
# for(t in 1:nsim){
#   plotSimmap(trees.sm[[t]], fsize=0.1, lwd=0.5, colors=cols)
#   add.simmap.legend(colors = cols, prompt=FALSE, x=0.95*par()$usr[1],y=0.1*par()$usr[3])
# }
# dev.off()

# 
# pdf('10fitness_migr_vs_within_grandchildren.pdf')
# for (t in 1:nsim){
#   tre.sm <- trees.sm[[t]]
# 
#   #identify the migration events
#   listClades <- list.clades(tre.sm)
#   membersClades <- getMembers(listClades, tre.sm)
# 
#   rates_migr_np2 <- c()
#   rates_within_np2 <- c()
#   rates_migr_np1 <- c()
#   rates_within_np1 <- c()
#   rates_migr_n <- c()
#   rates_within_n <- c()
#   # fitness_migr <- c()
#   dates_migr_np2 <- c()
#   # fitness_within <- c()
#   dates_within_np2 <- c()
#   dates_migr_n <- c()
#   dates_within_n <- c()
#   dates_migr_np1 <- c()
#   dates_within_np1 <- c()
#   count <- 0
#   GS <-c(getStates(tre.sm, type='tips'), getStates(tre.sm, type='nodes')) #region of each tip/node in the order : tips and then nodes
#   print(length(listClades$Node))
#   for (n in listClades$Node){
#     FP <- Fitness[getParent(tre.tt,n)]
#     reg=GS[n]
#     DP <- meta_tree[getParent(tre.tt,n),]$Decimal_Date
#     nodes_clade <- membersClades[[nodenumb.to.lab(n)]]
#     fam <- gd_children(n)
#     nodes_migr_np2 <- strtoi(check_within(fam$gd_ch, reg, GS))
#     nodes_migr_np1 <- strtoi(check_within(fam$ch, reg, GS))
#     nodes_migr_n <- strtoi(check_within(fam$node, reg, GS))
# 
#     # fitness_migr <- c(fitness_migr, Fitness[nodes_migr]-FP)
#     rates_mnp2 <- (Fitness[nodes_migr_np2]-FP)
#     rates_migr_np2 <- c(rates_migr_np2, rates_mnp2)
#     dates_migr_np2 <- c(dates_migr_np2, tre.tt$edge.length[nodes_migr_np2])
# 
#     rates_mnp1 <- (Fitness[nodes_migr_np1]-FP)
#     rates_migr_np1 <- c(rates_migr_np1, rates_mnp1)
#     dates_migr_np1 <- c(dates_migr_np1, tre.tt$edge.length[nodes_migr_np1])
# 
#     rates_mn <- (Fitness[nodes_migr_n]-FP)
#     rates_migr_n <- c(rates_migr_n, rates_mn)
#     dates_migr_n <- c(dates_migr_n, tre.tt$edge.length[nodes_migr_n])
# 
#     if(length(nodes_clade)>20){
#       for(gc in fam$gd_ch){
#         children <- tre.tt$edge[which(tre.tt$edge[,1]==gc),2]
#         for(c in children){
#           f <- gd_children(node=c)
#           if(length(f$gd_ch)==0) next()
#           if(length(check_within(unlist(f), reg, GS))==length(f)){ #if all the n, n+1 & n+2 are still within the same clade
# 
#             rates_wnp2 <- (Fitness[strtoi(f$gd_ch)]-Fitness[gc])
#             rates_within_np2 <- c(rates_within_np2, rates_wnp2)
#             dates_within_np2 <- c(dates_within_np2, tre.tt$edge.length[strtoi(f$gd_ch)])
# 
#             rates_wnp1 <- (Fitness[strtoi(f$ch)]-Fitness[gc])
#             rates_within_np1 <- c(rates_within_np1, rates_wnp1)
#             dates_within_np1 <- c(dates_within_np1, tre.tt$edge.length[strtoi(f$ch)])
# 
#             rates_wn <- (Fitness[strtoi(f$node)]-Fitness[gc])
#             rates_within_n <- c(rates_within_n, rates_wn)
#             dates_within_n <- c(dates_within_n, tre.tt$edge.length[strtoi(f$node)])
# 
#           }
#         }
#       }
#     }
#   }
#   plot.ecdf(rates_migr_np2, xlim=c(-15,10), ylim=c(0,1), col='green', main=paste('fitness evolution rate, within a region vs during a migration, n+2, tree',t ))
#   par(new=TRUE)
#   plot.ecdf(rates_within_np2, xlim=c(-15,10),ylim=c(0,1), col='red', main='')
#   plot.ecdf(dates_migr_np2, col='green', xlim=c(0,0.8), main='distribution of the branch lengths, n+2')
#   par(new=TRUE)
#   plot.ecdf(dates_within_np2, col='red', xlim=c(0,0.8), main='')
# 
#   plot.ecdf(rates_migr_np1, xlim=c(-15,10), ylim=c(0,1), col='green', main=paste('fitness evolution rate, within a region vs during a migration, n+1, tree',t ))
#   par(new=TRUE)
#   plot.ecdf(rates_within_np1, xlim=c(-15,10),ylim=c(0,1), col='red', main='')
#   plot.ecdf(dates_migr_np1, col='green', xlim=c(0,0.8), main='distribution of the branch lengths, n+1')
#   par(new=TRUE)
#   plot.ecdf(dates_within_np1, col='red', xlim=c(0,0.8), main='')
# 
#   plot.ecdf(rates_migr_n, xlim=c(-15,10), ylim=c(0,1), col='green', main=paste('fitness evolution rate, within a region vs during a migration, n, tree',t ))
#   par(new=TRUE)
#   plot.ecdf(rates_within_n, xlim=c(-15,10),ylim=c(0,1), col='red', main='')
#   plot.ecdf(dates_migr_n, col='green', xlim=c(0,0.8), main='distribution of the branch lengths, n')
#   par(new=TRUE)
#   plot.ecdf(dates_within_n, col='red', xlim=c(0,0.8), main='')
#   # plot.ecdf(rates_migr[!rates_migr==0], xlim=c(-120,70), ylim=c(0,1), col='green')
#   # par(new=TRUE)
#   # plot.ecdf(rates_within[!rates_within==0], xlim=c(-120,70),ylim=c(0,1), col='red')
# }
# 
# dev.off()
# 


pdf('10kaks_diff_migr_vs_within.pdf')
for (t in 1:nsim){
  tre.sm <- trees.sm[[t]]
  
  #identify the migration events
  listClades <- list.clades(tre.sm)
  print(length(listClades$Node))
  membersClades <- getMembers(listClades, tre.sm)
  
  rates_migr_np2 <- c()
  rates_within_np2 <- c()
  dates_migr_np2 <- c()
  dates_within_np2 <- c()

  GS <-c(getStates(tre.sm, type='tips'), getStates(tre.sm, type='nodes')) #region of each tip/node in the order : tips and then nodes
  
  for (n in listClades$Node){
    reg=GS[n]
    DP <- meta_tree[getParent(tre.tt,n),]$Decimal_Date
    nodes_clade <- membersClades[[nodenumb.to.lab(n)]]
    fam <- gd_children(n)
    nodes_migr_np2 <- strtoi(check_within(fam$gd_ch, reg, GS))
    
    ka_ks <- sapply(nodes_migr_np2, function(node) kaks(DNAbin2Aln(Alignment[c(node, strtoi(n))])) )
    ks <- unlist(ka_ks['ks',])
    ka <- unlist(ka_ks['ka',])[which(ks>=0)]
    ks <- ks[ks>=0]
    rates_migr_np2 <- c(rates_migr_np2, ka-ks)
    dates_migr_np2 <- c(dates_migr_np2, tre.tt$edge.length[nodes_migr_np2])
    
    if(length(nodes_clade)>20){
      for(gc in fam$gd_ch){
        children <- tre.tt$edge[which(tre.tt$edge[,1]==gc),2]
        for(c in children){
          f <- gd_children(node=c)
          if(length(f$gd_ch)==0) next()
          if(length(check_within(unlist(f), reg, GS))==length(f)){ #if all the n, n+1 & n+2 are still within the same clade
            
            ka_ks <- sapply(strtoi(f$gd_ch), function(node) kaks(DNAbin2Aln(Alignment[c(node, gc) ])))
            ks <- unlist(ka_ks['ks',])
            ka <- unlist(ka_ks['ka',])[which(ks>=0)]
            ks <- ks[ks>=0]
            rates_within_np2 <- c(rates_within_np2, ka-ks)
            dates_within_np2 <- c(dates_within_np2, tre.tt$edge.length[strtoi(f$gd_ch)])
          }
        }
      }
    }
  }
  # rates_migr_np2[is.infinite(rates_migr_np2)] <- 3
  # rates_within_np2[is.infinite(rates_within_np2)] <- 3
  # rates_migr_np2 <- rates_migr_np2[!c(is.na(rates_migr_np2)|rates_migr_np2<0)]
  # dates_migr_np2 <- dates_migr_np2[!c(is.na(rates_migr_np2)|rates_migr_np2<0)]
  # rates_migr_np2 <- rates_migr_np2[dates_migr_np2>0.02]
  # dates_migr_np2 <- dates_migr_np2[dates_migr_np2>0.02]
  
  # rates_within_np2 <- rates_within_np2[!c(is.na(rates_within_np2)|rates_within_np2<0)]
  # dates_within_np2 <- dates_within_np2[!c(is.na(rates_within_np2)|rates_within_np2<0)]
  # rates_within_np2 <- rates_within_np2[dates_within_np2>0.02]
  # dates_within_np2 <- dates_within_np2[dates_within_np2>0.02]
  # 
  
  print(paste('tree',t))
  print(length(rates_migr_np2[rates_migr_np2>0])/length(rates_migr_np2))
  print(length(rates_within_np2[rates_within_np2>0])/length(rates_within_np2))
  
  plot.ecdf(rates_migr_np2, xlim=c(-0.015,0.002), ylim=c(0,1), col='green', main=paste('fitness evolution rate, within a region vs during a migration, n+2, tree',t ))
  par(new=TRUE)
  plot.ecdf(rates_within_np2, xlim=c(-0.015,0.002),ylim=c(0,1), col='red', main='')
  plot.ecdf(dates_migr_np2, col='green', xlim=c(0,0.8), ylim=c(0,1), main='distribution of the branch lengths, n+2')
  par(new=TRUE)
  plot.ecdf(dates_within_np2, col='red', xlim=c(0,0.8), ylim=c(0,1), main='')

  # plot.ecdf(rates_migr[!rates_migr==0], xlim=c(-120,70), ylim=c(0,1), col='green')
  # par(new=TRUE)
  # plot.ecdf(rates_within[!rates_within==0], xlim=c(-120,70),ylim=c(0,1), col='red')
}

dev.off()

