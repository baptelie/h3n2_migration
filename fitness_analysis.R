library(phytools); library(ape); library(Biostrings); library(seqinr);
library(lubridate); library(magrittr); library(dplyr); library(stringr);

pattern='world_13-19'

setwd(paste("~/Documents/Phylogeny/", pattern, '/ML_tree_1',sep=""))

#Run treetime and load the tree
system(paste('/anaconda3/bin/treetime --aln nt_', pattern, '_subsamp.fasta --tree ', "tree_",pattern, '_subsamp --dates dates.csv --clock-rate 0.004 --outdir treetime_004', sep=''), wait=TRUE)
tre.tt <- ape::read.nexus("./treetime_004/timetree.nexus")
# tre.tt <- drop.tip(tre.tt, c('EPI_ISL_195549','EPI_ISL_197247'))

#Load and process metadata
meta = read.csv(paste('data_',pattern,'_subsamp.csv', sep=''))

sts <- decimal_date(mdy(meta$Collection_Date))
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
Alignment <- readDNAStringSet (paste("./treetime_004", "/ancestral_sequences.fasta", sep=""))
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
nsim=100

start=Sys.time()
trees.sm100 <- make.simmap(dst, region, model='ARD', nsim=nsim, Q='mcmc', burnin=10000, samplefreq=200, pi='estimated')
end=Sys.time()
print(end-start)

write.simmap(trees.sm100, 'trees.sm100')

pdf('stochastic_mappings.pdf')
cols=setNames(brewer.pal(9,'Set1'),sort(unique(region)))
for (t in seq(1,100,10)){
  plotSimmap(trees.sm100[[t]], fsize=0.1, lwd=1, colors=cols)
  add.simmap.legend(prompt=FALSE, x=0.9*par()$usr[1], y=0.9*par()$usr[4], colors=cols)
  axisPhylo()
}
dev.off()

#Analyze each stochastic mapping
nreg=length(levels(region))
C <- rep(0,nreg^2) #count the nb of times in which the probability of the migrant fitness givewn the donor region fitness distribution is >10%
D <- rep(0,nreg^2) #count the nb of times where each pair of moigration can be analyzed
E <- matrix(0, nrow=nreg, ncol=nreg)
rownames(E) <- levels(region)
colnames(E) <- levels(region)
slope <- array(0,dim=c(nreg, nreg, 3), dimnames=list(levels(region), levels(region), c('inf','equal','sup')))
# VS_slopes <- sapply(unique(meta_tree$VaxStrain), function(VS) VS_slope(VS))

# setwd("./fitness_evol")
for (t in 1:nsim){
  cat(paste('tree #', t))
  cnt <- rep(0,nreg^2)
  tre.sm <- trees.sm100[[t]]
  
  #identify the migration events
  listClades <- list.clades(tre.sm)
  membersClades <- getMembers(listClades, tre.sm)
  dt.edges <- dates.edges()
  PM<-pairsMigr(listClades)
  # print(PM)
  # print(sum(PM))
  
  #analysis of the migrant fitness probability
  LPM <- listPairsMigr(listClades)
  x <- sapply(levels(region), function(r) sapply(levels(region), function(d) proba.obs(d,r,LPM)))
  cnt[!is.na(x)] <- 1
  x[is.na(x)] <- 0
  sup <- ( x>0.1)*rep(1, length(x))
  C <- C+sup
  D <- D+cnt

  # plot the fitness evolution after migration
  # pdf('fitness_evol_plots_sm1.1.pdf')
  # DonorVS <- DonorVS_slopes2()

  for(don in levels(region)) for(rec in levels(region)){
    if(don==rec) next()
    FE <- fitness_evol_vs_TimeSlice(don, rec, LPM)
    if(!is.na(FE[1])) {
      E[don, rec] <- E[don, rec]+1
      nbobs <- length(FE$slope)
      for(i in 1:nbobs){
        if (FE$slope[i]<0 & FE$CI95[2*i]<0) {slope[don,rec,'inf'] <- slope[don,rec,'inf'] + 1/nbobs
        } else if (FE$slope[i]>0 & FE$CI95[2*i-1]>0) {slope[don,rec,'sup'] <- slope[don,rec,'sup'] + 1/nbobs
        } else slope[don,rec,'equal'] <- slope[don,rec,'equal'] + 1/nbobs
      }
    }
  }
  # dev.off()
}

for(i in 1:3) slope[,,i] <- slope[,,i]/E
# round(C/D, digits=2)
library(plotly)
xax <- list(title='Recipient', titlefont= list(size=20), tickfont=list(size=15))
yax <- list(title='Donor', titlefont= list(size=20), tickfont=list(size=15))
f0 = list(size=25)

Ctrim <- C
for(i in 1:(nreg*nreg)) if(D[i] < 70) Ctrim[i]<-NA

p <- plotly::plot_ly(z=Ctrim/D, type='heatmap', x=levels(region), y=levels(region), colors = c('red','green'), zmin=0, zmax=1) %>%
  layout(title='Frequency among the stochastic mappings that the probability of the migrants fitnesses given the donor region fitness distribution is >5%', xaxis=xax, yaxis=yax)
print(p)

p <- plotly::plot_ly(z=matrix(D, nrow=length(levels(region))), type='heatmap', x=levels(region), y=levels(region), colors = c('red','green')) %>%
  layout(title='total number of observations', xaxis=xax, yaxis=yax)
print(p)
# 
# for(don in levels(region)) for(rec in levels(region)){
#   slope[don,rec,] <- slope[don,rec,]/sum(slope[don,rec,])
# }
# round(slope, digits=2)

Strim <- slope
for(i in 1:(nreg*nreg)) if(E[i] < 70) {Strim[,,1][i]<-NA; Strim[,,2][i]<-NA; Strim[,,3][i]<-NA}

p <- plotly::plot_ly(z=Strim[,,'inf'], zmin=0, zmax=1, type='heatmap', x=levels(region), y=levels(region), colors = c('red','green')) %>%
  layout(title='frequency of evolution significantly negative', xaxis=xax, yaxis=yax)
print(p)
p <- plotly::plot_ly(z=Strim[,,'sup'],  zmin=0, zmax=1, type='heatmap', x=levels(region), y=levels(region), colors = c('red','green')) %>%
  layout(title='Frequency of evolution significantly positive', xaxis=xax, yaxis=yax, titlefont=f0 )
print(p)
p <- plotly::plot_ly(z=Strim[,,'equal'],  zmin=0, zmax=1, type='heatmap', x=levels(region), y=levels(region), colors = c('red','green')) %>%
  layout(title='Frequency of evolution not different', xaxis=xax, yaxis=yax, titlefont=f0 )
print(p)

ratio <- Strim[,,'inf']/Strim[,,'sup']
ratio[is.infinite(ratio)]<- NA
ratio[Strim[,,'equal']>0.75] <- NA

palette <- colorRampPalette(c("red", 'black', "darkgreen", "green", "lightgreen"))

p <- plotly::plot_ly(z=log10(ratio), type='heatmap', zmin=-0.6, zmax=2.5, x=levels(region), y=levels(region), colors = palette(50)) %>%
  layout(title='Ration of frequency of evolution significantly negative over positive', xaxis=xax, yaxis=yax, titlefont=f0 )
print(p)

meta_tree$abs_nb <- 1:nrow(meta_tree)

pdf('VaccineStrains2.pdf')
for(VS in unique(meta_tree$VaxStrain)[-6]){
  nodes <- (length(tre.tt$tip.label):nrow(meta_tree)) [meta_tree[length(tre.tt$tip.label):nrow(meta_tree),]$VaxStrain==VS]
  
  plotTree(tre.tt, fsize=0.01, main=VS, lwd=0.1)
  nodelabels(node=nodes, pch=21)
}
dev.off()