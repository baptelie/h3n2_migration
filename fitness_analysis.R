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
nsim=10

start=Sys.time()
trees.sm <- make.simmap(dst, region, model='ARD', nsim=nsim, Q='mcmc', burnin=10000, samplefreq=200, pi='estimated')
end=Sys.time()
print(end-start)

write.simmap(trees.sm, 'trees.sm')

pdf('stochastic_mappings.pdf')
cols=setNames(brewer.pal(9,'Set1'),sort(unique(region)))
for (t in 1:nsim){
  plotSimmap(trees.sm[[t]], fsize=0.1, lwd=1, colors=cols)
  add.simmap.legend(prompt=FALSE, x=0.9*par()$usr[1], y=0.9*par()$usr[4], colors=cols)
  axisPhylo()
}
dev.off()

#Analyze each stochastic mapping
nreg=length(levels(region))
C <- rep(0,nreg^2) #count the nb of times in which the probability of the migrant fitness givewn the donor region fitness distribution is >10%
D <- rep(0,nreg^2) #count the nb of times where each pair of moigration can be analyzed
slope <- array(0,dim=c(nreg, nreg, 3), dimnames=list(levels(region), levels(region), c('inf','equal','sup')))

# setwd("./fitness_evol")
for (t in 1:nsim){
  cnt <- rep(0,nreg^2)
  tre.sm <- trees.sm[[t]]
  
  #identify the migration events
  listClades <- list.clades(tre.sm)
  membersClades <- getMembers(listClades, tre.sm)
  dt.edges <- dates.edges()
  PM<-pairsMigr(listClades)
  print(PM)
  print(sum(PM))
  
  #analysis of the migrant fitness probability
  LPM <- listPairsMigr(listClades)
  x <- sapply(levels(region), function(r) sapply(levels(region), function(d) proba.obs(d,r,LPM)))
  cnt[!is.na(x)] <- 1
  x[is.na(x)] <- 0
  sup <- (x>0.025 & x<0.975)*rep(1, length(x))
  C <- C+sup
  D <- D+cnt

  #plot the fitness evolution after migration
  for(don in levels(region)) for(rec in levels(region)){
    if(don==rec) next()
    
    FE <- fitness_evol_sis(don, rec, LPM)
    if(!is.na(FE[1])) {
      for(i in 1:length(FE$slope)){
        if (FE$slope[i]<0 & FE$CI95[2*i]<0) {slope[don,rec,'inf'] <- slope[don,rec,'inf'] + 1
        } else if (FE$slope[i]>0 & FE$CI95[2*i-1]>0) {slope[don,rec,'sup'] <- slope[don,rec,'sup'] + 1
        } else slope[don,rec,'equal'] <- slope[don,rec,'equal'] + 1
      }
      slope[don,rec,] <- slope[don,rec,] / length(FE$slope)
    }
  }
}
# round(C/D, digits=2)
library(plotly)
xax <- list(title='Recipient', titlefont= list(size=20), tickfont=list(size=15))
yax <- list(title='Donor', titlefont= list(size=20), tickfont=list(size=15))
f0 = list(size=25)

p <- plotly::plot_ly(z=C/D, type='heatmap', x=levels(region), y=levels(region), colors = c('red','green')) %>%
  layout(title='Frequency among the stochastic mappings that the probability of the migrants fitnesses given the donor region fitness distribution is >10%', xaxis=xax, yaxis=yax)
print(p)

p <- plotly::plot_ly(z=matrix(D, nrow=length(levels(region))), type='heatmap', x=levels(region), y=levels(region), colors = c('red','green')) %>%
  layout(title='total number of observations', xaxis=xax, yaxis=yax)
print(p)

for(don in levels(region)) for(rec in levels(region)){
  slope[don,rec,] <- slope[don,rec,]/sum(slope[don,rec,])
}
# round(slope, digits=2)

p <- plotly::plot_ly(z=slope[,,'inf'], zmin=0, zmax=1, type='heatmap', x=levels(region), y=levels(region), colors = c('red','green')) %>%
  layout(title='frequency of evolution significantly negative', xaxis=xax, yaxis=yax)
print(p)
p <- plotly::plot_ly(z=slope[,,'sup'],  zmin=0, zmax=1, type='heatmap', x=levels(region), y=levels(region), colors = c('red','green')) %>%
  layout(title='Frequency of evolution significantly positive', xaxis=xax, yaxis=yax, titlefont=f0 )
print(p)

for(VS in unique(meta_tree$VaxStrain)){
  nodes <- meta_tree[meta_tree$VaxStrain==VS,]
  D <- nodes$Decimal_Date-2010.523
  LM <- lm(nodes$Fitness~D+0)
  print(VS)
  print(LM)
  plot(D, nodes$Fitness)
  abline(LM)
}