library(phytools); library(ape); library(Biostrings); library(seqinr);
library(lubridate); library(magrittr); library(dplyr); library(stringr); library(RColorBrewer)

pattern='world_13-19'

setwd(paste("~/Documents/Phylogeny/", pattern, '.4/small',sep=""))

#Run treetime and load the tree
system(paste('/anaconda3/bin/treetime --aln nt_', pattern, '_subsamp.fasta --tree ', "tree_",pattern, '_5run.treefile --clock-rate 0.004 --dates dates_tt.csv --outdir treetime004', sep=''), wait=TRUE)
tre.tt <- ape::read.nexus("./treetime_small_v2/timetree.nexus")

#Load and process metadata
meta = read.csv(paste('data_',pattern,'_subsamp_small.csv', sep=''))

sts <- decimal_date(ymd(meta$Collection_Date))
names(sts) <- meta$Isolate_Id
sts <- sts[tre.tt$tip.label]

region <- meta$Region
names(region) <- meta$Isolate_Id

nb_per_month_reg <- matrix(0,nrow=nreg, ncol=(2019.25-2013.5-(2/12))*12)
rownames(nb_per_month_reg)<- regions
colnames(nb_per_month_reg)<- round(seq(from=2013.5+(2/12), to=2019.25-(1/12), by=1/12), digits=2)

for(r in unique(region)){
  sub <- meta[(meta$Region == r),]
  col<-0
  for(m in seq(from=2013.5+(3/12), to=2019.25, by=1/12)){
    col<-col+1
    sub_month <- sub[sub$Decimal_Date <m & sub$Decimal_Date >= m-(1/12), ]
    nb_per_month_reg[r,col]<- nrow(sub_month)
  }
}

write.csv(nb_per_month_reg, 'nb_seq_per_month_reg.csv')

meta_nodes <- data.frame(Isolate_Id = tre.tt$node.label)
tip.order <- c(1:length(tre.tt$tip.label))
names(tip.order) <- meta$Isolate_Id
tip.order <- tip.order[tre.tt$tip.label]
meta_tree <- bind_rows(meta[tip.order, ], meta_nodes) #build a new metadata folder with all computed informations from all nodes, in the order of the absolute node numbers
meta_tree$Decimal_Date <- timeNodes(tre.tt, meta)
row.names(meta_tree) <- meta_tree$Isolate_Id
meta_tree$Isolate_Id <- NULL

rm(meta)

#compute the fitness on each node
Alignment <- readDNAStringSet (paste("./treetime", "004/ancestral_sequences.fasta", sep=""))
AA_Prefs <- read.csv(file= '/Users/belie/Documents/Phylogeny/New\ folder/AA_prefs_avg.csv', header = TRUE)
Alignment <- Alignment[c(tre.tt$tip.label, tre.tt$node.label)]
Alignment_AA <- Biostrings::translate(Alignment)
Alignment_AA <- Alignment_AA[c(tre.tt$tip.label, tre.tt$node.label)] #reorder in the order of the absolute node number
list_mut_AA = apply(tre.tt$edge, 1, function(edge) mismatches(Alignment_AA[ edge[2] ], Alignment_AA[ edge[1] ]))
# list_mut_nt = apply(tre.tt$edge, 1, function(edge) mismatches(Alignment[ edge[2] ], Alignment[ edge[1] ]))
# nb_mut_ns = sapply(list_mut_AA, function(e) nrow(e))
# nb_mut_s = sapply(list_mut_nt, function(e) nrow(e)) - nb_mut_ns
# ka_ks = apply(tre.tt$edge, 1, function(edge) kaks(DNAbin2Aln(Alignment[c(edge[2], edge[1]) ]) ))
# kadks = sapply(ka_ks, function(k) k$ka/k$ks)
list_ME = lapply(list_mut_AA, function(ListMut) Sum_ME(ListMut))
Fitness <- FitnessNodes(tre.tt) #compute the fitness of each node/tip, with 0 for the root
meta_tree$Fitness <- Fitness


#identify the homogeneous antigenic clades
meta_tree <- define_HAC(meta_tree)

#Run the MCMC stochastic mapping
# Q <- mcmcMk(dst,region, model='ARD', ngen=2000, print=500, prop.var=0.0005)
# Q <- mcmcMk(dst,region, model='ARD', prop.var=0.005, print=500)


# dst_small<- tre.tt_small
# dst_small$edge.length[dst_small$edge.length==0]<-max(nodeHeights(tre.tt_small))*1e-6 # set zero-length branches to be 1/1000000 total tree length
# nsim=25

dst<- tre.tt
dst$edge.length[dst$edge.length==0]<-max(nodeHeights(tre.tt))*1e-6 # set zero-length branches to be 1/1000000 total tree length
nsim=48

fit <- fitMk(dst, region, model='ARD', pi='estimated')
fittedQ<-matrix(NA,length(fit$states),length(fit$states))
fittedQ[]<-c(0,fit$rates)[fit$index.matrix+1]
diag(fittedQ)<-0
diag(fittedQ)<--rowSums(fittedQ)
colnames(fittedQ)<-rownames(fittedQ)<-fit$states
## ready to run our analysis
library(parallel)
trees.sm<-mclapply(1:4,function(n,tree,x,fixedQ) make.simmap(tree,x,Q=fixedQ,nsim=10),
                   tree=dst,x=region,fixedQ=fittedQ,mc.cores=if(.Platform$OS.type=="windows") 1L else 4L)

trees.sm<-do.call(c,trees.sm)
if(!("multiSimmap"%in%class(trees))) class(trees)<-c("multiSimmap",class(trees))
trees

tmp <- c()
for(t in trees.sm) tmp <- c(tmp, t)
trees.sm <- tmp
rm(tmp)

start=Sys.time()
trees.sm <- make.simmap(dst, region, model='ARD', nsim=nsim, Q='mcmc', burnin=10000, samplefreq=200, pi='estimated')
end=Sys.time()
print(end-start)

write.simmap(trees.sm, 'trees.sm25')
trees.sm <- read.simmap(file='trees.sm25')
pdf('stochastic_mappings_2.pdf')
cols=setNames(brewer.pal(9,'Set1'),sort(unique(region)))
for (t in 1:nsim){
  plotSimmap(trees.sm[[t]], fsize=0.1, lwd=1, colors=cols)
  add.simmap.legend(prompt=FALSE, x=0.9*par()$usr[1], y=0.9*par()$usr[4], colors=cols)
  axisPhylo()
}
dev.off()

#Identify the migration events of each stochastic realization
Migr_per_tree <- sapply(trees.sm, function(t){
  getMembers_succes(t)
})

#Analyze each stochastic mapping
nreg=length(levels(region))

regions <- levels(region)

An <-mclapply(1:40, function(t) {
  #t = tree number to analyze
    C <- rep(0,nreg^2) #count the nb of times in which the probability of the migrant fitness givewn the donor region fitness distribution is >10%
    D <- rep(0,nreg^2) #count the nb of times where each pair of moigration can be analyzed
    E <- matrix(0, nrow=nreg, ncol=nreg, dimnames=list(regions, regions))
    slope <- array(0,dim=c(nreg, nreg, 3), dimnames=list(regions, regions, c('inf','equal','sup')))
    cnt <- rep(0,nreg^2)
    tre.sm <- trees.sm[[t]]
    
    membersClades <- Migr_per_tree[[t]] #migration events per region pair
    dt.edges <- dates.edges(tree=tre.sm)
    PM<-pairsMigr(t)
    GS <-c(getStates(tree, type='tips'), getStates(tree, type='nodes')) #region of each tip/node in the order : tips and then nodes
    
    #analysis of the migrant fitness probability
    LPM <- listPairsMigr(PM,t)
    # x <- sapply(regions, function(r) sapply(regions, function(d) proba.obs(d,r,LPM,PM, dt.edges)))
    # cnt[!is.na(x)] <- 1
    # x[is.na(x)] <- 0
    # sup <- ( x>0.1)*rep(1, length(x))
    # C <- C+sup
    # D <- D+cnt
    
    for(don in regions) for(rec in regions){
      if(don==rec) next()
      FE <- fitness_evol_vs_TimeSlice(don, rec, LPM, PM, dt.edges, membersClades)
      if(!is.na(FE[1])) {
        E[don, rec] <- E[don, rec]+1
        nbobs <- length(FE$slope)
        for(i in 1:nbobs){
          if (FE$slope[i]<0 & FE$dif[i]) {slope[don,rec,'inf'] <- slope[don,rec,'inf'] + 1/nbobs
          } else if (FE$slope[i]>0 & FE$dif[i]) {slope[don,rec,'sup'] <- slope[don,rec,'sup'] + 1/nbobs
          } else slope[don,rec,'equal'] <- slope[don,rec,'equal'] + 1/nbobs
        }
      }
    }
    return(list(C=C,D=D,E=E,slope=slope))
  },
  mc.cores=if(.Platform$OS.type=="windows") 1L else 4L)

C <- rep(0,nreg^2) #count the nb of times in which the probability of the migrant fitness givewn the donor region fitness distribution is >10%
D <- rep(0,nreg^2) #count the nb of times where each pair of moigration can be analyzed
E <- matrix(0, nrow=nreg, ncol=nreg, dimnames = list(regions, regions))
slope <- array(0,dim=c(nreg, nreg, 3), dimnames=list(regions, regions, c('inf','equal','sup')))

for(i in 1:nsim) {
  slope <- slope+An[[i]]$slope
  E <- E+ An[[i]]$E
  C <- C + An[[i]]$C
  D <- D + An[[i]]$D
}

# setwd("./fitness_evol")
for (t in 1:nsim){
  cat(paste('tree #', t,'\n'))
  cnt <- rep(0,nreg^2)
  tre.sm <- trees.sm[[t]]
  
  #identify the migration events
  membersClades <- Migr_per_tree[[t]]
  dt.edges <- dates.edges()
  PM<-pairsMigr(t)
  cat('total number of events:',sum(PM),'\n')
  flush.console()
  
  #analysis of the migrant fitness probability
  LPM <- listPairsMigr(PM,t)
  x <- sapply(levels(region), function(r) sapply(levels(region), function(d) proba.obs(d,r,LPM)))
  cnt[!is.na(x)] <- 1
  x[is.na(x)] <- 0
  sup <- ( x>0.1)*rep(1, length(x))
  C <- C+sup
  D <- D+cnt

  # plot the fitness evolution after migration
  # pdf('fitness_evol_plots_vsTimeSlice.pdf')
  # DonorVS <- DonorVS_slopes2()

  for(don in levels(region)) for(rec in levels(region)){
    if(don==rec) next()
    FE <- fitness_evol_vs_TimeSlice(don, rec, LPM)
    if(!is.na(FE[1])) {
      E[don, rec] <- E[don, rec]+1
      nbobs <- length(FE$slope)
      for(i in 1:nbobs){
        if (FE$slope[i]<0 & FE$dif[i]) {slope[don,rec,'inf'] <- slope[don,rec,'inf'] + 1/nbobs
        } else if (FE$slope[i]>0 & FE$dif[i]) {slope[don,rec,'sup'] <- slope[don,rec,'sup'] + 1/nbobs
        } else slope[don,rec,'equal'] <- slope[don,rec,'equal'] + 1/nbobs
      }
    }
  }
  # dev.off()
}

for(i in 1:3) slope[,,i] <- slope[,,i]/E
# round(C/D, digits=2)
library(plotly)
xax <- list(title='Recipient', titlefont= list(size=18), tickfont=list(size=15))
yax <- list(title='Donor', titlefont= list(size=18), tickfont=list(size=15))
f0 = list(size=20)

Ctrim <- C
for(i in 1:(nreg*nreg)) if(D[i] < 30) Ctrim[i]<-NA

p <- plotly::plot_ly(z=Ctrim/D, type='heatmap', x=levels(region), y=levels(region), colors = c('red','green'), zmin=0, zmax=1) %>%
  layout(title='Frequency among the stochastic mappings that the probability of the migrants fitnesses given the donor region fitness distribution is >5%', xaxis=xax, yaxis=yax)
print(p)

p <- plotly::plot_ly(z=matrix(E, nrow=length(levels(region))), type='heatmap', x=levels(region), y=levels(region), colors = c('red','green')) %>%
  layout(title='total number of observations', xaxis=xax, yaxis=yax)
print(p)
# 
# for(don in levels(region)) for(rec in levels(region)){
#   slope[don,rec,] <- slope[don,rec,]/sum(slope[don,rec,])
# }
# round(slope, digits=2)


Strim <- slope
for(i in 1:(nreg*nreg)) if(E[i] < 30) {Strim[i]<-NA}

p <- plotly::plot_ly(z=Strim[,,'inf'], zmin=0, zmax=1, type='heatmap', x=levels(region), y=levels(region), colors = c('red','green')) %>%
  layout(title='Frequency of evolution significantly negative', xaxis=xax, yaxis=yax, titlefont=f0)
print(p)
p <- plotly::plot_ly(z=Strim[,,'sup'],  zmin=0, zmax=1, type='heatmap', x=levels(region), y=levels(region), colors = c('red','green')) %>%
  layout(title='Frequency of evolution significantly positive', xaxis=xax, yaxis=yax, titlefont=f0 )
print(p)
p <- plotly::plot_ly(z=Strim[,,'equal'],  zmin=0, zmax=1, type='heatmap', x=levels(region), y=levels(region), colors = c('red','green')) %>%
  layout(title='Frequency of evolution not different', xaxis=xax, yaxis=yax, titlefont=f0 )
print(p)
p <- plotly::plot_ly(z=Strim[,,'inf']-Strim[,,'sup'], zmin=-0.4, zmax=0.4, type='heatmap', x=levels(region), y=levels(region), colors = c('red','black','green')) %>%
  layout(title='Frequency of evolution significantly negative', xaxis=xax, yaxis=yax, titlefont=f0)
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


cntPerReg <- sapply(unique(region), function(r) {c <- 0; for(i in region) if(i==r) c<-c+1; c})
S <- summary(trees.sm)
Sum <- apply(S$Tr,2,function(i) sum(i))
plot(cntPerReg,Sum)
abline(lm(Sum~cntPerReg))