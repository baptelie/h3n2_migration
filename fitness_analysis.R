library(phytools); library(ape); library(Biostrings); library(seqinr);
library(lubridate); library(magrittr); library(dplyr); library(stringr);

pattern='15-18'

setwd('/Users/belie/Documents/Phylogeny/World_12-18/world_15-18')

#Run treetime and load the tree
system(paste('/anaconda3/bin/treetime ancestral --aln nt_world_', pattern, '_subsamp.fasta --tree ', "MCCtmp.nwk", ' --outdir treetime_ancestral', sep=''), wait=TRUE)
tre.tt <- ape::read.nexus("./treetime_ancestral/annotated_tree.nexus")

#Load and process metadata
meta = read.csv(paste('data_world_',pattern,'_subsamp.csv', sep=''))

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

#compute the fitness on each node
Alignment <- readDNAStringSet (paste("./treetime_ancestral", "/ancestral_sequences.fasta", sep=""))
AA_Prefs <- read.csv(file= '/Users/belie/Documents/Phylogeny/New\ folder/AA_prefs_avg.csv', header = TRUE)
Alignment_AA <- Biostrings::translate(Alignment)
Alignment_AA <- Alignment_AA[c(tre.tt$tip.label, tre.tt$node.label)] #reorder in the order of the absolute node number
list_mut_AA = apply(tre.tt$edge, 1, function(edge) mismatches(Alignment_AA[ edge[2] ], Alignment_AA[ edge[1] ]))
list_ME = lapply(list_mut_AA, function(ListMut) Sum_ME(ListMut))
Fitness <- FitnessNodes(tre.tt) #compute the fitness of each node/tip, with 0 for the root
meta_tree$Fitness <- Fitness

#identify the homogeneous antigenic clades
meta_tree <- define_HAC(meta_tree)

#Run the MCMC stochastic mapping
dst<- tre.tt
dst$edge.length[dst$edge.length==0]<-max(nodeHeights(tre.tt))*1e-6 # set zero-length branches to be 1/1000000 total tree length
nsim=10
trees.sm <- make.simmap(dst, region, model='ARD', nsim=nsim, Q='mcmc')

#Analyze each stochastic mapping
C <- rep(0,64)
for (t in 1:nsim){
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
  sup <- (x>0.1)*rep(1, length(x))
  C <- C+sup
}

LM <- matrix(C, nrow=length(unique(region)), byrow=TRUE)
rownames(LM) <- levels(region)
colnames(LM) <- levels(region)
LM