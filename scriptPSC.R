### Fitness analysis
mismatches <- function(query, ref){
  #find the mutations between a reference sequence and a query sequence
  query <- strsplit(as.character(query), split='')[[1]]
  ref <- strsplit(as.character(ref), split='')[[1]]
  Pos <- which(ref!=query)
  mism <- matrix(c(ref[Pos], Pos, query[Pos]), nrow=length(Pos), ncol=3)
  mism
}

AA_to_col <- function (AA){
  # Given an amino acid, it gives the corresponding column on the data frame AA_Prefs
  listAA = list("A", "C", "D",	"E",	"F",	"G",	"H",	"I",	"K",	"L",	"M",	"N",	"P",	"Q",	"R",	"S",	"T",	"V",	"W",	"Y","*")
  for (i in 1:22){
    if (i==22) stop(paste("Problem of translation, ", AA, " does not exists")) #Problem in the translation
    if (listAA[[i]]==AA) return(i+1)
  }
}

Mut_effect <- function(Mut){
  # Computes the mutational effect of a given mutation, as explained in Bloom et al. 2018
  prefWT <- AA_Prefs[strtoi(Mut[[2]]) , AA_to_col(Mut[[1]]) ]
  prefMut <- AA_Prefs[strtoi(Mut[[2]]) , AA_to_col(Mut[[3]]) ]
  
  if (prefWT == 0 | prefMut == 0) return(0) #assign ME=0 when a stop codon occurs, need to be fixed
  log2(prefMut/prefWT)
}

Sum_ME <- function(ListMut){
  # Computes the sum of the mutational effects between two sequences, taking into account the cases when there is only silent mutations
  if (nrow(ListMut) == 0) return(0)
  ME = sum( apply(ListMut, 1, function(mut) Mut_effect(mut) ) )
  ME
}

FitnessNodes <- function(tree){
  #compute the fitness of each node, as the cummulative ME from the root of the given tree
  print('compute the fitness of each tip and node, as the cummulative ME from the root')
  fitness <- c(rep(NA, length(tree$tip.label)), 0, rep(NA, length(tree$node.label)-1 )) #initialize with NA for each node, except 0 for the root
  names(fitness) <- c(tree$tip.label, tree$node.label)
  for (i in (1:nrow(tree$edge)) ){
    cat('\rWork in progress ... ',round(i/nrow(tree$edge)),'%', sep='')
    flush.console()
    fitness[tree$edge[i,2] ] <- fitness[[tree$edge[i,1] ]] + list_ME[[i ]]
  }
  fitness
}

### Analysis tools of the time aligned tree
timeNodes <- function(tree, meta){
  # List the infered time of each node and tips of a given tree
  tmrca <- max(sts[tree$tip.label])-max(nodeHeights(tree))
  dates <- nodeHeights(tree)[,2] + tmrca
  dates <- c(dates, tmrca) # add the age of the root node
  name.edge <- c(tree$tip.label, tree$node.label)[c(tree$edge[,2], length(tree$tip.label)+1)] #names in the order of tree$edge + the root node
  names(dates) <- name.edge
  
  dates[ c(tree$tip.label, tree$node.label)]
}

extract.clade.simmap2 <- function(tree, node){
  #corrects extract.clade.simmap errors
  clade <- extract.clade.simmap(tree, node=node)
  clade2 <- extract.clade(phy=tree, node=node)
  clade$node.label <- clade2$node.label
  clade$Nnode <- clade2$Nnode
  clade
}

dates.edges <- function(tree=tre.sm, meta=meta_tree){
  #creates a matrix giving for each edge the dates spent on each region
  DR <- array(dim = c(nrow(tree$edge),length(unique(region)),2))
  colnames(DR) <- levels(region)
  rownames(DR) <- sapply(tree$edge[,2], function(e) nodenumb.to.lab(e)) #each row=one edge, named by the child node (unique)
  rw=0
  for(maps in tree$maps){
    rw=rw+1
    d <- meta_tree$Decimal_Date[tree$edge[rw,1]] #take the date of the beginning of the edge
    for(reg in 1:length(maps)){
      DR[rw,names(maps[reg]),1] <- d
      d <- d+maps[reg]
      DR[rw,names(maps[reg]),2] <- d
    }
  }
  DR
}

DNAbin2Aln <- function(list_seq){
  L= length(list_seq[[1]])
  M = M=matrix(c(strsplit(as.character(list_seq[1]), split='')[[1]], strsplit(as.character(list_seq[2]), split='')[[1]] ), ncol=L, byrow=TRUE)
  rownames(M)=names(list_seq)
  ape::as.alignment(M)
}

nodelab.to.numb <- function(nodelab, tree=tre.tt){
  if(grepl('NODE',nodelab)) return(length(tree$tip.label) + which(tree$node.label==nodelab))
  if(grepl('ISL', nodelab)) return(which(tree$tip.label==nodelab))
  print(paste('error in nodelab.to.numb', nodelab, 'is not known', sep=' '))
}

nodenumb.to.lab <- function(nodenumb, tree=tre.tt){
  if(grepl('NODE',nodenumb) | grepl('EPI',nodenumb)) return(nodenumb) #if it is already a nodelab
  if(strtoi(nodenumb)<=length(tree$tip.label)) return(tree$tip.label[strtoi(nodenumb)])
  tree$node.label[strtoi(nodenumb)-length(tree$tip.label)]
}

### Identification of the clades of interest
nodes.sameReg <- function(reg, node, GS, tree, elts=node){
  #count the number of nodes and tips continuously in the same region
  children <- tree$edge[which(tree$edge[,1]==node),2]
  w <- children[which(GS[children]==reg)]
  elts <- c(elts, w)
  if (length(w>0)){
    for (i in w) {
      elts <- nodes.sameReg(reg, i, GS, tree, elts)
    }
  }
  return(elts)
}

getMembers <- function(listClades, tree){
  #lists all the nodes from each clade of interest continuously in the same region
  GS <-c(getStates(tree, type='tips'), getStates(tree, type='nodes')) #region of each tip/node in the order : tips and then nodes
  M <- apply(listClades, 1, function(c) sapply(nodes.sameReg(reg=c[2], node=c[3], GS, tree), function(n) nodenumb.to.lab(n) ) )
  names(M) <- sapply(listClades$Node, function(n) nodenumb.to.lab(n))
  M
}

length.clades <- function(tree, GS){
  #computes the number of nodes in the same region after an identified migration event
  edge.chg <- apply(tree$edge, 1, function(edge) {GS[edge[1]]!=GS[edge[2]] & edge[2]>tree$Nnode+1} ) #true/false in function of whether the edge changes of region, in the order of tree$edge
  node.chg <- tree$edge[,2][edge.chg] #root of each clade in the recipient region
  GSe <- GS[tree$edge[,2]] #reorder GS in function of the edge order
  reg.chg <- GSe[edge.chg] #recipient region
  print(length(node.chg))
  L <- lapply(1:length(reg.chg), function(edge) length(nodes.sameReg(reg.chg[edge], node.chg[edge], GS, tree))  )
  names(L) <- node.chg
  unlist(L)
}

list.clades <- function(tree){
  GS <- c(getStates(tree, type='tips'), getStates(tree, type='nodes'))
  LC <- length.clades(tree, GS)
  LC <- LC[LC>19] #keep clades with at least 20 continuous edges on the same region
  L <- names(LC)
  listRegParents <- unlist(lapply(L, function(node) GS[getParent(tree, node)] ))
  listRegReciep <- unlist(lapply(L, function(node) GS[node] ))
  
  result <- data.frame(Parent_Reg = listRegParents, Reciep_Reg = listRegReciep, Node = L, length = LC)
  result
}

pairsMigr <- function(listClades){
  donorRec_pairs <- matrix(data=0, nrow=length(unique(region)), ncol=length(unique(region)))
  rownames(donorRec_pairs) <- levels(region)
  colnames(donorRec_pairs) <- levels(region)
  
  for (r in 1:nrow(listClades)) donorRec_pairs[ toString(listClades[r,1]), toString(listClades[r,2]) ] <- donorRec_pairs[ toString(listClades[r,1]), toString(listClades[r,2]) ] + 1
  donorRec_pairs
}

listPairsMigr <- function(listClades){
  nbmax <- max(pairsMigr(listClades))
  donorRec_pairs <- array(data=NA, dim=c(length(unique(region)), length(unique(region)), nbmax) )
  rownames(donorRec_pairs) <- levels(region)
  colnames(donorRec_pairs) <- levels(region)
  
  for (r in 1:nrow(listClades)) {
    x <- donorRec_pairs[ toString(listClades[r,1]), toString(listClades[r,2]), ][!is.na(donorRec_pairs[ toString(listClades[r,1]), toString(listClades[r,2]), ])]
    donorRec_pairs[ toString(listClades[r,1]), toString(listClades[r,2]), ] <- c(x, nodenumb.to.lab(toString(listClades[r,3])), rep(NA, nbmax-1-length(x)) )  
  }
  donorRec_pairs
}

plot.clades <- function(presence.clades, threshold=0.6, listCladesUniques, listCladesTrees){
  cols<-setNames(palette(rainbow(18))[1:length(unique(region))],sort(unique(region))) #choose the colors for each region
  print("Saving the plot of the first stochastic mapped tree containing each clade in the folder clades_analysis")
  setwd('./clades_analysis')
  
  for (cladenb in (1:length(presence.clades[,1]))){
    if (length(presence.clades[cladenb,][!is.na(presence.clades[cladenb,])] ) < threshold*nsim) next()
    treenb <- 1
    j <- presence.clades[[cladenb,treenb]]
    while(is.na(j)){
      treenb <- treenb + 1
      j <- presence.clades[[cladenb,treenb]]
    }
    pdf(paste("clade ", listCladesUniques[cladenb], ".pdf", sep=''))
    simtree <- trees[[treenb]]
    listClades <- listCladesTrees[[treenb]]
    clade <- extract.clade.simmap2(simtree, node = getParent(tre.tt, strtoi(listClades[[j,3]])) )
    
    plotSimmap(clade,offset=0.5, fsize=min(0.4, 40/clade$Nnode), lwd=0.5, colors=cols)
    #nodelabels(text=extract.clade(tree, node = nodelab.to.numb(listClades[[i,3]], tree))$node.label, cex=min(0.5, 40/clade$Nnode), frame='none')
    #nodelabels(pie=to.matrix(getStates(clade, type='nodes'), sort(unique(region))), cex=0.2, piecol=cols)
    tiplabels(pie=to.matrix(region[clade$tip.label], sort(unique(region))), piecol=cols, cex=0.1)
    #axis(1,at=(0:3)*0.5, labels=c(2016,2016.5,2017,2017.5))
    add.simmap.legend(colors = cols, prompt=F, x=0.9*par()$usr[1], y=0.9*par()$usr[4] )
    dev.off()
  }
  setwd(paste("~/Documents/Phylogeny/world_",pattern, sep=''))
  
}

tips.keep <- function(HACnb, mode=c('booth', 'tips', 'nodes')){
  nodes.keep = c()
  tips.keep=c()
  nodes.rm <- which(nodes.chg%in% getDescendants(tre.tt, nodes.chg[HACnb]))
  if(mode=='booth'|mode=='tips') tips.keep <- HACs[[HACnb]]$tip.label [!HACs[[HACnb]]$tip.label %in% unique(unlist(sapply(HACs[nodes.rm], function(c) c$tip.label)))]
  if(mode=='booth'|mode=='nodes') nodes.keep <- HACs[[HACnb]]$node.label [!HACs[[HACnb]]$node.label %in% unique(unlist(sapply(HACs[nodes.rm], function(c) c$node.label)))]
  if (length(tips.keep)>100 | length(nodes.keep)>100) return(c(tips.keep, nodes.keep))
}

## Analysis of the fitness distribution
test.date.inRange <- function(start, end, date, tol){
  if(is.na(start)) return(FALSE)
  if({start>date-tol/2 & start<date+tol/2} | {end>date-tol/2 & end<date+tol/2} | {start<date-tol/2 & end>date+tol/2}) return(TRUE)
  FALSE
}

fitness.distr <- function(VaxStrain, reg, ageRoot, tol=0.125){
  # give a list of the fitnesses of the donor region around the same period as ageRoot, in the same vaccine strain
  edges.in <- apply(dt.edges[,reg,], 1, function(e) test.date.inRange(e[1],e[2],ageRoot, tol) )
  edges.in <- edges.in[rownames(meta_tree)]; edges.in[length(tre.tt$tip.label)+1] <- FALSE #reorder in the order of meta_tree, set the tree root to false
  nodes <- rownames(meta_tree)[grepl(VaxStrain,meta_tree$VaxStrain) & edges.in]
  Fitness[nodes]
}

time.slices <- function(HACs.tips, tree, file){
  wd <- toString(getwd())
  setwd(file)
  nodes <- 1
  for (nodesHAC in HACs.tips){
    for (t in seq(from=2013, to=2018, length.out=(3*5)+1)){
      print(t)
      pdf(paste(nodes,"_",format(round(t, 1), nsmall = 1), ".pdf", sep=''))
      c = 1
      fitness <- fitness.distr(tree, reg = "North_America", ageRoot = t)
      fitness <- fitness[names(fitness)%in%nodesHAC]
      if(length(fitness)>1){
        min <- min(c(fitness, -5)) -5
        max <- max(c(fitness, 0)) +5
        main=paste("Empirical CDF around ",format(date_decimal(t), format="%b %Y"), sep="")
        plot.ecdf(fitness, xlim=c(min,max), col=palette()[1], main=main, xlab="Fitness (measured as the cummulative mutational effect)", pch=0, verticals=TRUE, cex=0.4)
        legend("topleft", c("North_America", "SE_Asia", "Oceania", "Europe", "Eastern_Asia"), fill=c(palette()[1], palette()[2], palette()[3], palette()[4], palette()[5]), bty='n')
        
        for (reg in c("North_America", "SE_Asia", "Oceania", "Europe", "Eastern_Asia")){
          col = palette()[c]
          c <- c+1
          par(new=TRUE)
          fitness <- fitness.distr(tree, reg = reg, ageRoot = t)
          fitness <- fitness[names(fitness)%in%nodesHAC]
          if(length(fitness>0)) plot.ecdf(fitness, xlim=c(min,max), col=col, main='', xlab="", pch=0, verticals=TRUE)
        }
      }
      dev.off()
    }
    nodes <- nodes+1
  }
  dev.off()
  setwd(wd) #go back to the actual working directory
}

likelihood.obs <- function(Donor, Rec, LPM, listClades){
  clades <- LPM[toString(Donor), toString(Rec), ][!is.na(LPM[toString(Donor), toString(Rec), ])]
  L = 1
  for (c in clades){
    fitnessRoot <- Fitness[c]
    fitnessDon <- fitness.distr(meta_tree[c,]$VaxStrain, Donor, meta_tree[c,]$Decimal_Date)
    if(length(fitnessDon)<10) next()
    L = L * (1-ecdf(fitnessDon)(fitnessRoot))
  }
  log(L)
}

likelihood.sim <- function(Donor, Rec, LPM, nsim){
  clades <- LPM[toString(Donor), toString(Rec), ][!is.na(LPM[toString(Donor), toString(Rec), ])]
  Lmatrix = matrix(data=0, nrow=length(clades), ncol=nsim)
  rownames(Lmatrix)=clades
  for (c in clades){
    fitnessDon <- fitness.distr(meta_tree[c,]$VaxStrain, Donor, meta_tree[c,]$Decimal_Date)
    if(length(fitnessDon)<10) next()
    R <- sapply(1:nsim, function(i) sample.int(length(fitnessDon),1))
    e <- ecdf(fitnessDon)
    Lmatrix[c,] <- log(1-e(fitnessDon[R]))
  }
  apply(Lmatrix, 2, function(x) sum(x))
}

proba.obs <- function(Donor, Rec, LPM){
  if(PM[Donor,Rec]>=1){
    L <- likelihood.obs(Donor,Rec,LPM, listClades)
    sim <- likelihood.sim(Donor,Rec, LPM, 100)
    #hist(sim, breaks=50, main=c(Donor, Rec))
    e <- ecdf(sim)
    return(e(L))
  }
  return(NA)
}

define_HAC <- function(meta_tree){
  meta_tree$VaxStrain <- rep(NA, length(Alignment_AA))
  
  for(i in 1:length(Alignment_AA)){
    s <- as.vector(Alignment_AA[[i]])
    if(s[175]=='F') {meta_tree$VaxStrain[i] <- 'TX12'
    } else if(s[175]=='S') {meta_tree$VaxStrain[i] <- 'SW13'
    } else if(s[175]=='Y' & s[137]=='K' & s[187]=='K') {meta_tree$VaxStrain[i] <- 'SG16'
    } else if(s[137]=='N' & s[147]=='K') {meta_tree$VaxStrain[i] <- 'SW17'
    } else if(s[175]=='Y')meta_tree$VaxStrain[i] <- 'HK14'
  }
  meta_tree
  #listClades$HACnb <- sapply(listClades$Node, function(n) meta_tree$VaxStrain[strtoi(n)])
}

define_HAC2 <- function(meta_tree){
  meta_tree$VaxStrain <- rep(NA, length(Alignment_AA))
  NodesSG16 <- c(getDescendants(tre.tt, 15071), 15071)
  NodesSW17 <- c(getDescendants(tre.tt, 12879), 12879)
  NodesHK14 <- c(getDescendants(tre.tt, 10828), 10828)
  NodesHK14 <- NodesHK14[! NodesHK14 %in% c(NodesSG16, NodesSW17)]
  NodesSW13 <- c(getDescendants(tre.tt, 9953),9953)
  NodesTX12 <- (1:length(Alignment_AA))[! (1:length(Alignment_AA)) %in% c(NodesSG16, NodesSW17, NodesHK14, NodesSW13)]
  
  for (i in NodesSG16) meta_tree$VaxStrain[i] <- 'SG16'
  for (i in NodesSW17) meta_tree$VaxStrain[i] <- 'SW17'
  for (i in NodesHK14) meta_tree$VaxStrain[i] <- 'HK14'
  for (i in NodesSW13) meta_tree$VaxStrain[i] <- 'SW13'
  for (i in NodesTX12) meta_tree$VaxStrain[i] <- 'TX12'
  meta_tree
}

fitness_migr <- function(Donor, Rec, LPM){
  if(PM[Donor,Rec]<1) return(NA)
  clades <- LPM[toString(Donor), toString(Rec), ][!is.na(LPM[toString(Donor), toString(Rec), ])]
  fit <- c()
  t <- c()
  for (c in clades){
    tips <- grep('EPI',membersClades[[c]], value=TRUE)
    fit <- c(fit,Fitness[tips]-Fitness[c])
    t <- c(t,meta_tree[tips,]$Decimal_Date - meta_tree[c,]$Decimal_Date)
  }
  list(fitness=fit, time=t)
}

fitness_random <- function(nsim=100){
  #first, list the nodes with at least 15 children
  size.nodes <- function(node, tree, size=node){
    #like getDescendants, but stops at 15
    children <- tree$edge[which(tree$edge[,1]==node),2]
    w <- which(children>=length(tree$tip.label))
    size <- c(size,children)
    if (length(w>0) & length(size)<15) for (i in w) {
      size <- size.nodes(children[i], tree, size)
    }
    return(size)
  }
  
  nodesNb <- (length(tre.tt$tip.label)+1):(length(tre.tt$tip.label)+tre.tt$Nnode+1)
  Lnodes <- sapply(nodesNb, function(n) length(size.nodes(n,tre.tt)) )
  nodesOfI <- Lnodes[Lnodes>14]
  
  nodes <- runif(nsim)
  chain_length <- rnorm(nsim,40,10)
}


fitness_evol <- function(Donor, Rec, LPM, plot=FALSE){
  if(PM[Donor,Rec]<1) return(NA) #if no migration event for this pair of migration is identified
  clades <- LPM[toString(Donor), toString(Rec), ][!is.na(LPM[toString(Donor), toString(Rec), ])] #list the clades corresponding to this pair of migration
  slope_m <- c(); r2_m <- c(); CI95_m <- c() #initialise the output values
  slope_s <- c(); r2_s <- c(); CI95_s <- c()
  GS <- c(getStates(tre.sm, type='tips'), getStates(tre.sm, type='nodes'))
  
  for (c in clades){
    #compute the evolution in the migrating clade
    node_P <- getParent(tre.tt,nodelab.to.numb(c))
    tips_m <- grep('EPI',membersClades[[c]], value=TRUE)
    if(length(tips_m)<2) next()
    fit_m <- Fitness[tips_m]-Fitness[node_P]
    t_m <- meta_tree[tips_m,]$Decimal_Date - meta_tree[node_P,]$Decimal_Date
    reg_m <- lm(fit_m~t_m+0)
    if(plot){
      plot(t,fit,main=paste('migration from',Donor,'to',Rec,c), sub=paste('r^2:',round(summary(reg)$r.squared, digits=2),'95%CI:', unlist(round(confint(reg), digits=2))) )
      abline(reg)
    }
    slope_m <- c(slope_m, reg_m$coefficients)
    r2_m <- c(r2_m, summary(reg_m)$r.squared)
    CI95_m <- c(CI95_m, confint(reg_m) )
  }
  if(length(slope_m)==0) return(NA)
  list(slope_m=slope_m, r2_m=r2_m, CI95_m=CI95_m)
}

VS_slope <- function(VS, plot=FALSE){
  # #find the strain root to set it to 0
  # nodes <- meta_tree[meta_tree$VaxStrain==VS,]#meta data rows concerning the vaccine strain nodes
  # rootD <- min(nodes$Decimal_Date) #the oldest node within this vaccine strain
  # root <- rownames(nodes[nodes$Decimal_Date==rootD,]) #the corresponding node name
  if(VS=='SG16') {root <- 15071
  } else if(VS=='SW17') {root <- 12879
  } else if(VS=='HK14') {root <- 10828
  } else if(VS=='SW13') { root <- 9953
  } else if(VS=='TX12') { root <- length(tre.tt$tip.label)+1
  } else cat(paste('error in VS-slope:', VS, 'unknown'))
  
  #keep only the tips and recalibrate them given the root
  tips <- meta_tree[c(meta_tree[1:length(tre.tt$tip.label),]$VaxStrain==VS, rep(FALSE, length(tre.tt$node.label))),] #keep only the tips
  t <- tips$Decimal_Date - meta_tree$Decimal_Date[root]
  fit <- tips$Fitness - meta_tree$Fitness[root]
  
  #make the linear regression
  LM <- lm(fit~t+0)
  if (plot) { 
    newx = seq(min(t),max(t),by=0.1)
    conf_interval <- predict(LM, interval="confidence", level=0.99)
    plot(t,fit, main=VS )
    abline(LM, col='lightblue')
    matlines(t, conf_interval[,2:3], col = "blue", lty=2)
  }
  
  list(Coef=LM$coefficients, CI2.5=confint(LM)[1], CI97.5=confint(LM)[2])
}

fitness_evol_vs_VS <- function(Donor, Rec, LPM, VS_slopes, plot=FALSE){
  # Compare the slope of fitness evolution within the migrant clade vs the global 
  # evolution in the same vaccine strain
  
  if(PM[Donor,Rec]<1) return(NA) #if no migration event for this pair of migration is identified
  clades <- LPM[toString(Donor), toString(Rec), ][!is.na(LPM[toString(Donor), toString(Rec), ])] #list the clades corresponding to this pair of migration
  slope <- c(); CI95 <- c() #initialise the output values
  GS <- c(getStates(tre.sm, type='tips'), getStates(tre.sm, type='nodes'))
  
  for (c in clades){
    #compute the evolution in the migrating clade
    node_P <- getParent(tre.tt,nodelab.to.numb(c))
    tips_m <- grep('EPI',membersClades[[c]], value=TRUE)
    fit_m <- Fitness[tips_m]-Fitness[node_P]
    t_m <- meta_tree[tips_m,]$Decimal_Date - meta_tree[node_P,]$Decimal_Date
    if(length(tips_m)<5) next()
    reg_m <- lm(fit_m~t_m+0)
    if (plot) {
      plot(t_m, fit_m, main=paste('fitness evolution during migration from', Donor,'to', Rec), sub=paste('95%CI:', unlist(round(confint(reg_m), digits=2))[1:2]) )
      abline(reg_m)
    }
    
    VS <- meta_tree[c,]$VaxStrain
    slope <- c(slope, reg_m$coefficients - VS_slopes[1,VS][[1]])
    CI95 <- c(CI95, confint(reg_m) - unlist(VS_slopes[2:3,VS]) )
  }
  if(is.null(slope)) return(NA)
  list(slope=slope, CI95=CI95)
}

edges_same_time <- function(clade){
  # return the list of edges (identified by the child name) observed at the same 
  #time as the node of the clade of interest
  date <- meta_tree[clade,'Decimal_Date']
  nodes.bool <- apply(dt.edges,c(1,2), function(e) {
    if(!anyNA(e)){
      if(e[1]<date & e[2]>date) return(TRUE)}
  })
  nodes <- rownames(dt.edges)[apply(nodes.bool,1,
                                    function(r) {
                                      if(!is.null(unlist(r))) return(TRUE)
                                      FALSE
                                    }          
  )]
  nodes[!nodes==clade]
  #keep only the edges between nodes, and give the parent of the edge
  unique(unlist(sapply(nodes, function(n){
    if(grepl('NODE',n)) getParent(tre.tt, nodelab.to.numb(n))
  })))
}

Descendants_timeSlice <- function(node, SliceEnd, tree=tre.tt, elts=node){
  children <- tree$edge[which(tree$edge[,1]==node),2]
  w <- children[which(meta_tree$Decimal_Date[children]<SliceEnd)]
  elts <- c(elts, w)
  if (length(w>0)){
    for (i in w) {
      elts <- Descendants_timeSlice(i, SliceEnd, tree, elts)
    }
  }
  return(elts)
}


fitness_evol_vs_TimeSlice <- function(Donor, Rec, LPM, plot=FALSE){
  # Compare the slope of fitness evolution within the migrant clade vs the global 
  # evolution at the same time period (time slice)
  
  if(PM[Donor,Rec]<1) return(NA) #if no migration event for this pair of migration is identified
  clades <- LPM[toString(Donor), toString(Rec), ][!is.na(LPM[toString(Donor), toString(Rec), ])] #list the clades corresponding to this pair of migration
  slope <- c(); CI95 <- c() #initialise the output values
  
  for (c in clades){
    #compute the evolution in the migrating clade
    node_P <- getParent(tre.tt,nodelab.to.numb(c))
    tips_m <- grep('EPI',membersClades[[c]], value=TRUE)
    fit_m <- Fitness[tips_m]-Fitness[node_P]
    t_m <- meta_tree[tips_m,]$Decimal_Date - meta_tree[node_P,]$Decimal_Date
    if(length(tips_m)<5) next()
    reg_m <- lm(fit_m~t_m+0)
    if (plot) {
      plot(t_m, fit_m, main=paste('fitness evolution during migration from', Donor,'to', Rec))
      abline(reg_m)
    }
    
    #compute the global evolution in the same time slice
    end=max(t_m) + meta_tree[node_P,]$Decimal_Date
    edges <- edges_same_time(c)
    if(length(edges)<2) return(NA)
    tips_c <- sapply(edges, function(e){
      tips <- Descendants_timeSlice(e, end)
      tips <- tips[tips<=length(tre.tt$tip.label)] #keep only the tips
      if(length(tips)>5) return(c(e,tips)) #add the root of the clade for subsequent analysis
    })
    fit_c <- c()
    t_c <- c()
    for(e in tips_c){
      if(!is.null(e)){
        fit_c <- c(fit_c, Fitness[e[-1]] - Fitness[e[1]])
        t_c <- c(t_c, meta_tree$Decimal_Date[e[-1]] -  meta_tree$Decimal_Date[e[1]])
      }
    }
    if(is.null(t_c)) return(NA)
    reg_c <- lm(fit_c~t_c+0)
    if (plot) {
      plot(t_c, fit_c, main='fitness evolution in the sale time slice')
      abline(reg_c)
    }
    
    slope <- c(slope, reg_m$coefficients - reg_c$coefficients)
    CI95 <- c(CI95, confint(reg_m) - confint(reg_c) )
  }
  if(is.null(slope)) return(NA)
  list(slope=slope, CI95=CI95)
}

library(phytools); library(ape); library(Biostrings); library(seqinr);
library(lubridate); library(magrittr); library(dplyr); library(stringr);

tre.tt <- ape::read.nexus("timetree.nexus")
print(tre.tt)
meta = read.csv('data_world_13-19_sdubsamp.csv')
region <- meta$Region
names(region) <- meta$Isolate_Id

dst<- tre.tt
dst$edge.length[dst$edge.length==0]<-max(nodeHeights(tre.tt))*1e-6 # set zero-length branches to be 1/1000000 total tree length
nsim=25

start=Sys.time()
trees.sm <- make.simmap(dst, region, model='ARD', nsim=nsim, Q='mcmc', burnin=10000, samplefreq=200, pi='estimated')
end=Sys.time()
print(end-start)

write.simmap(trees.sm, 'trees.sm25')
print('wrote the 25 stochastic mappings as trees.sm25')