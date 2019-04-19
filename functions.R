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
    d <- meta$Decimal_Date[rw]
    for(reg in 1:length(maps)){
      DR[rw,names(maps[reg]),1] <- d
      d <- d+maps[reg]
      DR[rw,names(maps[reg]),2] <- d
    }
  }
  DR
}

select.date <- function(start, stop, dates.reg){
  # selects the nodes and tips from a region within a given period of time
  select_reg <- c()
  for (i in 1: nrow(dates.reg)){
    if( !is.na(dates.reg[i,1]) ){
      if (dates.reg[i,1]<start & dates.reg[i,2]>stop) {
        select_reg <- c(select_reg, names(dates.reg[i,1]))
      } else if (dates.reg[i,1]>start & dates.reg[i,1]<stop){
        select_reg <- c(select_reg, names(dates.reg[i,1]))
      } else if (dates.reg[i,2]>start & dates.reg[i,2]<stop){
        select_reg <- c(select_reg, names(dates.reg[i,1]))
      }
    }
  }
  select_reg
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
test.date <- function(start, end, date, tol){
  if(is.na(start)) return(FALSE)
  if({start>date-tol/2 & start<date+tol/2} | {end>date-tol/2 & end<date+tol/2} | {start<date-tol/2 & end>date+tol/2}) return(TRUE)
  FALSE
}

fitness.distr <- function(VaxStrain, reg, ageRoot, tol=0.125){
  edges.in <- apply(dt.edges[,reg,], 1, function(e) {test.date(e[1],e[2],ageRoot, tol)})
  nodes <- rownames(meta_tree)[grepl(VaxStrain,meta_tree$VaxStrain[-(length(tre.tt$tip.label)+1)]) & edges.in]
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


P <- function(fitnessClade, fitnessObs){
  #Compute the probability of a fitness value given the distribution within the clade
  d <- stats::density(fitnessClade, from=-40, to=10, n=501)
  rk <- (round(fitnessObs, digits = 1)+40)*10+1
  d$y[rk]
}

likelihood.obs <- function(Donor, Rec, LPM, listClades){
  clades <- LPM[toString(Donor), toString(Rec), ][!is.na(LPM[toString(Donor), toString(Rec), ])]
  L = 1
  for (c in clades){
    fitnessRoot <- Fitness[c]
    fitnessClade <- fitness.distr(meta_tree[c,]$VaxStrain, Donor, meta_tree[c,]$Decimal_Date)
    if(length(fitnessClade)<10) next()
    L = L * P(fitnessClade, fitnessRoot)
  }
  log(L)
}

likelihood.sim <- function(Donor, Rec, LPM, nsim){
  clades <- LPM[toString(Donor), toString(Rec), ][!is.na(LPM[toString(Donor), toString(Rec), ])]
  Lmatrix = matrix(data=0, nrow=length(clades), ncol=nsim)
  rownames(Lmatrix)=clades
  for (c in clades){
    fitnessClade <- fitness.distr(meta_tree[c,]$VaxStrain, Donor, meta_tree[c,]$Decimal_Date)
    if(length(fitnessClade)<10) next()
    R <- runif(nsim)
    e <- ecdf(fitnessClade)
    X <- seq(min(fitnessClade), max(fitnessClade), length.out=100)
    Y <- e(X)
    fsim <- sapply(R, function(r) X[which.min((Y-r)^2)])
    Lmatrix[c,] <- sapply(fsim, function(x) log(P(fitnessClade, x)))
  }
  apply(Lmatrix, 2, function(x) sum(x))
}

proba.obs <- function(Donor, Rec, LPM){
  if(PM[Donor,Rec]>=1){
    L <- likelihood.obs(Donor,Rec,LPM, listClades)
    sim <- likelihood.sim(Donor,Rec, LPM, 1000)
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


fitness_sister <- function(Donor, Rec, LPM){
  if(PM[Donor,Rec]<1) return(NA) #if no migration event for this pair of migration is identified
  clades <- LPM[toString(Donor), toString(Rec), ][!is.na(LPM[toString(Donor), toString(Rec), ])] #list the clades corresponding to this pair of migration
  slope <- c(); r2 <- c(); CI95 <- c() #initialise the output values
  GS <- c(getStates(tre.sm, type='tips'), getStates(tre.sm, type='nodes'))
  
  for (c in clades){
    node_P <- getParent(tre.tt,nodelab.to.numb(c))
    nodes <- nodes.sameReg(Donor, node_P, GS, tre.sm)
    tips <- nodes[nodes<=length(tre.sm$tip.label)]
    if(length(tips)<2) next()
    fit <- Fitness[tips]-Fitness[node_P]
    t <- meta_tree[tips,]$Decimal_Date - meta_tree[node_P,]$Decimal_Date
    reg <- lm(fit~t+0)
    
    slope <- c(slope, reg$coefficients)
    r2 <- c(r2, summary(reg)$r.squared)
    CI95 <- c(CI95, confint(reg) )
  }
  if(is.null(slope)) return(NA)
  list(slope=slope, r2=r2, CI95=CI95)
}


fitness_evol_sis <- function(Donor, Rec, LPM, plot=FALSE){
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
    if(length(tips_m)<2) next()
    reg_m <- lm(fit_m~t_m+0)
    
    #compute the evolution in the sister clade
    nodes <- nodes.sameReg(Donor, node_P, GS, tre.sm)
    tips_s <- nodes[nodes<=length(tre.sm$tip.label)]
    if(length(tips_s)<2) next()
    fit_s <- Fitness[tips_s]-Fitness[node_P]
    t_s <- meta_tree[tips_s,]$Decimal_Date - meta_tree[node_P,]$Decimal_Date
    reg_s <- lm(fit_s~t_s+0)
    
    slope <- c(slope, reg_m$coefficients-reg_s$coefficients)
    CI95 <- c(CI95, confint(reg_m)-confint(reg_s) )
  }
  if(is.null(slope)) return(NA)
  list(slope=slope, CI95=CI95)
}

fitness_evol_global <- function(Donor, Rec, LPM, plot=FALSE){
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
    reg_m <- lm(fit_m~t_m+0)
    
    #compute the global evolution during that period of time in the same vax strain clade
    VaxStrain <- meta_tree[c,]$VaxStrain
    tips_vax_strain <- rownames(meta_tree[1:length(tre.tt$tip.label),])[meta_tree[1:length(tre.tt$tip.label),]$VaxStrain==VaxStrain]
    tips_s <- tips_vax_strain[meta_tree[tips_vax_strain,]$Decimal_Date> meta_tree[node_P,]$Decimal_Date & meta_tree[tips_vax_strain,]$Decimal_Date < max(meta_tree[tips_m,]$Decimal_Date)]
    if(length(tips_s)<3) next()
    root_VS <- rownames(meta_tree[meta_tree$VaxStrain==VaxStrain,])[which(meta_tree[meta_tree$VaxStrain==VaxStrain,]$Decimal_Date==min(meta_tree[meta_tree$VaxStrain==VaxStrain,]$Decimal_Date))]
    fit_s <- Fitness[tips_s] - Fitness[root_VS]
    t_s <- meta_tree[tips_s,]$Decimal_Date - meta_tree[root_VS,]$Decimal_Date
    reg_s <- lm(fit_s~t_s+0)
    
    slope <- c(slope, reg_m$coefficients-reg_s$coefficients)
    CI95 <- c(CI95, confint(reg_m)-confint(reg_s) )
  }
  if(is.null(slope)) return(NA)
  list(slope=slope, CI95=CI95)
}
