### Fitness analysis
mismatches <- function(query, ref){
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
nameParent <- function(children){
  parent <- getParent(tre.tt, nodelab.to.numb(children, tre.tt))
  tre.tt$node.label[parent - length(tre.tt$tip.label)]
}

renameNodes <- function(tree){
  lab <- c(tree$tip.label,'NODE_0000000', rep(NA, tree$Nnode-1))
  while(anyNA(lab)){
    for(r in (1:length(tree$edge[,1]))){
      if(!is.na(lab[tree$edge[r,2]]) & is.na(lab[tree$edge[r,1]])) {
        lab[tree$edge[r,1]] <- nameParent(lab[tree$edge[r,2]])
      }
    }
  }
  node.label <- sort(lab[-(1:length(tree$tip.label))])
  tree$node.label <- node.label
  tree
}

timeNodes <- function(tree){
  # List the infered time of each node and tips of a given tree
  tmrca <- max(sts[tree$tip.label])-max(nodeHeights(tree))
  dates <- nodeHeights(tree)[,2] + tmrca
  dates <- c(dates, tmrca) # add the age of the root node
  name.edge <- c(tree$tip.label, tree$node.label)[c(tree$edge[,2], length(tree$tip.label)+1)] #names in the order of tree$edge + the root node
  names(dates) <- name.edge
  dates <- dates[ c(tree$tip.label, tree$node.label)]
  
  dates
}

extract.clade.simmap2 <- function(tree, node){
  #corrects extract.clade.simmap errors
  clade <- extract.clade.simmap(tree, node=node)
  clade2 <- extract.clade(phy=tree, node=node)
  clade$node.label <- clade2$node.label
  clade$Nnode <- clade2$Nnode
  clade
}

date.edge.reg <- function(reg, tree){
  #gives the dates of each edges of the tree spent within one region
  dates.nodes <- timeNodes(tree)
  edge.reg <- matrix(nrow=length(tree$maps), ncol=2)
  
  for (i in 1:length(tree$maps) ){
    for (map in 1:length(tree$maps[[i]]) ){
      
      if (names(tree$maps[[i]][map]) == reg ){
        if (map==1) {edge.reg[i,1] <- dates.nodes [ tree$edge[i,1] ]
        } else edge.reg[i,1] <- dates.nodes [ tree$edge[i,1] ] + sum(sapply(1:(map-1), function(r) tree$maps[[i]][[r]]))
        edge.reg[i,2] <- edge.reg[i,1] + tree$maps[[i]][[map]]
      }
    }
  }
  rownames(edge.reg) <- c(tree$tip.label, tree$node.label)[tree$edge[,2]]
  edge.reg
}

dates.edges <- function(tree=tre.ms, meta=meta_tree){
  DR <- array(dim = c(nrow(tree$edge),length(unique(region)),2))
  colnames(DR) <- levels(region)
  rownames(DR) <- sapply(tree$edge[,2], function(e) nodenumb.to.lab(e))
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

edge.in.interval <- function(){
  if (dates.reg[i,1]<start & dates.reg[i,2]>stop) {
    select_reg <- c(select_reg, names(dates.reg[i,1]))
  } else if (dates.reg[i,1]>start & dates.reg[i,1]<stop){
    select_reg <- c(select_reg, names(dates.reg[i,1]))
  } else if (dates.reg[i,2]>start & dates.reg[i,2]<stop){
    select_reg <- c(select_reg, names(dates.reg[i,1]))
  }
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

nodelab.to.numb <- function(nodelab, tree=tre.ms){
  if(grepl('NODE',nodelab)) return(length(tree$tip.label) + which(tree$node.label==nodelab))
  if(grepl('ISL', nodelab)) return(which(tree$tip.label==nodelab))
  print(paste('error in nodelab.to.numb', nodelab, 'is not known', sep=' '))
}

nodenumb.to.lab <- function(nodenumb, tree=tre.ms){
  if(grepl('NODE',nodenumb) | grepl('EPI',nodenumb)) return(nodenumb) #if it is already a nodelab
  if(strtoi(nodenumb)<=length(tree$tip.label)) return(tree$tip.label[nodenumb])
  tree$node.label[strtoi(nodenumb)-length(tree$tip.label)]
}

### Identification of the clades of interest
nodes.sameReg <- function(reg, node, GS, tree, elts=node){
  #count the number of nodes and tips continuously in the same region
  children <- tree$edge[which(tree$edge[,1]==node),2]
  w <- names(which(GS[children]==reg))
  elts <- c(elts, w)
  if (length(w>0)){
    for (i in w) {
      elts <- c(elts, nodes.sameReg(reg, i, GS, tree, elts))
    }
  }
  unique(elts)
}

getMembers <- function(listClades){
  GS <-c(getStates(tre.ms, type='tips'), getStates(tre.ms, type='nodes')) #region of each tip/node in the order : tips and then nodes
  M <- apply(listClades, 1, function(c) sapply(nodes.sameReg(reg=c[2], node=c[3], GS, tre.ms), function(n) {print(n);nodenumb.to.lab(n)}) )
  names(M) <- sapply(listClades$Node, function(n) nodenumb.to.lab(n))
  M
}

length.clades <- function(tree, GS){
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
  LC <- LC[LC>20] #keep clades with at least 20 continuous edges on the same region
  L <- names(LC)
  listRegParents <- unlist(lapply(L, function(node) GS[getParent(tree, node)] ))
  listRegReciep <- unlist(lapply(L, function(node) GS[node] ))
  
  result <- data.frame(Parent_Reg = listRegParents, Reciep_Reg = listRegReciep, Node = L, length = LC)
  result
}

test.Clade.Unique <- function(list.unique, i, j){
  if (list.unique[i]==list.unique[j]) return(j)
}

inventory.clades <- function(listCladesTrees, listCladesUniques){
  print('listing the presence/absence of the nodes of interest on each simulated tree ...')
  presence.clade <- matrix(ncol=nsim, nrow=length(listCladesUniques))
  regRec <- matrix(ncol=nsim, nrow=length(listCladesUniques))
  regOri <- matrix(ncol=nsim, nrow=length(listCladesUniques))
  rownames(presence.clade) <- listCladesUniques
  
  row = 0
  for (clade.unique in listCladesUniques){
    row <- row+1
    col <- 0
    for (listCladesSimtree in listCladesTrees){
      col <- col+1
      n <- 0
      for (node in listCladesSimtree$Node){
        n <- n+1
        if (clade.unique == node){
          presence.clade[row, col] <- n
          regRec[row, col] <- listCladesSimtree$Reciep_Reg[n]
          regOri[row, col] <- listCladesSimtree$Parent_Reg[n]
          break()
        }
      }
    }
    countPerRegRec <- sapply(unique(region), function(reg) sum( regRec[row, ] == reg ) )
    names(countPerRegRec) <- unique(region)
    regRecMaj <- names(which(countPerRegRec == max(countPerRegRec)))[1]
    presence.clade[row, ][!grepl(regRecMaj, regRec[row, ])] <- NA
    
    countPerRegOri <- sapply(unique(region), function(reg) sum( regOri[row, ] == reg ) )
    names(countPerRegOri) <- unique(region)
    regOriMaj <- names(which(countPerRegRec == max(countPerRegRec)))[1]
    presence.clade[row, ][!grepl(regOriMaj, regOri[row, ])] <- NA
  }
  presence.clade
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
    clade <- extract.clade.simmap2(simtree, node = getParent(tre.ms, strtoi(listClades[[j,3]])) )
    
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
  nodes.rm <- which(nodes.chg%in% getDescendants(tre.ms, nodes.chg[HACnb]))
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
  nodes <- meta_tree$Isolate_Id[grepl(VaxStrain,meta_tree$VaxStrain[-(length(tre.ms$tip.label)+1)]) & edges.in]
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
    cladenb <- which(listClades$Node==nodelab.to.numb(c))
    fitnessClade <- fitness.distr(listClades[cladenb,5], Donor, dates.nodes[c])
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
    cladenb <- which(listClades$Node==nodelab.to.numb(c))
    fitnessClade <- fitness.distr(listClades[cladenb,5], Donor, dates.nodes[c])
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
  if(PM[Donor,Rec]>1){
    L <- likelihood.obs(Donor,Rec,LPM, listClades)
    sim <- likelihood.sim(Donor,Rec, LPM, 1000)
    hist(sim, breaks=50, main=c(Donor, Rec))
    e <- ecdf(sim)
    return(e(L))
  }
  return(NA)
}

define_HAC <- function(){
  meta_tree$VaxStrain <- rep(NA, length(Alignment_AA))
  
  for(i in 1:length(Alignment_AA)){
    s <- as.vector(Alignment_AA[[i]])
    if(s[175]=='F') {meta_tree$VaxStrain[i] <- 'TX12'
    } else if(s[175]=='S') {meta_tree$VaxStrain[i] <- 'SW13'
    } else if(s[175]=='Y' & s[137]=='K' & s[187]=='K') {meta_tree$VaxStrain[i] <- 'SG16'
    } else if(s[137]=='N' & s[147]=='K') {meta_tree$VaxStrain[i] <- 'SW17'
    } else if(s[175]=='Y')meta_tree$VaxStrain[i] <- 'HK14'
  }
  
  listClades$HACnb <- sapply(listClades$Node, function(n) meta_tree$VaxStrain[strtoi(n)])
}  