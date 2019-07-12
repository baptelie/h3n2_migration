### Fitness analysis -----------------------------------------------------------
mismatches <- function(query, ref){
  #find the mutations between a reference sequence and a query sequence (two DNAStringSet)
  query <- strsplit(as.character(query), split='')[[1]]
  ref <- strsplit(as.character(ref), split='')[[1]]
  Pos <- which(ref!=query)
  mism <- matrix(c(ref[Pos], Pos, query[Pos]), nrow=length(Pos), ncol=3)
  mism
}

AA_to_col <- function (AA){
  # Given an amino acid, it gives the corresponding column on the data frame AA_Prefs of Lee et al. 2018 paper
  listAA = list("A", "C", "D",	"E",	"F",	"G",	"H",	"I",	"K",	"L",	"M",	"N",	"P",	"Q",	"R",	"S",	"T",	"V",	"W",	"Y","*")
  for (i in 1:22){
    if (i==22) stop(paste("Problem of translation, ", AA, " does not exists")) #Problem in the translation
    if (listAA[[i]]==AA) return(i+1)
  }
}

mut_effect <- function(Mut){
  # Computes the mutational effect of a given mutation, as explained in Lee et al. 2018
  prefWT <- AA_Prefs[strtoi(Mut[[2]]) , AA_to_col(Mut[[1]]) ]
  prefMut <- AA_Prefs[strtoi(Mut[[2]]) , AA_to_col(Mut[[3]]) ]
  
  if (prefWT == 0 | prefMut == 0) return(0) #assign ME=0 when a stop codon occurs, need to be fixed
  log2(prefMut/prefWT)
}

sum_ME <- function(ListMut){
  # Computes the sum of the mutational effects between two sequences, taking into account the cases when there is only silent mutations
  if (nrow(ListMut) == 0) return(0)
  ME = sum( apply(ListMut, 1, function(mut) mut_effect(mut) ) )
  ME
}

fitness_nodes <- function(tree){
  #compute the fitness of each node, as the cummulative ME from the root of the given tree
  cat('compute the fitness of each tip and node, as the cummulative ME from the root \n')
  fitness <- c(rep(NA, length(tree$tip.label)), 0, rep(NA, length(tree$node.label)-1 )) #initialize with NA for each node, except 0 for the root
  names(fitness) <- c(tree$tip.label, tree$node.label)
  for (i in (1:nrow(tree$edge)) ){
    fitness[tree$edge[i,2] ] <- fitness[[tree$edge[i,1] ]] + list_ME[[i ]]
  }
  fitness
}

### Analysis tools of the time aligned tree ----------------
time_nodes <- function(tree=tre.tt, meta=meta){
  # List the infered time of each node and tips of a given tree
  tmrca <- max(meta$Decimal_Date) - max(nodeHeights(tree)) #be careful to remove the outliers
  dates <- nodeHeights(tree)[,2] + tmrca
  dates <- c(dates, tmrca) # add the age of the root node
  name.edge <- c(tree$tip.label, tree$node.label)[c(tree$edge[,2], length(tree$tip.label)+1)] #names in the order of tree$edge + the root node
  names(dates) <- name.edge
  
  dates[ c(tree$tip.label, tree$node.label)]
}

dates_edges <- function(tree=tre.sm, meta=meta_tree){
  #creates a matrix giving for each edge the dates spent on each region
  regions <- colnames(tree$mapped.edge)
  DR <- array(dim = c(nrow(tree$edge),length(regions),2))
  colnames(DR) <- regions
  rownames(DR) <- sapply(tree$edge[,2], function(e) nodenumb_to_lab(e)) #each row=one edge, named by the child node (unique)
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

nodelab_to_numb <- function(nodelab, tree=tre.tt){
  #converts a node label into its absolute numbering (i.e order of tips + nodes)
  if(grepl('NODE',nodelab)) return(length(tree$tip.label) + which(tree$node.label==nodelab))
  if(grepl('ISL', nodelab)) return(which(tree$tip.label==nodelab))
  print(paste('error in nodelab.to.numb', nodelab, 'is not known', sep=' '))
}

nodenumb_to_lab <- function(nodenumb, tree=tre.tt){
  #converts an absolute node number (i.e order of tips + nodes) into its name as appearing on the alignment or on the tree
  if(grepl('NODE',nodenumb) | grepl('EPI',nodenumb)) return(nodenumb) #if it is already a nodelab
  if(strtoi(nodenumb)<=length(tree$tip.label)) return(tree$tip.label[strtoi(nodenumb)])
  tree$node.label[strtoi(nodenumb)-length(tree$tip.label)]
}

### Identification of the clades of interest ---------------------------
define_ag_clusters <- function(meta_tree, Alignment){
  meta_tree$VaxStrain <- rep('x', length(Alignment))
  
  for(i in 1:length(Alignment)){
    s <- as.vector(Alignment[[i]])
    if(s[175]=='F') {meta_tree$VaxStrain[i] <- 'TX12'
    } else if(s[175]=='S' & s[154]=='S') {meta_tree$VaxStrain[i] <- 'SW13'
    # } else if(s[175]=='Y' & s[137]=='K' & s[187]=='K') {meta_tree$VaxStrain[i] <- 'SG16'
    # } else if(s[137]=='N' & s[147]=='K') {meta_tree$VaxStrain[i] <- 'SW17'
    } else if(s[175]=='Y')meta_tree$VaxStrain[i] <- 'HK14'
  }
  meta_tree
}

nodes_same_reg <- function(reg, node, GS, tree, elts=node){
  #count the number of nodes and tips continuously in the same region
  children <- tree$edge[which(tree$edge[,1]==node),2]
  w <- children[which(GS[children]==reg)]
  elts <- c(elts, w)
  if (length(w>0)){
    for (i in w) {
      elts <- nodes_same_reg(reg, node=i, GS, tree, elts)
    }
  }
  return(elts)
}

get_members_all <- function(tree){
  #lists all the nodes from each clade of interest continuously in the same region
  GS <-c(getStates(tree, type='tips'), getStates(tree, type='nodes')) #region of each tip/node in the order : tips and then nodes
  
  edge.chg <- apply(tree$edge, 1, function(edge) {
    GS[edge[1]]!=GS[edge[2]] & edge[2]>length(tree$tip.label)
    }) #true/false in function of whether the edge changes of region, in the order of tree$edge
  
  node.chg <- tree$edge[,2][edge.chg] #root of each clade in the recipient region
  GSe <- GS[tree$edge[,2]] #reorder GS in function of the edge order
  chg <- data.frame(node=node.chg, reg=GSe[edge.chg]) #recipient region
  
  apply(chg, 1, function(c) nodes_same_reg(reg=c[[2]], node=c[[1]], GS, tree) )
}

get_members_succes <- function(tree, trunk.tree, meta){
  #filter for the clades staying a long enough time on the same recipient region
  GM_all <- get_members_all(tree)
  nb_tips <- length(tree$tip.label)
  lengthtips.GM_all <- sapply(GM_all, function(gm) sum(strtoi(gm)<=nb_tips) )
  GM_succ <- GM_all[ lengthtips.GM_all>=6 & !names(GM_all) %in% trunk.tree ]
  duration <- sapply( GM_succ, function(gm) max(meta$Decimal_Date[strtoi(gm)]) - min(meta$Decimal_Date[strtoi(gm)]) )
  e.length <- tree$edge.length
  names(e.length) <- tree$edge[,2]
  max.edge.length <- sapply(GM_succ, function(gm) max(e.length[gm[-1]] ) )
  GM_succ <- GM_succ[duration>=0.25 & max.edge.length<1]
  names(GM_succ) <- sapply(names(GM_succ), function(n) nodenumb_to_lab(n))
  GM_succ
}

list_clades <- function(tree, trenb){
  # summarise the migration clades of interest, given a stochastic mapped tree
  GS <- c(getStates(tree, type='tips'), getStates(tree, type='nodes'))
  L <- strtoi(sapply(Migr_per_tree[[trenb]], function(m) m[[1]]))
  listRegParents <- unlist(lapply(L, function(node) GS[getParent(tree, node)] ))
  listRegReciep <- unlist(lapply(L, function(node) GS[node] ))
  
  data.frame(Parent_Reg = listRegParents, Reciep_Reg = listRegReciep, Node = L)
}

pairs_migr <- function(tree, trenb){
  listClades<- list_clades(tree, trenb)
  regions <- colnames(tree$mapped.edge)
  donorRec_pairs <- matrix(data=0, nrow=length(regions), ncol=length(regions), dimnames=list(regions, regions))
  
  for (r in 1:nrow(listClades)) {
    donorRec_pairs[ toString(listClades[r,1]), toString(listClades[r,2]) ] <- donorRec_pairs[ toString(listClades[r,1]), toString(listClades[r,2]) ] + 1
  }
  donorRec_pairs
}

list_pairs_migr <- function(PM, tree, trenb){
  regions <- colnames(tree$mapped.edge)
  nbmax <- max(PM)
  listClades<- list_clades(tree, trenb)
  donorRec_pairs <- array(data=NA, dim=c(length(regions), length(regions), nbmax) )
  rownames(donorRec_pairs) <- colnames(donorRec_pairs) <- regions
  
  for (r in 1:nrow(listClades)) {
    x <- donorRec_pairs[ toString(listClades[r,1]), toString(listClades[r,2]), ][!is.na(donorRec_pairs[ toString(listClades[r,1]), toString(listClades[r,2]), ])]
    donorRec_pairs[ toString(listClades[r,1]), toString(listClades[r,2]), ] <- c(x, nodenumb_to_lab(toString(listClades[r,3])), rep(NA, nbmax-1-length(x)) )  
  }
  donorRec_pairs
}


trunk <- function(tree, meta=meta_tree){
  ntips <- length(tree$tip.label)
  tips_last_year <- tree$tip.label[meta$Decimal_Date[1:ntips]>(max(meta$Decimal_Date)-1)] #identify the tips sampled during the last year
  tips <- tips_last_year[sample.int(length(tips_last_year), 300)] #randomly take 300 of them
  
  getAllParents <- function(node,tree=tre.tt){
    Parents <- c()
    p <- nodelab_to_numb(node)
    while(p != length(tree$tip.label)+1){
      p <- getParent(tree,p)
      Parents <- c(Parents,p)
    }
    Parents
  }
  
  listParents <- lapply(tips, getAllParents)
  f <- table(unlist(listParents))
  names(f)[f >= 150] #keep the parent nodes shared by at least half of all randomly picked tips
}

trunk_location <- function(simtrees, treeTopo){
  trunkNodes <- trunk(treeTopo)
  locNodes <- lapply(simtrees, function(t){
    loc <- getStates(t, type='nodes')
    loc[trunkNodes]-length(treeTopo$tip.label)
  })
  datesNodes<-meta_tree$Decimal_Date[trunkNodes]
  locFreq <- matrix(0, nrow=length(unique(region)), ncol=round((max(datesNodes)-min(datesNodes))*12)+1 )
  rownames(locFreq)<- unique(region)
  for(i in 1:length(trunkNodes)){
    d <- (datesNodes[[i]]-min(datesNodes)) %/% (1/12)
    for(lt in 1:length(locNodes)){
      locFreq[locNodes[[lt]][[i]],d] <- locFreq[locNodes[[lt]][[i]],d] + 1
    }
  }
  locFreq
}

nb_per_month_reg <- function(){
  regions <- colnames(tree$mapped.edge)
  cnt <- matrix(0,nrow=nreg, ncol=(2019.25-2013.5-(2/12))*12)
  rownames(cnt)<- regions
  colnames(cnt)<- round(seq(from=2013.5+(2/12), to=2019.25-(1/12), by=1/12), digits=2)
  
  for(r in unique(region)){
    sub <- meta[(meta$Region == r),]
    col<-0
    for(m in seq(from=2013.5+(3/12), to=2019.25, by=1/12)){
      col<-col+1
      sub_month <- sub[sub$Decimal_Date <m & sub$Decimal_Date >= m-(1/12), ]
      cnt[r,col]<- nrow(sub_month)
    }
  }
  
  write.csv(cnt, 'nb_seq_per_month_reg.csv')
}

## Analysis of the fitness distribution ----------------------------------------
edges_same_time <- function(clade, Reg, GS, trnk, tol_before, tol_after, meta=meta_tree, refTree=tre.tt){
  # return the list of edges observed at the same 
  #time as the rootnode of the clade of interest, in the same vaccine clade, in a given region

  date <- meta[clade,'Decimal_Date']
  VaxCluster <- meta[clade,'VaxStrain']
  
  nodes <- apply(refTree$edge, 1, function(e){
    start <- meta$Decimal_Date[e[1]]
    stop <- meta$Decimal_Date[e[2]]
    if(! GS[e[2]] %in% Reg | meta$VaxStrain[e[2]] != VaxCluster | e[2] %in% trnk) return(NULL)
    if(start<date-tol_before & stop>date-tol_before & stop<date+tol_after) return(e[2]) #easiest case: one node falls within the time interval (ie one month before or after the date)
    
    if(start < date-tol_before & stop > date+tol_after & GS[e[1]] %in% Reg & meta$VaxStrain[e[1]] == VaxCluster & !e[1] %in% trnk  & stop-start < 0.5) {#if the edge is through the interval, and entirely in the good region & good vaccine cluster, return the closest node of the two
      return(e[1])
      # if(date-start < stop-date) return(e[1])
      # if(date-start > stop-date) return(e[2])
    }
  })
  unlist(nodes)
}

fitness_distr_donors <- function(Donor, Rec, LPM, dt.edges, meta, GS, trnk, refTree=tre.tt){
  clades <- LPM[Donor, toString(Rec), ][!is.na(LPM[Donor, toString(Rec), ])]
  F_distr <- lapply(clades, function(c){
    edges <- edges_same_time(clade = c, Reg = GS[getParent(refTree, node = nodelab_to_numb(c))], GS, trnk, tol_before = 1/12, tol_after = 0)
    meta$Fitness[edges]
    })
  names(F_distr) <- clades
  F_distr
}

likelihood_obs <- function(fit.distr.don, meta, plot=FALSE){
  # compute the likelihood of the observed migration events, given a pair of donor and recipient regions
  L =1
  for (c in names(fit.distr.don)){
    fitnessRoot <- meta[c,]$Fitness
    fitnessDon <- fit.distr.don[c][[1]]
    if(length(fitnessDon)<10) next()
    
    if(plot){
      hist(fitnessDon, breaks = round( (max(fitnessDon) - min(fitnessDon)) / 2 ) )
      norm.d <- rnorm(500, mean = mean(fitnessDon), sd = sd(fitnessDon))
      lines(density(norm.d)$x, density(norm.d)$y*50 )
    }
    # p = max(1-ecdf(fitnessDon)(fitnessRoot),0.01/length(fitnessDon))
    p = 1-pnorm(fitnessRoot, mean = mean(fitnessDon), sd = sd(fitnessDon))
    L = L * p
  }
  if(!exists('fitnessDon')) return(NA)
  log(L)
}

likelihood_sim <- function(fit.distr.don, nsim){
  # repeats nsim times a random choice of a migrant virus
  Lmatrix = matrix(data=0, nrow=length(fit.distr.don), ncol=nsim)
  rownames(Lmatrix) = names(fit.distr.don)
  for (c in names(fit.distr.don)){
    fitnessDon <- fit.distr.don[c][[1]]
    if(length(fitnessDon)<10) next()
    R <- sample.int(n = length(fitnessDon), size = nsim, replace= TRUE)
    # e <- ecdf(fitnessDon)
    # p <- 1-e(fitnessDon[R])
    # p[p<0.01/length(fitnessDon)] <- 0.01/length(fitnessDon)
    # Lmatrix[c,] <- log(p)
    p <- 1-pnorm(fitnessDon[R], mean = mean(fitnessDon), sd= sd(fitnessDon))
    Lmatrix[c,] <- log(p)
  }
  apply(Lmatrix, 2, function(x) sum(x))
}

proba_obs <- function(Rec, LPM, dt.edges, GS, trnk, meta = meta_tree){
  Trop <- c("China","SE_Asia","Southern_Asia")
  if(sum(!is.na(LPM[Trop,Rec,])) >= 1) { #if at least one migration event identified
    fit.distr.don <- lapply(Trop, function(d) fitness_distr_donors(d, Rec, LPM, dt.edges, meta = meta, GS, trnk))
    fit.distr.don <- unlist(fit.distr.don, recursive = FALSE)
    L <- likelihood_obs(fit.distr.don, meta = meta)
    sim <- likelihood_sim(fit.distr.don, nsim = 2000)
    e <- ecdf(sim)
    return(data.frame(P=e(L), bootstrap.distr=sim, lik = L))
  }
  return(list(P=NA, bootstrap.distr=NA, lik=NA))
}

summary_proba_obs <- function(Rec, simmaps, meta, Migr_per_tree, trnk, ncores){
  nsim=length(simmaps)
  distr <- mclapply(seq_along(simmaps), function(t){
    tre.sm <- simmaps[[t]]
    membersClades <- Migr_per_tree[[t]] #migration events per region pair
    dt.edges <- dates_edges(tree=tre.sm)
    PM <- pairs_migr(tre.sm, t)
    GS <- c(getStates(tre.sm, type='tips'), getStates(tre.sm, type='nodes')) #region of each tip/node in the order : tips and then nodes
    LPM <- list_pairs_migr(PM,tre.sm,t)

    proba_obs(Rec,LPM, dt.edges, GS, trnk, meta)
  }, mc.cores=if(.Platform$OS.type=="windows") 1L else ncores)

  distr.tbl <- as.data.frame(do.call(rbind, distr))
  
  ggp <- ggplot(distr.tbl) +
    # geom_density (aes(lik.distri, fill=tree),
    geom_density (aes(bootstrap.distr),
                  fill="#00AFBB",
                  alpha=0.3,
                  color=NA)
  # scale_fill_manual(values=rep("#00AFBB",nsim))
  
  ggp +
    # geom_area(aes(x = ifelse(bootstrap.distr<mean(unique(distr.tbl$lik)), bootstrap.distr, 0)),
    #           y = ..density..,
    #           fill="#00AFBB",
    #           alpha=0.6) +
    geom_segment(
      aes(x = mean(unique(distr.tbl$lik)), xend= mean(unique(distr.tbl$lik)), y=0, yend=layer_scales(ggp)$y$range$range[[2]] * 0.8, alpha=0.2),
      linetype = 'dashed',
      colour='gray40',
      size=0.5,
      alpha = 0.2) +
    # geom_histogram(aes(lik.obs, y=..density../10),
    #                binwidth=1,
    #                alpha=0.3) +
    # geom_density(
    #   aes(lik.obs, y=..density../5),
    #   color=NA,
    #   fill='gray40',
    #   adjust=2,
    #   alpha = 0.5) +
    annotate('text', x= mean(unique(distr.tbl$lik)), y = layer_scales(ggp)$y$range$range[[2]] * 0.83, label=paste('mean P:',round(mean(unique(distr.tbl$P)),2),'+/-',round(qt(0.975,df=nsim-1)*sd(unique(distr.tbl$P))/sqrt(nsim),2)), size=4, color='gray40') +
    labs(y="frequency fo the bootstrap", x='log likelihood of the migrant viruses',
         caption = paste('D. Likelihood distribution of the random migrant virus scenario \nand the actual likelihood of the observed migrant viruses from tropical Asian regions to',gsub('_',' ',Rec))) +
    theme(plot.caption = element_text(hjust=0.5, size=18, face='bold'),
          axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          legend.position = 'none')

# 
#   
#   for(t in 1:nsim){
#     tre.sm <- simmaps[[t]]
#     membersClades <- Migr_per_tree[[t]] #migration events per region pair
#     dt.edges <- dates_edges(tree=tre.sm)
#     PM <- pairs_migr(tre.sm, t)
#     GS <- c(getStates(tre.sm, type='tips'), getStates(tre.sm, type='nodes')) #region of each tip/node in the order : tips and then nodes
#     LPM <- list_pairs_migr(PM,tre.sm,t)
#     
#     x <- proba_obs(Rec,LPM, dt.edges, GS, trnk, meta)
#     distr <- rbind(distr, data.table( lik.distri = x$bootstrap.distr, tree=paste('tree',t), lik.obs=x$lik, p=x$P) )
#   }
#   P <- distr$p[!duplicated(distr$tree)]
#   ggp <- ggplot(distr) +
#     # geom_density (aes(lik.distri, fill=tree),
#     geom_density (aes(lik.distri),
#                   fill="#00AFBB",
#                   # alpha=1/nsim,
#                   color=NA)
#     # scale_fill_manual(values=rep("#00AFBB",nsim))
# 
#   ggp +
#       geom_segment(
#         aes(x = mean(distr$lik.obs), xend= mean(distr$lik.obs), y=0, yend=layer_scales(ggp)$y$range$range[[2]] * 0.8, alpha=0.2),
#         linetype = 'dashed',
#         colour='gray40',
#         size=0.5,
#         alpha = 0.2) +
#     # geom_histogram(aes(lik.obs, y=..density../10),
#     #                binwidth=1,
#     #                alpha=0.3) +
#     # geom_density(
#     #   aes(lik.obs, y=..density../5),
#     #   color=NA,
#     #   fill='gray40',
#     #   adjust=2,
#     #   alpha = 0.5) +
#     annotate('text', x=mean(distr$lik.obs), y = layer_scales(ggp)$y$range$range[[2]] * 0.83, label=paste('mean P:',round(mean(P),2),'+/-',round(qt(0.975,df=nsim-1)*sd(P)/sqrt(nsim),2)), size=4, color='gray40') +
#     labs(y="frequency fo the bootstrap", x='log likelihood of the migrant viruses',
#          caption = paste('D. Likelihood distribution of the random migrant virus scenario \nand the actual likelihood of the observed migrant viruses from tropical Asian regions to',gsub('_',' ',Rec))) +
#     theme(plot.caption = element_text(hjust=0.5, size=18, face='bold'),
#           axis.title = element_text(size=15),
#           axis.text = element_text(size=13),
#           legend.position = 'none')
}

descendants_time_slice <- function(node, SliceEnd, GS, tree=tre.tt, elts=node, meta=meta_tree){
  Trop <- c("China","SE_Asia","Southern_Asia")
  children <- tree$edge[which(tree$edge[,1]==node),2]
  w <- children[meta$Decimal_Date[children]<=SliceEnd & GS[children] %in% Trop]
  elts <- c(elts, w)
  if (length(w)>0){
    for (i in w) {
      elts <- descendants_time_slice(i, SliceEnd, GS, elts = elts)
    }
  }
  return(elts)
}

fitness_evol <- function(Rec, LPM, PM, dt.edges, membersClades, GS, trnk, meta = meta_tree, refTree = tre.tt, plot=FALSE){
  # Compare the slope of fitness evolution within the migrant clade vs the global 
  # evolution at the same time period (time slice)
  
  Trop <- c("China","SE_Asia","Southern_Asia")
  clades <- LPM[Trop,Rec,][!is.na(LPM[Trop,Rec,])] #list the clades corresponding to this pair of migration
  slope <- data.table(inf=0, equal=0, sup=0)
  
  for (c in clades){
    #compute the evolution in the migrating clade
    node_P <- getParent(refTree,nodelab_to_numb(c))
    tips_m <- strtoi(membersClades[[c]][strtoi(membersClades[[c]])<=length(refTree$tip.label)]) #get the tips in the recipient region
    fit_m <- meta$Fitness[tips_m] - meta$Fitness[node_P]
    t_m <- meta[tips_m,]$Decimal_Date - meta[node_P,]$Decimal_Date
    
    #compute the global evolution in the same time slice
    edges <- edges_same_time(node_P, Trop, GS, trnk, tol_before = 0.08, tol_after = 0.08)
    edges <- edges[edges>length(refTree$tip.label)] #keep only the internal nodes
    
    if(length(edges)<6) next()
    tips_c <- sapply(edges, function(e){
      nodes <- descendants_time_slice(node = e,SliceEnd = max(meta[tips_m,]$Decimal_Date), GS)
      e.length <- refTree$edge.length
      names(e.length) <- refTree$edge[,2]
      max.edge.length <- max(e.length[nodes] )
      tips <- nodes[nodes<=length(refTree$tip.label)] #keep only the tips
      if(length(tips)<6) return(NULL)
      duration <- max(meta$Decimal_Date[strtoi(nodes)]) - min(meta$Decimal_Date[strtoi(nodes)])
      if(duration >= 0.25 & max.edge.length < 1) return(c(e,tips)) #add the root of the clade for subsequent analysis
    })
    
    fit_c <- t_c <- tips.an <- c()
    nbCl <- 0
    for(e in tips_c){
      if(!is.null(e)){
        nbCl <- nbCl + 1
        tips.an <- c(tips.an, e[-1])
        fit_c <- c(fit_c, meta$Fitness[e[-1]] - meta$Fitness[e[1]])
        t_c <- c(t_c, meta_tree$Decimal_Date[e[-1]] -  meta_tree$Decimal_Date[e[1]])
      }
    }
    if(nbCl < 4 ) next()
    fit_c <- fit_c[!duplicated(tips.an)]
    t_c <- t_c[!duplicated(tips.an)]
    cat('C: ', fit_c,'\n')
    reg_c <- lm(fit_c~t_c+0)
    cat('M: ', fit_m,'\n')
    reg_m <- lm(fit_m~t_m+0)
    
    #ANCOVA test
    if(all(fit_c==0) & all(fit_m==0)){
      D=FALSE
    } else {
      FitEvol <- data.frame(clade = c(rep('migr',length(t_m)), rep('ctrl', length(t_c))), t=c(t_m,t_c), fit=c(fit_m, fit_c) )
      Ancova <- lm(fit~t*clade+0, data=FitEvol)
      if(summary(Ancova)$coefficients['t:clademigr','Pr(>|t|)']<0.5) D=TRUE else D=FALSE
    }
    
    if(D){
      if(reg_m$coefficients - reg_c$coefficients < 0) {slope$inf[1] <- slope$inf[1]+1
      } else if (reg_m$coefficients - reg_c$coefficients > 0) slope$sup[1] <- slope$sup[1]+1
    } else slope$equal[1] <- slope$equal[1]+1
    
    if (plot) {
      plot(t_m, fit_m, col=rgb(0.6,0,0.05),xlim=c(0,lgt), ylim=c(-15,5), pch=18,cex=1.2, sub=paste('slope statitically different ',D), xlab='time (year)', ylab='fitness')
      abline(reg_m, col=rgb(0.8,0,0.1))
      
      s <- sample.int(length(t_c),min(100, length(t_c)))
      par(new=TRUE)
      plot(t_c[s], fit_c[s], col=gray(0.3),xlim=c(0,lgt), ylim=c(-15,5), pch=15,cex=0.5, xlab='', ylab='')
      abline(reg_c, col=gray(0.7))
    }
  }
  slope/sum(slope)
}


fitness_evol2 <- function(Rec, LPM, PM, dt.edges, membersClades, GS, trnk, meta = meta_tree, refTree = tre.tt, plot=FALSE){
  # Compare the slope of fitness evolution within the migrant clade vs the global 
  # evolution at the same time period (time slice)
  
  Trop <- c("China","SE_Asia","Southern_Asia")
  clades <- LPM[Trop,Rec,][!is.na(LPM[Trop,Rec,])] #list the clades corresponding to this pair of migration
  slope <- data.table(inf=0, equal=0, sup=0)
  
  for (c in clades){
    #compute the evolution in the migrating clade
    node_P <- getParent(refTree,nodelab_to_numb(c))
    tips_m <- strtoi(membersClades[[c]][strtoi(membersClades[[c]])<=length(refTree$tip.label)]) #get the tips in the recipient region
    fit_m <- meta$Fitness[tips_m] - meta$Fitness[node_P]
    t_m <- meta[tips_m,]$Decimal_Date - meta[node_P,]$Decimal_Date
    reg_m <- lm(fit_m~t_m+0)
    lgt=max(t_m)
    
    #compute the global evolution in the same time slice
    edges <- edges_same_time(node_P, Trop, GS, trnk, tol_before = 0.125, tol_after = 0.125)
    edges <- edges[edges>length(refTree$tip.label)] #keep only the internal nodes
    
    if(length(edges)<6) next()
    compare <- sapply(edges, function(e){
      nodes <- descendants_time_slice(node = e,SliceEnd = max(meta[tips_m,]$Decimal_Date), GS)
      e.length <- refTree$edge.length
      names(e.length) <- refTree$edge[,2]
      max.edge.length <- max(e.length[nodes] )
      tips <- nodes[nodes<=length(refTree$tip.label)] #keep only the tips
      if(length(tips)<3) return(list(inf=NA, equal = NA, sup = NA))
      duration <- max(meta$Decimal_Date[strtoi(tips)]) - meta$Decimal_Date[e]
      if(duration < 0.2 | max.edge.length > 1) return(list(inf=NA, equal = NA, sup = NA))
      
      fit_c <- meta$Fitness[tips] - meta$Fitness[e]
      t_c <- meta_tree$Decimal_Date[tips] -  meta_tree$Decimal_Date[e]
      reg_c <- lm(fit_c~t_c+0)
      
      #ANCOVA test
      if(all(fit_c==0) & all(fit_m==0)) return(list(inf=0, equal = 1, sup = 0))
      
      FitEvol <- data.frame(clade = c(rep('migr',length(t_m)), rep('ctrl', length(t_c))), t=c(t_m,t_c), fit=c(fit_m, fit_c) )
      Ancova <- lm(fit~t*clade+0, data=FitEvol)
      if(summary(Ancova)$coefficients['t:clademigr','Pr(>|t|)']>0.05) return(list(inf=0, equal = 1, sup = 0))
      
      if(reg_m$coefficients - reg_c$coefficients < 0) return(list(inf=1, equal = 0, sup = 0))
      if (reg_m$coefficients - reg_c$coefficients > 0) return(list(inf=0, equal = 0, sup = 1))
    })
    
    nb_obs <- sum(unlist(compare), na.rm=TRUE)
    if(nb_obs < 6) next()
    comp.sum <- data.table(inf = unlist(compare[1,]) %>% sum(., na.rm=TRUE)/nb_obs,
                           equal = unlist(compare[2,]) %>% sum(., na.rm=TRUE)/nb_obs,
                           sup = unlist(compare[3,]) %>% sum(., na.rm=TRUE)/nb_obs)
    slope <- slope + comp.sum
    
    if (plot) {
      plot(t_m, fit_m, col=rgb(0.6,0,0.05),xlim=c(0,lgt), ylim=c(-15,5), pch=18,cex=1.2, sub=paste('slope statitically different ',D), xlab='time (year)', ylab='fitness')
      abline(reg_m, col=rgb(0.8,0,0.1))
      
      s <- sample.int(length(t_c),min(100, length(t_c)))
      par(new=TRUE)
      plot(t_c[s], fit_c[s], col=gray(0.3),xlim=c(0,lgt), ylim=c(-15,5), pch=15,cex=0.5, xlab='', ylab='')
      abline(reg_c, col=gray(0.7))
    }
  }
  slope/sum(slope)
}

fitness_evol3 <- function(Rec, LPM, PM, dt.edges, membersClades, GS, trnk, meta = meta_tree, refTree = tre.tt, plot=FALSE){
  # Compare the slope of fitness evolution within the migrant clade vs the global 
  # evolution at the same time period (time slice)
  
  Trop <- c("China","SE_Asia","Southern_Asia")
  clades <- LPM[Trop,Rec,][!is.na(LPM[Trop,Rec,])] #list the clades corresponding to this pair of migration
  # slope <- data.table(inf=0, equal=0, sup=0)
  reg_c <- reg_m <- c()
  
  for (c in clades){
    #compute the evolution in the migrating clade
    node_P <- getParent(refTree,nodelab_to_numb(c))
    tips_m <- strtoi(membersClades[[c]][strtoi(membersClades[[c]])<=length(refTree$tip.label)]) #get the tips in the recipient region
    fit_m <- meta$Fitness[tips_m] - meta$Fitness[node_P]
    t_m <- meta[tips_m,]$Decimal_Date - meta[node_P,]$Decimal_Date
    
    #compute the global evolution in the same time slice
    edges <- edges_same_time(node_P, Trop, GS, trnk, tol_before = 0.08, tol_after = 0.08)
    edges <- edges[edges>length(refTree$tip.label)] #keep only the internal nodes
    
    if(length(edges)<5) next()
    tips_c <- sapply(edges, function(e){
      nodes <- descendants_time_slice(node = e,SliceEnd = max(meta[tips_m,]$Decimal_Date), GS)
      e.length <- refTree$edge.length
      names(e.length) <- refTree$edge[,2]
      max.edge.length <- max(e.length[nodes] )
      tips <- nodes[nodes<=length(refTree$tip.label)] #keep only the tips
      if(length(tips)<3) return(NULL)
      duration <- max(meta$Decimal_Date[strtoi(nodes)]) - min(meta$Decimal_Date[strtoi(nodes)])
      if(length(tips)>5 & duration >= 0.25 & max.edge.length < 1) return(c(e,tips)) #add the root of the clade for subsequent analysis
    })
    
    fit_c <- t_c <- tips.an <- c()
    nbCl <- 0
    for(e in tips_c){
      if(!is.null(e)){
        nbCl <- nbCl + 1
        tips.an <- c(tips.an, e[-1])
        fit_c <- c(fit_c, meta$Fitness[e[-1]] - meta$Fitness[e[1]])
        t_c <- c(t_c, meta_tree$Decimal_Date[e[-1]] -  meta_tree$Decimal_Date[e[1]])
      }
    }
    if(nbCl < 4) next()
    fit_c <- fit_c[!duplicated(tips.an)]
    t_c <- t_c[!duplicated(tips.an)]
    reg_c <- c(reg_c, lm(fit_c~t_c+0)$coefficients)
    reg_m <- c(reg_m, lm(fit_m~t_m+0)$coefficients)
    
    
    # #ANCOVA test
    # if(all(fit_c==0) & all(fit_m==0)){
    #   D=FALSE
    # } else {
    #   FitEvol <- data.frame(clade = c(rep('migr',length(t_m)), rep('ctrl', length(t_c))), t=c(t_m,t_c), fit=c(fit_m, fit_c) )
    #   Ancova <- lm(fit~t*clade+0, data=FitEvol)
    #   if(summary(Ancova)$coefficients['t:clademigr','Pr(>|t|)']<0.05) D=TRUE else D=FALSE
    # }
    # 
    # if(D){
    #   if(reg_m$coefficients - reg_c$coefficients < 0) {slope$inf[1] <- slope$inf[1]+1
    #   } else if (reg_m$coefficients - reg_c$coefficients > 0) slope$sup[1] <- slope$sup[1]+1
    # } else slope$equal[1] <- slope$equal[1]+1
    # 
  }
  #t-test
  test <- t.test(x = reg_c, y = reg_m, paired = TRUE)
  
  if (plot) {
    data <- data.table(coef = c(reg_m, reg_c), cat = c(rep('migrant',length(reg_m)), rep('tropics', length(reg_c)) ) )
    
    # ggplot(data, aes(cat, coef)) +
    #   geom_boxplot(outlier.shape=NA) +
    #   geom_jitter (width = 0.2, size = 0.5) +
    #   labs(title = paste("Average proportion of the observed migration events from tropical Asian regions to",gsub('_',' ',Rec)),
    #        y='', caption = test$p.value)
    
    ggpaired(data, x = 'cat', y = 'coef',
             ylab = 'fitness evolution rate',
             color = 'cat', line.color = "gray", line.size = 0.4,
             palette = "jco")+
      stat_compare_means(paired = TRUE)
  }
}

summary_slopes <- function(Rec, simmaps, meta, ncores, trnk){
  nsim = length(simmaps)
  slope <- c()
  slope <- mclapply(seq_along(simmaps), function(t){
    tre.sm <- simmaps[[t]]
    membersClades <- Migr_per_tree[[t]] #migration events per region pair
    dt.edges <- dates_edges(tree=tre.sm)
    PM <- pairs_migr(tre.sm, t)
    GS <- c(getStates(tre.sm, type='tips'), getStates(tre.sm, type='nodes')) #region of each tip/node in the order : tips and then nodes
    LPM <- list_pairs_migr(PM,tre.sm,t)
    fitness_evol(Rec, LPM, PM, dt.edges, membersClades, GS, trnk)
  }, mc.cores=if(.Platform$OS.type=="windows") 1L else ncores)
  
  slope.tbl <- data.table(variable = rep(c('inf','equal','sup'), nsim), value = unlist(slope))
  
  ggplot(slope.tbl) +
    geom_boxplot(aes(variable, value), outlier.shape=NA) +
    geom_jitter(aes(variable, value), width = 0.2, size = 0.5) +
    
    
    
    # geom_bar(aes(variable, value), stat='summary', fun.y = 'mean') +
    labs(caption = paste("Average proportion of the observed migration events from tropical Asian regions to ",gsub('_',' ',Rec)),
         y='', x='slope difference between general tropical regions and the migration event') +
    theme(plot.caption = element_text(hjust=0.5, size=18, face='bold'),
          axis.title = element_text(size=15),
          axis.text = element_text(size=13),
          legend.position = 'none')
}

compare_slopes <- function(Rec, LPM, PM, dt.edges, membersClades, GS){
  regions <- colnames(PM)
  E <- matrix(0, nrow=nreg, ncol=nreg, dimnames=list(regions, regions))
  slope <- array(0,dim=c(nreg, nreg, 3), dimnames=list(regions, regions, c('inf','equal','sup')))
  
  for(don in Trop) {
    FE <- fitness_evol_vs_TimeSlice(don, Rec, LPM, PM, dt.edges, membersClades, GS)
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
  list(E=E,slope=slope)
}