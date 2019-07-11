library(phytools);

tre.tt <- ape::read.nexus("timetree.nexus")
cat('loaded the tree, with ', length(tre.tt$tip.label), ' tips and ', tre.tt$Nnode,' nodes\n')
meta = read.csv('data_12-19_subsamp.csv')
cat('loaded the metadata file')
region <- meta$Region
names(region) <- meta$Isolate_Id

dst<- tre.tt
dst$edge.length[dst$edge.length==0]<-max(nodeHeights(tre.tt))*1e-6 # set zero-length branches to be 1/1000000 total tree length

fittedQ.tmp <- read.csv(file='fittedQ_large.csv')
fittedQ <- data.matrix(fittedQ.tmp[,2:9])
rownames(fittedQ)<-colnames(fittedQ)
cat('loaded fittedQ', fittedQ,'\n')

library(parallel)
trees.sm<-mclapply(1:25,function(n,tree,x,fixedQ) make.simmap(tree,x,Q=fixedQ,nsim=2, pi='estimated'),
                   tree=dst,x=region,fixedQ=fittedQ,mc.cores=if(.Platform$OS.type=="windows") 1L else 25L)

tmp <- c()
for(t in trees.sm) tmp <- c(tmp, t)
trees.sm <- tmp
rm(tmp)

class(trees.sm)<-c("multiSimmap","multiPhylo")

write.simmap(trees.sm, 'trees.sm50')
cat('wrote the 50 stochastic mappings as trees.sm50')