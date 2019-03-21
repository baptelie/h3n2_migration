library(ape)

trees.test <- read.nexus("nt_world_15-18_subsamp2.trees")

(trees <- trees.test[48:57])
rm(trees.test)

nsim=length(trees)
setwd("./trees_analyzed")

lapply(1:nsim, function(t) write.tree(trees[[t]], file = paste("treeBeast_world_15-18_subsamp_", t,".nwk", sep='')))