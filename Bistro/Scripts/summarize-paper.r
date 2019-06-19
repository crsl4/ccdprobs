## R script to summarize bistro simulations for paper
## Claudia October 2018

library(ape)
library(phangorn)

data = "024"
rootname.bistro = paste0("bistro-sim-",data,"-nsites-1500-1000-1000")
truetree = read.tree(paste0("../Data/datasets/ExaBayes_ConsensusExtendedMajorityRuleNewick.",data,".cons"))
meanTreeBistro = read.tree(paste0("../Simulations/",rootname.bistro,".meanTree"))

## Distance of mean tree to true tree
dist.true = RF.dist(truetree, meanTreeBistro)

## MAP bistro with PP
df = read.table(paste0(rootname.bistro,"trprobs"), header=TRUE, sep=" ")
mapBistroTree = read.tree(text=df$tree[1])
mapBistroProb = df$prob[1]

## Distance of MAP bistro to true tree
dist.map = RF.dist(truetree, mapBistroTree)

## PP of true tree
for(i in 1:nrow(df)){
    tree = read.tree(text=df$tree[i])
    if(RF.dist(tree,truetree))
        break
}
ppTrueTree = df$prob[i]

## Topology Bistro ESS
lines = readLines(paste0(rootname.bistro,"trprobs"))
lines = readLines("mds.topPP")
ess0 = strsplit(lines[1],"=")[[1]][2]
ess = as.numeric(gsub(")", "", ess0))

## Missing: vstat summary
