## R script to simulate data with simData() in genetrees.R
## for ntaxa=3,4,...,12 to run bl
## Claudia June 2016

source("genetrees.R")
source("../../scripts/branch-length_lik.R")
seeds = c(1244, 41208, 9514, 118, 2367, 135, 28, 474, 818, 1916)
seed = seeds[1]
ntax = 4
nsites = 1500
alpha = 1000

## need to put the tree from cats dogs, and then take subsets for all number of taxa
set.seed(seed)
tree.text = "((1:0.11,2:0.078):0.03,3:0.098,4:0.091);"
tree = read.tree(text=tree.text)
Q = randomQ(4,rescale=TRUE)

write(paste("ntaxa =",ntax, "----------------"),file="../Examples/Simulations/settings.txt",append=TRUE)
write(paste("tree =",tree.text),file="../Examples/Simulations/settings.txt",append=TRUE)
write(paste("rates =",Q$r),file="../Examples/Simulations/settings.txt",append=TRUE)
write(paste("probs =",Q$r),file="../Examples/Simulations/settings.txt",append=TRUE)

x = simData(tree,nsites,Q,alpha)
