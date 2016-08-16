## R script to read the *.treeBL file and summarize branch lengths

library(ape)
library(phangorn)
trees = read.tree("../Examples/Artiodactyl/pars1k.treeBL")
length(trees)

tree1 = read.tree(text="(2,(((6,5),4),3),1);")
tree2 = read.tree(text="(2,((5,(6,4)),3),1);")

isTree1 = function(x){
    RF.dist(x,tree1)
}
isTree2 = function(x){
    RF.dist(x,tree2)
}

trees1 = trees[sapply(trees, isTree1) == 0]
length(trees1) # 365
trees2 = trees[sapply(trees, isTree2) == 0]
length(trees2) # 197

trees1[[1]]$edge
