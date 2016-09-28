## script to analyze different bistro comparisons
## Claudia September 2016

## ========================================================
## compare randQ, fixedQ, parsimony, lik, no weight
source("../../Scripts/readBistro.r")

data = computeEntropy("randQpars")
data = computeEntropy("randQlik")
data = computeEntropy("randQno-wt")
data = computeEntropy("fixedQpars")
data = computeEntropy("fixedQlik")

plotProb(data)

## for cats
data$mb = NA
data$mb[data$tree == "(1,(2,(((((7,8),9),10),12),11)),(3,((4,5),6)));"] = 0.774196
data$mb[data$tree == "(1,2,((3,((4,5),6)),(((((7,8),9),10),12),11)));"] = 0.150627
p2 <- ggplot(aes(x=tree,y=prob), data=data) +
    geom_point() + ggtitle("black-weightProb, red-bootstrap, blue-parsimony wt, green-lik wt") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_point(aes(y=count, col="red")) + guides(col=FALSE) +
    geom_point(aes(y=parsimonyWt, col="blue")) +
    geom_point(aes(y=loglikWt, col="green")) +
    geom_point(aes(y=mb), shape=8)
plot(p2)




bistro = readBistro("randQpars")
bistro = readBistro("randQlik")
bistro = readBistro("randQno-wt")
bistro = readBistro("fixedQpars")
bistro = readBistro("fixedQlik")

plotBistro(bistro)

library(ape)
tre=read.tree(text="(1,(2,(((((7,8),9),10),12),11)),(3,((4,5),6)));")
tre2 = read.tree(text="(1,2,((3,((4,5),6)),(((((7,8),9),10),12),11)));")
layout(matrix(c(1,2),1,2))
plot(tre)
plot(tre2)

## ========================================================
## compare the entropy and plots for artiodactyl, cats and
## simulated data with 6 and 2 taxa

source("../../Scripts/readBistro.r")
setwd("../Artiodactyl")
setwd("../cats-dogs")
setwd("../Simulations")

data = computeEntropy("randQ05")
data = computeEntropy("randQ20")
data = computeEntropy("randQ05-6")
data = computeEntropy("randQ20-6")
data = computeEntropy("randQ05-12")
data = computeEntropy("randQ20-12")

plotProb(data)

bistro = readBistro("randQ05")
bistro = readBistro("randQ20")
bistro = readBistro("randQ05-6")
bistro = readBistro("randQ20-6")
bistro = readBistro("randQ05-12")
bistro = readBistro("randQ20-12")

plotBistro(bistro)

## ========================================================
## studying the effect of different parsimony scores on
## artiodactyl and cats
## Conclusion: no apparent effect

source("../../Scripts/readBistro.r")

bistro = readBistro("randQ01")
bistro = readBistro("randQ03")
bistro = readBistro("randQ05")
bistro = readBistro("randQ07")
bistro = readBistro("randQ09")

setwd("../cats-dogs")

plotBistro(bistro)
