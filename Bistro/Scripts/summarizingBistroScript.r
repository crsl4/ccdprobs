## script to analyze different bistro comparisons
## Claudia September 2016

## ========================================================
## compare the entropy and plots for artiodactyl, cats and
## simulated data with 6 and 2 taxa

source("../../Scripts/readBistro.r")
library(ggplot2)
data = computeEntropy("randQ")
plotProb(data)


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
