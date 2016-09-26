## script to analyze different bistro comparisons
## Claudia September 2016

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
