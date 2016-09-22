## script with analysis with different parsimony scales
## Claudia September 2016

source("../../Scripts/readBistro.r")

bistro = readBistro("randQ01")
bistro = readBistro("randQ03")
bistro = readBistro("randQ05")
bistro = readBistro("randQ07")
bistro = readBistro("randQ09")

setwd("../cats-dogs")

plotBistro(bistro)
