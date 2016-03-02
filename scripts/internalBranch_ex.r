# example using functions in internalBranch.r
# Claudia March 2016

library(ape)
source("internalBranch.r")
library(ggplot2)

who = "birds"
d=read.dna("../datasets/birds4-clean.phy") #needs to be 4 taxa
seq1 = as.vector(unname(as.character(d[1,])))
seq2 = as.vector(unname(as.character(d[2,])))
seq3 = as.vector(unname(as.character(d[3,])))
seq4 = as.vector(unname(as.character(d[4,])))

## remove missing:
s1 <-seq1!="-"
s2 <- seq2!="-"
s3 <- seq3!="-"
s4 <- seq4!="-"
seq1 <- seq1[s1&s2&s3&s4]
seq2 <- seq2[s1&s2&s3&s4]
seq3 <- seq3[s1&s2&s3&s4]
seq4 <- seq4[s1&s2&s3&s4]

nsites= length(seq3)
d1x = 0.11
d2x = 0.078
d3y = 0.091
d4y = 0.098
dxy = 0.026

## need to run findMLE and see if it is close to 0.026, tweak the number of it
## make sure all intermediate functions do what they are supposed to do
## later need to worry about the variance: how to calculate E()?
