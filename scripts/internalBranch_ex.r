# example using functions in internalBranch.r
# Claudia March 2016

library(ape)
source("internalBranch.r")

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


t = findMLE(seq1,seq2,seq3,seq4,d1x,d2x,d3y,d4y)
## Newton-Raphson: does not work, negative t

t = findMLE2(seq1,seq2,seq3,seq4,d1x,d2x,d3y,d4y)
save(t,file="findMLE2.Rda")
plot(t$y,t$f)
abline(h=0)
t$f[t$f>0]
t$f[t$f<0]
t$f[968] ## 15.33
t$f[969] ## -6.5
t$y[968] ## 0.967
t$y[969] ## 0.968
## true value should be between 0.967, 0.968

exp(-4*dxy/3) ## 0.9659, so above good approx!
yest = t$y[968]
dxy_est=log(yest)*(-3/4)
dxy_est ## 0.02516
dxy ## 0.026

## later need to worry about the variance: how to calculate E()?
