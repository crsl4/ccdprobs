## r script to checl the nan problem in C++ code
## In Code/tree.C, when using dog, wolf and coyote, we get
## a positive diag in inverse hessian (should be negative)
## and Newton-Raphson does not converge
## Claudia May 2016

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')
library(ggplot2)
library(weights)
library(mvtnorm)

who = "cats"
dat=read.dna("catsdogs-nan.phy")
r=c(2.815,51.982,1.903,1.275,65.402,1.000)
p=c(0.2590,0.2379,0.1900,0.3131);
den = r[6]
r = r/den
Q = makeQ(r,p,4, rescale=TRUE)

seq1 = as.vector(unname(as.character(dat[1,])))
seq2 = as.vector(unname(as.character(dat[2,])))
seq3 = as.vector(unname(as.character(dat[3,])))

## remove missing:
s1 <-seq1!="-"
s2 <- seq2!="-"
s3 <- seq3!="-"
seq1 <- seq1[s1&s2&s3]
seq2 <- seq2[s1&s2&s3]
seq3 <- seq3[s1&s2&s3]

nsites = length(seq1)
if(length(seq2) != nsites && length(seq3) != nsites){
    stop("error in number of sites seq1,seq2,seq3")
}

seq1.dist = seqMatrix(seq1)
seq2.dist = seqMatrix(seq2)
seq3.dist = seqMatrix(seq3)

out12 = countsMatrix(seq1,seq2)
out13 = countsMatrix(seq1,seq3)
out23 = countsMatrix(seq2,seq3)

## starting point
eta = 1.0
jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)
t.lik12 = simulateBranchLength.norm(nsim=1, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
t.lik13 = simulateBranchLength.norm(nsim=1, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
t.lik23 = simulateBranchLength.norm(nsim=1, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)
d12 = t.lik12$t
d13 = t.lik13$t
d23 = t.lik23$t
d1x.t0 = (d12+d13-d23)/2
d2x.t0 = (d12+d23-d13)/2
d3x.t0 = (d13+d23-d12)/2
print(c(d1x.t0,d2x.t0,d3x.t0))
if(d1x.t0 < 0)
    d1x.t0 = 0.0001
if(d2x.t0 < 0)
    d2x.t0 = 0.0001
if(d3x.t0 < 0)
    d3x.t0 = 0.0001
print(c(d1x.t0,d2x.t0,d3x.t0))


tmle = findMLE3D(seq1.dist, seq2.dist, seq3.dist, Q, t0=c(d1x.t0,d2x.t0,d3x.t0), verbose=TRUE)
print(tmle$t)
print(tmle$obsInfo)

## CONCLUSION: not the same problem as in C++
