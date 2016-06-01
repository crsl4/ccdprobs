## R script to compare the loglik (and others) computation
## between R and C++ for gamma/beta and mvnormal approaches
## will read the logw.txt table from c++, and compute
## the desired quantities in R
## Claudia May 2016

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')
library(ggplot2)
library(weights)
library(mvtnorm)

## the seq data:
seed = 0646
set.seed(seed)
d1x0=0.11
d2x0=0.078
dxy0 = 0.03
d3y0 = 0.091
d4y0 = 0.098
eta = 0.5
nsites=1500
nuc <- c('a','c','g','t')
Q = randomQ(4,rescale=TRUE)
r=Q$r
p=Q$p
## simulate seqx
seqx = sample(nuc,size=nsites,prob=Q$p,replace=TRUE)
## simulate seq1
P = matrixExp(Q,d1x0)
seq1 = numeric(nsites)
for ( i in 1:nsites )
    seq1[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])
## simulate seq2
P = matrixExp(Q,d2x0)
seq2 = numeric(nsites)
for ( i in 1:nsites )
    seq2[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])
## simulate seqy
P = matrixExp(Q,dxy0)
seqy = numeric(nsites)
for ( i in 1:nsites )
    seqy[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])
## simulate seq3
P = matrixExp(Q,d3y0)
seq3 = numeric(nsites)
for ( i in 1:nsites )
    seq3[i] = sample(nuc,size=1,prob=P[which(nuc==seqy[i]),])
## simulate seq4
P = matrixExp(Q,d4y0)
seq4 = numeric(nsites)
for ( i in 1:nsites )
    seq4[i] = sample(nuc,size=1,prob=P[which(nuc==seqy[i]),])
seq1.dist = seqMatrix(seq1)
seq2.dist = seqMatrix(seq2)
seq3.dist = seqMatrix(seq3)
seq4.dist = seqMatrix(seq4)

## order: d2x,dxy,d3y,d4y,d1x
data = read.table("../BranchLengths/Code/test/logw646_norm_full.txt", sep=",", header=TRUE)
par3D = read.table("../BranchLengths/Code/test/par3D_646_norm_full.txt", sep=",", header=TRUE)
par2D = read.table("../BranchLengths/Code/test/par2D_646_norm_full.txt", sep=",", header=TRUE)

data = read.table("../BranchLengths/Code/test/logw646_norm_nonrand.txt", sep=",", header=TRUE)
par3D = read.table("../BranchLengths/Code/test/par3D_646_norm_nonrand.txt", sep=",", header=TRUE)
par2D = read.table("../BranchLengths/Code/test/par2D_646_norm_nonrand.txt", sep=",", header=TRUE)

##data = read.table("../BranchLengths/Code/test/logw646_gam_full.txt", sep=",", header=TRUE)

head(data)
summary(data)

mylogw = data$logweight - mean(data$logweight)
data$w = exp(mylogw)/sum(exp(mylogw))
data[data$w>0.01,]
length(data[data$w>0.01,]$w)
hist(data$w)
plot(1:length(data$w),cumsum(rev(sort(data$w))))
(1/sum(data$w^2))/1000

nrep = nrow(data)
loglikR = rep(NA,nrep)

for(r in 1:nrep){
    print(r)
    d2x=data$bl1[r]
    dxy=data$bl2[r]
    d3y=data$bl3[r]
    d4y=data$bl4[r]
    d1x=data$bl5[r]
    loglikR[r] = gtr.log.lik.all(d1x,d2x,dxy,d3y,d4y,seq1.dist, seq2.dist, seq3.dist, seq4.dist, Q)
}
data <- within(data,loglikR <- loglikR)
data <- within(data,logwR <- loglikR+logprior-logdens)
head(data)
summary(data)

mylogw = data$logwR - mean(data$logwR)
data$wR = exp(mylogw)/sum(exp(mylogw))
data[data$wR>0.01,]
length(data[data$wR>0.01,]$wR)
hist(data$wR)
plot(1:length(data$wR),cumsum(rev(sort(data$wR))))
(1/sum(data$wR^2))/nrep

plot(1:nrep,data$loglik)
plot(1:nrep, data$loglikR)
