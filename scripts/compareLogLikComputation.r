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



## ----------------------------------------------------------------
## comparing sample of 1000 between R and C++
## still differences in both mvnormal and gammabeta
## ----------------------------------------------------------------

## order: d2x,dxy,d3y,d4y,d1x
data = read.table("../BranchLengths/Code/test/logw646_gam_mle.txt", sep=",", header=TRUE)
par3D = read.table("../BranchLengths/Code/test/par3D_646_gam_mle.txt", sep=",", header=TRUE)
par2D = read.table("../BranchLengths/Code/test/par2D_646_gam_mle.txt", sep=",", header=TRUE)

data = read.table("../BranchLengths/Code/test/logw646_norm_mle.txt", sep=",", header=TRUE)
par3D = read.table("../BranchLengths/Code/test/par3D_646_norm_mle.txt", sep=",", header=TRUE)
par2D = read.table("../BranchLengths/Code/test/par2D_646_norm_mle.txt", sep=",", header=TRUE)

data = read.table("../BranchLengths/Code/test/logw646_gam_goodR.txt", sep=",", header=TRUE)
par3D = read.table("../BranchLengths/Code/test/par3D_646_gam_goodR.txt", sep=",", header=TRUE)
par2D = read.table("../BranchLengths/Code/test/par2D_646_gam_goodR.txt", sep=",", header=TRUE)

data = read.table("../BranchLengths/Code/test/logw646_norm_goodR.txt", sep=",", header=TRUE)
par3D = read.table("../BranchLengths/Code/test/par3D_646_norm_goodR.txt", sep=",", header=TRUE)
par2D = read.table("../BranchLengths/Code/test/par2D_646_norm_goodR.txt", sep=",", header=TRUE)

data = read.table("../BranchLengths/Code/test/logw646_norm_full.txt", sep=",", header=TRUE)
par3D = read.table("../BranchLengths/Code/test/par3D_646_norm_full.txt", sep=",", header=TRUE)
par2D = read.table("../BranchLengths/Code/test/par2D_646_norm_full.txt", sep=",", header=TRUE)

data = read.table("../BranchLengths/Code/test/logw646_norm_nonrand.txt", sep=",", header=TRUE)
par3D = read.table("../BranchLengths/Code/test/par3D_646_norm_nonrand.txt", sep=",", header=TRUE)
par2D = read.table("../BranchLengths/Code/test/par2D_646_norm_nonrand.txt", sep=",", header=TRUE)

data = read.table("../BranchLengths/Code/test/logw646_gam_full.txt", sep=",", header=TRUE)
par3D = read.table("../BranchLengths/Code/test/par3D_646_gam_full.txt", sep=",", header=TRUE)
par2D = read.table("../BranchLengths/Code/test/par2D_646_gam_full.txt", sep=",", header=TRUE)

data = read.table("../BranchLengths/Code/test/logw646_norm_tol.txt", sep=",", header=TRUE)
par3D = read.table("../BranchLengths/Code/test/par3D_646_norm_tol.txt", sep=",", header=TRUE)
par2D = read.table("../BranchLengths/Code/test/par2D_646_norm_tol.txt", sep=",", header=TRUE)

data = read.table("../BranchLengths/Code/test/logw.txt", sep=",", header=TRUE)

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

save(data,file="logw646_norm_full.Rda")
save(data,file="logw646_gam_full.Rda")

plot(data$loglik, data$loglikR)


load(file="simDataR_646.Rda")
load(file="simulations_gammabeta.Rda")
head(data)
my.logw.joint = data$logwv.joint - mean(data$logwv.joint)
data$w.joint = exp(my.logw.joint)/sum(exp(my.logw.joint))
data[data$w.joint>0.01,]
length(data[data$w.joint>0.01,]$w.joint)
hist(data$w.joint)
plot(1:length(data$w.joint),cumsum(rev(sort(data$w.joint))))
my.logw.cond = data$logwv.cond - mean(data$logwv.cond)
data$w.cond = exp(my.logw.cond)/sum(exp(my.logw.cond))
data[data$w.cond>0.01,]
length(data[data$w.cond>0.01,]$w.cond)
hist(data$w.cond)
plot(1:length(data$w.cond),cumsum(rev(sort(data$w.cond))))

(1/sum(data$w.joint^2))/nreps
(1/sum(data$w.cond^2))/nreps



## -------------------------------------------------
## comparing analysis all in R vs all in C++
## ------------------------------------------------
load(file="simDataR_646.Rda")
load(file="simulations_gammabeta.Rda")
dataR = data
load(file="logw646_norm_full.Rda")
load(file="logw646_gam_full.Rda")
dataC = data

hist(dataC$bl2) ## wider range
hist(dataR$dxy.cond)
summary(dataC$bl2)
summary(dataR$dxy.cond)

## -------------------------------------------
## comparing mean from gammas to mu
## -------------------------------------------

head(par3D)
par3D <- within(par3D, mugam1 <- a1/b1)
summary(par3D$mu1)
summary(par3D$mugam1)
par3D <- within(par3D, mugam2 <- a2/b2)
summary(par3D$mu2)
summary(par3D$mugam2)
par3D <- within(par3D, mugam3 <- a3/b3)
summary(par3D$mu3)
summary(par3D$mugam3)


## sample of 1 to study step by step ---------------------------------------------------------------
## compare with BL/Code/test/screen646_norm_1.txt
## different gradient and hessian computation
data = read.table("../BranchLengths/Code/test/logw646_norm_1.txt", sep=",", header=TRUE)
par3D = read.table("../BranchLengths/Code/test/par3D_646_norm_1.txt", sep=",", header=TRUE)
par2D = read.table("../BranchLengths/Code/test/par2D_646_norm_1.txt", sep=",", header=TRUE)

## We want to compare the Newton-Raphson (from screen646_norm_1.txt)
## Entering mleDistanceJoint
## lx,ly,lz: 0.0897102, 0.0262534, 0.125166
## sxy,sxz,syz: 0.115964, 0, 0
## Entering mleDistance2D
## Starting point in mleDistance2D
## mleDistance2D Newton-Raphson curr: 0.0897102  0.125166
## mleDistance2D Newton-Raphson gradient:  -2229.4 -496.485
## mleDistance2D Newton-Raphson inverse hessian:
## -1.44679e-05 -4.78424e-05
## -4.78424e-05 -0.000710913
## entered while
## mleDistance2D Newton-Raphson curr: 0.0897102  0.125166
## mleDistance2D Newton-Raphson gradient:  -2229.4 -496.485
## mleDistance2D Newton-Raphson inverse hessian:
## -1.44679e-05 -4.78424e-05
## -4.78424e-05 -0.000710913
## First delta: 0.0560078  0.459618
## New proposed: 0.0337024 -0.334451
## with sum: 0.115964
## ...
## Finally converged
## Gradient
## 8.88207e-06 0.000468723
## 2D mean: 0.0466979 0.0373409
## , cov matrix:
##  1.89557e-05 -1.40234e-06
## -1.40234e-06  4.61438e-05

## -------------------2D ------------------------------------
seqx.dist = sequenceDist(0.0782513, 0.112059,seq3.dist, seq4.dist, Q)
mu = findMLE2D(seq2.dist, seqx.dist, seq1.dist,Q, 0.115964, t0=c(0.08971,0.125166), verbose=TRUE)
## [1] "entering findMLE..."
## [1] 0.089710 0.125166
## [1] "gradient and hessian"
## [1] -157.90620  -78.70095
##            [,1]      [,2]
## [1,] -25201.942  1199.208
## [2,]   1199.208 -6550.252

hessian = matrix(c(-25201.942,  1199.208,1199.208, -6550.252), ncol=2)
solve(hessian)
##               [,1]          [,2]
## [1,] -4.002819e-05 -7.328287e-06
## [2,] -7.328287e-06 -1.540075e-04

mu$t ## 0.0826959, 0.11346942
-1*solve(mu$obsInfo)
##              [,1]         [,2]
## [1,] 4.352266e-05 2.821790e-06
## [2,] 2.821790e-06 1.249793e-04

d2 = simulateBranchLength.conditionalMultinorm(nsim=1,seq2.dist, seqx.dist, seq1.dist,Q,t0=c(0.08971, 0.125166),0.115964)
d2$t[1] ## 0.06916
d2$t[2] ## 0.11255
0.115964 - d2$t[1] ## 0.0468

## --------------------------- 3D ----------------------------------------------
## Entering mleDistance3D
## mleDistance3D Newton-Raphson curr: 0.0793232  0.113462  0.120932
## mleDistance3D Newton-Raphson gradient:  -45.688 -84.1692 -73.1214
## mleDistance3D Newton-Raphson inverse hessian:
## -9.85803e-05  1.18323e-05  9.08453e-06
##  1.18323e-05 -0.000148516  1.70883e-06
##  9.08453e-06  1.70883e-06 -0.000156222

mu = findMLE3D(seq3.dist, seq4.dist, seq2.dist,Q, t0=c(0.0793232,0.113462,0.120932), verbose=TRUE)
## [1] "entering findMLE..."
## [1] 0.0793232 0.1134620 0.1209320
## [1] -115.7679 -128.1574 -107.1479
##            [,1]       [,2]       [,3]
## [1,] -9206.0668 -1250.8563  -961.3468
## [2,] -1250.8563 -6653.9773  -414.5576
## [3,]  -961.3468  -414.5576 -6397.8609
