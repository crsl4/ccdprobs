## try regression method in Rherzky and Nei 1992
## Claudia April 2016

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')
library(ggplot2)
library(weights)
library(mvtnorm)

## ------------------
## Case (1,2)---(3,4)
## Simulating (d1x,d2x,dxy,d3y,d4y) jointly
## WEIGHTS: joint >> nj >> cond
## actually, nj not bad at all
## BRANCHES: nj are better centered around true values, but recall that
## the true likelihood was slightly shifted for 1500 sites
## and we want to match the likelihood. cond are very off.
## COV: cond cov matrix is very off, the variances in the diagonal
## should be similar, but they are not, and it has many zeros
who="(1,2)---(3,4)"
## d1x0=0.11
## d2x0=0.078
## dxy0 = 0.03
## d3y0 = 0.091
## d4y0 = 0.098
d1x0=0.1
d2x0=0.1
dxy0 = 0.1
d3y0 = 0.1
d4y0 = 0.1
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

## gamma density
out12 = countsMatrix(seq1,seq2)
out13 = countsMatrix(seq1,seq3)
out14 = countsMatrix(seq1,seq4)
out23 = countsMatrix(seq2,seq3)
out24 = countsMatrix(seq2,seq4)
out34 = countsMatrix(seq3,seq4)

nreps = 1000

logwv.joint = rep(0,nreps)
logl.joint = rep(0,nreps)
logdens.joint = rep(0,nreps)
d1x.joint=rep(0,nreps)
d2x.joint=rep(0,nreps)
d3y.joint=rep(0,nreps)
d4y.joint=rep(0,nreps)
dxy.joint=rep(0,nreps)
mat.joint <- c()
mean.joint <- c()

logwv.rn = rep(0,nreps)
logl.rn = rep(0,nreps)
logdens.rn = rep(0,nreps)
d1x.rn=rep(0,nreps)
d2x.rn=rep(0,nreps)
d3y.rn=rep(0,nreps)
d4y.rn=rep(0,nreps)
dxy.rn=rep(0,nreps)
mat.rn <- c()
mean.rn <- c()

for(nr in 1:nreps){
    print(nr)
    ## ----------------- joint ---------------------------------------------
    ## starting point BioNJ:
    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc14 = simulateBranchLength.jc(nsim=1,out14,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)
    jc24 = simulateBranchLength.jc(nsim=1,out24,eta=eta)
    jc34 = simulateBranchLength.jc(nsim=1,out34,eta=eta)
    t.lik12 = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
    t.lik13 = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
    t.lik14 = simulateBranchLength.lik(nsim=1, seq1.dist,seq4.dist,Q,t0=jc14$t,eta=eta)
    t.lik23 = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)
    t.lik24 = simulateBranchLength.lik(nsim=1, seq2.dist,seq4.dist,Q,t0=jc24$t,eta=eta)
    t.lik34 = simulateBranchLength.lik(nsim=1, seq3.dist,seq4.dist,Q,t0=jc34$t,eta=eta)
    d12 = t.lik12$t
    d13 = t.lik13$t
    d14 = t.lik14$t
    d23 = t.lik23$t
    d24 = t.lik24$t
    d34 = t.lik34$t
    v12 = d12
    v13 = d13
    v14 = d14
    v23 = d23
    v24 = d24
    v34 = d34
    s1 = d12+d13+d14
    s2 = d12+d23+d24
    d1x.t0 = (d12+1/2*(s1-s2))/2
    d2x.t0 = (d12+1/2*(s2-s1))/2
    lambda = 0.5 + ((v23-v13)+(v24-v14))/(2*2*v12)
    d3x = (lambda*d13+(1-lambda)*d23-lambda*d1x.t0-(1-lambda)*d2x.t0)
    d4x = (lambda*d14+(1-lambda)*d24-lambda*d1x.t0-(1-lambda)*d2x.t0)
    d3y.t0 = (d34+d3x-d4x)/2
    d4y.t0 = (d34+d4x-d3x)/2
    dxy.t0 = (d3x+d4x-d34)/2
    print(c(d1x.t0, d2x.t0, dxy.t0, d3y.t0, d4y.t0))

    d = simulateBranchLength.multinorm5D(nsim=1,seq1.dist, seq2.dist, seq3.dist,seq4.dist,Q,t0=c(d1x.t0, d2x.t0, d3y.t0,d4y.t0, dxy.t0), verbose=FALSE)
    d1x.joint[nr] = d$t[1]
    d2x.joint[nr] = d$t[2]
    d3y.joint[nr] = d$t[3]
    d4y.joint[nr] = d$t[4]
    dxy.joint[nr] = d$t[5]
    mat.joint[[nr]] = d$sigma
    mean.joint[[nr]] = d$mu

    if(d1x.joint[nr]<0 || d2x.joint[nr]<0 || d3y.joint[nr]<0 || dxy.joint[nr]<0 || d4y.joint[nr]<0){
        print("negative bl")
    } else{
        print("all positive")
        ##print(paste(d1x[nr],d2x[nr], dxy[nr], d3y[nr], d4y[nr]))
        logl.joint[nr] = gtr.log.lik.all(d1x.joint[nr],d2x.joint[nr],dxy.joint[nr],d3y.joint[nr],d4y.joint[nr],seq1.dist, seq2.dist, seq3.dist, seq4.dist, Q)
        logdens.joint[nr] = log(dmvnorm(d$t,mean=d$mu, sigma=d$sigma))
        logprior = logPriorExpDist(d1x.joint[nr], d2x.joint[nr], d3y.joint[nr], d4y.joint[nr], dxy.joint[nr], m=0.1)
        logwv.joint[nr] = logprior + logl.joint[nr] - logdens.joint[nr]
    }

    ## ----------------------- RN regression ------------------------------------------
    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc14 = simulateBranchLength.jc(nsim=1,out14,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)
    jc24 = simulateBranchLength.jc(nsim=1,out24,eta=eta)
    jc34 = simulateBranchLength.jc(nsim=1,out34,eta=eta)
    t.lik12 = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
    t.lik13 = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
    t.lik14 = simulateBranchLength.lik(nsim=1, seq1.dist,seq4.dist,Q,t0=jc14$t,eta=eta)
    t.lik23 = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)
    t.lik24 = simulateBranchLength.lik(nsim=1, seq2.dist,seq4.dist,Q,t0=jc24$t,eta=eta)
    t.lik34 = simulateBranchLength.lik(nsim=1, seq3.dist,seq4.dist,Q,t0=jc34$t,eta=eta)
    d12 = t.lik12$t
    d13 = t.lik13$t
    d14 = t.lik14$t
    d23 = t.lik23$t
    d24 = t.lik24$t
    d34 = t.lik34$t

    A = matrix(c(1,1,1,0,0,0,
        1,0,0,1,1,0,
        0,1,0,1,0,1,
        0,0,1,0,1,1,
        0,1,1,1,1,0), ncol=5)

    AtAinv = solve(t(A) %*% A)
    mean.rn[[nr]] <- AtAinv %*% t(A) %*% c(d12,d13,d14,d23,d24,d34)
    sigmahat <- sum((c(d12,d13,d14,d23,d24,d34) - A %*% mean.rn[[nr]])^2)
    mat.rn[[nr]] <- sigmahat * AtAinv

    d = rmvnorm(1, mean =mean.rn[[nr]], sigma=mat.rn[[nr]])
    d1x.rn[nr] = d[1]
    d2x.rn[nr] = d[2]
    d3y.rn[nr] = d[3]
    d4y.rn[nr] = d[4]
    dxy.rn[nr] = d[5]


    if(d1x.rn[nr]<0 || d2x.rn[nr]<0 || d3y.rn[nr]<0 || dxy.rn[nr]<0 || d4y.rn[nr]<0){
        print("negative bl")
    } else{
        print("all positive")
        logl.rn[nr] = gtr.log.lik.all(d1x.rn[nr],d2x.rn[nr],dxy.rn[nr],d3y.rn[nr],d4y.rn[nr],seq1.dist, seq2.dist, seq3.dist, seq4.dist, Q)
        logdens.rn[nr] = log(dmvnorm(d,mean=mean.rn[[nr]], sigma=mat.rn[[nr]]))
        logprior = logPriorExpDist(d1x.rn[nr], d2x.rn[nr], d3y.rn[nr], d4y.rn[nr], dxy.rn[nr], m=0.1)
        logwv.rn[nr] = logprior + logl.rn[nr] - logdens.rn[nr]
    }
}

data = data.frame(d1x.joint,d2x.joint,d3y.joint,d4y.joint,dxy.joint,logwv.joint, logl.joint, logdens.joint,
    d1x.rn,d2x.rn,d3y.rn,d4y.rn,dxy.rn,logwv.rn, logl.rn, logdens.rn)
head(data)
summary(data)
my.logw.joint = data$logwv.joint - mean(data$logwv.joint)
data$w.joint = exp(my.logw.joint)/sum(exp(my.logw.joint))
data[data$w.joint>0.01,]
length(data[data$w.joint>0.01,]$w.joint)
hist(data$w.joint)
plot(1:length(data$w.joint),cumsum(rev(sort(data$w.joint))))
data[data$logwv.rn==0,]
length(data[data$logwv.rn==0,]$logwv.rn)
data <- subset(data,logwv.rn!=0)
my.logw.rn = data$logwv.rn - mean(data$logwv.rn)
data$w.rn = exp(my.logw.rn)/sum(exp(my.logw.rn))
data[data$w.rn>0.01,]
length(data[data$w.rn>0.01,]$w.rn)
hist(data$w.rn)
plot(1:length(data$w.rn),cumsum(rev(sort(data$w.rn))))

## effective sample size:
(1/sum(data$w.joint^2))/nreps
(1/sum(data$w.rn^2))/nreps

save(data,file="data_RN.Rda")
save(mat.joint,mat.rn,mean.joint, mean.rn, file="mat_RN.Rda")

hist(d1x.joint)
abline(v=d1x0,col="red")
hist(d1x.rn)
abline(v=d1x0,col="red")

hist(dxy.joint)
abline(v=dxy0,col="red")
hist(data$dxy.rn)
abline(v=dxy0,col="red")

hist(d3y.joint)
abline(v=d3y0,col="red")
hist(d3y.rn)
abline(v=d3y0,col="red")

suma1=matrix(rep(0,25),ncol=5)
suma2=matrix(rep(0,25),ncol=5)
for(i in 1:nreps){
    suma1 = suma1+mat.joint[[i]]
    suma2 = suma2+mat.rn[[i]]
}
mean.mat.joint = suma1/nreps
mean.mat.rn = suma2/nreps

round(mean.mat.joint,5)
round(mean.mat.rn,5)
round(mean.mat.joint/max(mean.mat.joint),4)
round(mean.mat.rn/max(mean.mat.rn),4)


suma1=rep(0,5)
suma2=rep(0,5)
for(i in 1:nreps){
    suma1 = suma1+mean.joint[[i]]
    suma2 = suma2+mean.rn[[i]]
}

mean.mean.joint = suma1/nreps
mean.mean.rn = suma2/nreps

mean.mean.joint
mean.mean.rn


