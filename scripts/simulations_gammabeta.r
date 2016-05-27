## same as simulations_normal.r but with the new gamma/beta
## approach
## Claudia May 2016

## same as simulations.r for 4 taxa, but
## now using normal for the bl, and keeping the
## dependence from start (see ipad formulas)
## Claudia April 2016
## modified to save loglik separate from logw

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')
library(ggplot2)
library(weights)
library(mvtnorm)

## ------------------
## Case (1,2)---(3,4)
## Simulating (d1x,d2x,dxy,d3y,d4y) jointly
## WEIGHTS AND BRANCHES: work fine for cond!!
## for d??0 = 0.1 for all
## and for same bl as birds

seed = 0646 ## only for comparison with BranchLength/Code/test/bl
set.seed(seed) ## only for comparison with BranchLength/Code/test/bl

who="(1,2)---(3,4)"
d1x0=0.11
d2x0=0.078
dxy0 = 0.03
d3y0 = 0.091
d4y0 = 0.098
## d1x0=0.1
## d2x0=0.1
## dxy0 = 0.1
## d3y0 = 0.1
## d4y0 = 0.1
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
mat.joint <- c() ## better: vector("list",nreps)
mean.joint <- c()

logwv.cond = rep(0,nreps)
logl.cond = rep(0,nreps)
logdens.cond = rep(0,nreps)
d1x.cond=rep(0,nreps)
d2x.cond=rep(0,nreps)
d3y.cond=rep(0,nreps)
d4y.cond=rep(0,nreps)
dxy.cond=rep(0,nreps)
mat.cond <- c()
mean1.cond <- c()
mean2.cond <- c()


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

    ## ----------------------- conditional ------------------------------------------
    ## starting point:
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

    d = simulateBranchLength.gamma(nsim=1,seq1.dist, seq2.dist, seq3.dist,Q,t0=c(d1x.t0, d2x.t0, d3x.t0))
    d1x.cond[nr] = d$t[1]
    d2x.cond[nr] = d$t[2]
    d3x.cond = d$t[3]
    mean1.cond[[nr]] = d$mu

    seqx.dist = sequenceDist(d1x.cond[nr], d2x.cond[nr] ,seq1.dist, seq2.dist, Q)

    jcx3 = simulateBranchLength.jc(nsim=1,out34,eta=eta)
    jcx4 = simulateBranchLength.jc(nsim=1,out34,eta=eta)
    jc34 = simulateBranchLength.jc(nsim=1,out34,eta=eta)
    t.likx3 = simulateBranchLength.norm(nsim=1, seqx.dist,seq3.dist,Q,t0=jcx3$t,eta=eta)
    t.likx4 = simulateBranchLength.norm(nsim=1, seqx.dist,seq4.dist,Q,t0=jcx4$t,eta=eta)
    t.lik34 = simulateBranchLength.norm(nsim=1, seq3.dist,seq4.dist,Q,t0=jc34$t,eta=eta)
    d3x = t.likx3$t
    d4x = t.likx4$t
    d34 = t.lik34$t
    dxy.t0 = (d3x+d4x-d34)/2
    d3y.t0 = (d34+d3x-d4x)/2
    d4y.t0 = (d34+d4x-d3x)/2

    d2 = simulateBranchLength.conditionalGamma(nsim=1,seqx.dist, seq3.dist, seq4.dist,Q,t0=c(dxy.t0, d4y.t0),d3x.cond)
    dxy.cond[nr] = d2$t[1]
    d4y.cond[nr] = d2$t[2]
    d3y.cond[nr] = d3x.cond - dxy.cond[nr]
    mat.cond[[nr]]= matrix(c(d$sigma[1,1], d$sigma[1,2], d$sigma[1,3],0,0,
                d$sigma[1,2], d$sigma[2,2], d$sigma[2,3],0,0,
                d$sigma[1,3], d$sigma[2,3], d$sigma[3,3]+d2$sigma[1,1],(-1)*d2$sigma[2,1], (-1)*d2$sigma[1,1],
                0,0,(-1)*d2$sigma[2,1], d2$sigma[2,2], d2$sigma[1,2],
                0,0,(-1)*d2$sigma[1,1], d2$sigma[2,1], d2$sigma[1,1]),
                ncol=5)
    mean2.cond[[nr]] = d2$mu

    if(d1x.cond[nr]<0 || d2x.cond[nr]<0 || d3y.cond[nr]<0 || dxy.cond[nr]<0 || d4y.cond[nr]<0){
        print("negative bl")
    } else{
        print("all positive")
        ##print(paste(d1x[nr],d2x[nr], dxy[nr], d3y[nr], d4y[nr]))
        logl.cond[nr] = gtr.log.lik.all(d1x.cond[nr],d2x.cond[nr],dxy.cond[nr],d3y.cond[nr],d4y.cond[nr],seq1.dist, seq2.dist, seq3.dist, seq4.dist, Q)
        logdens.cond[nr] = d$logdens + d2$logdens
        logprior = logPriorExpDist(d1x.cond[nr], d2x.cond[nr], d3y.cond[nr], d4y.cond[nr], dxy.cond[nr], m=0.1)
        logwv.cond[nr] = logprior + logl.cond[nr] - logdens.cond[nr]
    }
}

data = data.frame(d1x.joint,d2x.joint,d3y.joint,d4y.joint,dxy.joint,logwv.joint, logl.joint, logdens.joint,
    d1x.cond,d2x.cond,d3y.cond,d4y.cond,dxy.cond,logwv.cond, logl.cond, logdens.cond,
    d1x.nj,d2x.nj,d3y.nj,d4y.nj,dxy.nj,logwv.nj, logl.nj, logdens.nj)
save(data,file="simulations_gammabeta.Rda")

if(FALSE){
head(data)
summary(data)
data[data$logwv.joint==0,]
data[data$logwv.cond==0,]

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

## effective sample size:
(1/sum(data$w.joint^2))/nreps
(1/sum(data$w.cond^2))/nreps


mat.joint[[1]]/max(mat.joint[[1]])
mat.cond[[1]]/max(mat.cond[[1]])


hist(d1x.joint)
abline(v=d1x0,col="red")
hist(d1x.cond)
abline(v=d1x0,col="red")

hist(dxy.joint)
abline(v=dxy0,col="red")
hist(data$dxy.cond)
abline(v=dxy0,col="red")

hist(d3y.joint)
abline(v=d3y0,col="red")
hist(d3y.cond)
abline(v=d3y0,col="red")

suma1=matrix(rep(0,25),ncol=5)
suma2=matrix(rep(0,25),ncol=5)
for(i in 1:nreps){
    suma1 = suma1+mat.joint[[i]]
    suma2 = suma2+mat.cond[[i]]
}
mean.mat.joint = suma1/nreps
mean.mat.cond = suma2/nreps

round(mean.mat.joint,5)
round(mean.mat.cond,5)
round(mean.mat.joint/max(mean.mat.joint),4)
round(mean.mat.cond/max(mean.mat.cond),4)

suma1=rep(0,5)
suma2=rep(0,3)
suma3=rep(0,2)
for(i in 1:nreps){
    suma1 = suma1+mean.joint[[i]]
    suma2 = suma2+mean1.cond[[i]]
    suma3 = suma3+mean2.cond[[i]]
}

mean.mean.joint = suma1/nreps
mean.mean1.cond = suma2/nreps
mean.mean2.cond = suma3/nreps

mean.mean.joint
mean.mean1.cond
mean.mean2.cond

p1 <- ggplot(data,aes(dxy.joint,d4y.joint))+geom_point(color="blue")+geom_point(aes(dxy.cond,d4y.cond),data=data,colour="red")+
    ggtitle("Blue=5D joint, Red= conditional")
plot(p1)
p2 <- ggplot(data,aes(d1x.joint,d2x.joint))+geom_point(color="blue")+geom_point(aes(d1x.cond,d2x.cond),data=data,colour="red")+
    ggtitle("Blue=5D joint, Red= conditional")
plot(p2)
}
