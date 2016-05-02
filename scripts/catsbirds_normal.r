## same as simulations_normal.r but with cats/birds data
## with topology uncertainty
## Claudia April 2016

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')
library(ggplot2)
library(weights)
library(mvtnorm)


who = "cats"
dat=read.dna("../datasets/4taxa-cats.phy") #needs to be 4 taxa
dat.tre=read.table("../datasets/4taxa-cats_ccdprobs.out", header=FALSE)
r=c(2.815,51.982,1.903,1.275,65.402,1.000)
p=c(0.2590,0.2379,0.1900,0.3131);
den = r[6]
r = r/den
Q = makeQ(r,p,4, rescale=TRUE)

who = "birds"
dat=read.dna("../datasets/birds4-clean.phy") #needs to be 4 taxa
dat.tre=read.table("../datasets/birds4-clean_ccdprobs.out", header=FALSE)
r = c(0.2463,0.1764,0.1231,0.0187,0.4185,0.0170)
p = c(0.2776,0.2937,0.1612,0.2675)
den = r[6]
r = r/den
Q = makeQ(r,p,4, rescale=TRUE)

nreps = 1000
trees = rep(NA,nreps)
tx1 = rep(NA,nreps)
tx2 = rep(NA,nreps)
tx3 = rep(NA,nreps)
tx4 = rep(NA,nreps)
eta = 0.5

logwv.joint = rep(0,nreps)
logl.joint = rep(0,nreps)
logdens.joint = rep(0,nreps)
d1x.joint=rep(0,nreps)
d2x.joint=rep(0,nreps)
d3y.joint=rep(0,nreps)
d4y.joint=rep(0,nreps)
dxy.joint=rep(0,nreps)
mat.joint <- vector("list",nreps)
mean.joint <- vector("list",nreps)

logwv.cond = rep(0,nreps)
logl.cond = rep(0,nreps)
logdens.cond = rep(0,nreps)
d1x.cond=rep(0,nreps)
d2x.cond=rep(0,nreps)
d3y.cond=rep(0,nreps)
d4y.cond=rep(0,nreps)
dxy.cond=rep(0,nreps)
mat.cond <- vector("list",nreps)
mean1.cond <- vector("list",nreps)
mean2.cond <- vector("list",nreps)

logwv.nj = rep(0,nreps)
logl.nj = rep(0,nreps)
logdens.nj = rep(0,nreps)
d1x.nj=rep(0,nreps)
d2x.nj=rep(0,nreps)
d3y.nj=rep(0,nreps)
d4y.nj=rep(0,nreps)
dxy.nj=rep(0,nreps)
mat.nj <- vector("list",nreps)

for(nr in 1:nreps){
    print(nr)
    t=sampleTopQuartet(dat.tre)
    tre = t$tre
    trees[nr] = write.tree(t$tre)
    leaves = tre$tip.label[tre$edge[which(tre$edge[,1]==sample(c(5,6),1) & tre$edge[,2] < 5),2]] #works only for unrooted quartet
    leaves = as.numeric(leaves)
    sis = leaves
    fth = setdiff(1:4,sis)
    print(paste("sis",sis))
    print(paste("fth",fth))
    tx1[nr] = sis[1]
    tx2[nr] = sis[2]
    tx3[nr] = fth[1]
    tx4[nr] = fth[2]
    seq1 = as.vector(unname(as.character(dat[sis[1],])))
    seq2 = as.vector(unname(as.character(dat[sis[2],])))
    seq3 = as.vector(unname(as.character(dat[fth[1],])))
    seq4 = as.vector(unname(as.character(dat[fth[2],])))
    # remove missing:
    s1 <-seq1!="-"
    s2 <- seq2!="-"
    s3 <- seq3!="-"
    s4 <- seq4!="-"
    seq1 <- seq1[s1&s2&s3&s4]
    seq2 <- seq2[s1&s2&s3&s4]
    seq3 <- seq3[s1&s2&s3&s4]
    seq4 <- seq4[s1&s2&s3&s4]

    nsites = length(seq1)
    if(length(seq2) != nsites && length(seq3) != nsites){
        stop("error in number of sites seq1,seq2,seq3")
    }

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

    d = simulateBranchLength.multinorm(nsim=1,seq1.dist, seq2.dist, seq3.dist,Q,t0=c(d1x.t0, d2x.t0, d3x.t0))
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

    d2 = simulateBranchLength.conditionalMultinorm(nsim=1,seqx.dist, seq3.dist, seq4.dist,Q,t0=c(dxy.t0, d4y.t0),d3x.cond)
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
        logdens.cond[nr] = log(dmvnorm(d$t,mean=d$mu, sigma=d$sigma))+log(dmvnorm(d2$t,mean=d2$mu,sigma=d2$sigma))
        logprior = logPriorExpDist(d1x.cond[nr], d2x.cond[nr], d3y.cond[nr], d4y.cond[nr], dxy.cond[nr], m=0.1)
        logwv.cond[nr] = logprior + logl.cond[nr] - logdens.cond[nr]
    }

    ## ------------------------------------- NJ -----------------------------------------------
    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc14 = simulateBranchLength.jc(nsim=1,out14,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)
    jc34 = simulateBranchLength.jc(nsim=1,out34,eta=eta)
    ## starting point for d3x,d4x
    t0=(max(c(jc12$t,jc13$t,jc23$t))+min(c(jc12$t,jc13$t,jc23$t)))/2

    t.lik12 = simulateBranchLength.norm(nsim=1, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
    t.lik13 = simulateBranchLength.norm(nsim=1, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
    t.lik14 = simulateBranchLength.norm(nsim=1, seq1.dist,seq4.dist,Q,t0=jc14$t,eta=eta)
    t.lik23 = simulateBranchLength.norm(nsim=1, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)
    t.lik34 = simulateBranchLength.norm(nsim=1, seq3.dist,seq4.dist,Q,t0=jc34$t,eta=eta)

    d12 = t.lik12$t
    d13 = t.lik13$t
    d14 = t.lik14$t
    d23 = t.lik23$t
    d34 = t.lik34$t

    d1x.nj[nr] = (d12+d13-d23)/2
    d2x.nj[nr] = (d12+d23-d13)/2
    dxy.nj[nr] = (d14-d34-d12+d23)/2
    d3y.nj[nr] = (d13-d14+d34)/2
    d4y.nj[nr] = (d14-d13+d34)/2
    A = matrix(c(1/2,1/2,-1/2,0,0,1/2,-1/2,0,1/2,-1/2,-1/2,1/2,1/2,0,0,0,0,1/2,-1/2,1/2,0,0,-1/2,1/2,1/2),ncol=5)
    S = diag(c(t.lik12$sigma^2, t.lik13$sigma^2, t.lik23$sigma^2, t.lik14$sigma^2, t.lik34$sigma^2))
    mat.nj[[nr]] = A %*% S %*% t(A)


    if(d1x.nj[nr]<0 || d2x.nj[nr]<0 || d3y.nj[nr]<0 || d4y.nj[nr]<0 || dxy.nj[nr]<0){
        print("negative bl")
    } else{
        print("all positive")
        ##print(paste(d1x[nr],d2x[nr], dxy[nr], d3y[nr], d4y[nr]))
        logl.nj[nr] = gtr.log.lik.all(d1x.nj[nr],d2x.nj[nr],dxy.nj[nr],d3y.nj[nr],d4y.nj[nr],seq1.dist, seq2.dist, seq3.dist, seq4.dist, Q)
        logdens.nj[nr] = logJointDensity.multinorm(d1x.nj[nr], d2x.nj[nr], dxy.nj[nr], d3y.nj[nr], d4y.nj[nr], t.lik12,t.lik13,t.lik23,t.lik14,t.lik34) #from multinorm
        logprior  = logPriorExpDist(d1x.nj[nr], d2x.nj[nr], d3y.nj[nr], d4y.nj[nr], dxy.nj[nr], m=0.1)
        logwv.nj[nr] = logprior +logl.nj[nr] - logdens.nj[nr]
    }
}

data = data.frame(trees,tx1,tx2,tx3,tx4,d1x.joint,d2x.joint,d3y.joint,d4y.joint,dxy.joint,logwv.joint, logl.joint, logdens.joint,
    d1x.cond,d2x.cond,d3y.cond,d4y.cond,dxy.cond,logwv.cond, logl.cond, logdens.cond,
    d1x.nj,d2x.nj,d3y.nj,d4y.nj,dxy.nj,logwv.nj, logl.nj, logdens.nj)
head(data)
summary(data)
data[data$logwv.joint==0,] ##0
data[data$logwv.cond==0,] ##0
data[data$logwv.nj==0,] ## too many!!
##data <- subset(data,logwv.nj!=0)

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


save(data,file=paste0("data_normal_",who,".Rda"))
save(mean.joint, mat.joint, mean1.cond, mean2.cond, mat.cond, file=paste0("meanmat_normal_",who,".Rda"))
load("data_normal_birds.Rda")

## effective sample size:
(1/sum(data$w.joint^2))/nreps
(1/sum(data$w.cond^2))/nreps

bigw = data[data$w.cond>0.01,]
subset(bigw,select=c(trees,tx1,tx2,tx3,tx4,d1x.cond, d2x.cond, d3y.cond, d4y.cond, dxy.cond, w.cond))
summary(bigw)

load("data_normal_birds.Rda")
summary(data)
data1 = data[data$trees=="((1,2),3,4);",]
data2 = data[data$trees=="((1,3),2,4);",]
data3 = data[data$trees=="(1,(2,3),4);",]

## to do: check summary statistics (weighted) and compare to MB for birds,
## see where there are shifts
## see weighted frequency for each tree also


## INTERNAL BL -----------------------------------------------------------------
## first tree (1,2)
## MB freq: 0.6819
sum(data1$w.joint) ## 0.5432
sum(data1$w.cond) ## 0.5172
w1.joint = data1$w.joint/sum(data1$w.joint)
w1.cond = data1$w.cond/sum(data1$w.cond)
## mean MB dxy: 0.0261
mm = weighted.mean(data1$dxy.joint,w1.joint) ## 0.0261
weighted.mean(data1$dxy.cond,w1.cond) ## 0.0251
## quartiles:
weighted.quantile(data1$dxy.joint,w1.joint,probs=0.25) ## 0.02164
weighted.quantile(data1$dxy.cond,w1.cond,probs=0.25) ## 0.02125
weighted.quantile(data1$dxy.joint,w1.joint,probs=0.75) ## 0.02998
weighted.quantile(data1$dxy.cond,w1.cond,probs=0.75) ## 0.02921

summary(data1$w.joint)
summary(data1$w.cond)

wtd.hist(data1$dxy.joint, weight=w1.joint)
abline(v=mm, col="red")
wtd.hist(data1$dxy.cond, weight=w1.cond)
abline(v=mm, col="red")

## second tree (1,3)
## MB freq: 0.0178
sum(data2$w.joint) ## 0.00839
sum(data2$w.cond) ## 0.0104
w2.joint = data2$w.joint/sum(data2$w.joint)
w2.cond = data2$w.cond/sum(data2$w.cond)
## mean MB dxy: 0.0213
mm = weighted.mean(data2$dxy.joint,w2.joint) ## 0.0212
weighted.mean(data2$dxy.cond,w2.cond) ## 0.01732
## quartiles:
weighted.quantile(data2$dxy.joint,w2.joint,probs=0.25) ## 0.0172
weighted.quantile(data2$dxy.cond,w2.cond,probs=0.25) ## 0.01472
weighted.quantile(data2$dxy.joint,w2.joint,probs=0.75) ## 0.02418
weighted.quantile(data2$dxy.cond,w2.cond,probs=0.75) ## 0.02066

summary(data2$w.joint)
summary(data2$w.cond)

wtd.hist(data2$dxy.joint, weight=w2.joint)
abline(v=mm, col="red")
wtd.hist(data2$dxy.cond, weight=w2.cond)
abline(v=mm, col="red")

## third tree (2,3)
## MB freq: 0.3003
sum(data3$w.joint) ## 0.4484
sum(data3$w.cond) ## 0.4723
w3.joint = data3$w.joint/sum(data3$w.joint)
w3.cond = data3$w.cond/sum(data3$w.cond)
## mean MB dxy: 0.0255
mm = weighted.mean(data3$dxy.joint,w3.joint) ## 0.0254
weighted.mean(data3$dxy.cond,w3.cond) ## 0.0257
## quartiles:
weighted.quantile(data3$dxy.joint,w3.joint,probs=0.25) ## 0.0214
weighted.quantile(data3$dxy.cond,w3.cond,probs=0.25) ## 0.02043
weighted.quantile(data3$dxy.joint,w3.joint,probs=0.75) ## 0.02937
weighted.quantile(data3$dxy.cond,w3.cond,probs=0.75) ## 0.03017

summary(data3$w.joint)
summary(data3$w.cond)

wtd.hist(data3$dxy.joint, weight=w3.joint)
abline(v=mm, col="red")
wtd.hist(data3$dxy.cond, weight=w3.cond)
abline(v=mm, col="red")


## EXTERNAL BL --------------------------------------------------------
## d1x does not mean the same across rows, need to see who is tx1
summary(data1$d1x.joint)
summary(data1$d1x.cond)

plot(data1$dxy.joint, data1$dxy.cond)

data1.1 = subset(data1,tx1 == 1)
summary(data1.1)
plot(data1.1$d1x.joint, data1.1$d1x.cond)
head(data1.1)
## mean MB d1x: 0.1108
w1.1.joint = data1.1$w.joint/sum(data1.1$w.joint)
w1.1.cond = data1.1$w.cond/sum(data1.1$w.cond)
weighted.mean(data1.1$d1x.joint,w1.1.joint) ## 0.1100
weighted.mean(data1.1$d1x.cond,w1.1.cond) ## 0.1105
## mean MB d2x: 0.0779
weighted.mean(data1.1$d2x.joint,w1.1.joint) ## 0.0778
weighted.mean(data1.1$d2x.cond,w1.1.cond) ## 0.0748
## mean MB d3y: 0.0909
weighted.mean(data1.1$d3y.joint,w1.1.joint) ## 0.0901
weighted.mean(data1.1$d3y.cond,w1.1.cond) ## 0.0883
## mean MB d4y: 0.0983
weighted.mean(data1.1$d4y.joint,w1.1.joint) ## 0.0988
weighted.mean(data1.1$d4y.cond,w1.1.cond) ## 0.0999
## mean MB dxy: 0.0261
weighted.mean(data1.1$dxy.joint,w1.1.joint) ## 0.0264
weighted.mean(data1.1$dxy.cond,w1.1.cond) ## 0.0236

data2.1 = subset(data2,tx1 == 1)
summary(data2.1)
plot(data2.1$d1x.joint, data2.1$d1x.cond)
head(data2.1)
## mean MB d1x: 0.1119
w2.1.joint = data2.1$w.joint/sum(data2.1$w.joint)
w2.1.cond = data2.1$w.cond/sum(data2.1$w.cond)
weighted.mean(data2.1$d1x.joint,w2.1.joint) ## 0.1109
weighted.mean(data2.1$d1x.cond,w2.1.cond) ## 0.1119
## mean MB d2x: 0.0909
weighted.mean(data2.1$d2x.joint,w2.1.joint) ## 0.0931
weighted.mean(data2.1$d2x.cond,w2.1.cond) ## 0.08833
## mean MB d3y: 0.0813
weighted.mean(data2.1$d3y.joint,w2.1.joint) ## 0.081
weighted.mean(data2.1$d3y.cond,w2.1.cond) ## 0.0723
## mean MB d4y: 0.103
weighted.mean(data2.1$d4y.joint,w2.1.joint) ## 0.1009
weighted.mean(data2.1$d4y.cond,w2.1.cond) ## 0.105
## mean MB dxy: 0.0213
weighted.mean(data2.1$dxy.joint,w2.1.joint) ## 0.0219
weighted.mean(data2.1$dxy.cond,w2.1.cond) ## 0.0202



data1.1 = subset(data1,tx1 == 3)
plot(data1.1$d1x.joint, data1.1$d1x.cond)
head(data1.1)
## mean MB d3y: 0.0909
w1.1.joint = data1.1$w.joint/sum(data1.1$w.joint)
w1.1.cond = data1.1$w.cond/sum(data1.1$w.cond)
weighted.mean(data1.1$d1x.joint,w1.1.joint) ## 0.0908
weighted.mean(data1.1$d1x.cond,w1.1.cond) ## 0.0930
## mean MB d4y: 0.0983
weighted.mean(data1.1$d2x.joint,w1.1.joint) ## 0.0981
weighted.mean(data1.1$d2x.cond,w1.1.cond) ## 0.1006
## mean MB d1x: 0.1108
weighted.mean(data1.1$d3y.joint,w1.1.joint) ## 0.1105
weighted.mean(data1.1$d3y.cond,w1.1.cond) ## 0.1011
## mean MB d2x: 0.0779
weighted.mean(data1.1$d4y.joint,w1.1.joint) ## 0.0786
weighted.mean(data1.1$d4y.cond,w1.1.cond) ## 0.0815
## mean MB dxy: 0.0261
weighted.mean(data1.1$dxy.joint,w1.1.joint) ## 0.0258
weighted.mean(data1.1$dxy.cond,w1.1.cond) ## 0.02705

## fill out a birds-comparison with this info


