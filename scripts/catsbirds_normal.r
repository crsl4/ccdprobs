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
    leaves = tre$tip.label[tre$edge[which(tre$edge[,1]==sample(c(5,6),1) & tre$edge[,2] < 5),2]] #works only for unrooted quartet
    leaves = as.numeric(leaves)
    sis = leaves
    fth = setdiff(1:4,sis)
    print(paste("sis",sis))
    print(paste("fth",fth))
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

data = data.frame(trees,d1x.joint,d2x.joint,d3y.joint,d4y.joint,dxy.joint,logwv.joint, logl.joint, logdens.joint,
    d1x.cond,d2x.cond,d3y.cond,d4y.cond,dxy.cond,logwv.cond, logl.cond, logdens.cond,
    d1x.nj,d2x.nj,d3y.nj,d4y.nj,dxy.nj,logwv.nj, logl.nj, logdens.nj)
head(data)
summary(data)
data[data$logwv.joint==0,]
data[data$logwv.cond==0,]
data[data$logwv.nj==0,]
data <- subset(data,logwv.nj!=0)

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
my.logw.nj = data$logwv.nj - mean(data$logwv.nj)
data$w.nj = exp(my.logw.nj)/sum(exp(my.logw.nj))
data[data$w.nj>0.01,]
length(data[data$w.nj>0.01,]$w.nj)
hist(data$w.nj)
plot(1:length(data$w.nj),cumsum(rev(sort(data$w.nj))))

save(data,file=paste0("data_normal_",who,".Rda"))

