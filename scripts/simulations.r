## r script to do simulations step by step, to understand how
## the likelihood is similar (or not) to the gamma dist of the BL,
## and how the weights behave
## Claudia March 2016

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')
library(ggplot2)
library(weights)

## how do weights behave??
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
nsites=150000
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
out23 = countsMatrix(seq2,seq3)
out34 = countsMatrix(seq3,seq4)

nreps = 1000
logwv = rep(0,nreps)
d1x=rep(0,nreps)
d2x=rep(0,nreps)
d3y=rep(0,nreps)
d4y=rep(0,nreps)
dxy=rep(0,nreps)

for(nr in 1:nreps){
    print(nr)
    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)
    jc34 = simulateBranchLength.jc(nsim=1,out34,eta=eta)
        ## starting point for d3x,d4x
    t0=(max(c(jc12$t,jc13$t,jc23$t))+min(c(jc12$t,jc13$t,jc23$t)))/2

    t.lik12 = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
    t.lik13 = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
    t.lik23 = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)

    d12 = t.lik12$t
    d13 = t.lik13$t
    d23 = t.lik23$t

    d1x[nr] = (d12+d13-d23)/2
    d2x[nr] = (d12+d23-d13)/2
    d3x = (d13+d23-d12)/2

    print(paste(d1x[nr], d2x[nr], d3x))

    ##seqx.dist = sequenceDist(d1x0,d2x0,seq1.dist,seq2.dist,Q) ## with true values d1x,d2x
    seqx.dist = sequenceDist(d1x[nr],d2x[nr],seq1.dist,seq2.dist,Q) ## with est values d1x,d2x

    ##t.lik3x = simulateBranchLength.lik(nsim=1, seqx.dist,seq3.dist,Q,t0=t0,eta=eta)
    t.lik4x = simulateBranchLength.lik(nsim=1, seqx.dist,seq4.dist,Q,t0=t0,eta=eta, verbose=FALSE)
    t.lik34 = simulateBranchLength.lik(nsim=1, seq3.dist,seq4.dist,Q,t0=jc34$t,eta=eta, verbose=FALSE)

    ##d3x = t.lik3x$t
    d4x = t.lik4x$t
    d34 = t.lik34$t

    d3y[nr] = (d34+d3x-d4x)/2
    d4y[nr] = (d34+d4x-d3x)/2
    dxy[nr] = (d3x+d4x-d34)/2


    if(d1x[nr]<0 || d2x[nr]<0 || d3y[nr]<0 || d4y[nr]<0 || dxy[nr]<0){
        print("negative bl")
    } else{
        print("all positive")
        print(paste(d1x[nr],d2x[nr], dxy[nr], d3y[nr], d4y[nr]))
        suma = gtr.log.lik.all(d1x[nr],d2x[nr],dxy[nr],d3y[nr],d4y[nr],seq1.dist, seq2.dist, seq3.dist, seq4.dist, Q)
        logdens = logJointDensity.lik(t.lik12, t.lik13, t.lik23, t.lik4x, t.lik34)
        logprior = logPriorExpDist(d1x[nr], d2x[nr], d3y[nr], d4y[nr], dxy[nr], m=0.1)
        print(paste('loglik', suma))
        print(paste('logdens', logdens))
        print(paste('logprior', logprior))
        logwv[nr] = logprior +suma - logdens
    }
}

data = data.frame(d1x,d2x,d3y,d4y,dxy,logwv)
head(data)
summary(data)
data[data$logwv==0,]
length(data[data$logwv==0,]$logwv)
data <- subset(data,logwv!=0)
my.logw = data$logwv - mean(data$logwv)
data$w = exp(my.logw)/sum(exp(my.logw))
data[data$w>0.01,]
length(data[data$w>0.01,]$w)
hist(data$w)
##save(data,file="data_simulations.Rda")

m.1x=weighted.mean(data$d1x,data$w)
m2.1x=weighted.mean(data$d1x^2,data$w)
v.1x=m2.1x-m.1x^2
m.1x
m.1x-2*sqrt(v.1x)
m.1x+2*sqrt(v.1x)
d1x0
weighted.quantile(data$d1x,data$w,probs=0.025)
weighted.quantile(data$d1x,data$w,probs=0.975)
plot(data$d1x,data$w, main="red=true, blue=weighted mean")
abline(v=d1x0, col="red")
abline(v=m.1x,col="blue")


m.2x=weighted.mean(data$d2x,data$w)
m2.2x=weighted.mean(data$d2x^2,data$w)
v.2x=m2.2x-m.2x^2
m.2x
m.2x-2*sqrt(v.2x)
m.2x+2*sqrt(v.2x)
d2x0
weighted.quantile(data$d2x,data$w,probs=0.025)
weighted.quantile(data$d2x,data$w,probs=0.975)
plot(data$d2x,data$w, main="red=true, blue=weighted mean")
abline(v=d2x0, col="red")
abline(v=m.2x,col="blue")


m.3y=weighted.mean(data$d3y,data$w)
m2.3y=weighted.mean(data$d3y^2,data$w)
v.3y=m2.3y-m.3y^2
m.3y
m.3y-2*sqrt(v.3y)
m.3y+2*sqrt(v.3y)
d3y0
weighted.quantile(data$d3y,data$w,probs=0.025)
weighted.quantile(data$d3y,data$w,probs=0.975)
plot(data$d3y,data$w, main="red=true, blue=weighted mean")
abline(v=d3y0, col="red")
abline(v=m.3y,col="blue")

m.4y=weighted.mean(data$d4y,data$w)
m2.4y=weighted.mean(data$d4y^2,data$w)
v.4y=m2.4y-m.4y^2
m.4y
m.4y-2*sqrt(v.4y)
m.4y+2*sqrt(v.4y)
d4y0
weighted.quantile(data$d4y,data$w,probs=0.025)
weighted.quantile(data$d4y,data$w,probs=0.975)
plot(data$d4y,data$w, main="red=true, blue=weighted mean")
abline(v=d4y0, col="red")
abline(v=m.4y,col="blue")


m.xy=weighted.mean(data$dxy,data$w)
m2.xy=weighted.mean(data$dxy^2,data$w)
v.xy=m2.xy-m.xy^2
m.xy
m.xy-2*sqrt(v.xy)
m.xy+2*sqrt(v.xy)
dxy0
weighted.quantile(data$dxy,data$w,probs=0.025)
weighted.quantile(data$dxy,data$w,probs=0.975)
plot(data$dxy,data$w, main="red=true, blue=weighted mean")
abline(v=dxy0, col="red")
abline(v=m.xy,col="blue")


## weighted histograms
wtd.hist(data$d1x,weight=data$w)
abline(v=d1x0, col="red")
abline(v=m.1x,col="blue")

wtd.hist(data$d2x,weight=data$w)
abline(v=d2x0, col="red")
abline(v=m.2x,col="blue")

wtd.hist(data$d3y,weight=data$w)
abline(v=d3y0, col="red")
abline(v=m.3y,col="blue")

wtd.hist(data$d4y,weight=data$w)
abline(v=d4y0, col="red")
abline(v=m.4y,col="blue")

wtd.hist(data$dxy,weight=data$w)
abline(v=dxy0, col="red")
abline(v=m.xy,col="blue")

## real histograms for nsites=150000 case
hist(data$d1x)
abline(v=d1x0, col="red")
abline(v=m.1x,col="blue")

hist(data$d2x)
abline(v=d2x0, col="red")
abline(v=m.2x,col="blue")

hist(data$d3y)
abline(v=d3y0, col="red")
abline(v=m.3y,col="blue")

hist(data$d4y)
abline(v=d4y0, col="red")
abline(v=m.4y,col="blue")

hist(data$dxy)
abline(v=dxy0, col="red")
abline(v=m.xy,col="blue")

## replicate Case (1,2)---(3,4)
## let's redo this case, but with fixed simulated data
## instead of drawing a new simulated data each time
who="(1,2)---(3,4)"
d1x0=0.11
d2x0=0.078
dxy0 = 0.03
d3y0 = 0.091
d4y0 = 0.098
d1x0=0.1
d2x0=0.1
dxy0 = 0.01
d3y0 = 0.1
d4y0 = 0.1
eta = 0.5
nsites=1500
nuc <- c('a','c','g','t')
nrep = 100

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

out12 = countsMatrix(seq1,seq2)
out13 = countsMatrix(seq1,seq3)
out23 = countsMatrix(seq2,seq3)
out34 = countsMatrix(seq3,seq4)

d1x=rep(0,nrep)
d2x=rep(0,nrep)
d3x=rep(0,nrep)
d3y=rep(0,nrep)
d4y=rep(0,nrep)
dxy=rep(0,nrep)

alpha12=rep(0,nrep)
beta12=rep(0,nrep)
alpha13=rep(0,nrep)
beta13=rep(0,nrep)
alpha23=rep(0,nrep)
beta23=rep(0,nrep)
alpha3x=rep(0,nrep)
beta3x=rep(0,nrep)
alpha4x=rep(0,nrep)
beta4x=rep(0,nrep)
alpha34=rep(0,nrep)
beta34=rep(0,nrep)


for(nr in 1:nrep){
    print(nr)
    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)
    jc34 = simulateBranchLength.jc(nsim=1,out34,eta=eta)
    ## starting point for d3x,d4x
    t0=(max(c(jc12$t,jc13$t,jc23$t))+min(c(jc12$t,jc13$t,jc23$t)))/2

    t.lik12 = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
    t.lik13 = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
    t.lik23 = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)

    alpha12[nr] = t.lik12$alpha
    beta12[nr] = t.lik12$beta
    alpha13[nr] = t.lik13$alpha
    beta13[nr] = t.lik13$beta
    alpha23[nr] = t.lik23$alpha
    beta23[nr] = t.lik23$beta

    d12 = t.lik12$t
    d13 = t.lik13$t
    d23 = t.lik23$t

    d1x[nr] = (d12+d13-d23)/2
    d2x[nr] = (d12+d23-d13)/2
    d3x[nr] = (d13+d23-d12)/2

    ##seqx.dist = sequenceDist(d1x0,d2x0,seq1.dist,seq2.dist,Q) ## with true values d1x,d2x
    seqx.dist = sequenceDist(d1x[nr],d2x[nr],seq1.dist,seq2.dist,Q) ## with est values d1x,d2x

    ##t.lik3x = simulateBranchLength.lik(nsim=1, seqx.dist,seq3.dist,Q,t0=t0,eta=eta)
    t.lik4x = simulateBranchLength.lik(nsim=1, seqx.dist,seq4.dist,Q,t0=t0,eta=eta)
    t.lik34 = simulateBranchLength.lik(nsim=1, seq3.dist,seq4.dist,Q,t0=jc34$t,eta=eta)

    ##d3x = t.lik3x$t
    d4x = t.lik4x$t
    d34 = t.lik34$t

    d3y[nr] = (d34+d3x[nr]-d4x)/2
    d4y[nr] = (d34+d4x-d3x[nr])/2
    dxy[nr] = (d3x[nr]+d4x-d34)/2

    alpha4x[nr] = t.lik4x$alpha
    beta4x[nr] = t.lik4x$beta
    alpha34[nr] = t.lik34$alpha
    beta34[nr] = t.lik34$beta
}

## todo: run this to compare case of fixed dataset
data=data.frame(d1x=d1x, d2x=d2x,d3x=d3x, d3y=d3y, d4y=d4y, dxy=dxy,alpha12=alpha12,beta12=beta12,alpha13=alpha13, beta13=beta13, alpha23=alpha23, beta23=beta23, mean12=alpha12/beta12, mean13=alpha13/beta13, mean23=alpha23/beta23, var12=alpha12/beta12^2, var13=alpha13/beta13^2, var23=alpha23/beta23^2, alpha4x=alpha4x, beta4x=beta4x, alpha34=alpha34, beta34=beta34, mean4x=alpha4x/beta4x, mean34=alpha34/beta34, var4x=alpha4x/beta4x^2, var34=alpha34/beta34^2)
head(data)


## plot d12,d13,d23
m12 = mean(data$mean12)
v12 = mean(data$var12)
m13 = mean(data$mean13)
v13 = mean(data$var13)
m23 = mean(data$mean23)
v23 = mean(data$var23)

w12 = rgamma(10000,m12^2/v12, m12/v12)
df12 = density(w12)
w13 = rgamma(10000,m13^2/v13, m13/v13)
df13 = density(w13)
w23 = rgamma(10000,m23^2/v23, m23/v23)
df23 = density(w23)

df = data.frame(x12=df12$x, y12=df12$y, x13=df13$x, y13=df13$y, x23=df23$x, y23=df23$y)

title = paste(who,'Gamma density with average mean&var, line true value,\n Blue=d12, Red=d13, Gold=d23,\n nsites', nsites, 'nrep', nrep)
p1 = ggplot(df,aes(x=x12,y=y12))+
    geom_line(colour="blue")+
    geom_vline(xintercept=d1x0+d2x0, colour="blue")+
    geom_line(aes(x=x13,y=y13), data=df, colour="red")+
    geom_vline(xintercept=d1x0+d3y0+dxy0, colour="red")+
    geom_line(aes(x=x23,y=y23), data=df, colour="gold")+
    geom_vline(xintercept=d2x0+d3y0+dxy0, colour="gold")+
    ggtitle(title)

plot(p1)

## d1x,d2x,d3x
title = paste(who,'Blue=d1x, Red=d2x, Green=d3x, nsites', nsites, 'nrep', nrep)
p1 = ggplot(data,aes(x=d1x))+
    geom_freqpoly(colour="blue")+
    geom_vline(xintercept=d1x0, colour="blue")+
    geom_vline(xintercept=mean(d1x), colour="blue", linetype="dashed")+
    geom_freqpoly(aes(x=d2x), data=data, colour="red")+
    geom_vline(xintercept=d2x0, colour="red")+
    geom_vline(xintercept=mean(d2x), colour="red", linetype="dashed")+
    geom_freqpoly(aes(x=d3x), data=data, colour="darkgreen")+
    geom_vline(xintercept=d3y0+dxy0, colour="darkgreen")+
    geom_vline(xintercept=mean(d3x), colour="darkgreen", linetype="dashed")+
    ggtitle(title)

plot(p1)

hist(data$d1x)
abline(v=d1x0, col="blue")
hist(data$d2x)
abline(v=d2x0, col="red")
hist(data$d3x)
abline(v=d3y0+dxy0, col="green")

## d4x, d34
m4x = mean(data$mean4x)
v4x = mean(data$var4x)
m34 = mean(data$mean34)
v34 = mean(data$var34)

w4x = rgamma(10000,m4x^2/v4x, m4x/v4x)
df4x = density(w4x)
w34 = rgamma(10000,m34^2/v34, m34/v34)
df34 = density(w34)

df = data.frame(x4x=df4x$x, y4x=df4x$y, x34=df34$x, y34=df34$y)

title = paste(who,'Gamma density with average mean&var, line true value,\n Red=d4x, Gold=d34,\n nsites', nsites, 'nrep', nrep)
p1 = ggplot(df,aes(x=x4x,y=y4x))+
    geom_line(colour="red")+
    geom_vline(xintercept=dxy0+d4y0, colour="red")+
    geom_vline(xintercept=m4x,colour="red", linetype="dashed")+
    geom_line(aes(x=x34,y=y34), data=df, colour="gold")+
    geom_vline(xintercept=d3y0+d4y0, colour="gold")+
    geom_vline(xintercept=m34,colour="gold", linetype="dashed")+
    ggtitle(title)

plot(p1)



## plot d3y,d4y,dxy
title = paste(who,'Blue=d3y, Red=d4y, Green=dxy, nsites', nsites, 'nrep', nrep)
p1 = ggplot(data,aes(x=d3y))+
    geom_freqpoly(colour="blue")+
    geom_vline(xintercept=d3y0, colour="blue")+
    geom_vline(xintercept=mean(d3y), colour="blue", linetype="dashed")+
    geom_freqpoly(aes(x=d4y), data=data, colour="red")+
    geom_vline(xintercept=d4y0, colour="red")+
    geom_vline(xintercept=mean(d4y), colour="red", linetype="dashed")+
    geom_freqpoly(aes(x=dxy), data=data, colour="darkgreen")+
    geom_vline(xintercept=dxy0, colour="darkgreen")+
    geom_vline(xintercept=mean(dxy), colour="darkgreen", linetype="dashed")+
    ggtitle(title)

plot(p1)

## replicate Case (1,2)---(3,4)
## conclusion: all branches seem to be correctly generated
## when we estimate seqx.dist with the true distances
## what happens when using the estimated distances?
## things look similar
who="(1,2)---(3,4)"
d1x0=0.11
d2x0=0.078
dxy0 = 0.03
d3y0 = 0.091
d4y0 = 0.098
d1x0=0.1
d2x0=0.1
dxy0 = 0.01
d3y0 = 0.1
d4y0 = 0.1
eta = 0.5
nsites=1500
nuc <- c('a','c','g','t')
nrep = 100

d1x=rep(0,nrep)
d2x=rep(0,nrep)
d3x=rep(0,nrep)
d3y=rep(0,nrep)
d4y=rep(0,nrep)
dxy=rep(0,nrep)

alpha12=rep(0,nrep)
beta12=rep(0,nrep)
alpha13=rep(0,nrep)
beta13=rep(0,nrep)
alpha23=rep(0,nrep)
beta23=rep(0,nrep)
alpha3x=rep(0,nrep)
beta3x=rep(0,nrep)
alpha4x=rep(0,nrep)
beta4x=rep(0,nrep)
alpha34=rep(0,nrep)
beta34=rep(0,nrep)


for(nr in 1:nrep){
    print(nr)
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
    out23 = countsMatrix(seq2,seq3)
    out34 = countsMatrix(seq3,seq4)
    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)
    jc34 = simulateBranchLength.jc(nsim=1,out34,eta=eta)
    ## starting point for d3x,d4x
    t0=(max(c(jc12$t,jc13$t,jc23$t))+min(c(jc12$t,jc13$t,jc23$t)))/2

    t.lik12 = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
    t.lik13 = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
    t.lik23 = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)

    alpha12[nr] = t.lik12$alpha
    beta12[nr] = t.lik12$beta
    alpha13[nr] = t.lik13$alpha
    beta13[nr] = t.lik13$beta
    alpha23[nr] = t.lik23$alpha
    beta23[nr] = t.lik23$beta

    d12 = t.lik12$t
    d13 = t.lik13$t
    d23 = t.lik23$t

    d1x[nr] = (d12+d13-d23)/2
    d2x[nr] = (d12+d23-d13)/2
    d3x[nr] = (d13+d23-d12)/2

    ##seqx.dist = sequenceDist(d1x0,d2x0,seq1.dist,seq2.dist,Q) ## with true values d1x,d2x
    seqx.dist = sequenceDist(d1x[nr],d2x[nr],seq1.dist,seq2.dist,Q) ## with est values d1x,d2x

    ##t.lik3x = simulateBranchLength.lik(nsim=1, seqx.dist,seq3.dist,Q,t0=t0,eta=eta)
    t.lik4x = simulateBranchLength.lik(nsim=1, seqx.dist,seq4.dist,Q,t0=t0,eta=eta)
    t.lik34 = simulateBranchLength.lik(nsim=1, seq3.dist,seq4.dist,Q,t0=jc34$t,eta=eta)

    ##d3x = t.lik3x$t
    d4x = t.lik4x$t
    d34 = t.lik34$t

    d3y[nr] = (d34+d3x[nr]-d4x)/2
    d4y[nr] = (d34+d4x-d3x[nr])/2
    dxy[nr] = (d3x[nr]+d4x-d34)/2

    alpha4x[nr] = t.lik4x$alpha
    beta4x[nr] = t.lik4x$beta
    alpha34[nr] = t.lik34$alpha
    beta34[nr] = t.lik34$beta
}

data=data.frame(d1x=d1x, d2x=d2x,d3x=d3x, d3y=d3y, d4y=d4y, dxy=dxy,alpha12=alpha12,beta12=beta12,alpha13=alpha13, beta13=beta13, alpha23=alpha23, beta23=beta23, mean12=alpha12/beta12, mean13=alpha13/beta13, mean23=alpha23/beta23, var12=alpha12/beta12^2, var13=alpha13/beta13^2, var23=alpha23/beta23^2, alpha4x=alpha4x, beta4x=beta4x, alpha34=alpha34, beta34=beta34, mean4x=alpha4x/beta4x, mean34=alpha34/beta34, var4x=alpha4x/beta4x^2, var34=alpha34/beta34^2)
head(data)


## plot d12,d13,d23
m12 = mean(data$mean12)
v12 = mean(data$var12)
m13 = mean(data$mean13)
v13 = mean(data$var13)
m23 = mean(data$mean23)
v23 = mean(data$var23)

w12 = rgamma(10000,m12^2/v12, m12/v12)
df12 = density(w12)
w13 = rgamma(10000,m13^2/v13, m13/v13)
df13 = density(w13)
w23 = rgamma(10000,m23^2/v23, m23/v23)
df23 = density(w23)

df = data.frame(x12=df12$x, y12=df12$y, x13=df13$x, y13=df13$y, x23=df23$x, y23=df23$y)

title = paste(who,'Gamma density with average mean&var, line true value,\n Blue=d12, Red=d13, Gold=d23,\n nsites', nsites, 'nrep', nrep)
p1 = ggplot(df,aes(x=x12,y=y12))+
    geom_line(colour="blue")+
    geom_vline(xintercept=d1x0+d2x0, colour="blue")+
    geom_line(aes(x=x13,y=y13), data=df, colour="red")+
    geom_vline(xintercept=d1x0+d3y0+dxy0, colour="red")+
    geom_line(aes(x=x23,y=y23), data=df, colour="gold")+
    geom_vline(xintercept=d2x0+d3y0+dxy0, colour="gold")+
    ggtitle(title)

plot(p1)

## d1x,d2x,d3x
title = paste(who,'Blue=d1x, Red=d2x, Green=d3x, nsites', nsites, 'nrep', nrep)
p1 = ggplot(data,aes(x=d1x))+
    geom_freqpoly(colour="blue")+
    geom_vline(xintercept=d1x0, colour="blue")+
    geom_vline(xintercept=mean(d1x), colour="blue", linetype="dashed")+
    geom_freqpoly(aes(x=d2x), data=data, colour="red")+
    geom_vline(xintercept=d2x0, colour="red")+
    geom_vline(xintercept=mean(d2x), colour="red", linetype="dashed")+
    geom_freqpoly(aes(x=d3x), data=data, colour="darkgreen")+
    geom_vline(xintercept=d3y0+dxy0, colour="darkgreen")+
    geom_vline(xintercept=mean(d3x), colour="darkgreen", linetype="dashed")+
    ggtitle(title)

plot(p1)

hist(data$d1x)
abline(v=d1x0, col="blue")
hist(data$d2x)
abline(v=d2x0, col="red")
hist(data$d3x)
abline(v=d3y0+dxy0, col="green")

## d4x, d34
m4x = mean(data$mean4x)
v4x = mean(data$var4x)
m34 = mean(data$mean34)
v34 = mean(data$var34)

w4x = rgamma(10000,m4x^2/v4x, m4x/v4x)
df4x = density(w4x)
w34 = rgamma(10000,m34^2/v34, m34/v34)
df34 = density(w34)

df = data.frame(x4x=df4x$x, y4x=df4x$y, x34=df34$x, y34=df34$y)

title = paste(who,'Gamma density with average mean&var, line true value,\n Red=d4x, Gold=d34,\n nsites', nsites, 'nrep', nrep)
p1 = ggplot(df,aes(x=x4x,y=y4x))+
    geom_line(colour="red")+
    geom_vline(xintercept=dxy0+d4y0, colour="red")+
    geom_vline(xintercept=m4x,colour="red", linetype="dashed")+
    geom_line(aes(x=x34,y=y34), data=df, colour="gold")+
    geom_vline(xintercept=d3y0+d4y0, colour="gold")+
    geom_vline(xintercept=m34,colour="gold", linetype="dashed")+
    ggtitle(title)

plot(p1)



## plot d3y,d4y,dxy
title = paste(who,'Blue=d3y, Red=d4y, Green=dxy, nsites', nsites, 'nrep', nrep)
p1 = ggplot(data,aes(x=d3y))+
    geom_freqpoly(colour="blue")+
    geom_vline(xintercept=d3y0, colour="blue")+
    geom_vline(xintercept=mean(d3y), colour="blue", linetype="dashed")+
    geom_freqpoly(aes(x=d4y), data=data, colour="red")+
    geom_vline(xintercept=d4y0, colour="red")+
    geom_vline(xintercept=mean(d4y), colour="red", linetype="dashed")+
    geom_freqpoly(aes(x=dxy), data=data, colour="darkgreen")+
    geom_vline(xintercept=dxy0, colour="darkgreen")+
    geom_vline(xintercept=mean(dxy), colour="darkgreen", linetype="dashed")+
    ggtitle(title)

plot(p1)


## replicate: Case x---(3,4)
## we want to know if from d3x,d4x,d34 we can get good estimated of d3y,d4y,dxy
## conclusion: d3y,d4y,dxy seem to be accurately computed
who="x---(3,4)"
d1x=0.1
d2x=0.1
dxy0 = 0.01
d3y0 = 0.1
d4y0 = 0.1
eta = 0.5
nsites=1500
nuc <- c('a','c','g','t')
nrep = 100
d3y=rep(0,nrep)
d4y=rep(0,nrep)
dxy=rep(0,nrep)

for(nr in 1:nrep){
    print(nr)
    Q = randomQ(4,rescale=TRUE)
    r=Q$r
    p=Q$p
    ## simulate seqx
    seqx = sample(nuc,size=nsites,prob=Q$p,replace=TRUE)

    ## simulate seq1
    P = matrixExp(Q,d1x)
    seq1 = numeric(nsites)
    for ( i in 1:nsites )
        seq1[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])
    ## simulate seq2
    P = matrixExp(Q,d2x)
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
    seqx.dist = sequenceDist(d1x,d2x,seq1.dist,seq2.dist,Q) ## with true values d1x,d2x
    seq3.dist = seqMatrix(seq3)
    seq4.dist = seqMatrix(seq4)


    ## gamma density
    out12 = countsMatrix(seq1,seq2)
    out13 = countsMatrix(seq1,seq3)
    out23 = countsMatrix(seq2,seq3)
    out34 = countsMatrix(seq3,seq4)
    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)
    jc34 = simulateBranchLength.jc(nsim=1,out34,eta=eta)
    ## starting point for d3x,d4x
    t0=(max(c(jc12$t,jc13$t,jc23$t))+min(c(jc12$t,jc13$t,jc23$t)))/2

    t.lik3x = simulateBranchLength.lik(nsim=1, seqx.dist,seq3.dist,Q,t0=t0,eta=eta)
    t.lik4x = simulateBranchLength.lik(nsim=1, seqx.dist,seq4.dist,Q,t0=t0,eta=eta)
    t.lik34 = simulateBranchLength.lik(nsim=1, seq3.dist,seq4.dist,Q,t0=jc34$t,eta=eta)

    d3x = t.lik3x$t
    d4x = t.lik4x$t
    d34 = t.lik34$t

    d3y[nr] = (d34+d3x-d4x)/2
    d4y[nr] = (d34+d4x-d3x)/2
    dxy[nr] = (d3x+d4x-d34)/2
}

data=data.frame(d3y=d3y, d4y=d4y,dxy=dxy)
head(data)

hist(d3y)
hist(d4y)
hist(dxy)

title = paste(who,'Blue=d3y, Red=d4y, Green=dxy, nsites', nsites, 'nrep', nrep)
p1 = ggplot(data,aes(x=d3y))+
    geom_freqpoly(colour="blue")+
    geom_vline(xintercept=d3y0, colour="blue")+
    geom_vline(xintercept=mean(d3y), colour="blue", linetype="dashed")+
    geom_freqpoly(aes(x=d4y), data=data, colour="red")+
    geom_vline(xintercept=d4y0, colour="red")+
    geom_vline(xintercept=mean(d4y), colour="red", linetype="dashed")+
    geom_freqpoly(aes(x=dxy), data=data, colour="darkgreen")+
    geom_vline(xintercept=dxy0, colour="darkgreen")+
    geom_vline(xintercept=mean(dxy), colour="darkgreen", linetype="dashed")+
    ggtitle(title)

plot(p1)

## replicate: Case x---(3,4)
## we will verify if we can also get good estimates of d3x,d4x,d34
## when we do not have a sequence in x, but a probability dist.
## recall that in the case (1,2)--3 we were able to generate d12,d13,d23
## close to the true values
## conclusion: d3x, d4x, d34 seem to be close to true values!
## that is, no bias for using seqx.dist estimated
who="x---(3,4)"
d1x=0.01
d2x=0.01
dxy = 0.01
d3y = 0.1
d4y = 0.1
eta = 0.5
nsites=1500
nuc <- c('a','c','g','t')
nrep = 100
alpha3x=rep(0,nrep)
beta3x=rep(0,nrep)
alpha4x=rep(0,nrep)
beta4x=rep(0,nrep)
alpha34=rep(0,nrep)
beta34=rep(0,nrep)

for(nr in 1:nrep){
    print(nr)
    Q = randomQ(4,rescale=TRUE)
    r=Q$r
    p=Q$p
    ## simulate seqx
    seqx = sample(nuc,size=nsites,prob=Q$p,replace=TRUE)

    ## simulate seq1
    P = matrixExp(Q,d1x)
    seq1 = numeric(nsites)
    for ( i in 1:nsites )
        seq1[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])
    ## simulate seq2
    P = matrixExp(Q,d2x)
    seq2 = numeric(nsites)
    for ( i in 1:nsites )
        seq2[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

    ## simulate seqy
    P = matrixExp(Q,dxy)
    seqy = numeric(nsites)
    for ( i in 1:nsites )
        seqy[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

    ## simulate seq3
    P = matrixExp(Q,d3y)
    seq3 = numeric(nsites)
    for ( i in 1:nsites )
        seq3[i] = sample(nuc,size=1,prob=P[which(nuc==seqy[i]),])
    ## simulate seq4
    P = matrixExp(Q,d4y)
    seq4 = numeric(nsites)
    for ( i in 1:nsites )
        seq4[i] = sample(nuc,size=1,prob=P[which(nuc==seqy[i]),])


    seq1.dist = seqMatrix(seq1)
    seq2.dist = seqMatrix(seq2)
    seqx.dist = sequenceDist(d1x,d2x,seq1.dist,seq2.dist,Q) ## with true values d1x,d2x
    seq3.dist = seqMatrix(seq3)
    seq4.dist = seqMatrix(seq4)


    ## gamma density
    out12 = countsMatrix(seq1,seq2)
    out13 = countsMatrix(seq1,seq3)
    out23 = countsMatrix(seq2,seq3)
    out34 = countsMatrix(seq3,seq4)
    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)
    jc34 = simulateBranchLength.jc(nsim=1,out34,eta=eta)
    ## starting point for d3x,d4x
    t0=(max(c(jc12$t,jc13$t,jc23$t))+min(c(jc12$t,jc13$t,jc23$t)))/2

    t.lik3x = simulateBranchLength.lik(nsim=1, seqx.dist,seq3.dist,Q,t0=t0,eta=eta)
    t.lik4x = simulateBranchLength.lik(nsim=1, seqx.dist,seq4.dist,Q,t0=t0,eta=eta)
    t.lik34 = simulateBranchLength.lik(nsim=1, seq3.dist,seq4.dist,Q,t0=jc34$t,eta=eta)

    alpha3x[nr] = t.lik3x$alpha
    beta3x[nr] = t.lik3x$beta
    alpha4x[nr] = t.lik4x$alpha
    beta4x[nr] = t.lik4x$beta
    alpha34[nr] = t.lik34$alpha
    beta34[nr] = t.lik34$beta
}

data=data.frame(alpha3x=alpha3x,beta3x=beta3x,alpha4x=alpha4x, beta4x=beta4x, alpha34=alpha34, beta34=beta34, mean3x=alpha3x/beta3x, mean4x=alpha4x/beta4x, mean34=alpha34/beta34, var3x=alpha3x/beta3x^2, var4x=alpha4x/beta4x^2, var34=alpha34/beta34^2)
head(data)

m3x = mean(data$mean3x)
v3x = mean(data$var3x)
m4x = mean(data$mean4x)
v4x = mean(data$var4x)
m34 = mean(data$mean34)
v34 = mean(data$var34)

w3x = rgamma(10000,m3x^2/v3x, m3x/v3x)
df3x = density(w3x)
w4x = rgamma(10000,m4x^2/v4x, m4x/v4x)
df4x = density(w4x)
w34 = rgamma(10000,m34^2/v34, m34/v34)
df34 = density(w34)

df = data.frame(x3x=df3x$x, y3x=df3x$y, x4x=df4x$x, y4x=df4x$y, x34=df34$x, y34=df34$y)

title = paste(who,'Gamma density with average mean&var, line true value,\n Blue=d3x, Red=d4x, Gold=d34,\n nsites', nsites, 'nrep', nrep, 'd12', d1x+d2x)
p1 = ggplot(df,aes(x=x3x,y=y3x))+
    geom_line(colour="blue")+
    geom_vline(xintercept=dxy+d3y, colour="blue")+
    geom_vline(xintercept=m3x,colour="blue", linetype="dashed")+
    geom_line(aes(x=x4x,y=y4x), data=df, colour="red")+
    geom_vline(xintercept=dxy+d4y, colour="red")+
    geom_vline(xintercept=m4x,colour="red", linetype="dashed")+
    geom_line(aes(x=x34,y=y34), data=df, colour="gold")+
    geom_vline(xintercept=d3y+d4y, colour="gold")+
    geom_vline(xintercept=m34,colour="gold", linetype="dashed")+
    ggtitle(title)

plot(p1)






## Case (1,2)---3
## importance weights:
## they look good, they weighted histograms are weighted statistics
## do cover the real values
who="Case (1,2)---3"
d1x0=0.1
d2x0=0.1
d3x0=0.1
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
## simulate seq3
P = matrixExp(Q,d3x0)
seq3 = numeric(nsites)
for ( i in 1:nsites )
    seq3[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

seq1.dist = seqMatrix(seq1)
seq2.dist = seqMatrix(seq2)
seq3.dist = seqMatrix(seq3)

## gamma density
out12 = countsMatrix(seq1,seq2)
out13 = countsMatrix(seq1,seq3)
out23 = countsMatrix(seq2,seq3)

nreps = 1000
logwv = rep(0,nreps)
d1x = rep(0,nreps)
d2x = rep(0,nreps)
d3x = rep(0,nreps)
for(nr in 1:nreps){
    print(nr)
    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)

    t.lik12 = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
    t.lik13 = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
    t.lik23 = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)

    d12 = t.lik12$t
    d13 = t.lik13$t
    d23 = t.lik23$t

    d1x[nr] = (d12+d13-d23)/2
    d2x[nr] = (d12+d23-d13)/2
    d3x[nr] = (d13+d23-d12)/2

    if(d1x[nr]<0 || d2x[nr]<0 || d3x[nr]<0){
        print("negative bl")
    } else{
        P1 = matrixExp(Q,d1x[nr])
        P2 = matrixExp(Q,d2x[nr])
        P3 = matrixExp(Q,d3x[nr])
        suma = 0
        for(ns in 1:nsites){
            l1 = P1 %*% seq1.dist[,ns]
            l2 = P2 %*% seq2.dist[,ns]
            l3 = P3 %*% seq3.dist[,ns]
            suma = suma + log(sum(Q$p * l1 * l2 * l3))
        }
        logdens = (t.lik12$alpha-1)*log(d12)-t.lik12$beta*d12 + (t.lik13$alpha-1)*log(d13)-t.lik13$beta*d13 + (t.lik23$alpha-1)*log(d23)-t.lik23$beta*d23
        logwv[nr] = -10*(d1x[nr]+d2x[nr]+d3x[nr]) +suma - logdens
    }
}

data = data.frame(d1x,d2x,d3x,logwv)
head(data)
summary(data)
data[data$logwv==0,]
length(data[data$logwv==0,]$logwv)
data <- subset(data,logwv!=0)
my.logw = data$logwv - mean(data$logwv)
data$w = exp(my.logw)/sum(exp(my.logw))
data[data$w>0.01,]
hist(data$w)
save(data,file="data1m.Rda")

m.1x=weighted.mean(data$d1x,data$w)
m2.1x=weighted.mean(data$d1x^2,data$w)
v.1x=m2.1x-m.1x^2
m.1x
m.1x-2*sqrt(v.1x)
m.1x+2*sqrt(v.1x)
d1x0
weighted.quantile(data$d1x,data$w,probs=0.025)
weighted.quantile(data$d1x,data$w,probs=0.975)
plot(data$d1x,data$w, main="red=true, blue=weighted mean")
abline(v=d1x0, col="red")
abline(v=m.1x,col="blue")


m.2x=weighted.mean(data$d2x,data$w)
m2.2x=weighted.mean(data$d2x^2,data$w)
v.2x=m2.2x-m.2x^2
m.2x
m.2x-2*sqrt(v.2x)
m.2x+2*sqrt(v.2x)
d2x0
weighted.quantile(data$d2x,data$w,probs=0.025)
weighted.quantile(data$d2x,data$w,probs=0.975)
plot(data$d2x,data$w, main="red=true, blue=weighted mean")
abline(v=d2x0, col="red")
abline(v=m.2x,col="blue")


m.3x=weighted.mean(data$d3x,data$w)
m2.3x=weighted.mean(data$d3x^2,data$w)
v.3x=m2.3x-m.3x^2
m.3x
m.3x-2*sqrt(v.3x)
m.3x+2*sqrt(v.3x)
d3x0
weighted.quantile(data$d3x,data$w,probs=0.025)
weighted.quantile(data$d3x,data$w,probs=0.975)
plot(data$d3x,data$w, main="red=true, blue=weighted mean")
abline(v=d3x0, col="red")
abline(v=m.3x,col="blue")


## weighted histograms
wtd.hist(data$d1x,weight=data$w)
abline(v=d1x0, col="red")
abline(v=m.1x,col="blue")

wtd.hist(data$d2x,weight=data$w)
abline(v=d2x0, col="red")
abline(v=m.2x,col="blue")

wtd.hist(data$d3x,weight=data$w)
abline(v=d3x0, col="red")
abline(v=m.3x,col="blue")




## to do: do similar replicate case for 3 branches, without the likelihood at this moment,
## to see if the gammas are close in mean/var to the true values of 12,13,23
## intuition says that we are estimating the sum d1x+d2x correctly (because we can
## estimate 12 correctly in the 1---2 case
## replicate: Case (1,2)---3
## conclusion: simulated d12,d13,d23 from gamma
## are close to true values on average
d1x=0.01
d2x=0.01
d3x=0.15
eta = 0.5
nsites=1500
nuc <- c('a','c','g','t')
nrep = 100
alpha12=rep(0,nrep)
beta12=rep(0,nrep)
alpha13=rep(0,nrep)
beta13=rep(0,nrep)
alpha23=rep(0,nrep)
beta23=rep(0,nrep)

for(nr in 1:nrep){
    print(nr)
    Q = randomQ(4,rescale=TRUE)
    r=Q$r
    p=Q$p
    ## simulate seqx
    seqx = sample(nuc,size=nsites,prob=Q$p,replace=TRUE)

    ## simulate seq1
    P = matrixExp(Q,d1x)
    seq1 = numeric(nsites)
    for ( i in 1:nsites )
        seq1[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])
    ## simulate seq2
    P = matrixExp(Q,d2x)
    seq2 = numeric(nsites)
    for ( i in 1:nsites )
        seq2[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])
    ## simulate seq3
    P = matrixExp(Q,d3x)
    seq3 = numeric(nsites)
    for ( i in 1:nsites )
        seq3[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

    seq1.dist = seqMatrix(seq1)
    seq2.dist = seqMatrix(seq2)
    seq3.dist = seqMatrix(seq3)

    ## gamma density
    out12 = countsMatrix(seq1,seq2)
    out13 = countsMatrix(seq1,seq3)
    out23 = countsMatrix(seq2,seq3)

    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)

    t.lik12 = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
    t.lik13 = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
    t.lik23 = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)

    alpha12[nr] = t.lik12$alpha
    beta12[nr] = t.lik12$beta
    alpha13[nr] = t.lik13$alpha
    beta13[nr] = t.lik13$beta
    alpha23[nr] = t.lik23$alpha
    beta23[nr] = t.lik23$beta
}

data=data.frame(alpha12=alpha12,beta12=beta12,alpha13=alpha13, beta13=beta13, alpha23=alpha23, beta23=beta23, mean12=alpha12/beta12, mean13=alpha13/beta13, mean23=alpha23/beta23, var12=alpha12/beta12^2, var13=alpha13/beta13^2, var23=alpha23/beta23^2)
head(data)

m12 = mean(data$mean12)
v12 = mean(data$var12)
m13 = mean(data$mean13)
v13 = mean(data$var13)
m23 = mean(data$mean23)
v23 = mean(data$var23)

w12 = rgamma(10000,m12^2/v12, m12/v12)
df12 = density(w12)
w13 = rgamma(10000,m13^2/v13, m13/v13)
df13 = density(w13)
w23 = rgamma(10000,m23^2/v23, m23/v23)
df23 = density(w23)

df = data.frame(x12=df12$x, y12=df12$y, x13=df13$x, y13=df13$y, x23=df23$x, y23=df23$y)

title = paste(who,'Gamma density with average mean&var, line true value,\n Blue=d12, Red=d13, Gold=d23,\n nsites', nsites, 'nrep', nrep)
p1 = ggplot(df,aes(x=x12,y=y12))+
    geom_line(colour="blue")+
    geom_vline(xintercept=d1x+d2x, colour="blue")+
    geom_line(aes(x=x13,y=y13), data=df, colour="red")+
    geom_vline(xintercept=d1x+d3x, colour="red")+
    geom_line(aes(x=x23,y=y23), data=df, colour="gold")+
    geom_vline(xintercept=d2x+d3x, colour="gold")+
    ggtitle(title)

plot(p1)


## replicate: Case (1,2)---3
## but now fixing the dataset
who="(1,2)---3"
d1x=0.1
d2x=0.1
d3x=0.1
eta = 0.5
nsites=1500
nuc <- c('a','c','g','t')
nrep = 100

Q = randomQ(4,rescale=TRUE)
r=Q$r
p=Q$p
## simulate seqx
seqx = sample(nuc,size=nsites,prob=Q$p,replace=TRUE)

## simulate seq1
P = matrixExp(Q,d1x)
seq1 = numeric(nsites)
for ( i in 1:nsites )
    seq1[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])
## simulate seq2
P = matrixExp(Q,d2x)
seq2 = numeric(nsites)
for ( i in 1:nsites )
    seq2[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])
## simulate seq3
P = matrixExp(Q,d3x)
seq3 = numeric(nsites)
for ( i in 1:nsites )
    seq3[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

seq1.dist = seqMatrix(seq1)
seq2.dist = seqMatrix(seq2)
seq3.dist = seqMatrix(seq3)

## gamma density
out12 = countsMatrix(seq1,seq2)
out13 = countsMatrix(seq1,seq3)
out23 = countsMatrix(seq2,seq3)


alpha12=rep(0,nrep)
beta12=rep(0,nrep)
alpha13=rep(0,nrep)
beta13=rep(0,nrep)
alpha23=rep(0,nrep)
beta23=rep(0,nrep)

for(nr in 1:nrep){
    print(nr)

    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)

    t.lik12 = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
    t.lik13 = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
    t.lik23 = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)

    alpha12[nr] = t.lik12$alpha
    beta12[nr] = t.lik12$beta
    alpha13[nr] = t.lik13$alpha
    beta13[nr] = t.lik13$beta
    alpha23[nr] = t.lik23$alpha
    beta23[nr] = t.lik23$beta
}

data=data.frame(alpha12=alpha12,beta12=beta12,alpha13=alpha13, beta13=beta13, alpha23=alpha23, beta23=beta23, mean12=alpha12/beta12, mean13=alpha13/beta13, mean23=alpha23/beta23, var12=alpha12/beta12^2, var13=alpha13/beta13^2, var23=alpha23/beta23^2)
head(data)

m12 = mean(data$mean12)
v12 = mean(data$var12)
m13 = mean(data$mean13)
v13 = mean(data$var13)
m23 = mean(data$mean23)
v23 = mean(data$var23)

w12 = rgamma(10000,m12^2/v12, m12/v12)
df12 = density(w12)
w13 = rgamma(10000,m13^2/v13, m13/v13)
df13 = density(w13)
w23 = rgamma(10000,m23^2/v23, m23/v23)
df23 = density(w23)

df = data.frame(x12=df12$x, y12=df12$y, x13=df13$x, y13=df13$y, x23=df23$x, y23=df23$y)

title = paste(who,'Gamma density with average mean&var, line true value,\n Blue=d12, Red=d13, Gold=d23,\n nsites', nsites, 'nrep', nrep)
p1 = ggplot(df,aes(x=x12,y=y12))+
    geom_line(colour="blue")+
    geom_vline(xintercept=d1x+d2x, colour="blue")+
    geom_line(aes(x=x13,y=y13), data=df, colour="red")+
    geom_vline(xintercept=d1x+d3x, colour="red")+
    geom_line(aes(x=x23,y=y23), data=df, colour="gold")+
    geom_vline(xintercept=d2x+d3x, colour="gold")+
    ggtitle(title)

plot(p1)


## replicate d1x,d2x,d3x from gamma and compare to true values
## conclusion: sampled d12,d13,d23 compute good estimates of d1x,d2x,d3x
who="(1,2)--3"
d1x0=0.01
d2x0=0.01
d3x0=0.7
eta = 0.5
nsites=1500
nuc <- c('a','c','g','t')
nrep = 100
d1x=rep(0,nrep)
d2x=rep(0,nrep)
d3x=rep(0,nrep)


for(nr in 1:nrep){
    print(nr)
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
    ## simulate seq3
    P = matrixExp(Q,d3x0)
    seq3 = numeric(nsites)
    for ( i in 1:nsites )
        seq3[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

    seq1.dist = seqMatrix(seq1)
    seq2.dist = seqMatrix(seq2)
    seq3.dist = seqMatrix(seq3)

    ## gamma density
    out12 = countsMatrix(seq1,seq2)
    out13 = countsMatrix(seq1,seq3)
    out23 = countsMatrix(seq2,seq3)

    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)

    t.lik12 = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
    t.lik13 = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
    t.lik23 = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)

    d12 = t.lik12$t
    d13 = t.lik13$t
    d23 = t.lik23$t

    d1x[nr] = (d12+d13-d23)/2
    d2x[nr] = (d12+d23-d13)/2
    d3x[nr] = (d13+d23-d12)/2
}

data=data.frame(d1x=d1x, d2x=d2x,d3x=d3x)
head(data)

hist(d1x)
abline(v=d1x0,col="blue")
hist(d2x)
abline(v=d2x0,col="red")
hist(d3x)
abline(v=d3x0,col="green")

title = paste(who,'Blue=d1x, Red=d2x, Green=d3x, nsites', nsites, 'nrep', nrep)
p1 = ggplot(data,aes(x=d1x))+
    geom_freqpoly(colour="blue")+
    geom_vline(xintercept=d1x0, colour="blue")+
    geom_vline(xintercept=mean(d1x), colour="blue", linetype="dashed")+
    geom_freqpoly(aes(x=d2x), data=data, colour="red")+
    geom_vline(xintercept=d2x0, colour="red")+
    geom_vline(xintercept=mean(d2x), colour="red", linetype="dashed")+
    geom_freqpoly(aes(x=d3x), data=data, colour="darkgreen")+
    geom_vline(xintercept=d3x0, colour="darkgreen")+
    geom_vline(xintercept=mean(d3x), colour="darkgreen", linetype="dashed")+
    ggtitle(title)

plot(p1)

## replicate d1x,d2x,d3x from gamma and compare to true values
## but now let's fix the dataset
## conclusion: this does not look as great, there are bad datasets and good
who="(1,2)--3"
d1x0=0.1
d2x0=0.1
d3x0=0.1
eta = 0.5
nsites=1500
nuc <- c('a','c','g','t')
nrep = 100

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
## simulate seq3
P = matrixExp(Q,d3x0)
seq3 = numeric(nsites)
for ( i in 1:nsites )
    seq3[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

seq1.dist = seqMatrix(seq1)
seq2.dist = seqMatrix(seq2)
seq3.dist = seqMatrix(seq3)

## gamma density
out12 = countsMatrix(seq1,seq2)
out13 = countsMatrix(seq1,seq3)
out23 = countsMatrix(seq2,seq3)

d1x=rep(0,nrep)
d2x=rep(0,nrep)
d3x=rep(0,nrep)


for(nr in 1:nrep){
    print(nr)

    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)

    t.lik12 = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
    t.lik13 = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
    t.lik23 = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)

    d12 = t.lik12$t
    d13 = t.lik13$t
    d23 = t.lik23$t

    d1x[nr] = (d12+d13-d23)/2
    d2x[nr] = (d12+d23-d13)/2
    d3x[nr] = (d13+d23-d12)/2
}

data=data.frame(d1x=d1x, d2x=d2x,d3x=d3x)
head(data)

hist(d1x)
hist(d2x)
hist(d3x)

title = paste(who,'Blue=d1x, Red=d2x, Green=d3x, nsites', nsites, 'nrep', nrep)
p1 = ggplot(data,aes(x=d1x))+
    geom_freqpoly(colour="blue")+
    geom_vline(xintercept=d1x0, colour="blue")+
    geom_vline(xintercept=mean(d1x), colour="blue", linetype="dashed")+
    geom_freqpoly(aes(x=d2x), data=data, colour="red")+
    geom_vline(xintercept=d2x0, colour="red")+
    geom_vline(xintercept=mean(d2x), colour="red", linetype="dashed")+
    geom_freqpoly(aes(x=d3x), data=data, colour="darkgreen")+
    geom_vline(xintercept=d3x0, colour="darkgreen")+
    geom_vline(xintercept=mean(d3x), colour="darkgreen", linetype="dashed")+
    ggtitle(title)

plot(p1)

## Case 1-----2
## need to study the effect ot t0, nsites, eta
## eta: 0.5 covers better the lik
## nsties, as expected, gives lower variance
## aroung t_MLE as nsites increases
## also, both lik and gamma centered in t0 (good)
who="Case 1---2"
t0=0.05
eta = 0.5
nsites=15000

nuc <- c('a','c','g','t')
Q = randomQ(4,rescale=TRUE)
r=Q$r
p=Q$p
## simulate seq1
seq1 = sample(nuc,size=nsites,prob=Q$p,replace=TRUE)
## simulate seq2
P = matrixExp(Q,t0)
seq2 = numeric(nsites)
for ( i in 1:nsites )
    seq2[i] = sample(nuc,size=1,prob=P[which(nuc==seq1[i]),])

seq1.dist = seqMatrix(seq1)
seq2.dist = seqMatrix(seq2)
out12 = countsMatrix(seq1,seq2)

delta=0.001
t=seq(t0/2,2*t0,by=delta)
ll=rep(0, length(t))
for(i in 1:length(t)){
    l = loglik(seq1.dist,seq2.dist,Q,t[i])
    ll[i] = l$ll
}

logl = ll - max(ll)
y2 = exp(logl)
y2 = y2/sum(y2)
y2 = y2 / delta
df.gtr = data.frame(x=t,y=y2)

nsim=10000
jc = simulateBranchLength.jc(nsim=1,out12,eta=eta)
t.lik = simulateBranchLength.lik(nsim, seq1.dist,seq2.dist,Q,t0=jc$t,eta=eta)
d.lik = density(t.lik$t)
df.lik = data.frame(x=d.lik$x,y=d.lik$y)


title = paste(who,'Blue = GTR, Red = Gamma,\n', 'eta', eta, 'nsites', nsites, 'Black line=true t, Gold line=MLE')
p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue") +
    geom_line(aes(x=x,y=y),data=df.lik,color="red",linetype="dashed") +
    xlab('branch length') +
    ylab('densities') +
    ggtitle(title)+
    geom_vline(xintercept = t0)+
    geom_vline(xintercept = t.lik$alpha/t.lik$beta, colour="gold", linetype="dashed")

plot(p1)

## weights
## works fine except for the lik bias on short bl
who="Case 1---2"
t0=0.15
eta = 0.5
nsites=15000

nuc <- c('a','c','g','t')
Q = randomQ(4,rescale=TRUE)
r=Q$r
p=Q$p
## simulate seq1
seq1 = sample(nuc,size=nsites,prob=Q$p,replace=TRUE)
## simulate seq2
P = matrixExp(Q,t0)
seq2 = numeric(nsites)
for ( i in 1:nsites )
    seq2[i] = sample(nuc,size=1,prob=P[which(nuc==seq1[i]),])
seq1.dist = seqMatrix(seq1)
seq2.dist = seqMatrix(seq2)
out12 = countsMatrix(seq1,seq2)

nreps = 1000
logwv = rep(0,nreps)
bl = rep(0,nreps)
for(i in 1:nreps){
    print(i)
    jc = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    t.lik = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=jc$t,eta=eta)
    bl[i] = t.lik$t
    l = loglik(seq1.dist,seq2.dist,Q,t.lik$t)
    logwv[i] = -10*t.lik$t+l$ll-((t.lik$alpha-1)*log(t.lik$t)-t.lik$beta*t.lik$t)
}

data = data.frame(bl,logwv)
head(data)
summary(data)
my.logw = logwv - mean(logwv)
data$w = exp(my.logw)/sum(exp(my.logw))
data[data$w>0.01,]
hist(data$w)

m=weighted.mean(data$bl,data$w)
m2=weighted.mean(data$bl^2,data$w)
v=m2-m^2
m-2*sqrt(v)
m+2*sqrt(v)
m
plot(data$bl,data$w, main="red=true, blue=weighted mean")
abline(v=t0, col="red")
abline(v=m,col="blue")


## Comparison between 1--x--2 and 1----2
## conclusion: these two cases behave similarly, there are
## no differences that stand out.
d1x=0.025
d2x=0.025
eta = 0.5
nsites=1500000

nuc <- c('a','c','g','t')
Q = randomQ(4,rescale=TRUE)
r=Q$r
p=Q$p
## simulate seqx
seqx = sample(nuc,size=nsites,prob=Q$p,replace=TRUE)

## simulate seq1
P = matrixExp(Q,d1x)
seq1 = numeric(nsites)
for ( i in 1:nsites )
    seq1[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])
## simulate seq2
P = matrixExp(Q,d2x)
seq2 = numeric(nsites)
for ( i in 1:nsites )
    seq2[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])


## simulate seq11
seq11 = sample(nuc,size=nsites,prob=Q$p,replace=TRUE)

## simulate seq22
P = matrixExp(Q,d1x+d2x)
seq22 = numeric(nsites)
for ( i in 1:nsites )
    seq22[i] = sample(nuc,size=1,prob=P[which(nuc==seq11[i]),])

seq1.dist = seqMatrix(seq1)
seq2.dist = seqMatrix(seq2)
seq11.dist = seqMatrix(seq11)
seq22.dist = seqMatrix(seq22)

out12 = countsMatrix(seq1,seq2)
out1122 = countsMatrix(seq11,seq22)

jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
jc1122 = simulateBranchLength.jc(nsim=1,out1122,eta=eta)
print(jc12$t)
print(jc1122$t)

nsim=10000
t.lik12 = simulateBranchLength.lik(nsim, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
t.lik1122 = simulateBranchLength.lik(nsim, seq11.dist,seq22.dist,Q,t0=jc1122$t,eta=eta)
print(t.lik12$t[1])
print(t.lik1122$t[1])

d.lik12 = density(t.lik12$t)
d.lik1122 = density(t.lik1122$t)

df.lik = data.frame(x12=d.lik12$x,y12=d.lik12$y,x1122=d.lik1122$x,y1122=d.lik1122$y)
head(df.lik)

title = paste('gamma density \n Blue = 1--x--2, Red = 1---2, nsites', nsites)
p1 = ggplot(df.lik, aes(x=x12,y=y12)) +
    geom_line(color="blue")+
    geom_vline(xintercept = d1x+d2x, colour="blue")+
    geom_vline(xintercept = t.lik12$alpha/t.lik12$beta, colour="blue", linetype="dashed")+
    geom_line(aes(x=x1122,y=y1122),data=df.lik,colour="red")+
    geom_vline(xintercept = d1x+d2x, colour="red")+
    geom_vline(xintercept = t.lik1122$alpha/t.lik1122$beta, colour="red", linetype="dashed")+
    ggtitle(title)

plot(p1)


## replicate the case of 1----2
## to see if on average MLE is close to true value
who="Case 1---2"
t0=0.1
eta = 0.5
nsites=1500
nuc <- c('a','c','g','t')
nrep = 100
alpha=rep(0,nrep)
beta=rep(0,nrep)

for(nr in 1:nrep){
    print(nr)
    Q = randomQ(4,rescale=TRUE)
    r=Q$r
    p=Q$p
    ## simulate seq1
    seq1 = sample(nuc,size=nsites,prob=Q$p,replace=TRUE)
    ## simulate seq2
    P = matrixExp(Q,t0)
    seq2 = numeric(nsites)
    for ( i in 1:nsites )
        seq2[i] = sample(nuc,size=1,prob=P[which(nuc==seq1[i]),])

    seq1.dist = seqMatrix(seq1)
    seq2.dist = seqMatrix(seq2)
    out12 = countsMatrix(seq1,seq2)
    jc = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    t.lik = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=jc$t,eta=eta)
    alpha[nr] = t.lik$alpha
    beta[nr] = t.lik$beta
}

data = data.frame(alpha=alpha, beta=beta, mean=alpha/beta, var=alpha/beta^2)
head(data)

m = mean(data$mean)
v = mean(data$var)

w = rgamma(10000,m^2/v, m/v)
d = density(w)
df = data.frame(x=d$x, y=d$y)

title = paste(who,'Gamma density with average mean&var, line true value\n nsites', nsites, 'nrep', nrep)
p1 = ggplot(df,aes(x=x,y=y))+
    geom_line(colour="blue")+
    geom_vline(xintercept=t0)+
    ggtitle(title)

plot(p1)

## -----------------------------------------------------------------------------------------------------
## case (1,2)--3
## computation and plot of the likelihood,
## not done successfully
delta=0.001
t1=seq(d1x/2,2*d1x,by=delta)
t2=seq(d2x/2,2*d2x,by=delta)
t3=seq(d3x/2,2*d3x,by=delta)
ll=rep(0, length(t1)*length(t2)*length(t3))
x12=rep(0, length(t1)*length(t2)*length(t3))
x13=rep(0, length(t1)*length(t2)*length(t3))
x23=rep(0, length(t1)*length(t2)*length(t3))
ind = 1
for(i in 1:length(t1)){
    for(j in 1:length(t2)){
        for(k in 1:length(t3)){
            print(ind)
            P1 = matrixExp(Q,t1[i])
            P2 = matrixExp(Q,t2[j])
            P3 = matrixExp(Q,t3[k])
            suma = 0
            for(ns in 1:nsites){
                l1 = P1 %*% seq1.dist[,ns]
                l2 = P2 %*% seq2.dist[,ns]
                l3 = P3 %*% seq3.dist[,ns]
                suma = suma + log(sum(Q$p * l1 * l2 * l3))
            }
            ll[ind] = suma
            x12[ind] = t1[i]
            x13[ind] = t2[j]
            x23[ind] = t3[k]
            ind = ind +1
        }
    }
}


df.gtr = data.frame(x12=x12,x13=x13,x23=x23,y=ll)
head(df.gtr)

## compare this behavior of d12 with simulating from d12 (instead of d1x,d2x)
## check for d1x,,d2x,d3x now
## aqui voy, tengo q ver como graficar las cosas bien, tengo q normalizar y
## parece qe p1 y p2 ya son comparable, pero no estan centradas donde deben!
df12 = subset(df.gtr,x13==d1x+d3x & x23==d2x+d3x)
head(df12)

df12=within(df12,y <- exp(y))
df12=within(df12,y <- y/sum(y))
df12=within(df12,y <- y/delta)


p2 = ggplot(df12, aes(x=x12,y=y))+
    geom_line(color="red")
plot(p2)

xlab('branch length') +
    ylab('densities') +
    ggtitle(title)+
    geom_vline(xintercept = t0)+
    geom_vline(xintercept = t.lik$alpha/t.lik$beta, colour="gold", linetype="dashed")

plot(p2)
