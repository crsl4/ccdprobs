# r script to test the performance of normal distribution
# Claudia February 2016

library(ape)
source("branch-length.R")

simulateBranchLength.normal = function(nsim,x,eta=0.9, verbose=FALSE) {
    n = sum(x)
    prop.ag = (x[1,3] + x[3,1]) / n
    prop.ct = (x[2,4] + x[4,2]) / n
    prop.tv = (x[1,2] + x[1,4] + x[2,1] + x[2,3] + x[3,2] + x[3,4] + x[4,1] + x[4,3]) / n
    if(verbose)
        print(paste("prop.ag",prop.ag,"prop.ct",prop.ct,"prop.tv",prop.tv))
    p.est = (apply(x,1,sum) + apply(x,2,sum)) / (2*n)
    p.a = p.est[1]
    p.c = p.est[2]
    p.g = p.est[3]
    p.t = p.est[4]
    p.r = sum(p.est[c(1,3)])
    p.y = sum(p.est[c(2,4)])
    if(verbose)
        print(paste("p.a",p.a,"p.c",p.c,"p.g",p.g,"p.t",p.t,"p.r",p.r,"p.y",p.y))
    numer1 = 2*p.a*p.g*p.r
    denom1 = numer1 - p.r^2*prop.ag - p.a*p.g*prop.tv
    c1 = numer1 / denom1
    numer2 = 2*p.c*p.t*p.y
    denom2 = numer2 - p.y^2*prop.ct - p.c*p.t*prop.tv
    c2 = numer2 / denom2
    c3 = (2*p.a^2*p.g^2) / (p.r * denom1) +
         (2*p.c^2*p.t^2) / (p.y * denom2) +
         (p.r^2 * (p.c^2 + p.t^2) + p.y^2 * (p.a^2 + p.g^2) ) / (2*p.r^2*p.y^2 - p.r*p.y*prop.tv)
    mu = -2 * ( (p.a*p.g/p.r) * log(1 - p.r*prop.ag/(2*p.a*p.g) - prop.tv/(2*p.r) ) +
                (p.c*p.t/p.y) * log(1 - p.y*prop.ct/(2*p.c*p.t) - prop.tv/(2*p.y) ) +
                (p.r*p.y - p.a*p.g*p.y/p.r - p.c*p.t*p.r/p.y) * log(1 - prop.tv/(2*p.r*p.y)) )
    v = (1/ eta) * ((c1^2*prop.ag + c2^2*prop.ct + c3^2*prop.tv) - (c1*prop.ag + c2*prop.ct + c3*prop.tv)^2)/n
    if(verbose)
        print(paste("mu",mu,"v",v))
    w = rnorm(nsim,mu,sqrt(v))
    return( list(t=w,mu=mu,sigma=sqrt(v)) )
}


## 1) Compare in x-----y, one branch length
nsites = 1500
branch.length = 0.15
delta=0.001
s = seq(delta,2*branch.length,delta)
Q = randomQ(4,rescale=TRUE)
print(round(Q$Q,4))
print(round(Q$p,4))
print(min(diag(Q$Q)))
x = simulateSequenceSummary(nsites,Q,branch.length)
eta = 0.5
nsim=10000

log.like = gtr.log.like(x,s,Q)
log.like = log.like - max(log.like)
y = exp(log.like)
y = y / sum(y)
y = y / delta
df.true = data.frame(s,y)

Qlist.gtr = optim.gtr(x,Q$r)
log.like2 = gtr.log.like(x,s,Qlist.gtr$Q)
log.like2 = log.like2 - max(log.like2)
y2 = exp(log.like2)
y2 = y2/sum(y2)
y2 = y2 / delta
df.gtr = data.frame(s,y2)


jc = simulateBranchLength.jc(nsim,x,eta)
d.jc = density(jc)
df.jc = data.frame(x=d.jc$x,y=d.jc$y)

tn = simulateBranchLength.tn(nsim,x,eta)
d.tn = density(tn$t)
df.tn = data.frame(x=d.tn$x,y=d.tn$y)

norm = simulateBranchLength.normal(nsim,x,eta)
d.n = density(norm$t)
df.n = data.frame(x=d.n$x,y=d.n$y)

library(ggplot2)
p1 = ggplot(df.true, aes(x=s,y=y)) +
    geom_line(color="blue") +
    geom_line(aes(x=x,y=y),data=df.jc,color="red",linetype="dashed") +
    geom_line(aes(x=x,y=y),data=df.tn,color="darkgreen",linetype="dashed") +
    geom_line(aes(x=s,y=y2),data=df.gtr,color="gold") +
    geom_line(aes(x=x,y=y),data=df.n,color="black", linetype="dashed") +
    xlab('branch length') +
    ylab('densities') +
    ggtitle('Blue = TrueQ, Red = JC, Green = TN, Gold = GTR, Black = Normal')
plot(p1)


## 2) joint density for 4-taxon tree
## birds data
library(ape)
library(ggplot2)
source("4taxa_cats_dogs_functions.r")
source('branch-length.r')

# loglik ((1,2)x,(3,4)y)
gtr.log.lik.all = function(d1x,d2x,dxy,d3y,d4y,seq1.dist,seq2.dist, seq3.dist,seq4.dist,Q){
    suma = 0
    for(s in 1:nsites){
        lik12 = siteLik(d1x,d2x,seq1.dist[,s],seq2.dist[,s],Q$Q)
        lik34 = siteLik(d3y,d4y,seq3.dist[,s],seq4.dist[,s],Q$Q)
        L = lik12 %*% t(lik34)
        Pxy = matrixExp(Q$Q,dxy)
        L2 = L*Pxy
        Lik = Q$Q$p * L2
        suma = suma+log(sum(Lik))
    }
    return ( suma )
}


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

##need seq1.dist, seq2.dist, seq3.dist which are matrices
n=4
seq1.dist = matrix(0,n,nsites)
seq2.dist = matrix(0,n,nsites)
seq3.dist = matrix(0,n,nsites)
seq4.dist = matrix(0,n,nsites)
for(i in 1:nsites){
    if(seq1[i] == 'a'){
        seq1.dist[1,i] = 1
    } else if(seq1[i] == 'c'){
        seq1.dist[2,i] = 1
    } else if(seq1[i] == 'g'){
        seq1.dist[3,i] = 1
    } else if(seq1[i] == 't'){
        seq1.dist[4,i] = 1
    }
    if(seq2[i] == 'a'){
        seq2.dist[1,i] = 1
    } else if(seq2[i] == 'c'){
        seq2.dist[2,i] = 1
    } else if(seq2[i] == 'g'){
        seq2.dist[3,i] = 1
    } else if(seq2[i] == 't'){
        seq2.dist[4,i] = 1
    }
    if(seq3[i] == 'a'){
        seq3.dist[1,i] = 1
    } else if(seq3[i] == 'c'){
        seq3.dist[2,i] = 1
    } else if(seq3[i] == 'g'){
        seq3.dist[3,i] = 1
    } else if(seq3[i] == 't'){
        seq3.dist[4,i] = 1
    }
    if(seq4[i] == 'a'){
        seq4.dist[1,i] = 1
    } else if(seq4[i] == 'c'){
        seq4.dist[2,i] = 1
    } else if(seq4[i] == 'g'){
        seq4.dist[3,i] = 1
    } else if(seq4[i] == 't'){
        seq4.dist[4,i] = 1
    }
}


nuc = c("a","c","g","t")
verbose=TRUE
eta=0.5
out12 = matrix(0,n,n) # distance between 1 and 2
for(i in 1:n)
    for(j in 1:n)
        out12[i,j] = sum(seq1==nuc[i] & seq2==nuc[j])
if(verbose)
    print(out12)
checkMatCounts(out12)

r = rep(1,6)
Q = optim.gtr(out12,r)

## fixed in true value
d1x = 0.11
d2x = 0.078
d3y = 0.091
d4y = 0.098
dxy = 0.026
delta = 0.01
branch.length=0.1

v.dxy=seq(delta, 2*branch.length, by=delta)
lik = gtr.log.lik.all(d1x,d2x,dxy,d3y,d4y,seq1.dist,seq2.dist, seq3.dist,seq4.dist,Q)

loglik = rep(0,length(v.dxy))
for(i in 1:length(v.dxy)){
    loglik[i] = gtr.log.lik.all(d1x,d2x,v.dxy[i],d3y,d4y,seq1.dist,seq2.dist, seq3.dist,seq4.dist,Q)
}
loglik = loglik - max(loglik)
y2 = exp(loglik)
y2 = y2/sum(y2)
y2 = y2 / delta
df.gtr = data.frame(x=v.dxy,y=y2)


p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue") +
    geom_vline(xintercept=dxy)+
    xlab('branch length') +
    ylab('densities')

plot(p1)


## TN estimate
d12.tn = simulateBranchLength.tn(nsim=1, out12, eta=eta)
d12=d12.tn$t
if(verbose)
    print(d12)

out13 = matrix(0,n,n)     # distance between 1 and 3
for(i in 1:n)
    for(j in 1:n)
        out13[i,j] = sum(seq1==nuc[i] & seq3==nuc[j])
if(verbose)
    print(out13)
checkMatCounts(out13)
d13.tn = simulateBranchLength.tn(nsim=1, out13, eta=eta)
d13 = d13.tn$t
if(verbose)
    print(d13)

out23 = matrix(0,n,n)     # distance between 2 and 3
for(i in 1:n)
    for(j in 1:n)
        out23[i,j] = sum(seq2==nuc[i] & seq3==nuc[j])
if(verbose)
    print(out23)
checkMatCounts(out23)
d23.tn = simulateBranchLength.tn(nsim=1, out23, eta=eta)
d23 = d23.tn$t
if(verbose)
    print(d23)


d1x = 0.11
d2x = 0.078
seqdist = sequenceDist(d1x,d2x, seq1.dist, seq2.dist,Q$Q)
seqx = sampleSeq(seqdist)
if(verbose)
    print(seqx)

out34 = matrix(0,n,n) # distance between 3 and 4
for(i in 1:n)
    for(j in 1:n)
        out34[i,j] = sum(seq3==nuc[i] & seq4==nuc[j])
if(verbose)
    print(out34)
checkMatCounts(out34)
d34.tn = simulateBranchLength.tn(nsim=1, out34, eta=eta)
d34 = d34.tn$t
if(verbose)
    print(d34)

out3x = matrix(0,n,n) # distance between 3 and x
for(i in 1:nsites){
    out3x = out3x + seq3.dist[,i]%*%t(seqx[,i])
}
if(verbose)
    print(out3x)
checkMatCounts(out3x)
d3x.tn = simulateBranchLength.tn(nsim=1, out3x, eta=eta)
d3x = d3x.tn$t
if(verbose)
    print(d3x)

out4x = matrix(0,n,n) # distance between 4 and x
for(i in 1:nsites){
    out4x = out4x + seq4.dist[,i]%*%t(seqx[,i])
}
if(verbose)
    print(out4x)
checkMatCounts(out4x)
d4x.tn = simulateBranchLength.tn(nsim=1, out4x, eta=eta)
d4x = d4x.tn$t
if(verbose)
    print(d4x)

density = rep(0,length(v.dxy))
for(i in 1:length(v.dxy)){
    density[i] = logJointDensity.tn(d1x,d2x,d3y,d4y,v.dxy[i],d12.tn,d13.tn,d23.tn,d3x.tn,d4x.tn,d34.tn)
}
density = density - max(density)
y = exp(density)
y = y/sum(y)
y = y / delta
df.tn = data.frame(x=v.dxy,y=y)

## Normal density
d12.n = simulateBranchLength.normal(nsim=1, out12, eta=eta)
d12=d12.n$t
if(verbose)
    print(d12)
d13.n = simulateBranchLength.normal(nsim=1, out13, eta=eta)
d13 = d13.n$t
if(verbose)
    print(d13)
d23.n = simulateBranchLength.normal(nsim=1, out23, eta=eta)
d23 = d23.n$t
if(verbose)
    print(d23)
d34.n = simulateBranchLength.normal(nsim=1, out34, eta=eta)
d34 = d34.n$t
if(verbose)
    print(d34)
d3x.n = simulateBranchLength.normal(nsim=1, out3x, eta=eta)
d3x = d3x.n$t
if(verbose)
    print(d3x)
d4x.n = simulateBranchLength.normal(nsim=1, out4x, eta=eta)
d4x = d4x.n$t
if(verbose)
    print(d4x)

density = rep(0,length(v.dxy))
for(i in 1:length(v.dxy)){
    density[i] = logJointDensity.norm(d1x,d2x,d3y,d4y,v.dxy[i],d12.n,d13.n,d23.n,d3x.n,d4x.n,d34.n)
}
density = density - max(density)
y = exp(density)
y = y/sum(y)
y = y / delta
df.norm = data.frame(x=v.dxy,y=y)



p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue") +
    geom_line(aes(x=x,y=y),data=df.tn,color="red",linetype="dashed") +
    geom_line(aes(x=x,y=y),data=df.norm,color="darkgreen",linetype="dashed") +
    geom_vline(xintercept=dxy)+
    xlab('branch length') +
    ylab('densities') +
    ggtitle('Blue = GTR, Red = TN, Green = Normal')

plot(p1)

## simulated data
library(ape)
source('branch-length.r')
source('4taxa_simulated_functions.r')
source('4taxa_cats_dogs_functions.r')

nsites = 1500
t="((1,2),3,4);"
tre=read.tree(text=t)
branch.length = c(0.026,0.11,0.078,0.091,0.098) #bl as birds mb: dxy,d1x,d2x,d3y,d4y
tre$edge.length <- branch.length
nuc <- c('a','c','g','t')

trueQ = randomQ(4,rescale=TRUE)
print(round(trueQ$Q,4))
print(round(trueQ$p,4))

# simulate seqx
seqx = sample(nuc,size=nsites,prob=trueQ$p,replace=TRUE)

#simulate seqy
s=branch.length[1]
P = matrixExp(trueQ,s)
seqy = numeric(nsites)
for ( i in 1:nsites )
    seqy[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

# simulate seq1
s=branch.length[2]
P = matrixExp(trueQ,s)
seq1 = numeric(nsites)
for ( i in 1:nsites )
    seq1[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

# simulate seq2
s=branch.length[3]
P = matrixExp(trueQ,s)
seq2 = numeric(nsites)
for ( i in 1:nsites )
    seq2[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

# simulate seq3
s=branch.length[4]
P = matrixExp(trueQ,s)
seq3 = numeric(nsites)
for ( i in 1:nsites )
    seq3[i] = sample(nuc,size=1,prob=P[which(nuc==seqy[i]),])

# simulate seq4
s=branch.length[5]
P = matrixExp(trueQ,s)
seq4 = numeric(nsites)
for ( i in 1:nsites )
    seq4[i] = sample(nuc,size=1,prob=P[which(nuc==seqy[i]),])



# loglik ((1,2)x,(3,4)y)
gtr.log.lik.all = function(d1x,d2x,dxy,d3y,d4y,seq1.dist,seq2.dist, seq3.dist,seq4.dist,Q){
    suma = 0
    for(s in 1:nsites){
        lik12 = siteLik(d1x,d2x,seq1.dist[,s],seq2.dist[,s],Q$Q)
        lik34 = siteLik(d3y,d4y,seq3.dist[,s],seq4.dist[,s],Q$Q)
        L = lik12 %*% t(lik34)
        Pxy = matrixExp(Q$Q,dxy)
        L2 = L*Pxy
        Lik = Q$Q$p * L2
        suma = suma+log(sum(Lik))
    }
    return ( suma )
}

##need seq1.dist, seq2.dist, seq3.dist which are matrices
n=4
seq1.dist = matrix(0,n,nsites)
seq2.dist = matrix(0,n,nsites)
seq3.dist = matrix(0,n,nsites)
seq4.dist = matrix(0,n,nsites)
for(i in 1:nsites){
    if(seq1[i] == 'a'){
        seq1.dist[1,i] = 1
    } else if(seq1[i] == 'c'){
        seq1.dist[2,i] = 1
    } else if(seq1[i] == 'g'){
        seq1.dist[3,i] = 1
    } else if(seq1[i] == 't'){
        seq1.dist[4,i] = 1
    }
    if(seq2[i] == 'a'){
        seq2.dist[1,i] = 1
    } else if(seq2[i] == 'c'){
        seq2.dist[2,i] = 1
    } else if(seq2[i] == 'g'){
        seq2.dist[3,i] = 1
    } else if(seq2[i] == 't'){
        seq2.dist[4,i] = 1
    }
    if(seq3[i] == 'a'){
        seq3.dist[1,i] = 1
    } else if(seq3[i] == 'c'){
        seq3.dist[2,i] = 1
    } else if(seq3[i] == 'g'){
        seq3.dist[3,i] = 1
    } else if(seq3[i] == 't'){
        seq3.dist[4,i] = 1
    }
    if(seq4[i] == 'a'){
        seq4.dist[1,i] = 1
    } else if(seq4[i] == 'c'){
        seq4.dist[2,i] = 1
    } else if(seq4[i] == 'g'){
        seq4.dist[3,i] = 1
    } else if(seq4[i] == 't'){
        seq4.dist[4,i] = 1
    }
}


nuc = c("a","c","g","t")
verbose=TRUE
eta=0.5
out12 = matrix(0,n,n) # distance between 1 and 2
for(i in 1:n)
    for(j in 1:n)
        out12[i,j] = sum(seq1==nuc[i] & seq2==nuc[j])
if(verbose)
    print(out12)
checkMatCounts(out12)

r = rep(1,6)
Q = optim.gtr(out12,r)

## fixed in true value
d1x = 0.11
d2x = 0.078
d3y = 0.091
d4y = 0.098
dxy = 0.026
delta = 0.01
branch.length=0.1

v.dxy=seq(delta, 2*branch.length, by=delta)
lik = gtr.log.lik.all(d1x,d2x,dxy,d3y,d4y,seq1.dist,seq2.dist, seq3.dist,seq4.dist,Q)

loglik = rep(0,length(v.dxy))
for(i in 1:length(v.dxy)){
    loglik[i] = gtr.log.lik.all(d1x,d2x,v.dxy[i],d3y,d4y,seq1.dist,seq2.dist, seq3.dist,seq4.dist,Q)
}
loglik = loglik - max(loglik)
y2 = exp(loglik)
y2 = y2/sum(y2)
y2 = y2 / delta
df.gtr = data.frame(x=v.dxy,y=y2)

library(ggplot2)
p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue") +
    geom_vline(xintercept=dxy)+
    xlab('branch length') +
    ylab('densities')

plot(p1)


## TN estimate
d12.tn = simulateBranchLength.tn(nsim=1, out12, eta=eta)
d12=d12.tn$t
if(verbose)
    print(d12)

out13 = matrix(0,n,n)     # distance between 1 and 3
for(i in 1:n)
    for(j in 1:n)
        out13[i,j] = sum(seq1==nuc[i] & seq3==nuc[j])
if(verbose)
    print(out13)
checkMatCounts(out13)
d13.tn = simulateBranchLength.tn(nsim=1, out13, eta=eta)
d13 = d13.tn$t
if(verbose)
    print(d13)

out23 = matrix(0,n,n)     # distance between 2 and 3
for(i in 1:n)
    for(j in 1:n)
        out23[i,j] = sum(seq2==nuc[i] & seq3==nuc[j])
if(verbose)
    print(out23)
checkMatCounts(out23)
d23.tn = simulateBranchLength.tn(nsim=1, out23, eta=eta)
d23 = d23.tn$t
if(verbose)
    print(d23)

## use true seqx
seqx.dist = matrix(0,n,nsites)
for(i in 1:nsites){
    if(seqx[i] == 'a'){
        seqx.dist[1,i] = 1
    } else if(seqx[i] == 'c'){
        seqx.dist[2,i] = 1
    } else if(seq3[i] == 'g'){
        seqx.dist[3,i] = 1
    } else if(seq3[i] == 't'){
        seqx.dist[4,i] = 1
    }
}
seqx=seqx.dist

## use estimated seqx
d1x = 0.11
d2x = 0.078
seqdist = sequenceDist(d1x,d2x, seq1.dist, seq2.dist,Q$Q)
seqx = sampleSeq(seqdist)

out34 = matrix(0,n,n) # distance between 3 and 4
for(i in 1:n)
    for(j in 1:n)
        out34[i,j] = sum(seq3==nuc[i] & seq4==nuc[j])
if(verbose)
    print(out34)
checkMatCounts(out34)
d34.tn = simulateBranchLength.tn(nsim=1, out34, eta=eta)
d34 = d34.tn$t
if(verbose)
    print(d34)

out3x = matrix(0,n,n) # distance between 3 and x
for(i in 1:nsites){
    out3x = out3x + seq3.dist[,i]%*%t(seqx[,i])
}
if(verbose)
    print(out3x)
checkMatCounts(out3x)
d3x.tn = simulateBranchLength.tn(nsim=1, out3x, eta=eta)
d3x = d3x.tn$t
if(verbose)
    print(d3x)

out4x = matrix(0,n,n) # distance between 4 and x
for(i in 1:nsites){
    out4x = out4x + seq4.dist[,i]%*%t(seqx[,i])
}
if(verbose)
    print(out4x)
checkMatCounts(out4x)
d4x.tn = simulateBranchLength.tn(nsim=1, out4x, eta=eta)
d4x = d4x.tn$t
if(verbose)
    print(d4x)

v.dxy2 = seq(-branch.length+delta, 2*branch.length-branch.length,by=delta)
## use v.dxy for estimated seq, ad v.dxy2 for true seq
density = rep(0,length(v.dxy))
for(i in 1:length(v.dxy)){
    density[i] = logJointDensity.tn(d1x,d2x,d3y,d4y,v.dxy[i],d12.tn,d13.tn,d23.tn,d3x.tn,d4x.tn,d34.tn)
}
density = density - max(density)
y = exp(density)
y = y/sum(y)
y = y / delta
df.tn = data.frame(x=v.dxy,y=y)

## Normal density
d12.n = simulateBranchLength.normal(nsim=1, out12, eta=eta)
d12=d12.n$t
if(verbose)
    print(d12)
d13.n = simulateBranchLength.normal(nsim=1, out13, eta=eta)
d13 = d13.n$t
if(verbose)
    print(d13)
d23.n = simulateBranchLength.normal(nsim=1, out23, eta=eta)
d23 = d23.n$t
if(verbose)
    print(d23)
d34.n = simulateBranchLength.normal(nsim=1, out34, eta=eta)
d34 = d34.n$t
if(verbose)
    print(d34)

d3x.n = simulateBranchLength.normal(nsim=1, out3x, eta=eta)
d3x = d3x.n$t
if(verbose)
    print(d3x)
d4x.n = simulateBranchLength.normal(nsim=1, out4x, eta=eta)
d4x = d4x.n$t
if(verbose)
    print(d4x)

density = rep(0,length(v.dxy))
## use v.dxy2 for true seqx, and v.dxy for estimated seqx
for(i in 1:length(v.dxy)){
    density[i] = logJointDensity.norm(d1x,d2x,d3y,d4y,v.dxy[i],d12.n,d13.n,d23.n,d3x.n,d4x.n,d34.n)
}
density = density - max(density)
y = exp(density)
y = y/sum(y)
y = y / delta
df.norm = data.frame(x=v.dxy,y=y)


library(ggplot2)
p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue") +
    geom_line(aes(x=x,y=y),data=df.tn,color="red",linetype="dashed") +
    geom_line(aes(x=x,y=y),data=df.norm,color="black", linetype="dashed") +
    geom_vline(xintercept=dxy)+
    xlab('branch length') +
    ylab('densities') +
    ggtitle('Blue = GTR, Red = TN, Black = Normal')

plot(p1)

