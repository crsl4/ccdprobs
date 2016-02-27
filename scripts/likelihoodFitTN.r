## r script to compare the TN approach to the likelihood
## on the internal branch length
## Claudia February 2016

## x-----y, bl=t
source('branch-length.r')
nsites=500
branch.length=0.15
eta.jc=0.5
eta.tn=0.5

p1 = doit(nsites=nsites,branch.length=branch.length,eta.jc=eta.jc, eta.tn=eta.tn)
plot(p1)

who = "birds"
d=read.dna("../datasets/birds4-clean.phy") #needs to be 4 taxa
seq1 = as.vector(unname(as.character(d[1,])))
seq2 = as.vector(unname(as.character(d[2,])))

## remove missing:
s1 <-seq1!="-"
s2 <- seq2!="-"
seq1 <- seq1[s1&s2&s3]
seq2 <- seq2[s1&s2&s3]
nsites=length(seq2)
nuc = c("a","c","g","t")
out12 = matrix(0,n,n) # counts between 1 and 2
for(i in 1:n)
    for(j in 1:n)
        out12[i,j] = sum(seq1==nuc[i] & seq2==nuc[j])
print(out12)
r = rep(1,6)
Q = optim.gtr(out12,r)

d1x = 0.11 #fixed in true value
d2x = 0.078 #fixed in true value

delta = 0.01
branch.length=0.15
d12=seq(delta, 2*branch.length, by=delta)
log.like2 = gtr.log.like(out12,d12,Q$Q)
log.like2 = log.like2 - max(log.like2)
y2 = exp(log.like2)
y2 = y2/sum(y2)
y2 = y2 / delta
df.gtr = data.frame(x=d12,y=y2)

p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue") +
    geom_vline(xintercept=d1x+d2x)+
    xlab('branch length') +
    ylab('densities')

plot(p1)

##Q = randomQ(4)
## branch.length=0.15
## nsites=1500
## x = simulateSequenceSummary(nsites,Q,branch.length)
## print(x)
## P = matrixExp(Q,branch.length)
## P
## diag(Q$p)
## diag(Q$p) %*% P
## log(diag(Q$p) %*% P)
## x * log(diag(Q$p) %*% P)


## ======================================================================================================================================
## case 2
## A>x----y, bl=t

library(ape)
library(ggplot2)
source("4taxa_cats_dogs_functions.r")
source('branch-length.r')

# loglik ((1,2)x,3)
gtr.log.lik.x = function(d1x,d2x,dxy,seq1.dist,seq2.dist, seq3.dist,Q){
    suma = 0
    for(site in 1:nsites){
        ##print(site)
        pA.given.x = siteLik(d1x,d2x,seq1.dist[,site], seq2.dist[,site], Q$Q)
        ##print(pA.given.x)
        Q$Q$p * pA.given.x
        col=which(seq3.dist[,site] == 1)
        Pxy = matrixExp(Q$Q,dxy)
        ##print(Pxy)
        ##print(Pxy[,col])
        ##print(log(sum(Q$Q$p * pA.given.x * Pxy[,col])))
        suma = suma + log(sum(Q$Q$p * pA.given.x * Pxy[,col]))
    }
    return ( suma )
}

who = "birds"
d=read.dna("../datasets/birds4-clean.phy") #needs to be 4 taxa
seq1 = as.vector(unname(as.character(d[1,])))
seq2 = as.vector(unname(as.character(d[2,])))
seq3 = as.vector(unname(as.character(d[3,])))

## remove missing:
s1 <-seq1!="-"
s2 <- seq2!="-"
s3 <- seq3!="-"
seq1 <- seq1[s1&s2&s3]
seq2 <- seq2[s1&s2&s3]
seq3 <- seq3[s1&s2&s3]

nsites= length(seq3)

##need seq1.dist, seq2.dist, seq3.dist which are matrices
n=4
seq1.dist = matrix(0,n,nsites)
seq2.dist = matrix(0,n,nsites)
seq3.dist = matrix(0,n,nsites)
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
}


nuc = c("a","c","g","t")
out12 = matrix(0,n,n) # counts between 1 and 2
for(i in 1:n)
    for(j in 1:n)
        out12[i,j] = sum(seq1==nuc[i] & seq2==nuc[j])
print(out12)

out13 = matrix(0,n,n) # counts between 1 and 3
for(i in 1:n)
    for(j in 1:n)
        out13[i,j] = sum(seq1==nuc[i] & seq3==nuc[j])
print(out13)

r = rep(1,6)
Q = optim.gtr(out12,r)

d1x = 0.11 #fixed in true value
d2x = 0.078 #fixed in true value
d3y = 0.091
dxy = 0.026
delta = 0.01
branch.length=0.15

v.dx3=seq(delta, 2*branch.length, by=delta)
lik = gtr.log.lik.x(d1x,d2x,v.dx3[1],seq1.dist,seq2.dist,seq3.dist,Q)

loglik = rep(0,length(v.dx3))
for(i in 1:length(v.dx3)){
    loglik[i] = gtr.log.lik.x(d1x,d2x,v.dx3[i],seq1.dist,seq2.dist,seq3.dist,Q)
}
loglik = loglik - max(loglik)
y2 = exp(loglik)
y2 = y2/sum(y2)
y2 = y2 / delta
df.gtr = data.frame(x=v.dx3,y=y2)


p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue") +
    geom_vline(xintercept=dxy+d3y)+
    xlab('branch length') +
    ylab('densities')

plot(p1)


## TN estimate
nsim = 10000
eta.tn = 0.5
tn = simulateBranchLength.tn(nsim,out13,eta.tn)
d.tn = density(tn$t)
df.tn = data.frame(x=d.tn$x-d1x,y=d.tn$y)

p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue") +
    geom_line(aes(x=x,y=y),data=df.tn,color="red",linetype="dashed") +
    geom_vline(xintercept=dxy+d3y)+
    xlab('branch length') +
    ylab('densities') +
    ggtitle('Blue = GTR, Red = TN')

plot(p1)

## see if TN is centered in 0.1, and the other one is centered in 0.02, as it should
