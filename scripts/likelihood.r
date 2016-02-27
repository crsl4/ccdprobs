library(ape)
source("branch-length.R")
library(ggplot2)

# --------------------------------- functions --------------------------------------------------------------------------
# 1-----x-----2
# d1x = distance from 1 to parent x, similarly d2x
# seq1.distj = jth column in seq1.dist matrix (for site j), similarly seq2.distj
# Q = estimated matrix of rates
# returns column of site likelihood for all 4 nucleotides
siteLik = function(d1x,d2x,seq1.distj,seq2.distj, Q, verbose=FALSE){
    P1 = matrixExp(Q,d1x)
    P2 = matrixExp(Q,d2x)
    lik = rep(0,4)
    for(i in 1:4){
        lik[i] = P1[i,]%*%seq1.distj * P2[i,]%*%seq2.distj
        if(verbose)
            print(lik)
    }
    return (lik)
}

# loglik ((1,2)x,3)
gtr.log.lik.x = function(d1x,d2x,dxy,seq1.dist,seq2.dist,seq3.dist,Q){
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


# ----------------------------------------------------- data -----------------------------------------------
who = "birds"
d=read.dna("../datasets/birds4-clean.phy") #needs to be 4 taxa

## 1) take three sequences: 1,2 are sisters
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

##need seq1.dist, seq2.dist, seq3.dist as matrices
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

## 2) Compute counts between seq1, seq2
nuc = c("a","c","g","t")
out12 = matrix(0,n,n) # counts between 1 and 2
for(i in 1:n)
    for(j in 1:n)
        out12[i,j] = sum(seq1==nuc[i] & seq2==nuc[j])
print(out12)

## 3) Estimate Q
r = rep(1,6)
Q = optim.gtr(out12,r)

## 4) Compute likelihood for vector of values for internal branch length: dx3
delta = 0.01
branch.length=0.15
dx3=seq(delta, 2*branch.length, by=delta)
d1x = 0.11 #fixed in true value
d2x = 0.078 #fixed in true value
loglik = rep(0,length(dx3))
for(i in 1:length(dx3)){
    loglik[i] = gtr.log.lik.x(d1x,d2x,dx3[i],seq1.dist,seq2.dist,seq3.dist,Q)
}
loglik = loglik - max(loglik)
y2 = exp(loglik)
y2 = y2/sum(y2)
y2 = y2 / delta
df.gtr = data.frame(x=dx3,y=y2)

# true values:
d3y = 0.091
dxy = 0.026

# plot likelihood with true value of dx3
p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue") +
    geom_vline(xintercept=dxy+d3y)+
    xlab('branch length') +
    ylab('densities')
plot(p1)


## 5) Compute a sequence of TN branch lengths for d13
out13 = matrix(0,n,n) # counts between 1 and 3
for(i in 1:n)
    for(j in 1:n)
        out13[i,j] = sum(seq1==nuc[i] & seq3==nuc[j])
print(out13)

nsim = 10000
eta.tn = 0.5
tn = simulateBranchLength.tn(nsim,out13,eta.tn)
d.tn = density(tn$t)
df.tn = data.frame(x=d.tn$x-d1x,y=d.tn$y) # substract d1x

p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue") +
    geom_line(aes(x=x,y=y),data=df.tn,color="red",linetype="dashed") +
    geom_vline(xintercept=dxy+d3y)+
    xlab('branch length') +
    ylab('densities') +
    ggtitle('Blue = GTR, Red = TN')

plot(p1)

pdf("plot.pdf")
dev.off()
