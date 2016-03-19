# example using functions in internalBranch.r
# Claudia March 2016

library(ape)
source("internalBranch.r")

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
d1x = 0.11
d2x = 0.078
d3y = 0.091
d4y = 0.098
dxy = 0.026


y = findMLE(seq1,seq2,seq3,seq4,d1x,d2x,d3y,d4y, y0=0.999) ## only works if y0>true value=0.96
est_dxy=-0.75*log(y[length(y)])
## t=exp(-4/3dxy)

est_var=1/obsInfo(seq1,seq2,seq3,seq4,d1x,d2x,d3y,d4y,y[length(y)])

mu=est_dxy
sigma=sqrt(est_var)
x=rnorm(10000,mu,sigma)
hist(x)

## ---------------------------------------------------
## compare likelihood to JC simulated branch length
library(ape)
library(ggplot2)
source("4taxa_functions.r")
source('branch-length.r')
source("internalBranch.r")

# loglik ((1,2)x,(3,4)y)
gtr.log.lik.all = function(d1x,d2x,dxy,d3y,d4y,seq1.dist,seq2.dist, seq3.dist,seq4.dist,Q){
    suma = 0
    for(s in 1:ncol(seq1.dist)){
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

out12 = countsMatrix(seq1,seq2)
r = rep(1,6)
Q = optim.gtr(out12,r)


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

eta = 0.2
nsim=10000
jc = simulateInternalBranchLength.jc.lik(nsim, seq1,seq2,seq3,seq4, d1x,d2x,d3y,d4y, eta=eta)
## density estimate
d.jc = density(jc)
## another data frame
df.jc = data.frame(x=d.jc$x,y=d.jc$y)

p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue") +
    geom_vline(xintercept=dxy)+
    xlab('branch length') + ylab('densities')+
    geom_line(aes(x=x,y=y),data=df.jc,color="red",linetype="dashed") +
    ggtitle('Blue = GTR, Red = JC, eta=0.2')

pdf("JC_lik.pdf")
plot(p1)
dev.off()

## when newton raphson was not working:
t = findMLE2(seq1,seq2,seq3,seq4,d1x,d2x,d3y,d4y) ## this takes forever!
#save(t,file="findMLE2.Rda")
pdf("plotLogLik.pdf")
plot(t$y,t$f)
abline(h=0)
dev.off()
t$f[t$f>0]
t$f[t$f<0]
t$f[968] ## 15.33
t$f[969] ## -6.5
t$y[968] ## 0.967
t$y[969] ## 0.968
## true value should be between 0.967, 0.968

## exp(-4*dxy/3) ## 0.9659, so above good approx!
## yest = t$y[968]
## dxy_est=log(yest)*(-3/4)
## dxy_est ## 0.02516
## dxy ## 0.026

