## r script to see how the plot of l(t) and l'(t)
## look like for different settings
## Claudia March 2016

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')

nsites = 1500
t="((1,2),3,4);"
tre=read.tree(text=t)
dxy0=0.03 ##small
dxy0=0.5 ##big
dxy0=0.003 ##smaller
branch.length = c(dxy0,0.11,0.078,0.091,0.098) #bl as birds mb: dxy,d1x,d2x,d3y,d4y
tre$edge.length <- branch.length
nuc <- c('a','c','g','t')

trueQ = randomQ(4,rescale=TRUE)
r=trueQ$r
p=trueQ$p
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

seq1.dist=seqMatrix(seq1)
seq2.dist=seqMatrix(seq2)
seq3.dist=seqMatrix(seq3)
seq4.dist=seqMatrix(seq4)


## -------------------------------------
## 1---2 (d12)
## -------------------------------------
dxy=branch.length[1]
d1x=branch.length[2]
d2x=branch.length[3]
d3y=branch.length[4]
d4y=branch.length[5]
r = rep(1,6)
title = paste("dxy",dxy0)

out12 = countsMatrix(seq1,seq2)
Q = optim.gtr(out12,r)
dist1 = seq1.dist
dist2 = seq2.dist
true_d=d1x+d2x
who="1---2"

out13= countsMatrix(seq1,seq3)
Q = optim.gtr(out13,r)
dist1 = seq1.dist
dist2 = seq3.dist
true_d=d1x+dxy+d3y
who="1---3"

out23= countsMatrix(seq2,seq3)
Q = optim.gtr(out23,r)
dist1 = seq2.dist
dist2 = seq3.dist
true_d=d2x+dxy+d3y
who="2---3"

t=seq(0.01,0.5,by=0.01)
ll=rep(0, length(t))
llprime=rep(0, length(t))
for(i in 1:length(t)){
    l = loglik(dist1,dist2,Q$Q,t[i])
    ll[i] = l$ll
    llprime[i] = l$ll_pr
}



pdf(paste0(who,"_lik_sim.pdf"))
layout(m=matrix(c(1,2),ncol=2))
plot(t,ll)
abline(v=true_d, col="red")
title(who)
plot(t,llprime)
abline(v=true_d, col="red")
abline(h=0, col="red")
title(title)
dev.off()

## -------------------------------------
## (1,2)x,(3,4)y
## -------------------------------------
nsites = 1500
t="((1,2),3,4);"
tre=read.tree(text=t)

dxy0=0.03 ##small
dxy0=0.5 ##big
dxy0=0.003 ##smaller

branch.length = c(dxy0,0.11,0.078,0.091,0.098) #bl as birds mb: dxy,d1x,d2x,d3y,d4y
tre$edge.length <- branch.length
nuc <- c('a','c','g','t')

trueQ = randomQ(4,rescale=TRUE)
r=trueQ$r
p=trueQ$p
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

seq1.dist=seqMatrix(seq1)
seq2.dist=seqMatrix(seq2)
seq3.dist=seqMatrix(seq3)
seq4.dist=seqMatrix(seq4)


dxy=branch.length[1]
d1x=branch.length[2]
d2x=branch.length[3]
d3y=branch.length[4]
d4y=branch.length[5]
Q = trueQ
title = paste("dxy",dxy0)


seqx.dist = sequenceDist(d1x,d2x,seq1.dist, seq2.dist, Q)
seqy.dist = sequenceDist(d3y,d4y,seq3.dist, seq4.dist, Q)

seqx.dist0 = seqMatrix(seqx)
seqy.dist0 = seqMatrix(seqy)


t=seq(0.0001,0.5,by=0.01)
ll=rep(0, length(t))
llprime=rep(0, length(t))
ll0=rep(0, length(t))
llprime0=rep(0, length(t))
for(i in 1:length(t)){
    l = loglik(seqx.dist,seqy.dist,Q,t[i])
    l0 = loglik(seqx.dist0,seqy.dist0,Q,t[i])
    ll[i] = l$ll
    llprime[i] = l$ll_pr
    ll0[i] = l0$ll
    llprime0[i] = l0$ll_pr
}

pdf("lik_plot_intBL.pdf")
layout(m=matrix(c(1,2,3,4),ncol=4))
plot(t,ll)
abline(v=dxy, col="red")
title("est seq")
plot(t,ll0,col="blue")
abline(v=dxy, col="red")
title("true seq")
plot(t,llprime)
abline(v=dxy, col="red")
abline(h=0, col="red")
plot(t,llprime0, col="blue")
abline(v=dxy, col="red")
abline(h=0, col="red")
title(title)
dev.off()


## ------------------------------------------
## lik for intBL for three topologies
## ------------------------------------------
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
seq1.dist=seqMatrix(seq1)
seq2.dist=seqMatrix(seq2)
seq3.dist=seqMatrix(seq3)
seq4.dist=seqMatrix(seq4)

r = c(0.2463,0.1764,0.1231,0.0187,0.4185,0.0170)
p = c(0.2776,0.2937,0.1612,0.2675)
den = r[6]
r = r/den
Q = makeQ(r,p,4, rescale=TRUE)

## tree 1 (1,2),3,4
dxy1 = 0.026
d1x = 0.11
d2x = 0.078
d3y = 0.091
d4y = 0.098
seqx.dist1 = sequenceDist(d1x,d2x,seq1.dist, seq2.dist, Q)
seqy.dist1 = sequenceDist(d3y,d4y,seq3.dist, seq4.dist, Q)

## tree 2 (1,3),2,4
dxy2=0.0213
d1x=0.112
d3x=0.0909
d2y=0.0813
d4y=0.103
seqx.dist2 = sequenceDist(d1x,d3x,seq1.dist, seq3.dist, Q)
seqy.dist2 = sequenceDist(d2y,d4y,seq2.dist, seq4.dist, Q)

## tree 3 (1,4),2,3
dxy3=0.0255
d1x=0.1085
d4x=0.0953
d2y=0.0816
d3y=0.0928
seqx.dist3 = sequenceDist(d1x,d4x,seq1.dist, seq4.dist, Q)
seqy.dist3 = sequenceDist(d2y,d3y,seq2.dist, seq3.dist, Q)

t=seq(0.01,0.5,by=0.01)
ll1=rep(0, length(t))
llprime1=rep(0, length(t))
ll2=rep(0, length(t))
llprime2=rep(0, length(t))
ll3=rep(0, length(t))
llprime3=rep(0, length(t))
for(i in 1:length(t)){
    l1 = loglik(seqx.dist1,seqy.dist1,Q,t[i])
    l2 = loglik(seqx.dist2,seqy.dist2,Q,t[i])
    l3 = loglik(seqx.dist3,seqy.dist3,Q,t[i])
    ll1[i] = l1$ll
    llprime1[i] = l1$ll_pr
    ll2[i] = l2$ll
    llprime2[i] = l2$ll_pr
    ll3[i] = l3$ll
    llprime3[i] = l3$ll_pr
}

pdf("lik_intBL_3trees.pdf")
layout(m=matrix(c(1,2,3,4,5,6),ncol=3))
plot(t,ll1)
abline(v=dxy1, col="red")
title("(1,2),3,4")
plot(t,llprime1)
abline(h=0, col="red")
abline(v=dxy1, col="red")
plot(t,ll2)
abline(v=dxy2, col="red")
title("(1,3),2,4")
plot(t,llprime2)
abline(h=0, col="red")
abline(v=dxy2, col="red")
plot(t,ll3)
abline(v=dxy3, col="red")
title("(1,4),2,3")
plot(t,llprime3)
abline(h=0, col="red")
abline(v=dxy3, col="red")
dev.off()


## ------------------------------------------------------------
## compare the gamma density to the 1---2, 1----3, 2----3
## cases for birds data
## ------------------------------------------------------------
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
seq1.dist=seqMatrix(seq1)
seq2.dist=seqMatrix(seq2)
seq3.dist=seqMatrix(seq3)
seq4.dist=seqMatrix(seq4)

r = rep(1,6)
dxy = 0.026
d1x = 0.11
d2x = 0.078
d3y = 0.091
d4y = 0.098

d1x=0.08141034
d2y=0.09337547
d3y=0.09561165
d4x=0.10845926
dxy=0.02568966
r = c(0.2463,0.1764,0.1231,0.0187,0.4185,0.0170)
p = c(0.2776,0.2937,0.1612,0.2675)
den = r[6]
r = r/den
Q = makeQ(r,p,4, rescale=TRUE)


out = countsMatrix(seq1,seq2)
Q = optim.gtr(out,r)
dist1 = seq1.dist
dist2 = seq2.dist
true_d=d1x+d2x
who="1---2"


out= countsMatrix(seq1,seq3)
Q = optim.gtr(out,r)
dist1 = seq1.dist
dist2 = seq3.dist
true_d=d1x+dxy+d3y
who="1---3"

out= countsMatrix(seq2,seq3)
Q = optim.gtr(out,r)
dist1 = seq2.dist
dist2 = seq3.dist
true_d=d2x+dxy+d3y
who="2---3"

dist1=seq2.dist
dist2=seq3.dist
out=countsMatrix(seq2,seq3)
true_d=d2y+d3y

delta=0.01
t=seq(delta,0.5,by=delta)
ll=rep(0, length(t))
for(i in 1:length(t)){
    l = loglik(dist1,dist2,Q,t[i])
    ll[i] = l$ll
}


logl = ll - max(ll)
y2 = exp(logl)
y2 = y2/sum(y2)
y2 = y2 / delta
df.gtr = data.frame(x=t,y=y2)

eta = 0.5
nsim=10000
jc = simulateBranchLength.jc(nsim=1,out,eta=eta)
t0=jc$t

t.lik = simulateBranchLength.lik(nsim, dist1,dist2,Q,t0,eta=eta)
d.lik = density(t.lik$t)
df.lik = data.frame(x=d.lik$x,y=d.lik$y)

library(ggplot2)
title = paste(who,'Blue = GTR, Red = Gamma')
p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue") +
    geom_line(aes(x=x,y=y),data=df.lik,color="red",linetype="dashed") +
    xlab('branch length') +
    ylab('densities') +
    ggtitle(title)+
    geom_vline(xintercept = true_d)

plot(p1)

pdf(paste0(who,"lik_birds.pdf"))
plot(p1)
dev.off()

