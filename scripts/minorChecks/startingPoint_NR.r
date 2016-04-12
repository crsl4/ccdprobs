## r script to study if JC (or TN) or which is a good
## starting point for the Newton-Raphson of findMLE
## Claudia March 2016


## one branch: 1------2
library(ape)
library(ggplot2)
source('branch-length_lik.r')

delta = 0.001
branch.length = 0.15
nsites = 1500
s = seq(delta,2*branch.length,delta)
trueQ = randomQ(4,rescale=TRUE)
print(round(trueQ$Q,4))
print(round(trueQ$p,4))
print(min(diag(trueQ$Q)))
nuc <- c('a','c','g','t')
# simulate seq1
seq1 = sample(nuc,size=nsites,prob=trueQ$p,replace=TRUE)

#simulate seq2
s=branch.length
P = matrixExp(trueQ,s)
seq2 = numeric(nsites)
for ( i in 1:nsites )
    seq2[i] = sample(nuc,size=1,prob=P[which(nuc==seq1[i]),])

x = countsMatrix(seq1,seq2)

Q = optim.gtr(x,trueQ$r)
print(Q$branch.length)
delta=0.001
s = seq(delta,2*branch.length,delta)
log.like2 = gtr.log.like(x,s,Q$Q)
log.like2 = log.like2 - max(log.like2)
y2 = exp(log.like2)
y2 = y2/sum(y2)
y2 = y2 / delta
df.gtr = data.frame(x=s,y=y2)

nsim=10000
seq1.dist = seqMatrix(seq1)
seq2.dist = seqMatrix(seq2)
w = simulateBranchLength.lik(nsim,seq1.dist,seq2.dist,Q,t0=0.1, eta=0.5)
d.w = density(w$t)
df.w = data.frame(x=d.w$x,y=d.w$y)

t.jc = simulateBranchLength.jc(nsim=1,x,eta=0.5)
w2 = simulateBranchLength.lik(nsim,seq1.dist,seq2.dist,Q,t0=t.jc$t, eta=0.5)
d.w2 = density(w2$t)
df.w2 = data.frame(x=d.w2$x,y=d.w2$y)

p1 = ggplot(df.gtr, aes(x=x,y=y)) +
    geom_line(color="blue")+
    geom_line(aes(x=x,y=y),data=df.w,color="red",linetype="dashed") +
    geom_line(aes(x=x,y=y),data=df.w2,color="darkgreen",linetype="dashed") +
    xlab('branch length') +
    ylab('densities') +
    ggtitle('Blue = GTR, Red = Simulated (t0=true), Green= Simulated (t0=JC), eta=0.5')

plot(p1)
