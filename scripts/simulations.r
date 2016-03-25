## r script to do simulations step by step, to understand how
## the likelihood is similar (or not) to the gamma dist of the BL,
## and how the weights behave
## Claudia March 2016

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')
library(ggplot2)


## Case (1,2)---3

who="Case (1,2)---3"
d1x=0.15
d2x=0.15
d3x=0.15
eta = 0.5
nsites=15000

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
## simulate seq3
P = matrixExp(Q,d3x)
seq3 = numeric(nsites)
for ( i in 1:nsites )
    seq3[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

seq1.dist = seqMatrix(seq1)
seq2.dist = seqMatrix(seq2)
seq3.dist = seqMatrix(seq3)


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

out12 = countsMatrix(seq1,seq2)
out13 = countsMatrix(seq1,seq3)
out23 = countsMatrix(seq2,seq3)

jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
jc13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
jc23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)

nsim=10000
t.lik12 = simulateBranchLength.lik(nsim, seq1.dist,seq2.dist,Q,t0=jc12$t,eta=eta)
t.lik13 = simulateBranchLength.lik(nsim, seq1.dist,seq3.dist,Q,t0=jc13$t,eta=eta)
t.lik23 = simulateBranchLength.lik(nsim, seq2.dist,seq3.dist,Q,t0=jc23$t,eta=eta)
d.lik12 = density(t.lik12$t)
d.lik13 = density(t.lik13$t)
d.lik23 = density(t.lik23$t)
df.lik = data.frame(x12=d.lik12$x,y12=d.lik12$y,x13=d.lik13$x,y13=d.lik13$y,x23=d.lik23$x,y23=d.lik23$y)
head(df.lik)

title = paste(who,'gamma density \n Blue = d12, Red = d13, Gold = d23', 'nsites', nsites)
p1 = ggplot(df.lik, aes(x=x12,y=y12)) +
    geom_line(color="blue")+
    geom_vline(xintercept = d1x+d2x, colour="blue")+
    geom_vline(xintercept = t.lik12$alpha/t.lik12$beta, colour="blue", linetype="dashed")+
    geom_line(aes(x=x13,y=y13),data=df.lik,colour="red")+
    geom_vline(xintercept = d1x+d3x, colour="red")+
    geom_vline(xintercept = t.lik13$alpha/t.lik13$beta, colour="red", linetype="dashed")+
    geom_line(aes(x=x23,y=y23),data=df.lik,colour="gold")+
    geom_vline(xintercept = d2x+d3x, colour="gold")+
    geom_vline(xintercept = t.lik23$alpha/t.lik23$beta, colour="gold", linetype="dashed")+
    ggtitle(title)

plot(p1)

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



## Case 1-----2
## need to study the effect ot t0, nsites, eta
## eta: 0.5 covers better the lik
## nsties, as expected, gives lower variance
## aroung t_MLE as nsites increases
## also, both lik and gamma centered in t0 (good)

who="Case 1---2"
t0=0.15
eta = 0.5
nsites=1500

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
