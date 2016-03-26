## r script to do simulations step by step, to understand how
## the likelihood is similar (or not) to the gamma dist of the BL,
## and how the weights behave
## Claudia March 2016

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')
library(ggplot2)
library(weights)

## Case (1,2)---3
## importance weights:
## they look good, they weighted histograms are weighted statistics
## do cover the real values
who="Case (1,2)---3"
d1x0=0.02
d2x0=0.1
d3x0=0.02
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
