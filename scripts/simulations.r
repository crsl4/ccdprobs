## r script to do simulations step by step, to understand how
## the likelihood is similar (or not) to the gamma dist of the BL,
## and how the weights behave
## Claudia March 2016

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')
library(ggplot2)

## Case 1-----2
## need to study the effect ot t0, nsites, eta
## eta: 0.5 covers better the lik
## nsties, as expected, gives lower variance
## aroung t_MLE as nsites increases
## also, both lik and gamma centered in t0 (good)

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
