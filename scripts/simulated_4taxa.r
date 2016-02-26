# r script to use a simulated data to compare bl summaries between mrbayes and importance sampling
# Claudia February 2016


## to do: modify below to simulate sequence at x, then at 1,2,y,3,4
## not uncertainty in topology, ccdprobs=1,0,0
## for a given tree, fix branch lengths, and simulate data
## then input data to compute branch lengths summaries, many many times
## with true sequence at x, and with estimated seq at x
## compare results
## also compare true Q with estimated Q

## this will let us know if the estimation of seqx is causing the problem

library(ape)
source('branch-length.r')
source('4taxa_simulated_functions.r')

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



b=sampleBLQuartetSim(seq1,seq2,seq3,seq4,seqx, verbose=TRUE,estSeqx=FALSE)
print(b)

#now, compute true likelihood
dxy=branch.length[1]
d1x=branch.length[2]
d2x=branch.length[3]
d3y=branch.length[4]
d4y=branch.length[5]


n=4
#need seq1.dist and seq2.dist which are matrices
seq1.dist = matrix(0,n,nsites)
seq2.dist = matrix(0,n,nsites)
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
}
#need seq3.dist and seq4.dist which are matrices
seq3.dist = matrix(0,n,nsites)
seq4.dist = matrix(0,n,nsites)
for(i in 1:nsites){
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

suma = 0
for(s in 1:nsites){
    lik12 = siteLik(d1x,d2x,seq1.dist[,s],seq2.dist[,s],trueQ)
    lik34 = siteLik(d3y,d4y,seq3.dist[,s],seq4.dist[,s],trueQ)
    L = lik12 %*% t(lik34)
    Pxy = matrixExp(trueQ,dxy)
    L2 = L*Pxy
    Lik = trueQ$p * L2
    suma = suma+log(sum(Lik))
}
print(suma)


logw = b$logprior+b$loglik-b$logdensity
print(logw)

nreps = 1000
logwv = rep(0,nreps)
branch.lengths = matrix(0,nreps,5)
err = 0
for(i in 1:nreps){
    b=try(sampleBLQuartetSim(seq1,seq2,seq3,seq4,seqx))
    if(class(b) == "try-error"){
        err = err + 1
    } else {
        branch.lengths[i,]=b$bl
        logw = b$logprior+b$loglik-b$logdensity
        #print(logw)
        logwv[i] = logw
    }
}
err #734
data = data.frame(logwv,branch.lengths)
who = rowSums(data) != 0
data = data[who,]
head(data)
summary(data)
my.logw = logwv[who] - mean(logwv[who])
data$w = exp(my.logw)/sum(exp(my.logw))
save(data,file="data_sim2.Rda")
#load("data_birds.Rda")

# with true seqx
head(data)
data[data$w>0.01,] #only one!

meanX1 = weighted.mean(data$X1,data$w)
meanX2 = weighted.mean(data$X2,data$w)
meanX3 = weighted.mean(data$X3,data$w)
meanX4 = weighted.mean(data$X4,data$w)
meanX5 = weighted.mean(data$X5,data$w)

q1.X1 = weighted.quantile(data$X1,data$w, probs=0.025)
q1.X2 = weighted.quantile(data$X2,data$w, probs=0.025)
q1.X3 = weighted.quantile(data$X3,data$w, probs=0.025)
q1.X4 = weighted.quantile(data$X4,data$w, probs=0.025)
q1.X5 = weighted.quantile(data$X5,data$w, probs=0.025)

q2.X1 = weighted.quantile(data$X1,data$w, probs=0.975)
q2.X2 = weighted.quantile(data$X2,data$w, probs=0.975)
q2.X3 = weighted.quantile(data$X3,data$w, probs=0.975)
q2.X4 = weighted.quantile(data$X4,data$w, probs=0.975)
q2.X5 = weighted.quantile(data$X5,data$w, probs=0.975)

branch.length
c(meanX1,meanX2,meanX3,meanX4,meanX5)
c(q1.X1, q1.X2, q1.X3, q1.X4, q1.X5)
c(q2.X1, q2.X2, q2.X3, q2.X4, q2.X5)

