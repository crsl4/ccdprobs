## r script to see if bigger internal BL makes a difference in
## negativeBL.r
## Claudia March 2016

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')

nsites = 1500
t="((1,2),3,4);"
tre=read.tree(text=t)
dxy0=0.03
dxy0=0.5
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


nreps = 1000
trees = rep(NA,nreps)
logwv = rep(0,nreps)
branch.lengths = matrix(0,nreps,6)
originalBL = matrix(0,nreps,5)
alpha = matrix(0,nreps,5)
beta = matrix(0,nreps,5)
seq = matrix(0,nreps,4)
den = r[6]
r = r/den
Q = makeQ(r,p,4, rescale=TRUE)
tree="((1,2),3,4);"
##tree="((1,3),2,4);"
##tree="((1,4),2,3);"
err=0
for(i in 1:nreps){
    print(i)
    tre=read.tree(text=tree)
    b=try(sampleBLQuartet_details_sim(seq1,seq2,seq3,seq4, trueT0=FALSE, Q=Q))
    if(class(b) == "try-error"){
        err = err+1
    } else{
        originalBL[i,] = c(b$d12$t, b$d23$t, b$d13$t, b$d4x$t, b$d34$t)
        alpha[i,] = c(b$d12$alpha, b$d23$alpha, b$d13$alpha, b$d4x$alpha, b$d34$alpha)
        beta[i,] = c(b$d12$beta, b$d23$beta, b$d13$beta, b$d4x$beta, b$d34$beta)
        branch.lengths[i,]=b$bl
        seq[i,] = b$seq
        trees[i] = write.tree(tre)
    }
}


data = data.frame(trees,seq,originalBL,alpha,beta,branch.lengths)
colnames(data) <- c("trees", "seq1","seq2","seq3","seq4","d12", "d23", "d13", "d4x", "d34",
                    "alpha12", "alpha23", "alpha13", "alpha4x", "alpha34",
                    "beta12", "beta23", "beta13", "beta4x", "beta34",
                    "d1x","d2x","d3x","d3y","d4y","dxy")
head(data)
data_big=data
err_big=err
data_small=data
err_small=err

save(data_big, err_big, data_small, err_small, file="simulatedNegBL.Rda")

## todo: finish this part, get to simplex, compare
## later, go to cov section, and plot histograms, are they different to birds?
## then, modify dxy to get different simplex and hist to compare


save(data,file=paste0("negBL_",who,"tree3.Rda"))
which(data$d1x<0)
which(data$d2x<0)
which(data$d3x<0)
which(data$d3y<0)
which(data$d4y<0)
which(data$dxy<0) #only dxy is negative as expected
negLines = which(data$dxy<0)
length(negLines) #407

library(ggplot2)
pl <- ggplot(data,aes(x=1:nreps,y=dxy))+
    geom_point(aes(color=trees))
plot(pl)


## --------------------------------------
## simplex between d3x,d4x,d34
## --------------------------------------
data <- within(data,sum <- d3x+d4x+d34)
data <- within(data,d3x.n <- d3x/sum)
data <- within(data,d4x.n <- d4x/sum)
data <- within(data,d34.n <- d34/sum)

#data <- data[negLines,]
library(ggtern)
newdat  <- data.frame(
    trees=data$trees,
    d3x = data$d3x.n,
    d4x = data$d4x.n,
    d34 = data$d34.n
)
pl3 <- ggtern(data=newdat,aes(d3x,d4x,d34)) +
    geom_point(aes(fill=trees),shape=21,size=1)+
    geom_Tline(Tintercept=0.5) + geom_Rline(Rintercept=0.5) + geom_Lline(Lintercept=0.5)
plot(pl3)

pdf(paste0("simplex_",who,"_tree3.pdf"))
plot(pl3)
dev.off()

## --------------------------------------
## simplex between d12,d23,d13
## --------------------------------------
data <- within(data,sum <- d12+d13+d23)
data <- within(data,d12.n <- d12/sum)
data <- within(data,d23.n <- d23/sum)
data <- within(data,d13.n <- d13/sum)
newdat2  <- data.frame(
    trees=data$trees,
    d12 = data$d12.n,
    d23 = data$d23.n,
    d13 = data$d13.n
)
pl4 <- ggtern(data=newdat2,aes(d12,d23,d13)) +
    geom_point(aes(fill=trees),shape=21,size=1)+
    geom_Tline(Tintercept=0.5) + geom_Rline(Rintercept=0.5) + geom_Lline(Lintercept=0.5)
plot(pl4)
