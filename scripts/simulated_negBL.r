## r script to see if bigger internal BL makes a difference in
## negativeBL.r
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
save(data,file="simulatedNegBL_smaller.Rda")

load("simulatedNegBL.Rda")
data=data_big
data=data_small

head(data)
which(data$d1x<0)
which(data$d2x<0)
which(data$d3x<0)
which(data$d3y<0)
which(data$d4y<0)
which(data$dxy<0)
negLines = which(data$dxy<0)
length(negLines) #407

library(ggplot2)
pl <- ggplot(data,aes(x=1:length(dxy),y=dxy))+
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

pdf(paste0("simplex_simSmaller.pdf"))
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

pdf(paste0("simplex_simSmaller2.pdf"))
plot(pl4)
dev.off()

## ------------------------------------------------------
## histograms BL
## ------------------------------------------------------
data <- within(data,var12 <- alpha12/beta12^2)
data <- within(data,var23 <- alpha23/beta23^2)
data <- within(data,var13 <- alpha13/beta13^2)
data <- within(data,var4x <- alpha4x/beta4x^2)
data <- within(data,var34 <- alpha34/beta34^2)
data <- within(data,var3x <- 0.25*(var23+var13+var12))

data <- within(data,cov1x.2x <- 0.25*(var12-var13-var23))
data <- within(data,cov2x.3x <- 0.25*(var23-var12-var13))
data <- within(data,cov1x.3x <- 0.25*(var13-var12-var23))
data <- within(data,cov3y.xy <- 0.25*(var34-var3x-var4x))
data <- within(data,cov4y.xy <- 0.25*(var4x-var34-var3x))
data <- within(data,cov3y.4y <- 0.25*(var3x-var34-var4x))
summary(data)

## for tree1
data1 <- data[!is.na(data$trees),]
summary(data1)

## birds tree1:
d1x = branch.length[2]
d2x = branch.length[3]
d3y = branch.length[4]
d4y = branch.length[5]
dxy = dxy0

d12=d1x+d2x
d13=d3y+dxy+d1x
d23=d3y+dxy+d2x
d34=d3y+d4y
d3x=d3y+dxy
d4x=d4y+dxy

dat=data1
## change file name and main title
pdf("simSmaller_histBL.pdf", height=6,width=15)
layout(m=matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow=2))
hist(dat$d12, main="Simulated (1,2),3,4")
abline(v=d12,col="red")
hist(dat$d13)
abline(v=d13,col="red")
hist(dat$d23)
abline(v=d23,col="red")
hist(dat$d34)
abline(v=d34,col="red")
hist(dat$d3x)
abline(v=d3x,col="red")
hist(dat$d4x)
abline(v=d4x,col="red")
hist(dat$d1x)
abline(v=d1x,col="red")
hist(dat$d2x)
abline(v=d2x,col="red")
hist(dat$d3y)
abline(v=d3y,col="red")
hist(dat$d4y)
abline(v=d4y,col="red")
hist(dat$dxy)
abline(v=dxy,col="red")
dev.off()
