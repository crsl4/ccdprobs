## r script to study the cases of negative BL
## Claudia March 2016
## birds:
## negBL.Rda: without choosing the first cherry with 1, estimating Q with out12, N-R old
## negBL2.Rda: choosing the first cherry with 1 always, estimating Q with out12, N-R old
## negBL3.Rda: without choosing the 1st cherry with 1, fixing Q to true matrix, better N-R

## then fix a tree as well: birds and cats
## then, compute pairwise cov and see if they are negative as expected
## what effect does eta have on the negative BL? (later what effect it has on the weights)

## ---------------------------------------------------
## study when do negative BL appear in the birds data
## sampling tree and BL (later sampling only BL)
## ---------------------------------------------------
library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')

who = "cats"
d=read.dna("../datasets/4taxa-cats.phy") #needs to be 4 taxa
dat.tre=read.table("../datasets/4taxa-cats_ccdprobs.out", header=FALSE)
r=c(2.815,51.982,1.903,1.275,65.402,1.000)
p=c(0.2590,0.2379,0.1900,0.3131);

who = "birds"
d=read.dna("../datasets/birds4-clean.phy") #needs to be 4 taxa
dat.tre=read.table("../datasets/birds4-clean_ccdprobs.out", header=FALSE)
r = c(0.2463,0.1764,0.1231,0.0187,0.4185,0.0170)
p = c(0.2776,0.2937,0.1612,0.2675)


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

for(i in 1:nreps){
    print(i)
    t=sampleTopQuartet(dat.tre)
    b=sampleBLQuartet_details(d,t$tre, trueT0=FALSE, verbose=FALSE, Q=Q)
    originalBL[i,] = c(b$d12$t, b$d23$t, b$d13$t, b$d4x$t, b$d34$t)
    alpha[i,] = c(b$d12$alpha, b$d23$alpha, b$d13$alpha, b$d4x$alpha, b$d34$alpha)
    beta[i,] = c(b$d12$beta, b$d23$beta, b$d13$beta, b$d4x$beta, b$d34$beta)
    branch.lengths[i,]=b$bl
    seq[i,] = b$seq
    trees[i] = write.tree(t$tre)
}

data = data.frame(trees,seq,originalBL,alpha,beta,branch.lengths)
colnames(data) <- c("trees", "seq1","seq2","seq3","seq4","d12", "d23", "d13", "d4x", "d34",
                    "alpha12", "alpha23", "alpha13", "alpha4x", "alpha34",
                    "beta12", "beta23", "beta13", "beta4x", "beta34",
                    "d1x","d2x","d3x","d3y","d4y","dxy")
head(data)
save(data,file=paste0("negBL_",who,".Rda"))

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

pl2 <- ggplot(data[negLines,], aes(trees))+
    geom_bar()
plot(pl2)


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

pdf("simplex_cats.pdf")
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


library(plot3D)
scatter3D(data$d3x,data$d4x,data$d34, bty="g",pch=18,
          col.var=as.integer(data$trees), col=c("#1B9E77", "#D95F02", "#7570B3"),
          pch = 18, ticktype = "detailed",
          colkey = list(at = c(2, 3, 4), side = 1,
              addlines = TRUE, length = 0.5, width = 0.5,
              labels = c("12","13","14")),
          phi=0, theta=0)



## ------------------------------
## plot by tree, no difference
## ------------------------------
load("negBL2.Rda")
data0 <- data
dataNeg <- data[negLines,]
head(dataNeg)
data1=dataNeg[dataNeg$trees=="((1,2),3,4);",]
data2=dataNeg[dataNeg$trees=="(1,(2,3),4);",]
data3=dataNeg[dataNeg$trees=="((1,3),2,4);",]

data <- data3
data <- within(data,sum <- d3x+d4x+d34)
data <- within(data,d3x.n <- d3x/sum)
data <- within(data,d4x.n <- d4x/sum)
data <- within(data,d34.n <- d34/sum)

pl.1<- ggtern(data=data,aes(d3x.n,d4x.n,d34.n)) +
    geom_point(aes(fill=trees),shape=21,size=1)+
    geom_Tline(Tintercept=0.5) + geom_Rline(Rintercept=0.5) + geom_Lline(Lintercept=0.5)
plot(pl.1)


## ---------------------------------------------------
## study when do negative BL appear in the birds data
## sampling only BL
## ---------------------------------------------------

library(ape)
source('branch-length_lik.r')
source('4taxa_functions.r')

who = "cats"
d=read.dna("../datasets/4taxa-cats.phy") #needs to be 4 taxa
dat.tre=read.table("../datasets/4taxa-cats_ccdprobs.out", header=FALSE)
r=c(2.815,51.982,1.903,1.275,65.402,1.000)
p=c(0.2590,0.2379,0.1900,0.3131);

who = "birds"
d=read.dna("../datasets/birds4-clean.phy") #needs to be 4 taxa
dat.tre=read.table("../datasets/birds4-clean_ccdprobs.out", header=FALSE)
r = c(0.2463,0.1764,0.1231,0.0187,0.4185,0.0170)
p = c(0.2776,0.2937,0.1612,0.2675)


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
tree="((1,3),2,4);"
tree="((1,4),2,3);"

for(i in 1:nreps){
    print(i)
    tre=read.tree(text=tree)
    b=sampleBLQuartet_details(d,tre, trueT0=FALSE, Q=Q)
    originalBL[i,] = c(b$d12$t, b$d23$t, b$d13$t, b$d4x$t, b$d34$t)
    alpha[i,] = c(b$d12$alpha, b$d23$alpha, b$d13$alpha, b$d4x$alpha, b$d34$alpha)
    beta[i,] = c(b$d12$beta, b$d23$beta, b$d13$beta, b$d4x$beta, b$d34$beta)
    branch.lengths[i,]=b$bl
    seq[i,] = b$seq
    trees[i] = write.tree(tre)
}

data = data.frame(trees,seq,originalBL,alpha,beta,branch.lengths)
colnames(data) <- c("trees", "seq1","seq2","seq3","seq4","d12", "d23", "d13", "d4x", "d34",
                    "alpha12", "alpha23", "alpha13", "alpha4x", "alpha34",
                    "beta12", "beta23", "beta13", "beta4x", "beta34",
                    "d1x","d2x","d3x","d3y","d4y","dxy")
head(data)
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


## ----------------------------------------------------------------
## covariances of branch lengths
## histograms of BL
## ---------------------------------------------------------------
load("negBL_birdstree1.Rda") # MLE tree birds
load("negBL_birdstree2.Rda") # tree birds
load("negBL_birdstree3.Rda") # NJ tree birds

load("negBL_catstree1.Rda") # MLE&NJ tree cats
load("negBL_catstree2.Rda") # tree cats
load("negBL_catstree3.Rda") # tree cats
head(data)

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
data1 <- subset(data,seq1+seq2==3) ## d12=dist 1,2
data2 <- subset(data,seq1+seq2==7) ## d12=dist 3,4
summary(data1) #to check that seq3 matches across
summary(data2)

## for tree2
data1 <- subset(data,seq1+seq2==4) ## d12=dist 1,3
data2 <- subset(data,seq1+seq2==6) ## d12=dist 2,4
summary(data1) #to check that seq3 matches across
summary(data2)

## for tree3
data1 <- subset(data,seq1==1) ## d12=dist 1,4
data2 <- subset(data,seq1==2) ## d12=dist 2,3
summary(data1) #to check that seq3 matches across
summary(data2)

## birds tree1:
d1x = 0.111
d2x = 0.078
d3y = 0.091
d4y = 0.098
dxy = 0.026
## birds tree2:
d1x = 0.112
d2x = 0.081
d3y = 0.091
d4y = 0.103
dxy = 0.0213
## birds tree3:
d1x = 0.110
d2x = 0.082
d3y = 0.093
d4y = 0.091
dxy = 0.0255

## cats tree1
d1x = 0.087
d2x = 0.0796
d3y = 0.0482
d4y = 0.0625
dxy = 0.0203
## cats tree2
d1x = 0.0937
d2x = 0.0853
d3y = 0.0435
d4y = 0.0563
dxy = 0.0193
## cats tree3
d1x = 0.105
d2x = 0.0881
d3y = 0.0503
d4y = 0.0617
dxy = 0.0186


d12=d1x+d2x
d13=d3y+dxy+d1x
d23=d3y+dxy+d2x
d34=d3y+d4y
d3x=d3y+dxy
d4x=d4y+dxy

dat=data1
## change file name and main title
pdf("cats_histBL_tree3.pdf", height=6,width=15)
layout(m=matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow=2))
hist(dat$d12, main="Cats (1,4),2,3")
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


## NJ estimate:
## library(ape)
## d=read.dna("~/Documents/phylo/projects/CladeCondProb/ccdprobs/datasets/birds4-clean.phy")
## ds=dist.dna(d)
## tr=nj(ds)
## plot(tr)
## tr$edge.length
## #0.087616636 0.082393027 0.003362989 0.099538427 0.092048942
## write.tree(tr)
## #"(owl:0.08761663563,penguin:0.08239302668,(kiwi:0.09953842733,parrot:0.09204894184):0.003362989414);"
