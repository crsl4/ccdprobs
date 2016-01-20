# R script to do bootstrap of sites and then neighbor joining to get a tree
# Claudia January 2016
# using 60_concat.in (from HGTinconsistency) as example
# R CMD BATCH bootstrap2NJ.r seqFile treeName
# generates tree file named treeName in datasets

if(length(args) == 0){
    file = "../datasets/60_concat.in"
    outfile = "../datasets/60_concat_bt.out"
}else{
    file = args[1]
    outfile = args[2]
}

# scenario 1: bootstrap + nj assuming JC69 model
library(ape)
data <- read.dna(file)
d <- dist.dna(data)
tr <- nj(d)
plot(tr)

# returns bootstrap datasets
bootstrapDNA <- function(x,B=100){
    # function based on ape boot.phylo
    boot.samp <- vector("list", B)
    y <- nc <- ncol(x)
    for (i in 1:B){
        index <- unlist(sample(y, replace = TRUE))
        boot.samp[[i]] <- x[, index]
    }
    ans <- boot.samp
    ans
}

# returns list of bootstrap trees obtained with nj
bootstrapNJTree <- function(x,model="JC69",B=100){
    # function based on ape boot.phylo
    boot.tree <- vector("list", B)
    y <- nc <- ncol(x)
    for (i in 1:B){
        boot.samp <- unlist(sample(y, replace = TRUE))
        boot.tree[[i]] <- nj(dist.dna(x[, boot.samp],model=model))
    }
    ans <- boot.tree
    ans
}


bd <- bootstrapDNA(data)
bt <- bootstrapNJTree(data)
plot(bt[[1]])


for(i in 1:100){
    write.tree(bt[[i]],file=outfile, append=TRUE)
}



# january 11-15 ------------------------------------------
if(false){
library(phangorn)
# read in data
data <- read.phyDat("60_concat.in", format="phylip", type="DNA")
# construct rooted tree with UPGMA or unrooted tree with NJ:
# Two options: dist.dna (ape) and dist.ml (phangorn)

# 1)
data2 <- as.DNAbin(data) #str(data) => class phyDat, and dist.dna needs DNAbin
# distance matrix using a model of DNA (ape)
dm2 <- dist.dna(data2, model="K80") #fixit: not sure how to choose GTR+G
# construct trees
treeUPGMA2 <- upgma(dm2)
treeNJ2 <- NJ(dm2)
# plot trees
layout(matrix(c(1,2),2,1),height=c(1,2))
par(mar=c(0,0,2,0)+0.1)
plot(treeUPGMA2, main="UPGMA")
plot(treeNJ2,"unrooted", main="NJ")

# 2)
dm <- dist.ml(data) #it says mostly for aminoacids in pdf, only F81 and JC69 for dna
treeUPGMA <- upgma(dm)
treeNJ <- NJ(dm)
# plot trees
layout(matrix(c(1,2),2,1),height=c(1,2))
par(mar=c(0,0,2,0)+0.1)
plot(treeUPGMA, main="UPGMA")
plot(treeNJ,"unrooted", main="NJ")


# MLE (for GTR+G)
fit <- pml(treeNJ,data=data)
fitGTR <- update(fit,k=4,inv=0.2) #GTR+G(4)+I
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,rearrangement = "NNI", control = pml.control(trace = 0))
# can we get the dist matrix from fitGTR? not clear in documentation
str(fitGTR)

# ---- bootstrap -----
# bootstrap.phyDat produces a list of bootstrapped dataset, but:
# 1) how to extract?
# 2) how to estimate trees the way we want?
bs <- bootstrap.phyDat(data,pratchet,bs=100) #fixt: pratchet=parsimony search to find best tree, NJ for neighbor joining
str(bs)
bs[[1]] #tree
str(bs[[1]]) #does not have the data

for(i in 1:100){
    write.tree(bs[[i]],file="trees.txt", append=TRUE)
}


# how to read and run with ccdprobs?
# need a new mbsum to read list of trees
}
