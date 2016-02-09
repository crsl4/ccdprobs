
library(ape)
source('branch-length.r')

d=read.dna("../datasets/4-taxon-cats.phy") #needs to be 4 taxa
tre=read.tree() #need tree for 4 taxa: run seq2ccdprobs, but then how to sample a tree? 4-taxon-cats_ccdprobs.out

# function to sample branch lengths for 4-taxon tree
sampleBL = function(d,tre,nsim=1, eta=0.5){
    thr = sample(1:4, size=3, replace=FALSE) # need to choose three taxa
    fth = setdiff(1:4,thr)
    seq1 = as.vector(unname(as.character(d[thr[1],])))
    seq2 = as.vector(unname(as.character(d[thr[2],])))
    seq3 = as.vector(unname(as.character(d[thr[3],])))
    seq4 = as.vector(unname(as.character(d[fth,])))
    nsites = length(seq1)
    if(length(seq2) != nsites && length(seq3) != nsites){
        error()
    }

    nuc = c("a","c","g","t")
    n = 4

    out = matrix(0,n,n) # distance between 1 and 2
    for(i in 1:n)
        for(j in 1:n)
            out[i,j] = sum(seq1==nuc[i] & seq2==nuc[j])
    print(out)
    # want to check if out is good for tn?
    # need to determine eta, and jc vs tn?
    d12 = simulateBranchLength.jc(nsim=nsim,out,eta=eta)
    #d12 = simulateBranchLength.tn(nsim=nsim, out, eta=eta)

    out = matrix(0,n,n)     # distance between 1 and 3
    for(i in 1:n)
        for(j in 1:n)
            out[i,j] = sum(seq1==nuc[i] & seq3==nuc[j])
    print(out)
    # want to check if out is good for tn?
    # need to determine eta, and jc vs tn?
    d13 = simulateBranchLength.jc(nsim=nsim,out,eta=eta)
    #d13 = simulateBranchLength.tn(nsim=nsim, out, eta=eta)

    out = matrix(0,n,n)     # distance between 2 and 3
    for(i in 1:n)
        for(j in 1:n)
            out[i,j] = sum(seq2==nuc[i] & seq3==nuc[j])
    print(out)
    # want to check if out is good for tn?
    # need to determine eta, and jc vs tn?
    d23 = simulateBranchLength.jc(nsim=nsim,out,eta=eta)
    #d23 = simulateBranchLength.tn(nsim=nsim, out, eta=eta)

    d1x = (d12+d13-d23)/2
    d2x = (d12+d23-d13)/2

    print(d1x)
    print(d2x)

    # need to update d1x, d2x in tre


    # need seq distribution at x: felsenstein algorithm
    seqx = matrix(rep(0.25,nsites*4),nrow=4) #temporary for next code

    # need to get d34,d3x,d4x to solve for all missing bl

    out = matrix(0,n,n) # distance between 3 and 4
    for(i in 1:n)
        for(j in 1:n)
            out[i,j] = sum(seq3==nuc[i] & seq4==nuc[j])
    print(out)
    # want to check if out is good for tn?
    # need to determine eta, and jc vs tn?
    d34 = simulateBranchLength.jc(nsim=nsim,out,eta=eta)
    #d34 = simulateBranchLength.tn(nsim=nsim, out, eta=eta)

    # need seq3.dist and seq4.dist which are matrices
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

    out = matrix(0,n,n) # distance between 3 and x
    for(i in 1:nsites){
        out = out + seq3.dist[,i]%*%t(seqx[,i])
    }
    print(out)
    # want to check if out is good for tn?
    # need to determine eta, and jc vs tn?
    d3x = simulateBranchLength.jc(nsim=nsim,out,eta=eta)
    #d3x = simulateBranchLength.tn(nsim=nsim, out, eta=eta)
    out = matrix(0,n,n) # distance between 4 and x
    for(i in 1:nsites){
        out = out + seq4.dist[,i]%*%t(seqx[,i])
    }
    print(out)
    # want to check if out is good for tn?
    # need to determine eta, and jc vs tn?
    d4x = simulateBranchLength.jc(nsim=nsim,out,eta=eta)
    #d4x = simulateBranchLength.tn(nsim=nsim, out, eta=eta)

    d3y = (d34+d3x-d4x)/2
    d4y = (d34+d4x-d3x)/2
    dxy = (d3x+d4x-d34)/2

                                        # need to update new bl in tre
}
