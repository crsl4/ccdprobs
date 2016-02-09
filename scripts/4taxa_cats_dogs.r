# list of doubts:
# 1) JC vs TN? eta=0.5
# 2) seq dist at x: how to estimate Q? which counts? the siteLik function is merely the lik, we want P(x|1,2)
# 3) we have ccdprobs, how to sample tree? from which output file?

library(ape)
source('branch-length.r')

d=read.dna("../datasets/4taxa-cats.phy") #needs to be 4 taxa
tre=read.tree("../datasets/4taxa-cats.tre") # fixit: how to sample a tree? from 4-taxon-cats_ccdprobs.out?
nsim=1
eta=0.5

# function to sample branch lengths for 4-taxon tree
sampleBL = function(d,tre,nsim=1, eta=0.5){
    thr = sample(1:4, size=3, replace=FALSE) # need to choose three taxa
    fth = setdiff(1:4,thr)
    seq1 = as.vector(unname(as.character(d[thr[1],]))) # fixit: need to make sure that 1,2 are sisters
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
    # want to check if out is good for tn? fixit
    # need to determine eta, and jc vs tn?
    d12 = simulateBranchLength.jc(nsim=nsim,out,eta=eta)
    d12 = simulateBranchLength.tn(nsim=nsim, out, eta=eta)
    print(d12)

    out = matrix(0,n,n)     # distance between 1 and 3
    for(i in 1:n)
        for(j in 1:n)
            out[i,j] = sum(seq1==nuc[i] & seq3==nuc[j])
    print(out)

    d13 = simulateBranchLength.jc(nsim=nsim,out,eta=eta)
    d13 = simulateBranchLength.tn(nsim=nsim, out, eta=eta)
    print(d13)

    out = matrix(0,n,n)     # distance between 2 and 3
    for(i in 1:n)
        for(j in 1:n)
            out[i,j] = sum(seq2==nuc[i] & seq3==nuc[j])
    print(out)

    d23 = simulateBranchLength.jc(nsim=nsim,out,eta=eta)
    d23 = simulateBranchLength.tn(nsim=nsim, out, eta=eta)
    print(d23)

    d1x = (d12+d13-d23)/2
    d2x = (d12+d23-d13)/2

    print(d1x)
    print(d2x)


    # need seq distribution at x: felsenstein algorithm
    #Q = optim.gtr(x,r) #fixit: x=matrix of counts, which counts when there are 4 sequences?
    Q = randomQ(4) #fixit: for the moment
    print(Q$Q)

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
    #print(seq3.dist)
    #print(seq4.dist)

    seqx = sequenceDist(d1x,d2x, seq1.dist, seq2.dist,Q) #fixit

    out = matrix(0,n,n) # distance between 3 and 4
    for(i in 1:n)
        for(j in 1:n)
            out[i,j] = sum(seq3==nuc[i] & seq4==nuc[j])
    print(out)

    d34 = simulateBranchLength.jc(nsim=nsim,out,eta=eta)
    d34 = simulateBranchLength.tn(nsim=nsim, out, eta=eta)
    print(d34)


    out = matrix(0,n,n) # distance between 3 and x
    for(i in 1:nsites){
        out = out + seq3.dist[,i]%*%t(seqx[,i])
    }
    print(out)

    d3x = simulateBranchLength.jc(nsim=nsim,out,eta=eta)
    d3x = simulateBranchLength.tn(nsim=nsim, out, eta=eta)
    print(d3x)

    out = matrix(0,n,n) # distance between 4 and x
    for(i in 1:nsites){
        out = out + seq4.dist[,i]%*%t(seqx[,i])
    }
    print(out)

    d4x = simulateBranchLength.jc(nsim=nsim,out,eta=eta)
    d4x = simulateBranchLength.tn(nsim=nsim, out, eta=eta)
    print(d4x)

    d3y = (d34+d3x-d4x)/2
    d4y = (d34+d4x-d3x)/2
    dxy = (d3x+d4x-d34)/2

    print(paste("d12",d12,"d13",d13,"d23",d23,"d34",d34))
    print(paste("d1x",d1x,"d2x",d2x,"d3y",d3y,"d4y",d4y,"dxy",dxy))

                                        # need to update new bl in tre
}

# d1x = distance from 1 to parent x, similarly d2x
# seq1.distj = jth column in seq1.dist matrix (for site j), similarly seq2.distj
# Q = matrix of rates
# returns column of site likelihood
siteLik = function(d1x,d2x,seq1.distj,seq2.distj, Q){
    P1 = matrixExp(Q,d1x)
    P2 = matrixExp(Q,d2x)
    lik = rep(0,4)
    for(i in 1:4){
        lik[i] = P1[i,]%*%seq1.distj * P2[i,]%*%seq2.distj
        print(lik)
    }
    return (lik)
}


# estimates the seq dist at x (parent of 1 and 2)
sequenceDist = function(d1x,d2x,seq1.dist,seq2.dist, Q){
    nsites = length(seq1)
    if(length(seq2) != nsites){
        error()
    }
    nuc = c('a','c','g','t')
    seqx = matrix(rep(0,nsites*4),nrow=4)

    for(i in 1:nsites){
        seqx[,i] = siteLik(d1x,d2x,seq1.dist[,i],seq2.dist[,i],Q)
    }
    print(seqx)
    # fixit: this is merely the likelihood, we want P(x|1,2)
}

