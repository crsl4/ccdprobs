# ===========================================================================================
# Functions
# ===========================================================================================

# 1-----x-----2
# d1x = distance from 1 to parent x, similarly d2x
# seq1.distj = jth column in seq1.dist matrix (for site j), similarly seq2.distj
# Q = estimated matrix of rates
# returns column of site likelihood
siteLik = function(d1x,d2x,seq1.distj,seq2.distj, Q, verbose=FALSE){
    P1 = matrixExp(Q,d1x)
    P2 = matrixExp(Q,d2x)
    lik = rep(0,4)
    for(i in 1:4){
        lik[i] = P1[i,]%*%seq1.distj * P2[i,]%*%seq2.distj
        if(verbose)
            print(lik)
    }
    return (lik)
}


# estimates the seq dist at x (parent of 1 and 2)
sequenceDist = function(d1x,d2x,seq1.dist,seq2.dist, Q, verbose=FALSE){
    nsites = length(seq1.dist[1,])
    if(length(seq2.dist[1,]) != nsites){
        stop("error in number of sites seq1,seq2")
    }
    nuc = c('a','c','g','t')
    seqx = matrix(rep(0,nsites*4),nrow=4)
    for(i in 1:nsites){
        seqx[,i] = siteLik(d1x,d2x,seq1.dist[,i],seq2.dist[,i],Q) * Q$p
        seqx[,i] = seqx[,i]/sum(seqx[,i])
    }
    if(verbose)
        print(seqx)
    return(seqx)
}

sampleSeq = function(seqx){
    seq = matrix(rep(0,ncol(seqx)*4),nrow=4)
    for(i in 1:ncol(seqx)){
        ind = sample(1:4,1,replace=FALSE,prob=seqx[,i])
        seq[ind,i] = 1
    }
    return (seq)
}

#warning if x (counts) with rows of zero (and also columns because error in density(tn)), also if few counts in one row/column
checkMatCounts = function(x){
    if((0 %in% rowSums(x)) || (0 %in% colSums(x)) || (any(rowSums(x)<5)) || (any(colSums(x)<5))){
        warning("the counts matrix rows/columns of zero, or very few counts per row/column, TN simulation of branch length could have an error")
    }
}

sampleTopQuartet = function(dat.tre, verbose=FALSE){
    u=runif(1,0,1)
    if(u<dat.tre$V2[1]){
        t = dat.tre$V1[1]
        prob=dat.tre$V2[1]
    } else if(u < dat.tre$V2[1]+dat.tre$V2[2]){
        t = dat.tre$V1[2]
        prob=dat.tre$V2[2]
    } else {
        t = dat.tre$V1[3]
        prob=dat.tre$V2[3]
    }
    tre = read.tree(text=as.character(t))
    if(verbose)
        print(write.tree(tre))
    return (list(tre=tre,prob=prob))
}

# function to sample branch lengths for 4-taxon tree
sampleBLQuartet = function(d,tre,eta=0.5, verbose=FALSE){
    sis = tre$tip.label[tre$edge[which(tre$edge[,1]==sample(c(5,6),1) & tre$edge[,2] < 5),2]] #works only for unrooted quartet
    fth = setdiff(1:4,sis)
    seq1 = as.vector(unname(as.character(d[as.numeric(sis[1]),])))
    seq2 = as.vector(unname(as.character(d[as.numeric(sis[2]),])))
    seq3 = as.vector(unname(as.character(d[fth[1],])))
    seq4 = as.vector(unname(as.character(d[fth[2],])))
    # remove missing:
    s1 <-seq1!="-"
    s2 <- seq2!="-"
    s3 <- seq3!="-"
    s4 <- seq4!="-"
    seq1 <- seq1[s1&s2&s3&s4]
    seq2 <- seq2[s1&s2&s3&s4]
    seq3 <- seq3[s1&s2&s3&s4]
    seq4 <- seq4[s1&s2&s3&s4]

    nsites = length(seq1)
    if(length(seq2) != nsites && length(seq3) != nsites){
        stop("error in number of sites seq1,seq2,seq3")
    }

    nuc = c("a","c","g","t")
    n = 4

    out12 = matrix(0,n,n) # distance between 1 and 2
    for(i in 1:n)
        for(j in 1:n)
            out12[i,j] = sum(seq1==nuc[i] & seq2==nuc[j])
    if(verbose)
        print(out12)
    checkMatCounts(out12)

    #d12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    d12.tn = simulateBranchLength.tn(nsim=1, out12, eta=eta)
    d12=d12.tn$t
    if(verbose)
        print(d12)

    out13 = matrix(0,n,n)     # distance between 1 and 3
    for(i in 1:n)
        for(j in 1:n)
            out13[i,j] = sum(seq1==nuc[i] & seq3==nuc[j])
    if(verbose)
        print(out13)
    checkMatCounts(out13)

    #d13 = simulateBranchLength.jc(nsim=1,out13,eta=eta)
    d13.tn = simulateBranchLength.tn(nsim=1, out13, eta=eta)
    d13 = d13.tn$t
    if(verbose)
        print(d13)

    out23 = matrix(0,n,n)     # distance between 2 and 3
    for(i in 1:n)
        for(j in 1:n)
            out23[i,j] = sum(seq2==nuc[i] & seq3==nuc[j])
    if(verbose)
        print(out23)
    checkMatCounts(out23)

    #d23 = simulateBranchLength.jc(nsim=1,out23,eta=eta)
    d23.tn = simulateBranchLength.tn(nsim=1, out23, eta=eta)
    d23 = d23.tn$t
    if(verbose)
        print(d23)

    d1x = (d12+d13-d23)/2
    d2x = (d12+d23-d13)/2

    if(verbose){
        print(d1x)
        print(d2x)
    }

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

    #need seq distribution at x: felsenstein algorithm
    r = rep(1,6)
    #r = c(2.815,51.982,1.903,1.275,65.402,1.000) #mrbayes
    Q = optim.gtr(out12,r)
    if(verbose){
        print(Q$Q$Q)
        print(Q$Q$p)
    }
    #seqx = sequenceDist(d1x,d2x, seq1.dist, seq2.dist,Q$Q)

    ## sample sequence at x
    seqdist = sequenceDist(d1x,d2x, seq1.dist, seq2.dist,Q$Q)
    seqx = sampleSeq(seqdist)
    if(verbose)
        print(seqx)

    out34 = matrix(0,n,n) # distance between 3 and 4
    for(i in 1:n)
        for(j in 1:n)
            out34[i,j] = sum(seq3==nuc[i] & seq4==nuc[j])
    if(verbose)
        print(out34)
    checkMatCounts(out34)

    #d34 = simulateBranchLength.jc(nsim=1,out34,eta=eta)
    d34.tn = simulateBranchLength.tn(nsim=1, out34, eta=eta)
    d34 = d34.tn$t
    if(verbose)
        print(d34)

    out3x = matrix(0,n,n) # distance between 3 and x
    for(i in 1:nsites){
        out3x = out3x + seq3.dist[,i]%*%t(seqx[,i])
    }
    if(verbose)
        print(out3x)
    checkMatCounts(out3x)

    #d3x = simulateBranchLength.jc(nsim=1,out3x,eta=eta)
    d3x.tn = simulateBranchLength.tn(nsim=1, out3x, eta=eta)
    d3x = d3x.tn$t
    if(verbose)
        print(d3x)

    out4x = matrix(0,n,n) # distance between 4 and x
    for(i in 1:nsites){
        out4x = out4x + seq4.dist[,i]%*%t(seqx[,i])
    }
    if(verbose)
        print(out4x)
    checkMatCounts(out4x)

    #d4x = simulateBranchLength.jc(nsim=1,out4x,eta=eta)
    d4x.tn = simulateBranchLength.tn(nsim=1, out4x, eta=eta)
    d4x = d4x.tn$t
    if(verbose)
        print(d4x)

    d3y = (d34+d3x-d4x)/2
    d4y = (d34+d4x-d3x)/2
    dxy = (d3x+d4x-d34)/2

    if(verbose){
        print(paste("d12",d12,"d13",d13,"d23",d23,"d34",d34))
        print(paste("d1x",d1x,"d2x",d2x,"d3y",d3y,"d4y",d4y,"dxy",dxy))
    }

    bl <- rep(0,5) #works only for unrooted quartet
    ed1x = which(tre$edge[,2] == which(tre$tip.label == sis[1]))
    ed2x = which(tre$edge[,2] == which(tre$tip.label == sis[2]))
    ed3y = which(tre$edge[,2] == which(tre$tip.label == as.character(fth[1])))
    ed4y = which(tre$edge[,2] == which(tre$tip.label == as.character(fth[2])))
    bl[ed1x] = d1x
    bl[ed2x] = d2x
    bl[ed3y] = d3y
    bl[ed4y] = d4y
    ind = which(bl==0)
    bl[ind] = dxy
    if(verbose)
        print(bl)

    # now, compute likelihood of quartet with bl
    suma = 0
    for(s in 1:nsites){
        lik12 = siteLik(d1x,d2x,seq1.dist[,s],seq2.dist[,s],Q$Q)
        lik34 = siteLik(d3y,d4y,seq3.dist[,s],seq4.dist[,s],Q$Q)
        L = lik12 %*% t(lik34)
        Pxy = matrixExp(Q$Q,dxy)
        L2 = L*Pxy
        Lik = Q$Q$p * L2
        suma = suma+log(sum(Lik))
    }
    if(verbose)
        print(suma)

    # we need to compute density(bl|top)
    dens = logJointDensity.tn(d1x,d2x,d3y,d4y,dxy,d12.tn,d13.tn,d23.tn,d3x.tn,d4x.tn,d34.tn)
    prior = logPriorExpDist(d1x,d2x,d3y,d4y,dxy,0.1) #could be other priors
    return (list(bl=bl, loglik=suma, logdensity=dens, logprior=prior))
}

# d12,d13,d23,d3x,d4x,d34 = simulateBranchLength.tn
logJointDensity.tn = function(d1x,d2x,d3y,d4y,dxy,d12.tn,d13.tn,d23.tn,d3x.tn,d4x.tn,d34.tn, verbose=FALSE){
    logd = (d12.tn$alpha-1)*log(d1x+d2x)-d12.tn$beta*(d1x+d2x)+(d34.tn$alpha-1)*log(d3y+d4y)-d34.tn$beta*(d3y+d4y)-d23.tn$beta*(d2x-d1x)+(d3x.tn$alpha-1)*log(d3y+dxy)-
        d3x.tn$beta*(d3y+dxy)+(d4x.tn$alpha-1)*log(d4y+dxy)-d4x.tn$beta*(d4y+dxy)
    if(verbose)
        print(logd)

    # fixit: integral shown finite, but R says it diverges
    ## integrand = function(x){
    ##     x^(d13.tn$alpha-1)*(x-d1x+d2x)^(d23.tn$alpha-1)*exp((-d13.tn$beta-d23.tn$beta)*x)
    ## }
    ## # also, beta~500, so exp almost zero for x>10
    ## I = integrate(integrand,lower=0,upper=10) # if upper=Inf, error
    ## logd = logd + log(I$value)

    return ( logd ) #fixit: ignoring integral
}

# m= mean
logPriorExpDist = function(d1x,d2x,d3y,d4y,dxy,m, verbose=FALSE){
    logp = (-1/m)*(d1x+d2x+d3y+d4y+dxy)
    if(verbose)
        print(logp)
    return ( logp )
}


weighted.quantile = function(x,w,probs=0.25){
    if(length(x) != length(w))
        stop("x and w should have same length")
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]/sum(w)
    y <- which(cumsum(w)<probs)
    if(length(y) == 0){
        return ( 0.0 )
    } else {
        return( x[y[length(y)]] )
    }
}

weighted.mean = function(x,w){
    if(length(x) != length(w))
        stop("x and w should have same length")
    w <- w/sum(w) # need to normalize
    return (x %*% w)
}


## NORMAL
## d12,d13,d23,d3x,d4x,d34 = simulateBranchLength.norm
logJointDensity.norm = function(d1x,d2x,d3y,d4y,dxy,d12.n,d13.n,d23.n,d3x.n,d4x.n,d34.n, verbose=FALSE){
    logd = (d1x+d2x-d12.n$mu)^2/d12.n$sigma^2 + (d3y+d4y-d34.n$mu)^2/d34.n$sigma^2 + (d3y+dxy-d3x.n$mu)^2/d3x.n$sigma^2 + (d4y+dxy-d4x.n$mu)^2/d4x.n$sigma^2 +
        (-d1x+d2x-d23.n$mu)^2/(d23.n$sigma^2) -
                ((d23.n$sigma^2*d13.n$mu+d13.n$sigma^2*d1x-d13.n$sigma^2*d2x+d13.n$sigma^2*d23.n$mu)^2/(d13.n$sigma^2*d23.n$sigma^2*(d23.n$sigma^2+d13.n$sigma^2)))
    logd = -0.5 * logd
    if(verbose)
        print(logd)
    return ( logd )
}
