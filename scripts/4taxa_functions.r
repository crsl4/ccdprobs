# ===========================================================================================
# Functions
# ===========================================================================================

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
sampleBLQuartet = function(d,tre,eta=0.5, verbose=FALSE, trueT0=FALSE, Q){
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

    seq1.dist = seqMatrix(seq1)
    seq2.dist = seqMatrix(seq2)
    seq3.dist = seqMatrix(seq3)
    seq4.dist = seqMatrix(seq4)

    out12 = countsMatrix(seq1,seq2)
    if(verbose)
        print(out12)
    ##r = rep(1,6)
    ##Q = optim.gtr(out12,r) ## fixit: estimating Q for 1,2 only and using everywhere
    ##if(verbose){
    ##    print(Q$Q$Q)
    ##    print(Q$Q$p)
    ##}
    jc = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    t0=jc$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        ##t0 = 0.17 #close to true value
        t0=0.14
        if(verbose)
            print(paste("true t0",t0))
    }
    d12.lik = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=t0,eta=eta)
    d12 = d12.lik$t
    if(verbose)
        print(d12)

    out13 = countsMatrix(seq1,seq3)
    if(verbose)
        print(out13)
    jc = simulateBranchLength.jc(nsim=1, out13, eta=eta)
    t0=jc$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        ##t0 = 0.2 #close to true value
        t0=0.14
        if(verbose)
            print(paste("true t0",t0))
    }
    d13.lik = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=t0,eta=eta)
    d13 = d13.lik$t
    if(verbose)
        print(d13)

    out23 = countsMatrix(seq2,seq3)
    if(verbose)
        print(out23)
    jc = simulateBranchLength.jc(nsim=1, out23, eta=eta)
    t0=jc$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        ##t0 = 0.18 #close to true value
        t0=0.14
        if(verbose)
            print(paste("true t0",t0))
    }
    d23.lik = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=t0,eta=eta)
    d23 = d23.lik$t
    if(verbose)
        print(d23)

    d1x = (d12+d13-d23)/2
    d2x = (d12+d23-d13)/2
    d3x = (d13+d23-d12)/2
    if(verbose){
        print(d1x)
        print(d2x)
        print(d3x)
    }

    seqx.dist = sequenceDist(d1x,d2x,seq1.dist,seq2.dist,Q)
    ##t0 = 0.11
    t0=0.1
    ##jc = simulateBranchLength.jc(nsim=1, out12, eta=eta) ## fixit
    ##t0=jc$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        t0 = 0.11 #close to true value
        if(verbose)
            print(paste("true t0",t0))
    }
    d4x.lik = simulateBranchLength.lik(nsim=1, seqx.dist,seq4.dist,Q,t0=t0,eta=eta)
    d4x = d4x.lik$t
    if(verbose)
        print(d4x)

    out34 = countsMatrix(seq3,seq4)
    if(verbose)
        print(out34)
    jc = simulateBranchLength.jc(nsim=1, out34, eta=eta)
    t0=jc$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        ##t0 = 0.16 #close to true value
        t0=0.1
        if(verbose)
            print(paste("true t0",t0))
    }
    d34.lik = simulateBranchLength.lik(nsim=1, seq3.dist,seq4.dist,Q,t0=t0,eta=eta)
    d34 = d34.lik$t
    if(verbose)
        print(d34)

    d3y = (d34+d3x-d4x)/2
    d4y = (d34+d4x-d3x)/2
    dxy = (d3x+d4x-d34)/2

    if(dxy<0 || d1x<0 || d2x<0 || d3y<0 || d4y<0)
        stop("dxy negative")

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

    ## now, compute likelihood of quartet with bl
    suma = gtr.log.lik.all(d1x,d2x,dxy,d3y,d4y,seq1.dist, seq2.dist, seq3.dist, seq4.dist, Q)

    # we need to compute density(bl|top)
    dens = logJointDensity.lik(d12.lik,d13.lik,d23.lik,d4x.lik,d34.lik)
    prior = logPriorExpDist(d1x,d2x,d3y,d4y,dxy,0.1) #could be other priors
    return (list(bl=bl, loglik=suma, logdensity=dens, logprior=prior))
}

## d12,d13,d23,d4x,d34 = simulateBranchLength.lik
logJointDensity.lik = function(d12.lik,d13.lik,d23.lik,d4x.lik,d34.lik){
    logd = (d12.lik$alpha-1)*log(d12.lik$t)-d12.lik$beta*d12.lik$t +
        (d13.lik$alpha-1)*log(d13.lik$t)-d13.lik$beta*d13.lik$t +
            (d23.lik$alpha-1)*log(d23.lik$t)-d23.lik$beta*d23.lik$t +
                (d4x.lik$alpha-1)*log(d4x.lik$t)-d4x.lik$beta*d4x.lik$t +
                    (d34.lik$alpha-1)*log(d34.lik$t)-d34.lik$beta*d34.lik$t
    return ( logd )
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

## JC beta
## d12,d13,d23,d3x,d4x,d34 = simulateBranchLength.jc
logJointDensity.jc = function(d1x,d2x,d3y,d4y,dxy,d12.jc,d13.jc,d23.jc,d3x.jc,d4x.jc,d34.jc, verbose=FALSE){
    logd = (-4/3)*(d12.jc$beta*(d1x+d2x)+d34.jc$beta*(d3y+d4y)+d23.jc$beta*(-d1x+d2x)+d3x.jc$beta*(d3y+dxy)+d4x.jc$beta*(d4y+dxy))+
        (d12.jc$alpha-1)*log(1-exp((-4/3)*(d1x+d2x)))+
            (d34.jc$alpha-1)*log(1-exp((-4/3)*(d3y+d4y)))+
                (d3x.jc$alpha-1)*log(1-exp((-4/3)*(d3y+dxy)))+
                    (d4x.jc$alpha-1)*log(1-exp((-4/3)*(d4y+dxy)))
    if(verbose)
        print(logd)
    return ( logd )
}


## for negativeBL.r
sampleBLQuartet_details= function(d,tre,eta=0.5, verbose=FALSE, trueT0=FALSE, Q=diag(4)){
    if(typeof(Q) == "list"){
        estQ = FALSE
    } else{
        estQ = TRUE
    }
    leaves = tre$tip.label[tre$edge[which(tre$edge[,1]==sample(c(5,6),1) & tre$edge[,2] < 5),2]] #works only for unrooted quartet
    leaves = as.numeric(leaves)
    ## if(leaves[1] == 1 || leaves[2] == 1){
    ##     sis = leaves
    ##     fth = setdiff(1:4,sis)
    ## } else{
    ##     fth = leaves
    ##     sis = setdiff(1:4,fth)
    ## }
    sis = leaves
    fth = setdiff(1:4,sis)
    seq1 = as.vector(unname(as.character(d[sis[1],])))
    seq2 = as.vector(unname(as.character(d[sis[2],])))
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

    seq1.dist = seqMatrix(seq1)
    seq2.dist = seqMatrix(seq2)
    seq3.dist = seqMatrix(seq3)
    seq4.dist = seqMatrix(seq4)

    out12 = countsMatrix(seq1,seq2)
    if(verbose)
        print(out12)
    if(estQ){
        r = rep(1,6)
        Q = optim.gtr(out12,r) ## fixit: estimating Q for 1,2 only and using everywhere
        if(verbose){
            print(Q$Q$Q)
            print(Q$Q$p)
        }
    }

    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    t0=jc12$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        ##t0 = 0.17 #close to true value
        t0=0.14
        if(verbose)
            print(paste("true t0",t0))
    }
    if(estQ){
        d12.lik = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q$Q,t0=t0,eta=eta)
    } else{
        d12.lik = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=t0,eta=eta)
    }
    d12 = d12.lik$t
    if(verbose)
        print(d12)

    out13 = countsMatrix(seq1,seq3)
    if(verbose)
        print(out13)
    jc13 = simulateBranchLength.jc(nsim=1, out13, eta=eta)
    t0=jc13$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        ##t0 = 0.2 #close to true value
        t0=0.14
        if(verbose)
            print(paste("true t0",t0))
    }
    if(estQ){
        d13.lik = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q$Q,t0=t0,eta=eta)
    } else{
        d13.lik = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=t0,eta=eta)
    }
    d13 = d13.lik$t
    if(verbose)
        print(d13)

    out23 = countsMatrix(seq2,seq3)
    if(verbose)
        print(out23)
    jc23 = simulateBranchLength.jc(nsim=1, out23, eta=eta)
    t0=jc23$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        ##t0 = 0.18 #close to true value
        t0=0.14
        if(verbose)
            print(paste("true t0",t0))
    }
    if(estQ){
        d23.lik = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q$Q,t0=t0,eta=eta)
    } else{
        d23.lik = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=t0,eta=eta)
    }
    d23 = d23.lik$t
    if(verbose)
        print(d23)

    d1x = (d12+d13-d23)/2
    d2x = (d12+d23-d13)/2
    d3x = (d13+d23-d12)/2
    if(verbose){
        print(d1x)
        print(d2x)
        print(d3x)
    }

    if(estQ){
        seqx.dist = sequenceDist(d1x,d2x,seq1.dist,seq2.dist,Q$Q)
    } else{
        seqx.dist = sequenceDist(d1x,d2x,seq1.dist,seq2.dist,Q)
    }
    t0=(max(c(jc12$t,jc13$t,jc23$t))+min(c(jc12$t,jc13$t,jc23$t)))/2
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        t0 = 0.11 #close to true value
        if(verbose)
            print(paste("true t0",t0))
    }
    if(estQ){
        d4x.lik = simulateBranchLength.lik(nsim=1, seqx.dist,seq4.dist,Q$Q,t0=t0,eta=eta)
    } else{
        d4x.lik = simulateBranchLength.lik(nsim=1, seqx.dist,seq4.dist,Q,t0=t0,eta=eta)
    }
    d4x = d4x.lik$t
    if(verbose)
        print(d4x)

    out34 = countsMatrix(seq3,seq4)
    if(verbose)
        print(out34)
    jc = simulateBranchLength.jc(nsim=1, out34, eta=eta)
    t0=jc$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        ##t0 = 0.16 #close to true value
        t0=0.1
        if(verbose)
            print(paste("true t0",t0))
    }
    if(estQ){
        d34.lik = simulateBranchLength.lik(nsim=1, seq3.dist,seq4.dist,Q$Q,t0=t0,eta=eta)
    } else{
        d34.lik = simulateBranchLength.lik(nsim=1, seq3.dist,seq4.dist,Q,t0=t0,eta=eta)
    }
    d34 = d34.lik$t
    if(verbose)
        print(d34)

    d3y = (d34+d3x-d4x)/2
    d4y = (d34+d4x-d3x)/2
    dxy = (d3x+d4x-d34)/2

    return (list(d12=d12.lik,d23=d23.lik,d13=d13.lik,d4x=d4x.lik,d34=d34.lik,
                 bl=c(d1x,d2x,d3x,d3y,d4y,dxy),
                 seq=c(as.numeric(sis[1]), as.numeric(sis[2]),fth[1],fth[2])))
}

## for simulated_negBL.r
sampleBLQuartet_details_sim= function(seq1,seq2,seq3,seq4,eta=0.5, verbose=FALSE, trueT0=FALSE, Q=diag(4)){
    if(typeof(Q) == "list"){
        estQ = FALSE
    } else{
        estQ = TRUE
    }

    nsites = length(seq1)
    if(length(seq2) != nsites && length(seq3) != nsites){
        stop("error in number of sites seq1,seq2,seq3")
    }

    nuc = c("a","c","g","t")
    n = 4

    seq1.dist = seqMatrix(seq1)
    seq2.dist = seqMatrix(seq2)
    seq3.dist = seqMatrix(seq3)
    seq4.dist = seqMatrix(seq4)

    out12 = countsMatrix(seq1,seq2)
    if(verbose)
        print(out12)
    if(estQ){
        r = rep(1,6)
        Q = optim.gtr(out12,r) ## fixit: estimating Q for 1,2 only and using everywhere
        if(verbose){
            print(Q$Q$Q)
            print(Q$Q$p)
        }
    }

    jc12 = simulateBranchLength.jc(nsim=1,out12,eta=eta)
    t0=jc12$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        ##t0 = 0.17 #close to true value
        t0=0.14
        if(verbose)
            print(paste("true t0",t0))
    }
    if(estQ){
        d12.lik = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q$Q,t0=t0,eta=eta)
    } else{
        d12.lik = simulateBranchLength.lik(nsim=1, seq1.dist,seq2.dist,Q,t0=t0,eta=eta)
    }
    d12 = d12.lik$t
    if(verbose)
        print(d12)

    out13 = countsMatrix(seq1,seq3)
    if(verbose)
        print(out13)
    jc13 = simulateBranchLength.jc(nsim=1, out13, eta=eta)
    t0=jc13$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        ##t0 = 0.2 #close to true value
        t0=0.14
        if(verbose)
            print(paste("true t0",t0))
    }
    if(estQ){
        d13.lik = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q$Q,t0=t0,eta=eta)
    } else{
        d13.lik = simulateBranchLength.lik(nsim=1, seq1.dist,seq3.dist,Q,t0=t0,eta=eta)
    }
    d13 = d13.lik$t
    if(verbose)
        print(d13)

    out23 = countsMatrix(seq2,seq3)
    if(verbose)
        print(out23)
    jc23 = simulateBranchLength.jc(nsim=1, out23, eta=eta)
    t0=jc23$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        ##t0 = 0.18 #close to true value
        t0=0.14
        if(verbose)
            print(paste("true t0",t0))
    }
    if(estQ){
        d23.lik = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q$Q,t0=t0,eta=eta)
    } else{
        d23.lik = simulateBranchLength.lik(nsim=1, seq2.dist,seq3.dist,Q,t0=t0,eta=eta)
    }
    d23 = d23.lik$t
    if(verbose)
        print(d23)

    d1x = (d12+d13-d23)/2
    d2x = (d12+d23-d13)/2
    d3x = (d13+d23-d12)/2
    if(verbose){
        print(d1x)
        print(d2x)
        print(d3x)
    }

    if(estQ){
        seqx.dist = sequenceDist(d1x,d2x,seq1.dist,seq2.dist,Q$Q)
    } else{
        seqx.dist = sequenceDist(d1x,d2x,seq1.dist,seq2.dist,Q)
    }
    t0=(max(c(jc12$t,jc13$t,jc23$t))+min(c(jc12$t,jc13$t,jc23$t)))/2
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        t0 = 0.11 #close to true value
        if(verbose)
            print(paste("true t0",t0))
    }
    if(estQ){
        d4x.lik = simulateBranchLength.lik(nsim=1, seqx.dist,seq4.dist,Q$Q,t0=t0,eta=eta)
    } else{
        d4x.lik = simulateBranchLength.lik(nsim=1, seqx.dist,seq4.dist,Q,t0=t0,eta=eta)
    }
    d4x = d4x.lik$t
    if(verbose)
        print(d4x)

    out34 = countsMatrix(seq3,seq4)
    if(verbose)
        print(out34)
    jc = simulateBranchLength.jc(nsim=1, out34, eta=eta)
    t0=jc$t
    if(verbose)
        print(paste("JC t0", t0))
    if(trueT0){
        ##t0 = 0.16 #close to true value
        t0=0.1
        if(verbose)
            print(paste("true t0",t0))
    }
    if(estQ){
        d34.lik = simulateBranchLength.lik(nsim=1, seq3.dist,seq4.dist,Q$Q,t0=t0,eta=eta)
    } else{
        d34.lik = simulateBranchLength.lik(nsim=1, seq3.dist,seq4.dist,Q,t0=t0,eta=eta)
    }
    d34 = d34.lik$t
    if(verbose)
        print(d34)

    d3y = (d34+d3x-d4x)/2
    d4y = (d34+d4x-d3x)/2
    dxy = (d3x+d4x-d34)/2

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

    ## now, compute likelihood of quartet with bl
    suma = gtr.log.lik.all(d1x,d2x,dxy,d3y,d4y,seq1.dist, seq2.dist, seq3.dist, seq4.dist, Q)

    # we need to compute density(bl|top)
    dens = logJointDensity.lik(d12.lik,d13.lik,d23.lik,d4x.lik,d34.lik)
    prior = logPriorExpDist(d1x,d2x,d3y,d4y,dxy,0.1) #could be other priors
    return (list(bl=bl, loglik=suma, logdensity=dens, logprior=prior))

    ## return (list(d12=d12.lik,d23=d23.lik,d13=d13.lik,d4x=d4x.lik,d34=d34.lik,
    ##              bl=c(d1x,d2x,d3x,d3y,d4y,dxy),
    ##              seq=c(1,2,3,4)))
}
