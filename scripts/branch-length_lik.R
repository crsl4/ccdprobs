## Bret Larget
## Janaury 29, 2016
## Code to test branch length proposal idea

## Code to make a Q matrix and its spectral decomposition
##   input is:
##     x which is the vector for parameters in lower triangle of rate matrix
##     p which is the stationary distribution
##     n which is the number of states
##     a flag which will rescale so mean number of transitions per unit time is one
##   output is a list with the Q matrix itself, matrices for the eigenvectors and their inverse, and a vector of the eigenvalues
##

## Modified to use likelihood for proposal
## Claudia March 2016

makeQ = function(r,p,n,rescale=FALSE,symmetric=TRUE) {
  Q = matrix(0,n,n)
  Q[row(Q) > col(Q)] = r
  Q = Q + t(Q)
  Q = Q %*% diag(p)
  diag(Q) = -apply(Q,1,sum)
  if ( rescale ) {
      mu = sum( -p*diag(Q) )
      Q = Q/mu
  }
  if ( symmetric ) {
      p.sqrt = sqrt(p)
      S = diag(p.sqrt) %*% Q %*% diag(1/p.sqrt)
      eig = eigen(S,symmetric=TRUE) ## Error in eigen(S, symmetric = TRUE) (from branch-length.R#27) : infinite or missing values in 'x'
      V = diag(1/p.sqrt) %*% eig$vectors
      lambda = eig$values
      Vinv = t(eig$vectors) %*% diag(p.sqrt)
  }
  else {
      eig = eigen(Q)
      V = eig$vectors
      lambda = eig$values
      Vinv = solve(V)
  }
  return( list(Q=Q,V=V,Vinv=Vinv,lambda=lambda,p=p,r=r) )
}


##
matrixExp = function(Q,s)   {
  return( Q$V %*% diag( exp(Q$lambda*s) ) %*% Q$Vinv )
}

## Make a random GTR Q
randomQ = function(n,rescale=TRUE) {
    p = runif(n); p = p/sum(p)
    r = runif((n^2-n)/2)
    Q = makeQ(r,p,n,rescale=rescale)
    return( Q )
}

## Make matrix of that random pairs of states, x from stationary distribution, y state after time s
simulateSequenceSummary = function(nsites,Q,s) {
    n = nrow(Q$Q)
    P = matrixExp(Q,s)
    x = sample(1:n,size=nsites,prob=Q$p,replace=TRUE)
    y = numeric(nsites)
    for ( i in 1:nsites )
        y[i] = sample(1:n,size=1,prob=P[x[i],])
    out = matrix(0,n,n)
    for ( i in 1:n )
        for ( j in 1:n )
            out[i,j] = sum( x==i & y==j )
    return( out )
}

## GTR log-likelihood function for optimization
## theta is log(branch.length,rac,rag,rat,rcg,rct) because rgt=1 is the standardization
## this function is made to use with optim
## p is stationary distribution (using observed counts, not truth)
## x is matrix of observed counts
logl.gtr = function(theta,p,x) {
    t0 = exp(theta[1])
    r = c(exp(theta[2:6]),1)
    Q = makeQ(r,p,n=4,rescale=TRUE,symmetric=TRUE)
    P = Q$V %*% diag( exp(Q$lambda*t0) ) %*% Q$Vinv
    logl = sum(x * log(diag(p) %*% P))
    return (logl)
}

## Optimize gtr parameters, return Q matrix
## cheat and use real r as starting values
optim.gtr = function(x,r0) {
    n = sum(x)
    p = (apply(x,1,sum) + apply(x,2,sum)) / (2*n) #observed counts for p
    ## use JC distance formula to get starting point for branch length; others will all start at 1 (log(1) = 0)
    jc.x = n - sum(diag(x))
    jc.d = -3 * log(1 - 4*jc.x/(3*n))/4
    r = r0/r0[6]
    theta0 = c(log(jc.d),log(r[-6]))
#    gtr.out = optim(par=theta0,fn=logl.gtr,p=p,x=x,method="L-BFGS-B",control=list(fnscale=-1))
    gtr.out = optim(par=theta0,fn=logl.gtr,p=p,x=x,control=list(fnscale=-1))
    bl.opt = exp(gtr.out$par[1])
    r.opt = c(exp(gtr.out$par[2:6]),1)
    Q.gtr = makeQ(r.opt,p,n=4,rescale=TRUE,symmetric=TRUE)
    return(list(Q=Q.gtr,branch.length=bl.opt))
}

## Simulate branch length using JC
simulateBranchLength.jc = function(nsim,x, eta=0.9) {
    n = sum(x)
    changes = n - sum(diag(x))
    theta = 4*changes/3
    w = -0.75 * log( rbeta(nsim,eta*(n - theta),eta*theta) )
    return( list(t=w, beta=eta*(n-theta), alpha=eta*theta) )
}

## Now try the Tamura-Nei way
simulateBranchLength.tn = function(nsim,x,eta=0.9, verbose=FALSE) {
    n = sum(x)
    prop.ag = (x[1,3] + x[3,1]) / n
    prop.ct = (x[2,4] + x[4,2]) / n
    prop.tv = (x[1,2] + x[1,4] + x[2,1] + x[2,3] + x[3,2] + x[3,4] + x[4,1] + x[4,3]) / n
    if(verbose)
        print(paste("prop.ag",prop.ag,"prop.ct",prop.ct,"prop.tv",prop.tv))
    p.est = (apply(x,1,sum) + apply(x,2,sum)) / (2*n)
    p.a = p.est[1]
    p.c = p.est[2]
    p.g = p.est[3]
    p.t = p.est[4]
    p.r = sum(p.est[c(1,3)])
    p.y = sum(p.est[c(2,4)])
    if(verbose)
        print(paste("p.a",p.a,"p.c",p.c,"p.g",p.g,"p.t",p.t,"p.r",p.r,"p.y",p.y))
    numer1 = 2*p.a*p.g*p.r
    denom1 = numer1 - p.r^2*prop.ag - p.a*p.g*prop.tv
    c1 = numer1 / denom1
    numer2 = 2*p.c*p.t*p.y
    denom2 = numer2 - p.y^2*prop.ct - p.c*p.t*prop.tv
    c2 = numer2 / denom2
    c3 = (2*p.a^2*p.g^2) / (p.r * denom1) +
         (2*p.c^2*p.t^2) / (p.y * denom2) +
         (p.r^2 * (p.c^2 + p.t^2) + p.y^2 * (p.a^2 + p.g^2) ) / (2*p.r^2*p.y^2 - p.r*p.y*prop.tv)
    mu = -2 * ( (p.a*p.g/p.r) * log(1 - p.r*prop.ag/(2*p.a*p.g) - prop.tv/(2*p.r) ) +
                (p.c*p.t/p.y) * log(1 - p.y*prop.ct/(2*p.c*p.t) - prop.tv/(2*p.y) ) +
                (p.r*p.y - p.a*p.g*p.y/p.r - p.c*p.t*p.r/p.y) * log(1 - prop.tv/(2*p.r*p.y)) )
    v = (1/ eta) * ((c1^2*prop.ag + c2^2*prop.ct + c3^2*prop.tv) - (c1*prop.ag + c2*prop.ct + c3*prop.tv)^2)/n
    if(verbose)
        print(paste("mu",mu,"v",v))
    w = rgamma(nsim,mu^2/v,mu/v)
    return( list(t=w,alpha=mu^2/v,beta=mu/v) )
}

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


## function to get matrix of counts
countsMatrix = function(seq1,seq2, verbose=FALSE){
    nuc = c("a","c","g","t")
    n = 4
    out12 = matrix(0,n,n) # distance between 1 and 2
    for(i in 1:n)
        for(j in 1:n)
            out12[i,j] = sum(seq1==nuc[i] & seq2==nuc[j])
    if(verbose)
        print(out12)
    return ( out12 )
}


## function to convert a sequence into a matrix 4*nsites with one 1 per column
seqMatrix = function(seq1){
    n=4
    nsites=length(seq1)
    seq1.dist = matrix(0,n,nsites)
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
    }
    return ( seq1.dist )
}

## --------------------------------------------------------------------------------------
## likelihood functions

## loglik ((1,2)x,(3,4)y)
## fixit: need to modify to use P(Ak|x) and P(Bk|y)
## instead of d1x,d2x,d3y,d4y,seq?.dist
## maybe create another function that calls this one
gtr.log.lik.all = function(d1x,d2x,dxy,d3y,d4y,seq1.dist,seq2.dist, seq3.dist,seq4.dist,Q){
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
    return ( suma )
}

## pa = P(Ak|x), column vector of size 4
## pb = P(Bk|y), column vector of size 4
## returns fk, fk_prime, fk_doubleprime
## fixit: matrix multiplication can be more efficient with eigenvector decomp
fk = function(pa,pb,Q,t){
    Pxy = matrixExp(Q$Q,t)
    A=diag(Q$Q$p * pa)
    B=diag(pb)
    fk = rep(1,4) * A * Pxy * B * t(rep(1,4))
    fk.pr = rep(1,4) * A * Q$Q* Pxy * B * t(rep(1,4))
    fk.doublepr = rep(1,4) * A * Q$Q * Q$Q * Pxy * B * t(rep(1,4))
    return ( list(fk=fk, fk_pr=fk.pr, fk_doublepr= fk.doublepr) )
}

## returns loglik, llprime, lldoubleprime
loglik = function(seq1.dist, seq2.dist, Q, t){
    if(ncol(seq1.dist) != ncol(seq2.dist))
        stop("wrong number of sites")
    suma = 0
    sum2 = 0
    suma3 = 0
    for(i in 1:ncol(seq1.dist)){
        f = fk(seq1.dist[,i],seq2.dist[,i],Q,t)
        suma = suma + log(f$fk)
        suma2 = suma2 + f$fk_pr/f$fk
        suma3 = suma3 + (f$fk*f$fk_doublepr - f$fk_pr*f$fk_pr)/(f$fk^2)
    }
    return ( list(ll=suma, ll_pr=suma2, ll_doublepr=suma3) )
}

## t0= starting point for Newton-Raphson
findMLE = function(seq1.dist, seq2.dist, Q, t0, tol=0.0001, Nmax=10000){
    t = t0
    error = 1
    i = 1
    while(error > tol & i < Nmax){
        f =loglik(seq1.dist, seq2.dist, Q, t[i])
        t[i+1] = t[i] - f$ll/f$ll_pr
        error = abs(t[i+1]-t[i])
        i = i+1
    }
    if(i>=Nmax)
        warning("Newton-Rapshon did not converge")
    return ( list(t=t[length(t)], obsInfo=f$ll_doublepr) )
}

simulateBranchLength.lik = function(nsim,seq1,seq2, Q){
    seq1.dist = seqMatrix(seq1)
    seq2.dist = seqMatrix(seq2)
    t0 = 0.999
    mu = findMLE(seq1.dist, seq2.dist, Q, t0)
    w = rgamma(nsim, mu$t^2*mu$obsInfo, mu$t*mu$obsInfo)
    return ( w )
}


