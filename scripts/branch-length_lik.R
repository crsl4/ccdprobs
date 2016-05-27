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

makeQscale = function(Q0,scale,p,rescale=FALSE,symmetric=TRUE) {
  Q = scale * Q0
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


## Calculate GTR log-likelihood at a sequence s for an observed matrix x and matrix Q
gtr.log.like = function(x,s,Q) {
    ns = length(s)
    log.like = numeric( ns )
    for ( i in 1:ns ) {
        P = matrixExp(Q,s[i])
        log.like[i] = sum( x * log(diag(Q$p) %*% P) )
    }
    return( log.like )
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
        seqx[,i] = siteLik(d1x,d2x,seq1.dist[,i],seq2.dist[,i],Q) #* Q$p
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
    for(s in 1:ncol(seq1.dist)){
        lik12 = siteLik(d1x,d2x,seq1.dist[,s],seq2.dist[,s],Q)
        lik34 = siteLik(d3y,d4y,seq3.dist[,s],seq4.dist[,s],Q)
        L = lik12 %*% t(lik34)
        Pxy = matrixExp(Q,dxy)
        L2 = L*Pxy
        Lik = Q$p * L2
        suma = suma+log(sum(Lik))
    }
    return ( suma )
}

## pa = P(Ak|x), column vector of size 4
## pb = P(Bk|y), column vector of size 4
## returns fk, fk_prime, fk_doubleprime
## fixit: matrix multiplication can be more efficient with eigenvector decomp
fk = function(pa,pb,Q,t){
    Pxy = matrixExp(Q,t)
    A=diag(Q$p * pa)
    B=diag(pb)
    fk = rep(1,4) %*% A %*% Pxy %*% B %*% rep(1,4)
    fk.pr = rep(1,4) %*% A %*% Q$Q %*% Pxy %*% B %*% rep(1,4)
    fk.doublepr = rep(1,4) %*% A %*% Q$Q %*% Q$Q %*% Pxy %*% B %*% rep(1,4)
    return ( list(fk=as.numeric(fk), fk_pr=as.numeric(fk.pr), fk_doublepr= as.numeric(fk.doublepr)) )
}

## returns loglik, llprime, lldoubleprime
loglik = function(seq1.dist, seq2.dist, Q, t){
    if(ncol(seq1.dist) != ncol(seq2.dist))
        stop("wrong number of sites")
    suma = 0
    suma2 = 0
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
findMLE = function(seq1.dist, seq2.dist, Q, t0=0.1, tol=0.0001, Nmax=10000, verbose=FALSE){
    if(verbose)
        print("entering findMLE...")
    t = rep(0,Nmax)
    if(t0<0)
        to=0.00001
    t[1] = t0 ## fixit: later do a binary search before choosing t0
    error = 1
    i = 1
    while(error > tol & i < Nmax){
        if(verbose)
            print(t[i])
        f =loglik(seq1.dist, seq2.dist, Q, t[i])
        gap = f$ll_pr/f$ll_doublepr
        t[i+1] = t[i] - gap
        while(t[i+1]<0){ #avoid negative BL
            gap = gap/2
            if(verbose)
                print("found negative candidate t[i+1]")
            t[i+1] = t[i] - gap
        }
        ## after finding a positive candidate:
        f2 =loglik(seq1.dist, seq2.dist, Q, t[i+1])
        if(verbose){
            print("f$ll, f2$ll")
            print(f$ll)
            print(f2$ll)
        }
        if(f$ll * f2$ll < 0){ #one positive, and one negative
            while(abs(f2$ll) > 2*abs(f$ll)){ #t[i+1] too far from root
                if(verbose)
                    print("t[i+1] positive, but too far from root")
                t[i+1] = t[i+1]/2
                f2 =loglik(seq1.dist, seq2.dist, Q, t[i+1])
            }
        }
        error = abs(t[i+1]-t[i])
        i = i+1
    }
    t = t[which(t>0)]
    if(i>=Nmax)
        warning("Newton-Rapshon did not converge")
    return ( list(t=t[length(t)], obsInfo=f2$ll_doublepr) )
}

simulateBranchLength.lik = function(nsim,seq1.dist,seq2.dist, Q, t0, eta=0.5, verbose=FALSE){
    mu = findMLE(seq1.dist, seq2.dist, Q, t0, verbose=verbose)
    w = rgamma(nsim, mu$t^2*(-mu$obsInfo)*eta, mu$t*(-mu$obsInfo)*eta)
    return ( list(t=w, alpha=mu$t^2*(-mu$obsInfo)*eta, beta=mu$t*(-mu$obsInfo)*eta) )
}

simulateBranchLength.norm = function(nsim,seq1.dist,seq2.dist, Q, t0, eta=0.5, verbose=FALSE){
    mu = findMLE(seq1.dist, seq2.dist, Q, t0, verbose=verbose)
    w = rnorm(nsim, mu$t, sqrt(-1/(eta*mu$obsInfo)))
    return ( list(t=w, mu=mu$t, sigma=sqrt(-1/(eta*mu$obsInfo))) )
}



simulateData = function(Q,branch.length, nsites, filename="simSeq.txt"){
    nuc <- c('a','c','g','t')
    ## simulate seqx
    seqx = sample(nuc,size=nsites,prob=Q$p,replace=TRUE)

    ## simulate seqy
    s=branch.length[1]
    P = matrixExp(Q,s)
    seqy = numeric(nsites)
    for ( i in 1:nsites )
        seqy[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

    ## simulate seq1
    s=branch.length[2]
    P = matrixExp(Q,s)
    seq1 = numeric(nsites)
    for ( i in 1:nsites )
        seq1[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

    ## simulate seq2
    s=branch.length[3]
    P = matrixExp(Q,s)
    seq2 = numeric(nsites)
    for ( i in 1:nsites )
        seq2[i] = sample(nuc,size=1,prob=P[which(nuc==seqx[i]),])

    ## simulate seq3
    s=branch.length[4]
    P = matrixExp(Q,s)
    seq3 = numeric(nsites)
    for ( i in 1:nsites )
        seq3[i] = sample(nuc,size=1,prob=P[which(nuc==seqy[i]),])

    ## simulate seq4
    s=branch.length[5]
    P = matrixExp(Q,s)
    seq4 = numeric(nsites)
    for ( i in 1:nsites )
        seq4[i] = sample(nuc,size=1,prob=P[which(nuc==seqy[i]),])

    l1 = paste("6",nsites)
    l2 = paste("1",paste0(seq1,collapse=""))
    l3 = paste("2",paste0(seq2,collapse=""))
    l4 = paste("3",paste0(seq3,collapse=""))
    l5 = paste("4",paste0(seq4,collapse=""))
    l6 = paste("x",paste0(seqx,collapse=""))
    l7 = paste("y",paste0(seqy,collapse=""))

    write(l1,file=filename)
    write(l2,file=filename, append=TRUE)
    write(l3,file=filename, append=TRUE)
    write(l4,file=filename, append=TRUE)
    write(l5,file=filename, append=TRUE)
    write(l6,file=filename, append=TRUE)
    write(l7,file=filename, append=TRUE)

    d = read.dna(filename)

    return ( d )
}




## -------------------------
## 3D likelihood

## p = P(Ai|xi) column vector of size 4
## Q = matrix
## returns column vector 4x1, with rows as y
## for S, S', S''
Sfn = function(t,p,Q){
    P = matrixExp(Q,t)
    S = P %*% diag(p) %*% rep(1,4)
    Spr = Q$Q %*% P %*% diag(p) %*% rep(1,4)
    Sdoublepr = Q$Q %*% Q$Q %*% P %*% diag(p) %*% rep(1,4)
    return (list(S=S, Spr=Spr, Sdoublepr=Sdoublepr))
}

## p1 = P(A1|x1), column vector of size 4
## p2 = P(A2|x2), column vector of size 4
## p3 = P(A3|x3), column vector of size 4
## returns fk, fk_prime1, fk_prime2, fk_prime3,
## fk_doubleprime11, fk_doubleprime12, fk_doubleprime13, fk_doubleprime22,
## fk_doubleprime23, fk_doubleprime33
## fixit: matrix multiplication can be more efficient with eigenvector decomp
fk3D = function(p1,p2,p3,Q,t1,t2,t3){
    S1fn = Sfn(t1,p1,Q)
    S2fn = Sfn(t2,p2,Q)
    S3fn = Sfn(t3,p3,Q)
    S1=diag(c(S1fn$S))
    S2=diag(c(S2fn$S))
    S3=diag(c(S3fn$S))
    S1pr = diag(c(S1fn$Spr))
    S2pr = diag(c(S2fn$Spr))
    S3pr = diag(c(S3fn$Spr))
    S1doublepr = diag(c(S1fn$Sdoublepr))
    S2doublepr = diag(c(S2fn$Sdoublepr))
    S3doublepr = diag(c(S3fn$Sdoublepr))
    fk = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% S3 %*% rep(1,4)
    fk.pr1 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2 %*% S3 %*% rep(1,4)
    fk.pr2 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2pr %*% S3 %*% rep(1,4)
    fk.pr3 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% S3pr %*% rep(1,4)
    fk.doublepr11 = rep(1,4) %*% diag(Q$p) %*% S1doublepr %*% S2 %*% S3 %*% rep(1,4)
    fk.doublepr22 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2doublepr %*% S3 %*% rep(1,4)
    fk.doublepr33 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% S3doublepr %*% rep(1,4)
    fk.doublepr12 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2pr %*% S3 %*% rep(1,4)
    fk.doublepr13 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2 %*% S3pr %*% rep(1,4)
    fk.doublepr23 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2pr %*% S3pr %*% rep(1,4)

    return ( list(fk=as.numeric(fk), fk_pr1=as.numeric(fk.pr1),fk_pr2=as.numeric(fk.pr2),fk_pr3=as.numeric(fk.pr3),
                  fk_doublepr11= as.numeric(fk.doublepr11),fk_doublepr12= as.numeric(fk.doublepr12),fk_doublepr13= as.numeric(fk.doublepr13),
                  fk_doublepr22= as.numeric(fk.doublepr22),fk_doublepr23= as.numeric(fk.doublepr23),fk_doublepr33= as.numeric(fk.doublepr33)))
}

## returns loglik, gradient (vector 3x1), hessian (matrix 3x3)
## t = vector 3x1 t1,t2,t3
loglik3D = function(seq1.dist, seq2.dist, seq3.dist,Q, t){
    if(ncol(seq1.dist) != ncol(seq2.dist))
        stop("wrong number of sites")
    logl = 0
    dll1 = 0
    dll2 = 0
    dll3 = 0
    d2ll_11 = 0
    d2ll_12 = 0
    d2ll_13 = 0
    d2ll_22 = 0
    d2ll_23 = 0
    d2ll_33 = 0
    for(i in 1:ncol(seq1.dist)){
        f = fk3D(seq1.dist[,i],seq2.dist[,i],seq3.dist[,i],Q,t[1],t[2],t[3])
        logl = logl + log(f$fk)
        dll1 = dll1 + f$fk_pr1/f$fk
        dll2 = dll2 + f$fk_pr2/f$fk
        dll3 = dll3 + f$fk_pr3/f$fk
        d2ll_11= d2ll_11 + (f$fk*f$fk_doublepr11 - f$fk_pr1*f$fk_pr1)/(f$fk^2)
        d2ll_12= d2ll_12 + (f$fk*f$fk_doublepr12 - f$fk_pr1*f$fk_pr2)/(f$fk^2)
        d2ll_13= d2ll_13 + (f$fk*f$fk_doublepr13 - f$fk_pr1*f$fk_pr3)/(f$fk^2)
        d2ll_22= d2ll_22 + (f$fk*f$fk_doublepr22 - f$fk_pr2*f$fk_pr2)/(f$fk^2)
        d2ll_23= d2ll_23 + (f$fk*f$fk_doublepr23 - f$fk_pr2*f$fk_pr3)/(f$fk^2)
        d2ll_33= d2ll_33 + (f$fk*f$fk_doublepr33 - f$fk_pr3*f$fk_pr3)/(f$fk^2)
    }
    return ( list(ll=logl, gradient=c(dll1,dll2,dll3), hessian=matrix(c(d2ll_11, d2ll_12, d2ll_13, d2ll_12, d2ll_22, d2ll_23, d2ll_13, d2ll_23, d2ll_33),ncol=3)) )
}

## t0= starting point for Newton-Raphson
## fixit: jumps could be farther from root, but no way to fix this in 3D (bret gave idea, need to code it)
findMLE3D = function(seq1.dist, seq2.dist, seq3.dist, Q, t0=c(0.1,0.1,0.1), tol=0.0001, Nmax=10000, verbose=FALSE){
    if(verbose)
        print("entering findMLE...")
    tnew = rep(0,3) # will not save all sequence
    if(any(t0<0))
        t0[which(t0<0)] = 0.0001
    told = t0 ## fixit: later do a binary search before choosing t0
    error = 1
    i = 1
    while(error > tol & i < Nmax){
        if(verbose)
            print(told)
        f =loglik3D(seq1.dist, seq2.dist, seq3.dist,Q, told)
        if(verbose){
            print(f$gradient)
            print(f$hessian)
        }
        gap = solve(f$hessian) %*% f$gradient
        tnew = told - gap
        while(any(tnew<0)){ #avoid negative BL
            gap = gap/2
            if(verbose)
                print("found negative candidate tnew")
            tnew = told - gap
        }
        ## after finding a positive candidate:
        f2 =loglik3D(seq1.dist, seq2.dist, seq3.dist, Q, tnew)
        ## if(verbose){
        ##     print("f$ll, f2$ll")
        ##     print(f$ll)
        ##     print(f2$ll)
        ## }
        ## if(f$ll * f2$ll < 0){ #one positive, and one negative
        ##     while(abs(f2$ll) > 2*abs(f$ll)){ #t[i+1] too far from root
        ##         if(verbose)
        ##             print("t[i+1] positive, but too far from root")
        ##         t[i+1] = t[i+1]/2
        ##         f2 =loglik(seq1.dist, seq2.dist, Q, t[i+1])
        ##     }
        ## }
        error = max(abs(tnew-told))
        if(verbose)
            print(paste("Error: ", error))
        i = i+1
        told = tnew
    }
    if(i>=Nmax)
        warning("Newton-Rapshon did not converge")
    print("Gradient at the end:")
    print(f2$gradient)
    return ( list(t=tnew, obsInfo=f2$hessian) )
}

## simulates d1x,d2x,d3x jointly
simulateBranchLength.multinorm = function(nsim,seq1.dist,seq2.dist, seq3.dist, Q, t0, verbose=FALSE){
    mu = findMLE3D(seq1.dist, seq2.dist, seq3.dist,Q, t0, verbose=verbose)
    Sigma = solve(mu$obsInfo)
    w = rmvnorm(nsim, mu$t, -Sigma)
    return ( list(t=w, mu=mu$t, sigma=-Sigma) )
}

## simulates d1x,d2x,d3x jointly
## with gamma/beta
simulateBranchLength.gamma = function(nsim,seq1.dist,seq2.dist, seq3.dist, Q, t0, verbose=FALSE){
    mu = findMLE3D(seq1.dist, seq2.dist, seq3.dist,Q, t0, verbose=verbose)
    Sigma = solve(mu$obsInfo)
    R = chol(-Sigma)
    L = t(R)
    ## simulate T1
    alpha1 = mu$t[1]^2/(L[1,1]*L[1,1])
    lambda1 = mu$t[1]/(L[1,1]*L[1,1])
    t1 = rgamma(1,shape=alpha1,rate=lambda1)
    ## simulate T2|T1
    z1 = (t1-mu$t[1])/L[1,1]
    num = mu$t[2] + L[2,1]*z1
    alpha2 = num^2/L[2,2]^2
    lambda2 = num/L[2,2]^2
    t2 = rgamma(1,shape=alpha2,rate=lambda2)
    ## simulate T3|T1,T2
    z2 = (t2-mu$t[2]-L[2,1]*z1)/L[2,2]
    num = mu$t[3] + L[3,1]*z1 + L[3,2]*z2
    alpha3 = num^2/L[3,3]^2
    lambda3 = num/L[3,3]^2
    t3 = rgamma(1,shape=alpha3,rate=lambda3)
    w = c(t1,t2,t3)
    ## logdensity
    logdens = (alpha1-1)*log(t1)-lambda1*t1+(alpha2-1)*log(t2)-lambda2*t2+(alpha3-1)*log(t3)-lambda3*t3
    return ( list(t=w, mu=mu$t, sigma=-Sigma, logdens=logdens) )
}

## ----------------------------------------------------
## 2D conditional likelihood

## p1 = P(A1|x1), column vector of size 4
## p2 = P(A2|x2), column vector of size 4
## input t1,t2 and sum t1+t3=s13
## returns fk, fk_prime1, fk_prime2, fk_prime3,
## fk_doubleprime11, fk_doubleprime12, fk_doubleprime13, fk_doubleprime22,
## fk_doubleprime23, fk_doubleprime33
## fixit: matrix multiplication can be more efficient with eigenvector decomp
fk2D = function(p1,p2,p3,Q,t1,t2,s13){
    if(s13-t1<0)
        stop("problem where sum is less than summand: t1+t3=s13")
    S1fn = Sfn(t1,p1,Q)
    S2fn = Sfn(t2,p2,Q)
    S3fn = Sfn(s13-t1,p3,Q)
    S1=diag(c(S1fn$S))
    S2=diag(c(S2fn$S))
    S3=diag(c(S3fn$S))
    S1pr = diag(c(S1fn$Spr))
    S2pr = diag(c(S2fn$Spr))
    S3pr = (-1)*diag(c(S3fn$Spr))
    S1doublepr = diag(c(S1fn$Sdoublepr))
    S2doublepr = diag(c(S2fn$Sdoublepr))
    S3doublepr = diag(c(S3fn$Sdoublepr))
    fk = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% S3 %*% rep(1,4)
    fk.pr1 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2 %*% S3 %*% rep(1,4)
    fk.pr2 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2pr %*% S3 %*% rep(1,4)
    fk.pr3 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% S3pr %*% rep(1,4)
    fk.doublepr11 = rep(1,4) %*% diag(Q$p) %*% S1doublepr %*% S2 %*% S3 %*% rep(1,4)
    fk.doublepr22 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2doublepr %*% S3 %*% rep(1,4)
    fk.doublepr33 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% S3doublepr %*% rep(1,4)
    fk.doublepr12 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2pr %*% S3 %*% rep(1,4)
    fk.doublepr13 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2 %*% S3pr %*% rep(1,4)
    fk.doublepr23 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2pr %*% S3pr %*% rep(1,4)

    return ( list(fk=as.numeric(fk), fk_pr1=as.numeric(fk.pr1),fk_pr2=as.numeric(fk.pr2),fk_pr3=as.numeric(fk.pr3),
                  fk_doublepr11= as.numeric(fk.doublepr11),fk_doublepr12= as.numeric(fk.doublepr12),fk_doublepr13= as.numeric(fk.doublepr13),
                  fk_doublepr22= as.numeric(fk.doublepr22),fk_doublepr23= as.numeric(fk.doublepr23),fk_doublepr33= as.numeric(fk.doublepr33)))
}


## returns loglik, gradient (vector 2x1), hessian (matrix 2x2)
## t = vector 2x1 t1,t2
loglik2D = function(seqx.dist, seq3.dist, seq4.dist,Q, t, d3x){
    if(ncol(seq3.dist) != ncol(seq4.dist))
        stop("wrong number of sites")
    logl = 0
    dll1 = 0
    dll2 = 0
    d2ll_11 = 0
    d2ll_12 = 0
    d2ll_22 = 0
    for(i in 1:ncol(seqx.dist)){
        f = fk2D(seqx.dist[,i],seq4.dist[,i],seq3.dist[,i],Q,t[1],t[2],d3x) ## t1=dxy, t2=d4y, t3=d3y
        logl = logl + log(f$fk)
        dll1 = dll1 + (f$fk_pr1+f$fk_pr3)/f$fk
        dll2 = dll2 + f$fk_pr2/f$fk
        d2ll_11= d2ll_11 + (f$fk*(f$fk_doublepr11+2*f$fk_doublepr13+f$fk_doublepr33) - (f$fk_pr1+f$fk_pr3)^2)/(f$fk^2)
        d2ll_12= d2ll_12 + (f$fk*(f$fk_doublepr12+f$fk_doublepr23) - f$fk_pr2*(f$fk_pr1+f$fk_pr3))/(f$fk^2)
        d2ll_22= d2ll_22 + (f$fk*f$fk_doublepr22 - f$fk_pr2*f$fk_pr2)/(f$fk^2)
    }
    return ( list(ll=logl, gradient=c(dll1,dll2), hessian=matrix(c(d2ll_11, d2ll_12, d2ll_12, d2ll_22),ncol=2)) )
}

## t0= starting point for Newton-Raphson
## fixit: jumps could be farther from root, but no way to fix this in 3D (bret gave idea, need to code it)
findMLE2D = function(seqx.dist, seq3.dist, seq4.dist, Q, d3x, t0=c(0.1,0.1), tol=0.0001, Nmax=10000, verbose=FALSE){
    if(verbose)
        print("entering findMLE...")
    tnew = rep(0,2) # will not save all sequence
    if(any(t0<0))
        t0[which(t0<0)] = 0.0001
    told = t0 ## fixit: later do a binary search before choosing t0
    error = 1
    i = 1
    while(error > tol & i < Nmax){
        if(verbose)
            print(told)
        f =loglik2D(seqx.dist, seq3.dist, seq4.dist,Q, told, d3x)
        gap = solve(f$hessian) %*% f$gradient
        tnew = told - gap
        while(any(tnew<0)){ #avoid negative BL
            gap = gap/2
            if(verbose)
                print("found negative candidate tnew")
            tnew = told - gap
        }
        ## after finding a positive candidate:
        f2 =loglik2D(seqx.dist, seq3.dist, seq4.dist, Q, tnew, d3x)
        error = max(abs(tnew-told))
        if(verbose)
            print(error)
        i = i+1
        told = tnew
    }
    if(i>=Nmax)
        warning("Newton-Rapshon did not converge")
    return ( list(t=tnew, obsInfo=f2$hessian) )
}


## simulates dxy,d4y because d3y is conditional on dxy+d3y=d3x
## need the sum (d3x) as input
simulateBranchLength.conditionalMultinorm = function(nsim,seqx.dist,seq3.dist, seq4.dist, Q, t0, d3x,verbose=FALSE){
    mu = findMLE2D(seqx.dist, seq3.dist, seq4.dist,Q, d3x, t0, verbose=verbose)
    Sigma = solve(mu$obsInfo)
    w = rmvnorm(nsim, mu$t, -Sigma)
    return ( list(t=w, mu=mu$t, sigma=-Sigma) )
}

## simulates dxy,d4y because d3y is conditional on dxy+d3y=d3x
## need the sum (d3x) as input
simulateBranchLength.conditionalGamma = function(nsim,seqx.dist,seq3.dist, seq4.dist, Q, t0, d3x,verbose=FALSE){
    mu = findMLE2D(seqx.dist, seq3.dist, seq4.dist,Q, d3x, t0, verbose=verbose)
    Sigma = solve(mu$obsInfo)
    R = chol(-Sigma)
    L = t(R)
    ## simulate beta
    alpha1 = mu$t[1]^2*(d3x-mu$t[1])/(d3x*L[1,1]^2) - mu$t[1]/d3x
    beta1 = mu$t[1]*(d3x-mu$t[1])^2/(d3x*L[1,1]^2) - (d3x-mu$t[1])/d3x
    b = rbeta(1,alpha1,beta1)
    t1 = b*d3x
    ## simulate T2|T1
    z1 = (t1-mu$t[1])/L[1,1]
    num = mu$t[2]+L[2,1]*z1
    alpha2 = num^2/L[2,2]^2
    lambda2 = num/L[2,2]^2
    t2 = rgamma(1,shape=alpha2,rate=lambda2)
    w = c(t1,t2)
    ## logdensity
    logdens = (alpha1-1)*log(t1)+(beta1-1)*log(d3x-t1)-(alpha1+beta1-1)*log(d3x)+(alpha2-1)*log(t2)-lambda2*t2
    return ( list(t=w, mu=mu$t, sigma=-Sigma, logdens = logdens) )
}


## ------------------------------------------------------------
## 5D likelihood

## p1 = P(A1|x1), column vector of size 4
## p2 = P(A2|x2), column vector of size 4
## p3 = P(A3|y3), column vector of size 4
## p4 = P(A4|y4), column vector of size 4
## order of bl: d1x,d2x,d3y,d4y,dxy
## returns fk, fk_prime1, fk_prime2, fk_prime3, fk_prime4,fk_prime5,
## fk_doubleprime11, fk_doubleprime12, fk_doubleprime13, fk_doubleprime22,
## fk_doubleprime23, fk_doubleprime33
## fixit: matrix multiplication can be more efficient with eigenvector decomp
fk5D = function(p1,p2,p3,p4,Q,t1,t2,t3,t4,t5){
    S1fn = Sfn(t1,p1,Q)
    S2fn = Sfn(t2,p2,Q)
    S3fn = Sfn(t3,p3,Q)
    S4fn = Sfn(t4,p4,Q)
    S1=diag(c(S1fn$S))
    S2=diag(c(S2fn$S))
    S3=diag(c(S3fn$S))
    S4=diag(c(S4fn$S))
    S1pr = diag(c(S1fn$Spr))
    S2pr = diag(c(S2fn$Spr))
    S3pr = diag(c(S3fn$Spr))
    S4pr = diag(c(S4fn$Spr))
    S1doublepr = diag(c(S1fn$Sdoublepr))
    S2doublepr = diag(c(S2fn$Sdoublepr))
    S3doublepr = diag(c(S3fn$Sdoublepr))
    S4doublepr = diag(c(S4fn$Sdoublepr))
    ##print("S4doublepr")
    ##print(S4doublepr)
    Pxy = matrixExp(Q,t5)
    fk = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% Pxy %*% S3 %*% S4 %*% rep(1,4)
    fk.pr1 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2 %*% Pxy %*% S3 %*% S4 %*% rep(1,4)
    fk.pr2 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2pr %*% Pxy %*% S3 %*% S4 %*% rep(1,4)
    fk.pr3 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% Pxy %*% S3pr %*% S4 %*% rep(1,4)
    fk.pr4 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% Pxy %*% S3 %*% S4pr %*% rep(1,4)
    fk.pr5 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% Q$Q %*% Pxy %*% S3 %*% S4 %*% rep(1,4)

    fk.doublepr11 = rep(1,4) %*% diag(Q$p) %*% S1doublepr %*% S2 %*% Pxy %*% S3 %*% S4 %*% rep(1,4)
    fk.doublepr22 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2doublepr %*% Pxy %*% S3 %*% S4 %*% rep(1,4)
    fk.doublepr33 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% Pxy %*% S3doublepr %*% S4 %*% rep(1,4)
    fk.doublepr44 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% Pxy %*% S3 %*% S4doublepr %*% rep(1,4)
    fk.doublepr55 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% Q$Q %*% Q$Q %*% Pxy %*% S3 %*% S4 %*% rep(1,4)

    ##print("fk.doublepr44")
    ##print(fk.doublepr44)

    fk.doublepr12 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2pr %*% Pxy %*% S3 %*% S4 %*% rep(1,4)
    fk.doublepr13 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2 %*% Pxy %*% S3pr %*% S4 %*% rep(1,4)
    fk.doublepr14 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2 %*% Pxy %*% S3 %*% S4pr %*% rep(1,4)
    fk.doublepr15 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2 %*% Q$Q %*% Pxy %*% S3 %*% S4 %*% rep(1,4)
    fk.doublepr23 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2pr %*% Pxy %*% S3pr %*% S4 %*% rep(1,4)
    fk.doublepr24 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2pr %*% Pxy %*% S3 %*% S4pr %*% rep(1,4)
    fk.doublepr25 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2pr %*% Q$Q %*% Pxy %*% S3 %*% S4 %*% rep(1,4)
    fk.doublepr34 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% Pxy %*% S3pr %*% S4pr %*% rep(1,4)
    fk.doublepr35 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% Q$Q %*% Pxy %*% S3pr %*% S4 %*% rep(1,4)
    fk.doublepr45 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% Q$Q %*% Pxy %*% S3 %*% S4pr %*% rep(1,4)

    return ( list(fk=as.numeric(fk), fk_pr1=as.numeric(fk.pr1),fk_pr2=as.numeric(fk.pr2),fk_pr3=as.numeric(fk.pr3), fk_pr4=as.numeric(fk.pr4), fk_pr5=as.numeric(fk.pr5),
                  fk_doublepr11= as.numeric(fk.doublepr11),fk_doublepr12= as.numeric(fk.doublepr12),fk_doublepr13= as.numeric(fk.doublepr13),
                  fk_doublepr14= as.numeric(fk.doublepr14),fk_doublepr15= as.numeric(fk.doublepr15),
                  fk_doublepr22= as.numeric(fk.doublepr22),fk_doublepr23= as.numeric(fk.doublepr23),fk_doublepr24= as.numeric(fk.doublepr24),fk_doublepr25= as.numeric(fk.doublepr25),
                  fk_doublepr33= as.numeric(fk.doublepr33), fk_doublepr34= as.numeric(fk.doublepr34), fk_doublepr35= as.numeric(fk.doublepr35),
                  fk_doublepr44= as.numeric(fk.doublepr44),fk_doublepr45= as.numeric(fk.doublepr45),fk_doublepr55= as.numeric(fk.doublepr55)
                  )
            )
}

## returns loglik, gradient (vector 5x1), hessian (matrix 5x5)
## t = vector 5x1 d1x,d2x,d3y,d4y,dxy
loglik5D = function(seq1.dist, seq2.dist, seq3.dist,seq4.dist,Q, t, verbose=FALSE){
    if(ncol(seq1.dist) != ncol(seq2.dist))
        stop("wrong number of sites")
    logl = 0
    dll1 = 0
    dll2 = 0
    dll3 = 0
    dll4 = 0
    dll5 = 0
    d2ll_11 = 0
    d2ll_12 = 0
    d2ll_13 = 0
    d2ll_14 = 0
    d2ll_15 = 0
    d2ll_22 = 0
    d2ll_23 = 0
    d2ll_24 = 0
    d2ll_25 = 0
    d2ll_33 = 0
    d2ll_34 = 0
    d2ll_35 = 0
    d2ll_44 = 0
    d2ll_45 = 0
    d2ll_55 = 0
    for(i in 1:ncol(seq1.dist)){
        f = fk5D(seq1.dist[,i],seq2.dist[,i],seq3.dist[,i],seq4.dist[,i],Q,t[1],t[2],t[3],t[4],t[5])
        logl = logl + log(f$fk)
        dll1 = dll1 + f$fk_pr1/f$fk
        dll2 = dll2 + f$fk_pr2/f$fk
        dll3 = dll3 + f$fk_pr3/f$fk
        dll4 = dll4 + f$fk_pr4/f$fk
        dll5 = dll5 + f$fk_pr5/f$fk
        d2ll_11= d2ll_11 + (f$fk*f$fk_doublepr11 - f$fk_pr1*f$fk_pr1)/(f$fk^2)
        d2ll_12= d2ll_12 + (f$fk*f$fk_doublepr12 - f$fk_pr1*f$fk_pr2)/(f$fk^2)
        d2ll_13= d2ll_13 + (f$fk*f$fk_doublepr13 - f$fk_pr1*f$fk_pr3)/(f$fk^2)
        d2ll_14= d2ll_14 + (f$fk*f$fk_doublepr14 - f$fk_pr1*f$fk_pr4)/(f$fk^2)
        d2ll_15= d2ll_15 + (f$fk*f$fk_doublepr15 - f$fk_pr1*f$fk_pr5)/(f$fk^2)
        d2ll_22= d2ll_22 + (f$fk*f$fk_doublepr22 - f$fk_pr2*f$fk_pr2)/(f$fk^2)
        d2ll_23= d2ll_23 + (f$fk*f$fk_doublepr23 - f$fk_pr2*f$fk_pr3)/(f$fk^2)
        d2ll_24= d2ll_24 + (f$fk*f$fk_doublepr24 - f$fk_pr2*f$fk_pr4)/(f$fk^2)
        d2ll_25= d2ll_25 + (f$fk*f$fk_doublepr25 - f$fk_pr2*f$fk_pr5)/(f$fk^2)
        d2ll_33= d2ll_33 + (f$fk*f$fk_doublepr33 - f$fk_pr3*f$fk_pr3)/(f$fk^2)
        d2ll_34= d2ll_34 + (f$fk*f$fk_doublepr34 - f$fk_pr3*f$fk_pr4)/(f$fk^2)
        d2ll_35= d2ll_35 + (f$fk*f$fk_doublepr35 - f$fk_pr3*f$fk_pr5)/(f$fk^2)
        d2ll_44= d2ll_44 + (f$fk*f$fk_doublepr44 - f$fk_pr4*f$fk_pr4)/(f$fk^2)
        d2ll_45= d2ll_45 + (f$fk*f$fk_doublepr45 - f$fk_pr4*f$fk_pr5)/(f$fk^2)
        d2ll_55= d2ll_55 + (f$fk*f$fk_doublepr55 - f$fk_pr5*f$fk_pr5)/(f$fk^2)
    }
    g = c(dll1,dll2,dll3,dll4,dll5)
    M = matrix(c(d2ll_11, d2ll_12, d2ll_13, d2ll_14,d2ll_15,
                      d2ll_12, d2ll_22, d2ll_23, d2ll_24, d2ll_25,
                      d2ll_13, d2ll_23, d2ll_33, d2ll_34, d2ll_35,
                      d2ll_14, d2ll_24, d2ll_34, d2ll_44, d2ll_45,
                      d2ll_15, d2ll_25, d2ll_35, d2ll_45, d2ll_55),ncol=5)
    if(verbose){
        print("gradient")
        print(g)
        print("vector for hessian")
        print("1")
        print(d2ll_11)
        print("2")
        print(d2ll_12)
        print("3")
        print(d2ll_13)
        print("4")
        print(d2ll_14)
        print("5")
        print(d2ll_15)
        print("6")
        print(d2ll_12)
        print("7")
        print(d2ll_22)
        print("8")
        print(d2ll_23)
        print("9")
        print(d2ll_24)
        print("10")
        print(d2ll_25)
        print("11")
        print(d2ll_13)
        print("12")
        print(d2ll_23)
        print("13")
        print(d2ll_33)
        print("14")
        print(d2ll_34)
        print("15")
        print(d2ll_35)
        print("16")
        print(d2ll_14)
        print("17")
        print(d2ll_24)
        print("18")
        print(d2ll_34)
        print("19")
        print(d2ll_44)
        print("20")
        print(d2ll_45)
        print("21")
        print(d2ll_15)
        print("22")
        print(d2ll_25)
        print("23")
        print(d2ll_35)
        print("24")
        print(d2ll_45)
        print("25")
        print(d2ll_55)
        print("hessian")
        print(M)
    }
    return ( list(ll=logl, gradient=g, hessian=M) )
}

## t0= starting point for Newton-Raphson
## fixit: jumps could be farther from root, but no way to fix this in 3D (bret gave idea, need to code it)
findMLE5D = function(seq1.dist, seq2.dist,seq3.dist, seq4.dist, Q, t0=rep(0.1,5), tol=0.0001, Nmax=10000, verbose=FALSE){
    if(verbose)
        print("entering findMLE...")
    tnew = rep(0,5) # will not save all sequence
    if(any(t0<0))
        t0[which(t0<0)] = 0.0001
    told = t0 ## fixit: later do a binary search before choosing t0
    error = 1
    i = 1
    while(error > tol & i < Nmax){
        if(verbose)
            print(told)
        f =loglik5D(seq1.dist, seq2.dist,seq3.dist, seq4.dist,Q, told, verbose=verbose)
        gap = solve(f$hessian) %*% f$gradient
        tnew = told - gap
        while(any(tnew<0)){ #avoid negative BL
            gap = gap/2
            if(verbose)
                print("found negative candidate tnew")
            tnew = told - gap
        }
        ## after finding a positive candidate:
        f2 =loglik5D(seq1.dist, seq2.dist,seq3.dist, seq4.dist, Q, tnew)
        error = max(abs(tnew-told))
        if(verbose)
            print(error)
        i = i+1
        told = tnew
    }
    if(i>=Nmax)
        warning("Newton-Rapshon did not converge")
    return ( list(t=tnew, obsInfo=f2$hessian) )
}

## simulates d1x,d2x,d3y,d4y,dxy jointly
simulateBranchLength.multinorm5D = function(nsim,seq1.dist,seq2.dist,seq3.dist, seq4.dist, Q, t0, verbose=FALSE){
    mu = findMLE5D(seq1.dist, seq2.dist,seq3.dist, seq4.dist,Q, t0, verbose=verbose)
    Sigma = solve(mu$obsInfo)
    if(verbose){
        print(mu$t)
        print(Sigma)
    }
    w = rmvnorm(nsim, mu$t, -Sigma)
    return ( list(t=w, mu=mu$t, sigma=-Sigma) )
}


## ----------------------------------------------------------------------------------------------
## 5 taxa case: conditional on two sums


## loglik ((1,2)x,(3,4)z,5)y
## fixit: need to modify to use P(Ak|x) and P(Bk|y)
## instead of d1x,d2x,dxy,dyz,d3z,d4z,d5y, seq?.dist
## maybe create another function that calls this one
gtr.log.lik.all.5taxa = function(d1x,d2x,dxy,dyz,d3z,d4z,d5y,seq1.dist,seq2.dist, seq3.dist,seq4.dist,seq5.dist,Q){
    suma = 0
    for(s in 1:ncol(seq1.dist)){
        lik12 = siteLik(d1x,d2x,seq1.dist[,s],seq2.dist[,s],Q)
        lik34 = siteLik(d3z,d4z,seq3.dist[,s],seq4.dist[,s],Q)
        Pxy = matrixExp(Q,dxy)
        Pyz = matrixExp(Q,dyz)
        P5y = matrixExp(Q,d5y)
        vx = Pxy %*% diag(lik12) %*% rep(1,4)
        vz = Pyz %*% diag(lik34) %*% rep(1,4)
        v5 = P5y %*% seq5.dist[,s]
        fk = rep(1,4) %*% diag(Q$p) %*% diag(c(vx)) %*% diag(c(vz)) %*% diag(c(v5)) %*% rep(1,4)
        suma = suma+log(fk)
    }
    return ( suma )
}


## p1 = P(A1|x1), column vector of size 4: seqx
## p2 = P(A2|x2), column vector of size 4: seqz
## p3 = P(A3|x3), column vector of size 4: seq5
## s1=d3x-d3z, s2=d5z, t=dyz
## returns fk, fk_prime1, fk_prime2, fk_prime3,
## fk_doubleprime11, fk_doubleprime12, fk_doubleprime13, fk_doubleprime22,
## fk_doubleprime23, fk_doubleprime33
## fixit: matrix multiplication can be more efficient with eigenvector decomp
fk1D = function(p1,p2,p3,Q,s1,s2,t){
    if(s1<t || s2<t)
        stop("error: sum smaller than summand")
    S1fn = Sfn(s1-t,p1,Q)
    S2fn = Sfn(t,p2,Q)
    S3fn = Sfn(s2-t,p3,Q)
    S1=diag(c(S1fn$S))
    S2=diag(c(S2fn$S))
    S3=diag(c(S3fn$S))
    S1pr = diag(-c(S1fn$Spr))
    S2pr = diag(c(S2fn$Spr))
    S3pr = diag(-c(S3fn$Spr))
    S1doublepr = diag(c(S1fn$Sdoublepr))
    S2doublepr = diag(c(S2fn$Sdoublepr))
    S3doublepr = diag(c(S3fn$Sdoublepr))
    fk = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% S3 %*% rep(1,4)
    fk.pr1 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2 %*% S3 %*% rep(1,4)
    fk.pr2 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2pr %*% S3 %*% rep(1,4)
    fk.pr3 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% S3pr %*% rep(1,4)
    fk.doublepr11 = rep(1,4) %*% diag(Q$p) %*% S1doublepr %*% S2 %*% S3 %*% rep(1,4)
    fk.doublepr22 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2doublepr %*% S3 %*% rep(1,4)
    fk.doublepr33 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2 %*% S3doublepr %*% rep(1,4)
    fk.doublepr12 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2pr %*% S3 %*% rep(1,4)
    fk.doublepr13 = rep(1,4) %*% diag(Q$p) %*% S1pr %*% S2 %*% S3pr %*% rep(1,4)
    fk.doublepr23 = rep(1,4) %*% diag(Q$p) %*% S1 %*% S2pr %*% S3pr %*% rep(1,4)

    return ( list(fk=as.numeric(fk), fk_pr1=as.numeric(fk.pr1),fk_pr2=as.numeric(fk.pr2),fk_pr3=as.numeric(fk.pr3),
                  fk_doublepr11= as.numeric(fk.doublepr11),fk_doublepr12= as.numeric(fk.doublepr12),fk_doublepr13= as.numeric(fk.doublepr13),
                  fk_doublepr22= as.numeric(fk.doublepr22),fk_doublepr23= as.numeric(fk.doublepr23),fk_doublepr33= as.numeric(fk.doublepr33)))
}

## returns loglik, llprime, lldoubleprime
## for the conditionl case of dyz
loglik1D = function(seqx.dist, seqz.dist, seq5.dist,Q, s1,s2,t){
    if(ncol(seqx.dist) != ncol(seqz.dist))
        stop("wrong number of sites")
    suma = 0
    suma2 = 0
    suma3 = 0
    for(i in 1:ncol(seqx.dist)){
        f = fk1D(seqx.dist[,i],seqz.dist[,i],seq5.dist[,i],Q,s1,s2,t)
        suma = suma + log(f$fk)
        suma2 = suma2 + (f$fk_pr1+f$fk_pr2+f$fk_pr3)/f$fk
        suma3 = suma3 + (f$fk*(f$fk_doublepr11 + 2*f$fk_doublepr12 + 2*f$fk_doublepr13 + 2*f$fk_doublepr23 + f$fk_doublepr22 + f$fk_doublepr33) - (f$fk_pr1 + f$fk_pr2+f$fk_pr3)^2)/(f$fk^2)
    }
    return ( list(ll=suma, ll_pr=suma2, ll_doublepr=suma3) )
}


## t0= starting point for Newton-Raphson
findMLE1D = function(seqx.dist, seqz.dist, seq5.dist, Q, s1,s2, t0=0.1, tol=0.0001, Nmax=10000, verbose=FALSE){
    if(verbose)
        print("entering findMLE...")
    t = rep(0,Nmax)
    if(t0<0)
        to=0.00001
    t[1] = t0 ## fixit: later do a binary search before choosing t0
    error = 1
    i = 1
    while(error > tol & i < Nmax){
        if(verbose)
            print(t[i])
        f =loglik1D(seqx.dist, seqz.dist, seq5.dist,Q, s1,s2,t[i])
        gap = f$ll_pr/f$ll_doublepr
        t[i+1] = t[i] - gap
        while(t[i+1]<0){ #avoid negative BL
            gap = gap/2
            if(verbose)
                print("found negative candidate t[i+1]")
            t[i+1] = t[i] - gap
        }
        ## after finding a positive candidate:
        f2 =loglik1D(seqx.dist, seqz.dist, seq5.dist,Q, s1,s2,t[i+1])
        if(verbose){
            print("f$ll, f2$ll")
            print(f$ll)
            print(f2$ll)
        }
        if(f$ll * f2$ll < 0){ #one positive, and one negative
            while(abs(f2$ll) > 2*abs(f$ll)){ #t[i+1] too far from root
                if(verbose)
                    print("t[i+1] positive, but too far from root")
                t[i+1] = t[i+1]/2
                f2 =loglik(seq1.dist, seq2.dist, Q, t[i+1])
            }
        }
        error = abs(t[i+1]-t[i])
        i = i+1
    }
    t = t[which(t>0)]
    if(i>=Nmax)
        warning("Newton-Rapshon did not converge")
    return ( list(t=t[length(t)], obsInfo=f2$ll_doublepr) )
}


## simulates dyz from
## need the sum (d3x) as input
simulateBranchLength.conditionalNorm = function(nsim,seqx.dist,seqz.dist, seq5.dist, Q, t0, s1,s2,verbose=FALSE){
    mu = findMLE1D(seqx.dist, seqz.dist, seq5.dist,Q, s1,s2, t0, verbose=verbose)
    w = rnorm(nsim, mu$t, -1/mu$obsInfo)
    return ( list(t=w, mu=mu$t, sigma=-1/mu$obsInfo) )
}


## aqui voy: tengo q hacer conditionalNorm, and later gtr.loglik.all.5taxa, then run and test
