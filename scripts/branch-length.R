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
    return( w )
}

## Now try the Tamura-Nei way
simulateBranchLength.tn = function(nsim,x,eta=0.9) {
    n = sum(x)
    prop.ag = (x[1,3] + x[3,1]) / n
    prop.ct = (x[2,4] + x[4,2]) / n
    prop.tv = (x[1,2] + x[1,4] + x[2,1] + x[2,3] + x[3,2] + x[3,4] + x[4,1] + x[4,3]) / n
    print(paste("prop.ag",prop.ag,"prop.ct",prop.ct,"prop.tv",prop.tv))
    p.est = (apply(x,1,sum) + apply(x,2,sum)) / (2*n)
    p.a = p.est[1]
    p.c = p.est[2]
    p.g = p.est[3]
    p.t = p.est[4]
    p.r = sum(p.est[c(1,3)])
    p.y = sum(p.est[c(2,4)])
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
    print(paste("mu",mu,"v",v))
    w = rgamma(nsim,mu^2/v,mu/v)
    return( w )
}

## Plot likelihood and density from simulated sample
comparePlot = function(x,s,Q,nsim=10000,eta.jc=0.5,eta.tn=0.9) {
    require(ggplot2)
    delta = s[2] - s[1] # (assumes that times s are a regular sequence)
    ## compute the likelihood using counts x and generator Q for each s
    log.like = gtr.log.like(x,s,Q)
    ## turn into a density estimate proportional to likelihood, so truncated uniform prior in t
    log.like = log.like - max(log.like) # going to rescale anyway, may exp() accurate
    y = exp(log.like)
    ## make y sum to one
    y = y / sum(y)
    ## density at each point is y divided by delta
    y = y / delta
    ## Make a data frame
    df.true = data.frame(s,y)
    ## GTR Q optimized over data
    Qlist.gtr = optim.gtr(x,Q$r) ## sometimes: Error in eigen(S, symmetric = TRUE) (from branch-length.R#27) : infinite or missing values in 'x'
    print(Qlist.gtr$branch.length)
    ## compute this density also
    log.like2 = gtr.log.like(x,s,Qlist.gtr$Q)
    log.like2 = log.like2 - max(log.like2)
    y2 = exp(log.like2)
    y2 = y2/sum(y2)
    y2 = y2 / delta
    df.gtr = data.frame(s,y2)

    ## now get a sample of JCish points to compare
    jc = simulateBranchLength.jc(nsim,x,eta.jc)
    ## density estimate
    d.jc = density(jc)
    ## another data frame
    df.jc = data.frame(x=d.jc$x,y=d.jc$y)
    ## TN estimate
    tn = simulateBranchLength.tn(nsim,x,eta.tn)
    d.tn = density(tn)
    df.tn = data.frame(x=d.tn$x,y=d.tn$y)
    ## plot it
    p1 = ggplot(df.true, aes(x=s,y=y)) +
        geom_line(color="blue") +
        geom_line(aes(x=x,y=y),data=df.jc,color="red",linetype="dashed") +
        geom_line(aes(x=x,y=y),data=df.tn,color="darkgreen",linetype="dashed") +
        geom_line(aes(x=s,y=y2),data=df.gtr,color="gold") +
        xlab('branch length') +
        ylab('densities') +
            ggtitle('Blue = TrueQ, Red = JC, Green = TN, Gold = GTR')
    return( p1 )
}


## doit
#Q.random = randomQ(4)
nsites = 500
branch.length = 0.15

doit = function(nsites, branch.length, eta.jc=0.5, eta.tn=0.8, nsim=10000, delta = 0.001) {
    s = seq(delta,2*branch.length,delta)
    Q = randomQ(4,rescale=TRUE)
    print(round(Q$Q,4))
    print(round(Q$p,4))
    print(min(diag(Q$Q)))
    x = simulateSequenceSummary(nsites,Q,branch.length)
    print(x)
    p1 = comparePlot(x,s,Q,nsim,eta.jc,eta.tn)
#    plot(p1)
    return (p1)
}

#for(i in 1:100){
#    doit(nsites,branch.length)}

 ##        [,1]    [,2]    [,3]    [,4]
## [1,] -4.8076  3.6641  0.9282  0.2153
## [2,]  0.0038 -0.6253  0.5871  0.0344
## [3,]  0.0040  2.4673 -4.2804  1.8091
## [4,]  0.0003  0.0393  0.4921 -0.5316
## [1] 0.0005 0.4731 0.1126 0.4139
## [1] -4.807612
##      [,1] [,2] [,3] [,4]
## [1,]    0    0    0    0
## [2,]    0  223   15    7
## [3,]    0   22   42    9
## [4,]    0    4   10  168
## Error in eigen(S, symmetric = TRUE) (from branch-length.R#27) : infinite or missing values in 'x'
## In addition: There were 14 warnings (use warnings() to see them)
## > warnings()
## Warning messages:
## 1: In log(diag(Q$p) %*% P) : NaNs produced
## ...

## NOTE: error when row of zeros in x (counts), p has a very small number 0.0005

#################
# after symmetric=FALSE in all, new error:
#Error in optim(par = theta0, fn = logl.gtr, p = p, x = x, control = list(fnscale = -1)) (from branch-length.R#106) :
#  function cannot be evaluated at initial parameters

##         [,1]    [,2]    [,3]    [,4]
## [1,] -1.1240  0.8402  0.2836  0.0002
## [2,]  0.5079 -0.9974  0.4854  0.0041
## [3,]  0.2355  0.6667 -0.9027  0.0005
## [4,]  0.0238  0.6653  0.0617 -0.7508
## [1] 0.2585 0.4276 0.3113 0.0026
## [1] -1.124033
##      [,1] [,2] [,3] [,4]
## [1,]  135    9    1    0
## [2,]    9  177   12    0
## [3,]    8   17  132    0
## [4,]    0    0    0    0

## NOTE: same error when row of zeros in x (counts), p has a very small number 0.0026
