## Bret Larget
## Janaury 29, 2016
## Code to test branch length proposal idea

## Code to make a Q matrix and its spectral decomposition
##   input is:
##     x which is the vector for parameters in lower triangle of rate matrix (r_AC, r_AG, r_AT, r_CG, r_CT, r_GT)
##     p which is the stationary distribution
##     n which is the number of states
##     a flag which will rescale so mean number of transitions per unit time is one
##   output is a list with the Q matrix itself, matrices for the eigenvectors and their inverse, and a vector of the eigenvalues
##

makeQ = function(x,p,n,rescale=FALSE) {
  Q = matrix(0,n,n)
  Q[row(Q) > col(Q)] = x
  Q = Q + t(Q)
  Q = Q %*% diag(p)
  diag(Q) = -apply(Q,1,sum)
  if ( rescale ) {
      mu = sum( -p*diag(Q) )
      Q = Q/mu
  }
  eig = eigen(Q)
  V = eig$vectors
  lambda = eig$values
  Vinv = solve(V)
  return( list(Q=Q,V=V,Vinv=Vinv,lambda=lambda,p=p) )
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

## Calculate GTR log-likelihood at a sequence s (times) for an observed matrix x and matrix Q
gtr.log.like = function(x,s,Q) {
    ns = length(s)
    log.like = numeric( ns )
    for ( i in 1:ns ) {
        P = matrixExp(Q,s[i])
        log.like[i] = sum( x * log(diag(Q$p) %*% P) )
    }
    return( log.like )
}

## Simulate branch length using JC
simulateBranchLength = function(nsim,x, eta=0.9) {
    n = sum(x)
    changes = n - sum(diag(x))
    theta = 4*changes/3
    w = -0.75 * log( rbeta(nsim,eta*(n - theta),eta*theta) )
    return( w )
}

## Plot likelihood and density from simulated sample
comparePlot = function(x,s,Q,nsim=10000,eta=0.5) {
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
    ## now get a sample of JCish points to compare
    jc = simulateBranchLength(nsim,x,eta)
    ## density estimate
    d = density(jc)
    ## another data frame
    df.est = data.frame(x=d$x,y=d$y)
    ## plot it
    p1 = ggplot(df.true, aes(x=s,y=y)) +
        geom_line(color="blue") +
        geom_line(aes(x=x,y=y),data=df.est,color="red",linetype="dashed") +
        xlab('branch length') +
        ylab('densities') +
            ggtitle('Blue = GTR Likelihood, Red = JC Approximation')
    return( p1 )
}

## doit
Q.random = randomQ(4)
nsites = 500
branch.length = 0.15

doit = function(nsites, branch.length, eta=0.5, nsim=10000, delta = 0.001) {
    s = seq(delta,2*branch.length,delta)
    Q = randomQ(4,rescale=TRUE)
    print(round(Q$Q,4))
    print(round(Q$p,4))
    x = simulateSequenceSummary(nsites,Q,branch.length)
    p1 = comparePlot(x,s,Q,nsim,eta)
    plot(p1)
}
