## we want to know how the covariance matrix

## we need to modify the functions below to only use t1,t2 as variables
## we don't care about the sampling, we only want to get the obs Info matrix
## for fixed values of t3

## simulates d1x,d2x,d3x jointly
simulateBranchLength.multinorm = function(nsim,seq1.dist,seq2.dist, seq3.dist, Q, t0, verbose=FALSE){
    mu = findMLE3D(seq1.dist, seq2.dist, seq3.dist,Q, t0, verbose=verbose)
    Sigma = solve(mu$obsInfo)
    w = rmvnorm(nsim, mu$t, -Sigma)
    return ( list(t=w, mu=mu$t, sigma=-Sigma) )
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
