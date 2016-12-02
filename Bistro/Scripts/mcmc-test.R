## mcmc-test.R
## Bret Larget
## November 18, 2016

## f(x) = dexp(x,1)
## proposal is y = r*x where r = exp(lambda*(U-0.5)) where lambda = 0.19 and U ~ Uniform(0,1)
## for a single edge, Jacobian is r, not r^3 !!!

mcmc = function(ngen=100000,lambda=0.19)
{
    curr = 1
    x = numeric(ngen)
    logh.curr = dexp(curr,rate=1,log=TRUE)
    accept = numeric(ngen)
    r = exp(lambda*(runif(ngen)-0.5))
    u = runif(ngen)
    for ( i in 1:ngen )
    {
        prop = r[i]*curr
        logh.prop = dexp(prop,rate=1,log=TRUE)
        p = exp( min( c(0, logh.prop - logh.curr + log(r[i])) ) )
        accept[i] = p
        if ( u[i] < p ) ## accept
        {
            curr = prop
            logh.curr = logh.prop
        }
        x[i] = curr
    }
    message(paste("Mean acceptance =",round(mean(accept),4)))
    return( data.frame(x=x,accept=accept,density=dexp(x,rate=1), r=r) )
}
