## Importance sampling example
## Dirichlet density target
## Use a wider scale for proposal

ddirichlet = function(x,alpha,log=FALSE)
{
  if ( log )
    return( sum( (alpha-1)*log(x) ) )
  else
    return( prod( x^(alpha-1) ) )
}

rdirichlet = function(alpha)
{
  x = rgamma(length(alpha),alpha,1)
  return(x/sum(x))
}

importance.sample = function(target,proposal,n=1000)
{
  logWt = numeric(n)
  logTarget = numeric(n)
  logProposal = numeric(n)
  for ( i in 1:n )
  {
    repeat {
      x = rdirichlet(proposal)
      logTarget[i] = ddirichlet(x,target,log=TRUE)
      logProposal[i] = ddirichlet(x,proposal,log=TRUE)
      logWt[i] = logTarget[i] - logProposal[i]
      if ( !is.na(logWt[i]) )
        break
    }
  }
  w = exp(logWt - max(logWt))
  w = w/sum(w)
  message(paste("ESS% =",round(100/sum(w^2)/n,2)))
  return( data.frame(logTarget,logProposal,logWt,w) )
}

#target = rgamma(10,2,1)
#target = target/sum(target) * 50
#proposal = proposal*0.8
#out = importance.sample(target,proposal)