## Bret Larget
## October 18
## Importance Sampling Test
## Use t(22) as proposal (10% extra variance)
## Target is multivariate independent standard normal

isamp = function(n,d,dof,bias=0)
{
  f = function(x,dof)
  {
    lognum = sum( dnorm(x,log=TRUE) )
    logden = sum( dt(x,dof,log=TRUE))
    return( c(lognum,logden) )
  }
  
  ## add a bit of bias
  s = matrix(rt(n*d,dof),n,d) + bias
  df = as.data.frame(t(apply(s,1,f,dof=dof))) ## unnormalized log densities
  names(df) = c("log.target","log.proposal")
  m = with(df, max(log.target-log.proposal))
  df$w = with(df, exp(log.target-log.proposal-m))
  df$w = with(df, w/sum(w))
  ess = with(df, 1 / sum( w^2 ))
  message(paste("ESS = ",round(ess),", %ESS = ",round(100*ess/n,2),"%",sep=""))

  require(ggplot2)
  require(viridis)
  p = ggplot(df,aes(x=log.target,y=log.proposal,color=w)) +
    geom_point() +
    scale_color_viridis()
  return(p)
}