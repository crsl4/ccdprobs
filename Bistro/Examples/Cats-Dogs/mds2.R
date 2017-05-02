## MDS plot to compare bistro run with MrBayes.

## Read in matrix
require(ggplot2)
require(viridis)
require(dplyr)

d.bistro = matrix(scan("bistro4---0-249.treeBL.distances"),101,101)
d.mb = matrix(scan("cats-dogs-2.nex.run1.tre.distances"),101,101)

mds.phylo = function(d,k=2,lambda=200)
{
  if ( !is.matrix(d) )
    stop("d needs to be a matrix")
  n = nrow(d)
  mds = cmdscale(d,k=k)
  if ( k >= 2 )
  {
    df = data.frame(x=mds[,1], y=mds[,2], w=exp(-lambda*d[,n]))
    ggplot(filter(df,w>0), aes(x=x,y=y,color=w)) +
      geom_point() +
      scale_color_viridis() +
      geom_point(data=df[n,],color="red") +
      theme_bw()
  }
}

meanDistance = function(d)
{
  return( mean(d[row(d)>col(d)]) )
}

mds.phylo(d.bistro,lambda=1)
mds.phylo(d.mb)
meanDistance(d.bistro)
meanDistance(d.mb)
