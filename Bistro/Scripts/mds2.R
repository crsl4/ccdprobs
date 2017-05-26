## MDS plot to compare bistro run with MrBayes.

## Read in matrix
require(ggplot2)
require(viridis)
require(dplyr)


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


mds.phylo.combined = function(d,top,k=2,lambda=200)
{
  if ( !is.matrix(d) )
      stop("d needs to be a matrix")
  if (nrow(d) != nrow(top) )
      stop("distance matrix and list of trees need to have the same dimension")
  n = nrow(d)
  mds = cmdscale(d,k=k)
  if ( k >= 2 )
  {
      names(top) = "top"
      temp.trees = as.character(top$top)
      tab = rev(sort(table(temp.trees)))
      if ( nrow(tab) > 6 ){
          top.trees = names(tab)[1:5]
          temp.trees[ !(top$top %in% top.trees) ] = "other"
      }
      df = data.frame(x=mds[,1], y=mds[,2], w=exp(-lambda*d[,n]),top=temp.trees,
          who = c(rep("bistro",(n-2)/2),rep("mb",(n-2)/2),"mean.bistro","mean.mb"))
      ggplot(filter(df,w>0), aes(x=x,y=y,color=top,shape=who)) +
          geom_point() +
#              scale_color_viridis() +
                  geom_point(data=df[n-1,],color="red") +
                      geom_point(data=df[n,],color="red") +
                          theme_bw()+
                              theme(legend.position="top")
  }
}

meanDistance = function(d)
{
  return( mean(d[row(d)>col(d)]) )
}


#d.bistro = matrix(scan("bistro4---0-249.treeBL.distances"),100,100)
#d.mb = matrix(scan("cats-dogs-2.nex.run1.tre.distances"),100,100)
d.combined = matrix(scan("combined.distances"),102,102)
d.trees = read.table("combined.listTrees",header=FALSE)

mds.phylo.combined(d.combined,d.trees)
pdf("combined024.pdf")
mds.phylo.combined(d.combined,d.trees)
dev.off()
#mds.phylo(d.bistro,lambda=200)
#mds.phylo(d.mb)
#meanDistance(d.bistro)
#meanDistance(d.mb)
