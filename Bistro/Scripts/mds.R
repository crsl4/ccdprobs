## Make mds plot of distance data

## Read in matrix
require(ggplot2)
require(viridis)
require(dplyr)

d = matrix(scan("mds.mat"),101,101)
out = cmdscale(d,k=2)

top = read.table("mds.bootstrap",header=FALSE)
names(top) = "top"
temp.trees = as.character(top$top)
tab = rev(sort(table(temp.trees)))
if ( nrow(tab) > 6 )
{
  top.trees = names(tab)[1:5]
  temp.trees[ !(top$top %in% top.trees) ] = "other"
}
trees = data.frame(top=c(temp.trees,"(1,(2,((3,((4,5),6)),(((((7,8),9),10),11),12))));"))
trees = trees %>%
  mutate(rank=6)
n = nrow(tab)
for ( i in 1:6 )
{
  if ( i < n )
  {
    trees$rank[trees$top==names(tab)[i]] = i
  }
}

df = data.frame(x=out[,1],y=out[,2],w=exp(-200*d[,101]),top=trees$top,rank=trees$rank)
df$top = reorder(df$top,df$rank)
df = df %>%
  mutate(w = w / sum(w[-101]))
df$w[101] = 0
ggplot(filter(df,w>0), aes(x=x,y=y,color=w,shape=top)) +
  geom_point() +
  scale_color_viridis() +
  geom_point(data=df[101,],color="red") +
  theme_bw()

df2 = data.frame(dist=d[,101],mds.dist=sqrt(df$x^2+df$y^2),top=df$top,rank=df$rank)
df2$top = reorder(df2$top,df2$rank)
ggplot(df2, aes(x=dist,y=mds.dist,shape=top)) +
  geom_point()
