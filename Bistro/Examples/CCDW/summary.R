#### Combined graphs, all code
#### Bret
#### 25 October 2016

## mrbayes
top1 = read.table("ccdw.run1.top")
top2 = read.table("ccdw.run2.top")
tre1 = read.table("ccdw.run1.tre")
tre2 = read.table("ccdw.run2.tre")

## remove burnin
top = c(as.character(top1[-(1:2001),]),as.character(top2[-(1:2001),]))
tre = c(as.character(tre1[-(1:2001),]),as.character(tre2[-(1:2001),]))

## best tree
keep = which(top=="(1,2,(3,4));")

## sample 1000 for better plots
keep = sort(sample(keep,1000))

tre.mb = tre[keep]

## remove big stuff
rm(top1,top2,tre1,tre2,top,tre,keep)

require(ape)
N.mb = length(tre.mb) ### 1000 now!!!

## get number of edges
tree = read.tree(text=tre.mb[1])
m.mb = matrix(0,N.mb,length(tree$edge.length))
colnames(m.mb) = paste("b",tree$edge[,2],tree$edge[,1],sep=".")

for ( i in 1:N.mb )
{
  tree = read.tree(text=tre.mb[i])
  m.mb[i,] = tree$edge.length
}

### now do bistro
source("../../Scripts/readBistro.r")

bistro = readBistro("test1")
data = readDataSort("test1")


keep = which(bistro$tree=="(1,2,(3,4));")
tre.bistro = data$V1[keep]
require(ape)
tree=read.tree(text=as.character(tre.bistro[1]))


N.bistro = length(tre.bistro)
## get number of edges
m.bistro = matrix(0,N.bistro,length(tree$edge.length))
colnames(m.bistro) = paste("b",tree$edge[,2],tree$edge[,1],sep=".")

## check same names!!!
##cbind(colnames(m.mb),colnames(m.bistro))
all( colnames(m.mb) == colnames(m.bistro) )

for ( i in 1:N.bistro )
{
  tree = read.tree(text=as.character(tre.bistro[i]))
  m.bistro[i,] = tree$edge.length
}

## combine data for combined plots
m = rbind(m.bistro,m.mb)
df = data.frame(m)
df$set = factor( c(rep("Bistro",N.bistro),rep("MrBayes",N.mb)) )

## make the plots
adjEdges = getAdjacentEdges(tree)

require(viridis)
require(ggplot2)
require(dplyr)

pdf("scatter-plots-test4.pdf")
vpal = viridis(2,end=0.8)
for( i in 1:nrow(adjEdges))
{
  i1 = adjEdges[i,1]
  i2 = adjEdges[i,2]
  median.bistro.1 = median( drop(as.matrix(filter(df,set=="Bistro") %>% select(i1))) )
  median.bistro.2 = median( drop(as.matrix(filter(df,set=="Bistro") %>% select(i2))) )
  median.mb.1 = median( drop(as.matrix(filter(df,set=="MrBayes") %>% select(i1))) )
  median.mb.2 = median( drop(as.matrix(filter(df,set=="MrBayes") %>% select(i2))) )

  plot(
    ggplot(df,aes(x=df[,adjEdges[i,1]],
                  y=df[,adjEdges[i,2]],
                  color=set)
           ) +
      geom_point(alpha=0.5) +
      scale_color_manual(values=vpal) +
      geom_vline(xintercept=median.bistro.1,color=vpal[1]) +
      geom_vline(xintercept=median.mb.1,color=vpal[2]) +
      geom_hline(yintercept=median.bistro.2,color=vpal[1]) +
      geom_hline(yintercept=median.mb.2,color=vpal[2]) +
   #   facet_grid(set ~ .) +
      ggtitle(paste(names(df)[adjEdges[i,1]],names(df)[adjEdges[i,2]])) +
      theme_bw()
  )
}
dev.off()

pdf("test1.pdf")
vpal = viridis(2,end=0.8)
for(i in 1:ncol(m))
{
    median.bistro = median( drop(as.matrix(filter(df,set=="Bistro") %>% select(i))) )
    median.mb = median( drop(as.matrix(filter(df,set=="MrBayes") %>% select(i))) )

    plot(ggplot(df,aes(x=df[,i],col=set))+geom_density() +
         scale_color_manual(values=vpal) +
         geom_vline(xintercept=median.bistro,color=vpal[1]) +
         geom_vline(xintercept=median.mb,color=vpal[2]) +
         ggtitle(paste(names(df)[i])) +
         theme_bw())

}
dev.off()

