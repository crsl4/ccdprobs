## Filter out burnin-in plus trees that do not have the most common topology

top1 = read.table("sim-cats-dogs.run1.top")
top2 = read.table("sim-cats-dogs.run2.top")
tre1 = read.table("sim-cats-dogs.run1.tre")
tre2 = read.table("sim-cats-dogs.run2.tre")

top = c(as.character(top1[-(1:2001),]),as.character(top2[-(1:2001),]))
tre = c(as.character(tre1[-(1:2001),]),as.character(tre2[-(1:2001),]))

keep = which(top=="((((((2,3),4),5),(6,(((7,(8,9)),10),11))),12),1);")

tre = tre[keep]

require(ape)
N = length(tre)
m = matrix(0,N,ncol=22) ## 22 works in this case!!!
for ( i in 1:N )
{
  tree = read.tree(text=tre[i])
  m[i,] = tree$edge.length
}

colnames(m) = paste(tree$edge[,2],tree$edge[,1],sep=":")

## get rid of first column which is a fake 0 length edge from root
m = m[,-1]

round(cor(m),1)

plot(tree,use.edge.length = FALSE)
nodelabels()
edgelabels()

