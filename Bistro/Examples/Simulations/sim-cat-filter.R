## Filter out burnin-in plus trees that do not have the most common topology

## after simulated data with different sizes on the internal branch length:
rootname = "sim-cats-short1"
rootname = "sim-cats-short2"
rootname = "sim-cats-short3"
top1 = read.table(paste0(rootname,"-run1.top"))
top2 = read.table(paste0(rootname,"-run2.top"))
tre1 = read.table(paste0(rootname,"-run1.tre"))
tre2 = read.table(paste0(rootname,"-run2.tre"))

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

library(corrplot)
corrplot(cor(m))

####################
library(ape)
truetree = read.tree(text = "(cheetah:0.11850973637208513101,((((snow_leopard:0.04020488777776567990,leopard:0.03846672365840048818):0.01445254156731264929,tiger:0.07079306712878565000):0.01190623639760595223,clouded_leopard:0.10461902411036745619):0.04344639957238824457,(red_fox:0.11974055940327851810,(((coyote:0.00840558068745050208,(gray_wolf:0.00206050882985083861,dog:0.00185256446369396789):0.03205946058703370433):0.02609285257533808938,dhole:0.07049077201732806275):0.13276609809571754406,raccoon_dog:0.15542990325076813662):0.07955504846187926027):0.79869116234474835103):0.03995629108638096977,cat:0.03751335233479641956):0.0;")
plot(root(truetree,node=14))
nodelabels()

## want to study the correlation matrix and density plots
## of bistro for the main tree
source("../../Scripts/readBistro.r")

bistro = readBistro("eta11-bl-root")
data = readDataSort("eta11-bl-root")

bistro = readBistro("eta11-bl")
data = readDataSort("eta11-bl")

keep = which(bistro$tree=="(1,((((2,3),4),5),(6,(((7,(8,9)),10),11))),12);")
head(data)
tre=data$V1[keep]
require(ape)

tree=read.tree(text=as.character(tre[1]))
plot(tree)
length(tree$edge.length)

N = length(tre)
m = matrix(0,N,ncol=21) ## 21 works in this case!!!
for ( i in 1:N )
{
  tree = read.tree(text=as.character(tre[i]))
  m[i,] = tree$edge.length
}

colnames(m) = paste(tree$edge[,2],tree$edge[,1],sep=":")
round(cor(m),1)

plot(tree,use.edge.length = FALSE)
nodelabels()

library(corrplot)
corrplot(cor(m))

pdf("bistro-bl-densities.pdf")
for(i in 1:ncol(m)){
    df = data.frame(x=m[,i])
    plot(ggplot(df,aes(x=x))+geom_density() + ggtitle(colnames(m)[i]))
}
dev.off()

## bivariate plots
require(ggplot2)
df = data.frame(x=m[,2],y=m[,21])
ggplot(df,aes(x=x,y=y)) + geom_point()


## adding info on which was the root
datarand = readData("eta11-bl-root")
treesrand = datarand$V1[keep]
rt=rep(NA,length(treesrand))
for(i in 1:length(treesrand)){
    tree = read.tree(text=as.character(treesrand[i]))
    rt[i] = isRootOnCats(tree)
}

df = data.frame(x=m[,2],y=m[,21], r = rt)
ggplot(df,aes(x=x,y=y, color=r)) + geom_point()


## want to study the correlation matrix and density plots
## of MrBayes branch lengths for the main tree
## after running badger:
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


## plot each branch length
library(ggplot2)
ggplot(as.data.frame(m),aes(x="15:14"))+geom_density()
colnames(m)

pdf("mb-bl-densities.pdf")
for(i in 1:ncol(m)){
    df = data.frame(x=m[,i])
    plot(ggplot(df,aes(x=x))+geom_density() + ggtitle(colnames(m)[i]))
}
dev.off()

ggplot(data.frame(x=m[,1],y=m[,2]), aes(x=x,y=y)) + geom_point()
