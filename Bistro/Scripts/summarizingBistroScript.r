## script to analyze different bistro comparisons
## Claudia September 2016

## after alpha<1 gone
source("../../Scripts/readBistro.r")
bistro = readBistro("fixedQlik2")
plotBistro(bistro)




bistro[bistro$logl+bistro$logPrior< -9400,] ## 728
## in fixedQlik---500-749.treeBL, line 228:
tree ="((5:0.0890146,((12:0.0324841,1:0.1280738):0.0326449,(((10:0.0837049,(7:0.1706488,(8:0.0018742,9:0.0000000):0.0251494):0.0000033):0.0050635,11:0.1718760):0.0637874,6:0.0853249):0.7941643):0.0479423):0.0232329,(2:0.0448230,3:0.0307667):0.0142521,4:0.0871662);"
## it seems that this is a bad tree, could it be the 0.0000?

source("../../Scripts/readBistro.r")
bistro = readBistro("fixedQlik")

data1= read.table("fixedQlik---0-249.treeBL")
data2= read.table("fixedQlik---250-499.treeBL")
data3= read.table("fixedQlik---500-749.treeBL")
data4= read.table("fixedQlik---750-999.treeBL")
data=rbind(data1,data2,data3,data4)

f = function(x){
    grepl("0.0000000",x)
}

bistro$zero = f(as.character(data$V1))

require(viridis)
temp.tree = as.character(bistro$tree)
tab = with(bistro, rev(sort(table(tree))))
if ( length(levels(bistro$tree)) > 6 ) {
    ##        tab = with(bistro, rev(sort(table(tree))))
    top.trees = names(tab)[1:5]
    temp.tree[ !(bistro$tree %in% top.trees) ] = "other"
}
bistro$Tree = factor(temp.tree)
rm(temp.tree)
bistro$Rank = 6
n = length(names(tab))
for ( i in 1:6 )
    {
        if ( i < n )
            bistro$Rank[bistro$Tree==names(tab)[i]] = i
    }
bistro$Tree = with( bistro, reorder(Tree,Rank) )
viridis.scale = viridis(n=length(levels(bistro$Tree)))

my.plot = ggplot(bistro,aes(x=logl+logPrior,y=logQ+logTop+logBL,color=w,shape=Tree)) +
    geom_point() +
    scale_color_viridis() +
    coord_fixed()+
    xlim(c(-9175,-9075)) +
    theme(legend.position="top")
plot(my.plot)


## ========================================================
## comparing the std of bistro and mb
mb = c(1.098572e-04,3.097000e-05,2.969937e-05,6.043099e-05,9.758691e-05,3.449542e-04,9.141917e-06,1.846828e-06,1.398818e-06,7.264603e-05,2.181405e-04,3.624530e-05,1.769607e-05,3.830197e-05,1.879410e-05,2.159238e-04,2.695235e-03,3.226004e-04,1.257593e-04,2.197222e-05,1.234506e-04)
mb=sqrt(mb)
## weighted:
bistro=c(0.009437,0.006876,0.005003,0.007381,0.009572,0.019396,0.00275,0.001217,0.001023,0.00804,0.01211,0.00608,0.004251,0.005228,0.004287,0.013104,0.04518,0.016327,0.009719,0.004497,0.010514)
## sampled
bistro=c(0.01188,0.006337,0.006122,0.008422,0.010742,0.038234,0.010255,0.005689,0.000961,0.01229,0.016588,0.013138,0.005532,0.007896,0.004136,0.017659,0.053969,0.019247,0.014317,0.012532,0.015986)

plot(mb)
points(bistro,col="red")

plot(mb,bistro,xlim=c(0,0.06), ylim=c(0,0.06))
lines(x=c(0,0.06),y=c(0,0.06))



## want to compare randQpars (cats) and randQpars-art
## why is art much better in ESS than cats?
## cloud is thinner! seems to be a problem with logBL
source("../../Scripts/readBistro.r")
bistro = readBistro("randQpars")
bistro = readBistro("randQlik")
bistro = readBistro("fixedQlik")
bistro = readBistro("fixedQlik-art")
bistro = readBistro("randQlik-art")

plotBistro(bistro)

plotBistro(subset(bistro,logl+logPrior>-9250 & logQ+logTop+logBL<200))

head(bistro)
temp.tree = as.character(bistro$tree)
tab = with(bistro, rev(sort(table(tree))))
if ( length(levels(bistro$tree)) > 6 ) {
    top.trees = names(tab)[1:5]
    temp.tree[ !(bistro$tree %in% top.trees) ] = "other"
}
bistro$Tree = factor(temp.tree)
rm(temp.tree)
bistro$Rank = 6
n = length(names(tab))
for ( i in 1:6 )
    {
        if ( i < n )
            bistro$Rank[bistro$Tree==names(tab)[i]] = i
    }
bistro$Tree = with( bistro, reorder(Tree,Rank) )
viridis.scale = viridis(n=length(levels(bistro$Tree)))

my.plot = ggplot(bistro,aes(x=logl+logPrior,y=logQ+logTop+logBL,color=w,shape=Tree))
my.plot = ggplot(bistro,aes(x=logl,y=logQ,color=w,shape=Tree))
my.plot = ggplot(bistro,aes(x=logl,y=logTop,color=w,shape=Tree))
my.plot = ggplot(bistro,aes(x=logl,y=logBL,color=w,shape=Tree))
my.plot = ggplot(bistro,aes(x=1:length(logBL),y=logBL,color=w,shape=Tree))
my.plot = ggplot(bistro,aes(x=1:length(logl),y=logl,color=w,shape=Tree))



my.plot = my.plot +geom_point(size=2) +
    scale_color_viridis() +
    ##coord_fixed() +
    theme(legend.position="top")
plot(my.plot)

## compare randQ, fixedQ, parsimony, lik, no weight
source("../../Scripts/readBistro.r")

bistro = readBistro("randQpars")
bistro = readBistro("randQlik")
bistro = readBistro("randQno-wt")
bistro = readBistro("fixedQpars")
bistro = readBistro("fixedQlik")
plotBistro(bistro)

bistro = readBistro("randQpars-art")
bistro = readBistro("randQlik-art")
bistro = readBistro("randQno-wt-art")
bistro = readBistro("fixedQpars-art")
bistro = readBistro("fixedQlik-art")
plotBistro(bistro)

data = computeEntropy("randQpars")
data = computeEntropy("randQpars-art")
plotProb(data)
head(data)

## ========================================================
## compare randQ, fixedQ, parsimony, lik, no weight
source("../../Scripts/readBistro.r")

data = computeEntropy("randQpars")
data = computeEntropy("randQlik")
data = computeEntropy("randQno-wt")
data = computeEntropy("fixedQpars")
data = computeEntropy("fixedQlik")

plotProb(data)

## for cats
data$mb = NA
data$mb[data$tree == "(1,(2,(((((7,8),9),10),12),11)),(3,((4,5),6)));"] = 0.774196
data$mb[data$tree == "(1,2,((3,((4,5),6)),(((((7,8),9),10),12),11)));"] = 0.150627
p2 <- ggplot(aes(x=tree,y=prob), data=data) +
    geom_point() + ggtitle("black-weightProb, red-bootstrap, blue-parsimony wt, green-lik wt") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_point(aes(y=count, col="red")) + guides(col=FALSE) +
    geom_point(aes(y=parsimonyWt, col="blue")) +
    geom_point(aes(y=loglikWt, col="green")) +
    geom_point(aes(y=mb), shape=8)
plot(p2)




bistro = readBistro("randQpars")
bistro = readBistro("randQlik")
bistro = readBistro("randQno-wt")
bistro = readBistro("fixedQpars")
bistro = readBistro("fixedQlik")

plotBistro(bistro)

library(ape)
tre=read.tree(text="(1,(2,(((((7,8),9),10),12),11)),(3,((4,5),6)));")
tre2 = read.tree(text="(1,2,((3,((4,5),6)),(((((7,8),9),10),12),11)));")
layout(matrix(c(1,2),1,2))
plot(tre)
plot(tre2)

## ========================================================
## compare the entropy and plots for artiodactyl, cats and
## simulated data with 6 and 2 taxa

source("../../Scripts/readBistro.r")
setwd("../Artiodactyl")
setwd("../cats-dogs")
setwd("../Simulations")

data = computeEntropy("randQ05")
data = computeEntropy("randQ20")
data = computeEntropy("randQ05-6")
data = computeEntropy("randQ20-6")
data = computeEntropy("randQ05-12")
data = computeEntropy("randQ20-12")

plotProb(data)

bistro = readBistro("randQ05")
bistro = readBistro("randQ20")
bistro = readBistro("randQ05-6")
bistro = readBistro("randQ20-6")
bistro = readBistro("randQ05-12")
bistro = readBistro("randQ20-12")

plotBistro(bistro)

## ========================================================
## studying the effect of different parsimony scores on
## artiodactyl and cats
## Conclusion: no apparent effect

source("../../Scripts/readBistro.r")

bistro = readBistro("randQ01")
bistro = readBistro("randQ03")
bistro = readBistro("randQ05")
bistro = readBistro("randQ07")
bistro = readBistro("randQ09")

setwd("../cats-dogs")

plotBistro(bistro)
