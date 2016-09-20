## read in the different .out files
## kluge here assumes 4 cores with 0-249, 250-499, 500-749, 750-999

readBistro = function (stem)
{
    bistro.1 = read.table(paste(stem,"---0-249.out",sep=""),header=TRUE)
    bistro.2 = read.table(paste(stem,"---250-499.out",sep=""),header=TRUE)
    bistro.3 = read.table(paste(stem,"---500-749.out",sep=""),header=TRUE)
    bistro.4 = read.table(paste(stem,"---750-999.out",sep=""),header=TRUE)
    bistro = rbind(bistro.1,bistro.2,bistro.3,bistro.4)
    rm(bistro.1,bistro.2,bistro.3,bistro.4)
    bistro$w = with(bistro, exp(logWt - max(logWt)))
    bistro$w = with(bistro, w/sum(w))
    bistro$logQ = with(bistro, logl + logPrior - logTop - logProp - logWt)
    require(ggplot2)
    require(viridis)
    viridis.scale = viridis(n=length(levels(bistro$tree)))
##  my.plot = ggplot(bistro,aes(x=logl+logPrior,y=logQ+logTop+logProp,color=w,shape=tree)) + geom_point() + scale_color_viridis() + coord_fixed()
    my.plot = ggplot(bistro,aes(x=logl+logPrior,y=logQ+logTop+logProp,color=w)) +
        geom_point() +
        scale_color_viridis() +
##        guides(color=FALSE) +
        coord_fixed()
    plot(my.plot)
    return(bistro)
}
