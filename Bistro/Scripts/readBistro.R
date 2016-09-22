## read in the different .out files

readBistro = function(stem)
{
    files = list.files(pattern=paste0("^",stem,"---.*\\.out"))
    bistro = data.frame()
    for ( f in files ) {
        temp = read.table(f,header=TRUE)
        bistro = rbind(bistro,temp)
        rm(temp)
    }
    bistro$w = with(bistro, exp(logWt - max(logWt)))
    bistro$w = with(bistro, w/sum(w))
    return ( bistro )
}

plotBistro = function(bistro) {
    require(ggplot2)
    require(viridis)
    temp.tree = as.character(bistro$tree)
    if ( length(levels(bistro$tree)) > 6 ) {
        tab = with(bistro, rev(sort(table(tree))))
        top.trees = names(tab)[1:5]
        temp.tree[ !(bistro$tree %in% top.trees) ] = "other"
    }
    bistro$Tree = factor(temp.tree)
    rm(temp.tree)
    viridis.scale = viridis(n=length(levels(bistro$Tree)))
    my.plot = ggplot(bistro,aes(x=logl+logPrior,y=logQ+logTop+logBL,color=w,shape=Tree)) + geom_point() + scale_color_viridis() + coord_fixed()
    plot(my.plot)
    return(invisible(my.plot))
}
