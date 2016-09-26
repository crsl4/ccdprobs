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

computeEntropy = function(stem){
    require(entropy)
    file1 = read.table(paste0(stem,".topPP"), header=TRUE)
    file2 = read.table(paste0(stem,".topCounts"), header=FALSE)
    colnames(file2) <- c("tree", "bootstrapCount", "parsimonyWeight", "parsimonyScore", "difference")
    file2 <- within(file2,bootstrapCount <- bootstrapCount/sum(bootstrapCount))
    file2 <- within(file2,parsimonyWeight <- parsimonyWeight/sum(parsimonyWeight))
    entropyBootstrap = entropy.empirical(file2$bootstrapCount)
    entropyParsimony = entropy.empirical(file2$parsimonyWeight)
    entropyWeight = entropy.empirical(file1$prob)
    print(paste("Entropy bootstrap sample:", entropyBootstrap))
    print(paste("Entropy bootstrap sample after parsimony weight:", entropyParsimony))
    print(paste("Entropy weighted sample:", entropyWeight))
    m=merge(file1,file2,by="tree", all=TRUE)
    return (m)
}


plotBistro = function(bistro) {
    require(ggplot2)
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
    my.plot = ggplot(bistro,aes(x=logl+logPrior,y=logQ+logTop+logBL,color=w,shape=Tree)) + geom_point() + scale_color_viridis() + coord_fixed()
    plot(my.plot)
    return(invisible(my.plot))
}

plotProb = function(data){
    require(ggplot2)
    p2 <- ggplot(aes(x=tree,y=prob), data=data) +
    geom_point() + ylab("black-weightProb, red-bootstrap, blue-parsimony weight") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_point(aes(y=bootstrapCount, col="red")) + guides(col=FALSE) +
    geom_point(aes(y=parsimonyWeight, col="blue"))
    plot(p2)
    return(invisible(p2))
}
