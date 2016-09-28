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
    file2 = read.table(paste0(stem,".topCounts"), header=TRUE)
    file2 <- within(file2,count <- count/sum(count))
    file2 <- within(file2,parsimonyWt <- parsimonyWt/sum(parsimonyWt))
    file2 <- within(file2,loglikWt <- loglikWt/sum(loglikWt))
    entropyBootstrap = entropy.empirical(file2$count)
    entropyParsimony = entropy.empirical(file2$parsimonyWt)
    entropyLik = entropy.empirical(file2$loglikWt)
    entropyWeight = entropy.empirical(file1$prob)
    print(paste("Entropy bootstrap sample:", entropyBootstrap))
    print(paste("Entropy bootstrap sample after parsimony weight:", entropyParsimony))
    print(paste("Entropy bootstrap sample after lik weight:", entropyLik))
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
        geom_point() + ggtitle("black-weightProb, red-bootstrap, blue-parsimony wt, green-lik wt") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        geom_point(aes(y=count, col="red")) + guides(col=FALSE) +
        geom_point(aes(y=parsimonyWt, col="blue")) +
        geom_point(aes(y=loglikWt, col="green"))
    plot(p2)
    return(invisible(p2))
}
