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

readData = function(stem)
{
    files = list.files(pattern=paste0("^",stem,"---.*\\.treeBL$"))
    bistro = data.frame()
    for ( f in files ) {
        temp = read.table(f,header=FALSE)
        bistro = rbind(bistro,temp)
        rm(temp)
    }
    return ( bistro )
}

readDataSort = function(stem)
{
    files = list.files(pattern=paste0("^",stem,"---.*\\.treeBLsort"))
    bistro = data.frame()
    for ( f in files ) {
        temp = read.table(f,header=FALSE)
        bistro = rbind(bistro,temp)
        rm(temp)
    }
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


plotBistro = function(bistro,filterData=FALSE) {
    require(ggplot2)
    require(viridis)
    require(dplyr)
    if(filterData){
        bistro = bistro %>% filter(logWt > max(logWt)-20)
    }
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
    for ( i in 1:6 ){
        if ( i < n )
            bistro$Rank[bistro$Tree==names(tab)[i]] = i
    }
    bistro$Tree = with( bistro, reorder(Tree,Rank) )
##    viridis.scale = viridis(n=length(levels(bistro$Tree)))
    my.plot = ggplot(bistro,aes(x=logl+logPrior,y=logQ+logTop+logBL,color=w,shape=Tree)) +
        geom_point() +
            scale_color_viridis() +
                    theme(legend.position="top")
    if(filterData){
        my.plot = my.plot + coord_fixed()
    }
    lim = range(bistro$logQ+bistro$logTop+bistro$logBL)
    my.plot = my.plot +
        geom_abline(intercept=seq(lim[1], lim[2]),
                    slope=1,
                    colour="white")
    plot(my.plot)
    return(invisible(my.plot))
}

plotBistro2 = function(bistro,filterData=FALSE) {
  require(dplyr)
  require(ggplot2)
  require(viridis)
  if(filterData){
      bistro = bistro %>% filter(logWt > max(logWt)-20)
  }
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
  ##    viridis.scale = viridis(n=length(levels(bistro$Tree)))
  bistro = bistro %>% mutate(logTarget = logl + logPrior, logProposal = logBL + logTop + logQ)
  my.plot = ggplot(bistro,aes(x=logProposal,y=logTarget-logProposal,color=w,shape=Tree)) +
    geom_point() +
    scale_color_viridis() +
    theme(legend.position="top")
  if(filterData){
      my.plot = my.plot + coord_fixed()
  }
  lim = range(bistro$logQ+bistro$logTop+bistro$logBL)
  print(lim)
  my.plot = my.plot +
      geom_abline(intercept=seq(lim[1], lim[2]),
                  slope=1,
                  colour="white")
  plot(my.plot)
  return(invisible(my.plot))
}

plotProb = function(data){
    require(ggplot2)
    p2 <- ggplot(aes(x=tree,y=prob), data=data) +
    geom_point() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("black-weightProb, red-bootstrap, blue-parsimony wt, green-lik wt") +
    geom_point(aes(x=tree,y=count), col="red",data=data) + guides(col=FALSE) +
    geom_point(aes(x=tree,y=parsimonyWt), col="blue",data=data) +
    geom_point(aes(x=tree,y=loglikWt), col="green", data=data)
    plot(p2)
    return(invisible(p2))
}


## to identify the bad trees, whether they have a bad root
getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

getRoot <- function(tree){
    totalN = tree$Nnode + length(tree$tip.label)
    al = seq(1,totalN,by=1)
    root = 0
    for(i in 1:totalN){
        desc = getDescendants(tree, node=i)
        if(length(desc) > 0 & length(desc) == length(al[-i]))
            if(all(sort(desc) == al[-i]))
                root = i
    }
    return(root)
}

## need to modify this to not use the max, but how to find the nodes objects?
isRootGood <- function(tree){
    root = getRoot(tree)
    ind = which(tree$edge[,1] == root)
    other = tree$edge[ind,2]
    taxa = tree$tip.label
    cats = c(1,2,3,4,5,12)
    dogs = c(6,7,8,9,10,11)

    res = FALSE
    for(i in 1:3){
        des = getDescendants(tree,other[i])
        desc = taxa[des]
        truedesc = desc[!is.na(desc)]
        sorttruedesc = sort(as.numeric(truedesc))
        if(length(sorttruedesc) == 6){
            res = all(sorttruedesc == cats) | all(sorttruedesc == dogs)
            break
        }
    }

    return(res)
}

isRootOnCats <- function(tree){
    root = getRoot(tree)
    ind = which(tree$edge[,1] == root)
    other = tree$edge[ind,2]
    taxa = tree$tip.label
    cats = c(1,2,3,4,5,12)
    dogs = c(6,7,8,9,10,11)

    res = FALSE
    for(i in 1:3){
        des = getDescendants(tree,other[i])
        desc = taxa[des]
        truedesc = desc[!is.na(desc)]
        sorttruedesc = sort(as.numeric(truedesc))
        if(length(sorttruedesc) == 6){
            res = all(sorttruedesc == cats)
            break
        }
    }

    return(res)
}

## Function to find all pairs of adjacent edges
getAdjacentEdges = function(tree)
{
  require(ape)
  edges = tree$edge
  internal = unique(edges[,1])
  n = length(internal)
  out = matrix(0,3*n,2)
  curr = 1
  for( i in internal )
  {
    keep = which(edges[,1]==i | edges[,2]==i)
    if ( length(keep) !=3 )
    {
      print(edges[keep,])
      stop("Something screwed up")
    }
    out[curr,] = c(keep[1],keep[2])
    curr = curr+1
    out[curr,] = c(keep[1],keep[3])
    curr = curr+1
    out[curr,] = c(keep[2],keep[3])
    curr = curr+1
  }
  return( out )
}


plotConvergence = function(stem){
    require(ggplot2)
    mcmc_out = read.table(paste0(stem,".mcmc.par"),header=FALSE)
    names(mcmc_out) = c("logl",paste0("pi",1:4),paste0("s",1:6))
    n = length(mcmc_out$logl)

    print(ggplot(mcmc_out,aes(x=1:n,y=logl)) + geom_point())

    print(ggplot(mcmc_out,aes(x=1:n)) +
        geom_point(aes(y=pi1),col="blue") +
        geom_point(aes(y=pi2),col="red") +
        geom_point(aes(y=pi3),col="green") +
        geom_point(aes(y=pi4),col="yellow"))

    print(ggplot(mcmc_out,aes(x=1:n)) +
        geom_point(aes(y=s1),col="blue") +
        geom_point(aes(y=s2),col="red") +
        geom_point(aes(y=s3),col="green") +
        geom_point(aes(y=s4),col="yellow") +
        geom_point(aes(y=s5),col="orange") +
        geom_point(aes(y=s6),col="purple"))
}

plotLengthFrequency = function(stem){
    tab1 = read.table(paste0(stem,".bootstrapCladeBL"),header=TRUE)
    tab2 = read.table(paste0(stem,".vstat"),header=TRUE)
    myfun = function(x){
        a = grepl("-",x)
        b = grepl(",",x)
        return(a|b)
    }
    ind1 = sapply(tab1$clade,FUN=myfun)
    ind2 = sapply(tab2$clade,FUN=myfun)
    tab1 = tab1[ind1,]
    tab2 = tab2[ind2,]
    print(ggplot(tab1,aes(x=weighted.mean.BL,y=weight)) + geom_point() +
        geom_point(data=tab2,color="red") + ggtitle("black=bootstrap, red=IS"))
}


plotBistroPDF = function(stem){
    pdf(paste0(stem,".pdf"))
    plotBistro(stem)
    plotConvergence(stem)
    dev.off()
}
