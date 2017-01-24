## stem is the rootname for bistro
## mb=TRUE will compare to the ccdw runs of MrBayes
## o.w. it will compare to mcmc1.tre, mcmc1.par from bmcmc

## artiodactyl "(1,2,(3,(4,(5,6))));"
## cats-dogs "(1,(2,(((((3,4),5),6),7),8)),(9,((10,11),12)));"
## bistro writes it like this: "(1,(2,(((((7,8),9),10),12),11)),(3,((4,5),6)));"
## primates: "(1,2,((((((3,4),5),6),7),(((8,9),10),11)),12));"

summaryBistro = function(stem, mb=FALSE, besttree="(1,2,(3,4));", bmcmc="mcmc1"){
    if(mb)
        {
            ## mrbayes
            top1 = read.table("ccdw.run1.top")
            top2 = read.table("ccdw.run2.top")
            tre1 = read.table("ccdw.run1.tre")
            tre2 = read.table("ccdw.run2.tre")

            ## remove burnin
            top = c(as.character(top1[-(1:2001),]),as.character(top2[-(1:2001),]))
            tre = c(as.character(tre1[-(1:2001),]),as.character(tre2[-(1:2001),]))
            ## best tree
            keep = which(top==besttree)

            ## sample 1000 for better plots
            keep = sort(sample(keep,1000))
            tre.mb = tre[keep]

            ## remove big stuff
            rm(top1,top2,tre1,tre2,top,tre,keep)

            other = "MrBayes"
        }
    else
        {
            ## New using bmcmc
            tre = read.table(paste0(bmcmc,".tre"))
            n = nrow(tre)
            burn = n/11
            ## remove burnin, first 1000 saved values
            tre = as.character(tre[-(1:burn),])
            if( n > 1000 )
                tre.mb = tre[sample(1:length(tre),1000,replace=FALSE)]
            else
                tre.mb = tre
            other = "MCMC"
        }

    require(ape)
    N.mb = length(tre.mb) ### 1000 now!!!

    tree = read.tree(text=tre.mb[1])
    ## get number of edges
    m.mb = matrix(0,N.mb,length(tree$edge.length))
    colnames(m.mb) = paste("b",tree$edge[,2],tree$edge[,1],sep=".")

    for ( i in 1:N.mb )
        {
            tree = read.tree(text=tre.mb[i])
            m.mb[i,] = tree$edge.length
        }

    ## now do bistro
    source("../../Scripts/readBistro.r")
    bistro = readBistro(stem)
    data = readDataSort(stem)
    pdf(paste0(stem,other,"-cloud.pdf"))
    plotBistro(bistro)
    dev.off()

    ## we need to treat the cats-dogs differently because canonical sort in bistro does not
    ## write it in the same way as in mcmc
    if(besttree == "(1,(2,(((((3,4),5),6),7),8)),(9,((10,11),12)));"){
        keep = which(bistro$tree=="(1,(2,(((((7,8),9),10),12),11)),(3,((4,5),6)));")
        tre.bistro = data$V1[keep]
    }else{
        keep = which(bistro$tree==besttree)
        tre.bistro = data$V1[keep]
        tree=read.tree(text=as.character(tre.bistro[1]))
    }


    N.bistro = length(tre.bistro)
    ## get number of edges
    m.bistro = matrix(0,N.bistro,length(tree$edge[,1]))
    colnames(m.bistro) = paste("b",tree$edge[,2],tree$edge[,1],sep=".")

    ## check same names!!!
    all( colnames(m.mb) == colnames(m.bistro) )

    for ( i in 1:N.bistro )
        {
            tree = read.tree(text=as.character(tre.bistro[i]))
            m.bistro[i,] = tree$edge.length
        }

    require(corrplot)
    pdf(paste0(stem,other,"-corr.pdf"))
    corrplot(cor(m.bistro))
    corrplot(cor(m.mb))
    dev.off()

    ## combine data for combined plots
    m = rbind(m.bistro,m.mb)
    df = data.frame(m)
    df$set = factor( c(rep("Bistro",N.bistro),rep(other,N.mb)) )

    ## make the plots
    adjEdges = getAdjacentEdges(tree)

    require(viridis)
    require(ggplot2)
    require(dplyr)

    pdf(paste0(stem,other,"-scatter.pdf"))
    vpal = viridis(2,end=0.8)
    for( i in 1:nrow(adjEdges))
        {
            i1 = adjEdges[i,1]
            i2 = adjEdges[i,2]
            median.bistro.1 = median( drop(as.matrix(filter(df,set=="Bistro") %>% select(i1))) )
            median.bistro.2 = median( drop(as.matrix(filter(df,set=="Bistro") %>% select(i2))) )
            median.mb.1 = median( drop(as.matrix(filter(df,set==other) %>% select(i1))) )
            median.mb.2 = median( drop(as.matrix(filter(df,set==other) %>% select(i2))) )

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

    pdf(paste0(stem,other,"-density.pdf"))
    vpal = viridis(2,end=0.8)
    for(i in 1:ncol(m))
        {
            median.bistro = mean( drop(as.matrix(filter(df,set=="Bistro") %>% select(i))) )
            median.mb = mean( drop(as.matrix(filter(df,set==other) %>% select(i))) )

            plot(ggplot(df,aes(x=df[,i],col=set))+geom_density() +
                 scale_color_manual(values=vpal) +
                 geom_vline(xintercept=median.bistro,color=vpal[1]) +
                 geom_vline(xintercept=median.mb,color=vpal[2]) +
                 ggtitle(paste(names(df)[i])) +
                 theme_bw())

        }
    dev.off()

########################################################################################
    if(mb)
        {
            foo1 = read.table("ccdw.run1.p", header=TRUE)
            foo2 = read.table("ccdw.run2.p", header=TRUE)
            foo = rbind(foo1[-(1:2001),],foo2[-(1:2001),])
            foo$sAC = with(foo, r.A...C. * pi.A. * pi.C.)
            foo$sAG = with(foo, r.A...G. * pi.A. * pi.G.)
            foo$sAT = with(foo, r.A...T. * pi.A. * pi.T.)
            foo$sCG = with(foo, r.C...G. * pi.C. * pi.G.)
            foo$sCT = with(foo, r.C...T. * pi.C. * pi.T.)
            foo$sGT = with(foo, r.G...T. * pi.G. * pi.T.)
            s = with(foo,sAC+sAG+sAT+sCG+sCT+sGT)
            foo$sAC = foo$sAC/s
            foo$sAG = foo$sAG/s
            foo$sAT = foo$sAT/s
            foo$sCG = foo$sCG/s
            foo$sCT = foo$sCT/s
            foo$sGT = foo$sGT/s
            N.mb2 = length(foo$Gen)
            foo2 = subset(foo,select=c("pi.A.","pi.C.","pi.G.","pi.T.","sAC","sAG","sAT","sCG","sCT","sGT"))
            names(foo2) = c("pi1","pi2","pi3","pi4","s1","s2","s3","s4","s5","s6")
        }
    else
        {
            foo2 = read.table(paste0(bmcmc,".par"))
            foo2 = foo2[-(1:burn),-1]
            if ( n > 1000 )
                foo2 = foo2[sample(1:nrow(foo2),1000,replace=FALSE),]
            names(foo2) = c("pi1","pi2","pi3","pi4","s1","s2","s3","s4","s5","s6")
            N.mb2 = nrow(foo2)
        }
    bistro2 = subset(bistro,select=c("pi1","pi2","pi3","pi4","s1","s2","s3","s4","s5","s6"))
    N.bistro2 = nrow(bistro2)

    df2 = rbind(bistro2,foo2)
    df2$set = factor( c(rep("Bistro",N.bistro2),rep(other,N.mb2)) )


    pdf(paste0(stem,other,"-rates-density.pdf"))
    vpal = viridis(2,end=0.8)
    for(i in 1:(ncol(df2)-1))
        {
            median.bistro = mean( drop(as.matrix(filter(df2,set=="Bistro") %>% select(i))) )
            median.mb = mean( drop(as.matrix(filter(df2,set==other) %>% select(i))) )
            print(median.mb)
            plot(ggplot(df2,aes(x=df2[,i],col=set))+geom_density() +
                 scale_color_manual(values=vpal) +
                 geom_vline(xintercept=median.bistro,color=vpal[1]) +
                 geom_vline(xintercept=median.mb,color=vpal[2]) +
                 ggtitle(paste(names(df2)[i])) +
                 theme_bw())
        }
    dev.off()
}

## same as summaryBistro, but plots the true values as well
## using true bl and q=pi,s values of artiodactyl now
summaryBistroSim = function(stem, mb=FALSE, besttree="(1,2,(3,(4,(5,6))));", bmcmc="mcmc1", bl=c(1.460221e-01,1.214757e-01,1.072255e-01,8.941711e-02,8.610459e-02,7.642583e-02,1.846791e-02,5.146787e-02,2.783294e-02), q=c(0.2870113,0.2715106,0.1457293,0.2957489,0.12041579,0.07811855,0.42186604,0.09645864,0.24588493,0.03725606)){
    if(mb)
        {
            ## mrbayes
            top1 = read.table("ccdw.run1.top")
            top2 = read.table("ccdw.run2.top")
            tre1 = read.table("ccdw.run1.tre")
            tre2 = read.table("ccdw.run2.tre")

            ## remove burnin
            top = c(as.character(top1[-(1:2001),]),as.character(top2[-(1:2001),]))
            tre = c(as.character(tre1[-(1:2001),]),as.character(tre2[-(1:2001),]))
            ## best tree
            keep = which(top==besttree)

            ## sample 1000 for better plots
            keep = sort(sample(keep,1000))
            tre.mb = tre[keep]

            ## remove big stuff
            rm(top1,top2,tre1,tre2,top,tre,keep)

            other = "MrBayes"
        }
    else
        {
            ## New using bmcmc
            tre = read.table(paste0(bmcmc,".tre"))
            n = nrow(tre)
            burn = n/11
            ## remove burnin, first 1000 saved values
            tre = as.character(tre[-(1:burn),])
            if( n > 1000 )
                tre.mb = tre[sample(1:length(tre),1000,replace=FALSE)]
            else
                tre.mb = tre
            other = "MCMC"
        }

    require(ape)
    N.mb = length(tre.mb) ### 1000 now!!!

    tree = read.tree(text=tre.mb[1])
    ## get number of edges
    m.mb = matrix(0,N.mb,length(tree$edge.length))
    colnames(m.mb) = paste("b",tree$edge[,2],tree$edge[,1],sep=".")

    for ( i in 1:N.mb )
        {
            tree = read.tree(text=tre.mb[i])
            m.mb[i,] = tree$edge.length
        }

    ## now do bistro
    source("../../Scripts/readBistro.r")
    bistro = readBistro(stem)
    data = readDataSort(stem)
    pdf(paste0(stem,other,"-cloud.pdf"))
    plotBistro(bistro)
    dev.off()

    ## we need to treat the cats-dogs differently because canonical sort in bistro does not
    ## write it in the same way as in mcmc
    if(besttree == "(1,(2,(((((3,4),5),6),7),8)),(9,((10,11),12)));"){
        keep = which(bistro$tree=="(1,(2,(((((7,8),9),10),12),11)),(3,((4,5),6)));")
        tre.bistro = data$V1[keep]
    }else{
        keep = which(bistro$tree==besttree)
        tre.bistro = data$V1[keep]
        tree=read.tree(text=as.character(tre.bistro[1]))
    }


    N.bistro = length(tre.bistro)
    ## get number of edges
    m.bistro = matrix(0,N.bistro,length(tree$edge[,1]))
    colnames(m.bistro) = paste("b",tree$edge[,2],tree$edge[,1],sep=".")

    ## check same names!!!
    all( colnames(m.mb) == colnames(m.bistro) )

    for ( i in 1:N.bistro )
        {
            tree = read.tree(text=as.character(tre.bistro[i]))
            m.bistro[i,] = tree$edge.length
        }

    require(corrplot)
    pdf(paste0(stem,other,"-corr.pdf"))
    corrplot(cor(m.bistro))
    corrplot(cor(m.mb))
    dev.off()

    ## combine data for combined plots
    m = rbind(m.bistro,m.mb)
    df = data.frame(m)
    df$set = factor( c(rep("Bistro",N.bistro),rep(other,N.mb)) )

    ## make the plots
    adjEdges = getAdjacentEdges(tree)

    require(viridis)
    require(ggplot2)
    require(dplyr)

    pdf(paste0(stem,other,"-scatter.pdf"))
    vpal = viridis(2,end=0.8)
    for( i in 1:nrow(adjEdges))
        {
            i1 = adjEdges[i,1]
            i2 = adjEdges[i,2]
            median.bistro.1 = median( drop(as.matrix(filter(df,set=="Bistro") %>% select(i1))) )
            median.bistro.2 = median( drop(as.matrix(filter(df,set=="Bistro") %>% select(i2))) )
            median.mb.1 = median( drop(as.matrix(filter(df,set==other) %>% select(i1))) )
            median.mb.2 = median( drop(as.matrix(filter(df,set==other) %>% select(i2))) )

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

    pdf(paste0(stem,other,"-density.pdf"))
    vpal = viridis(2,end=0.8)
    if(ncol(m) != length(bl))
        stop("input bl vector with true values does not match the number of bl")
    for(i in 1:ncol(m))
        {
            median.bistro = mean( drop(as.matrix(filter(df,set=="Bistro") %>% select(i))) )
            median.mb = mean( drop(as.matrix(filter(df,set==other) %>% select(i))) )
            trueBL = bl[i]

            plot(ggplot(df,aes(x=df[,i],col=set))+geom_density() +
                 scale_color_manual(values=vpal) +
                 geom_vline(xintercept=median.bistro,color=vpal[1]) +
                 geom_vline(xintercept=median.mb,color=vpal[2]) +
                 geom_vline(xintercept=trueBL,color="black") +
                 ggtitle(paste(names(df)[i])) +
                 theme_bw())

        }
    dev.off()

########################################################################################
    if(mb)
        {
            foo1 = read.table("ccdw.run1.p", header=TRUE)
            foo2 = read.table("ccdw.run2.p", header=TRUE)
            foo = rbind(foo1[-(1:2001),],foo2[-(1:2001),])
            foo$sAC = with(foo, r.A...C. * pi.A. * pi.C.)
            foo$sAG = with(foo, r.A...G. * pi.A. * pi.G.)
            foo$sAT = with(foo, r.A...T. * pi.A. * pi.T.)
            foo$sCG = with(foo, r.C...G. * pi.C. * pi.G.)
            foo$sCT = with(foo, r.C...T. * pi.C. * pi.T.)
            foo$sGT = with(foo, r.G...T. * pi.G. * pi.T.)
            s = with(foo,sAC+sAG+sAT+sCG+sCT+sGT)
            foo$sAC = foo$sAC/s
            foo$sAG = foo$sAG/s
            foo$sAT = foo$sAT/s
            foo$sCG = foo$sCG/s
            foo$sCT = foo$sCT/s
            foo$sGT = foo$sGT/s
            N.mb2 = length(foo$Gen)
            foo2 = subset(foo,select=c("pi.A.","pi.C.","pi.G.","pi.T.","sAC","sAG","sAT","sCG","sCT","sGT"))
            names(foo2) = c("pi1","pi2","pi3","pi4","s1","s2","s3","s4","s5","s6")
        }
    else
        {
            foo2 = read.table(paste0(bmcmc,".par"))
            foo2 = foo2[-(1:burn),-1]
            if ( n > 1000 )
                foo2 = foo2[sample(1:nrow(foo2),1000,replace=FALSE),]
            names(foo2) = c("pi1","pi2","pi3","pi4","s1","s2","s3","s4","s5","s6")
            N.mb2 = nrow(foo2)
        }
    bistro2 = subset(bistro,select=c("pi1","pi2","pi3","pi4","s1","s2","s3","s4","s5","s6"))
    N.bistro2 = nrow(bistro2)

    df2 = rbind(bistro2,foo2)
    df2$set = factor( c(rep("Bistro",N.bistro2),rep(other,N.mb2)) )


    pdf(paste0(stem,other,"-rates-density.pdf"))
    vpal = viridis(2,end=0.8)
    if(length(q) != 10)
        stop("input vector q with true values does not have 10 elements: 4 pi, 6 s")
    for(i in 1:(ncol(df2)-1))
        {
            median.bistro = mean( drop(as.matrix(filter(df2,set=="Bistro") %>% select(i))) )
            median.mb = mean( drop(as.matrix(filter(df2,set==other) %>% select(i))) )
            trueQ = q[i]
            print(median.mb)
            plot(ggplot(df2,aes(x=df2[,i],col=set))+geom_density() +
                 scale_color_manual(values=vpal) +
                 geom_vline(xintercept=median.bistro,color=vpal[1]) +
                 geom_vline(xintercept=median.mb,color=vpal[2]) +
                 geom_vline(xintercept=trueQ,color="black") +
                 ggtitle(paste(names(df2)[i])) +
                 theme_bw())
        }
    dev.off()
}
