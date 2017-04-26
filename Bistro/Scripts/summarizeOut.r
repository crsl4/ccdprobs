## r script to summarize .out files
## based on readBistro.R by Bret Larget
## will read all root---*.out files and
## create a pstat file to compare with mrbayes
## need to be in the folder with files
## Claudia September 2016

##bistroroot = "goo"

library(SDMTools)
library(reldist)

summarizeOut = function(bistroroot){
    files = list.files(pattern=paste0("^",bistroroot,"---.*\\.out"))
    bistro = data.frame()
    for ( f in files ) {
        temp = read.table(f,header=TRUE)
        bistro = rbind(bistro,temp)
        rm(temp)
    }
    bistro$w = with(bistro, exp(logWt - max(logWt)))
    bistro$w = with(bistro, w/sum(w))
    str(bistro)

    ## transform s to mrbayes rates
    bistro$r1 <- with(bistro, s1/(pi1*pi2))
    bistro$r2 <- with(bistro, s2/(pi1*pi3))
    bistro$r3 <- with(bistro, s3/(pi1*pi4))
    bistro$r4 <- with(bistro, s4/(pi2*pi3))
    bistro$r5 <- with(bistro, s5/(pi2*pi4))
    bistro$r6 <- with(bistro, s6/(pi3*pi4))
    sum <- bistro$r1+bistro$r2+bistro$r3+bistro$r4+bistro$r5+bistro$r6
    bistro$r1 <- bistro$r1/sum
    bistro$r2 <- bistro$r2/sum
    bistro$r3 <- bistro$r3/sum
    bistro$r4 <- bistro$r4/sum
    bistro$r5 <- bistro$r5/sum
    bistro$r6 <- bistro$r6/sum

    means=rep(0,10)
    means[1] <- wt.mean(bistro$r1,bistro$w)
    means[2] <- wt.mean(bistro$r2,bistro$w)
    means[3] <- wt.mean(bistro$r3,bistro$w)
    means[4] <- wt.mean(bistro$r4,bistro$w)
    means[5] <- wt.mean(bistro$r5,bistro$w)
    means[6] <- wt.mean(bistro$r6,bistro$w)
    means[7] <- wt.mean(bistro$pi1,bistro$w)
    means[8] <- wt.mean(bistro$pi2,bistro$w)
    means[9] <- wt.mean(bistro$pi3,bistro$w)
    means[10] <- wt.mean(bistro$pi4,bistro$w)

    vars=rep(0,10)
    vars[1] <- wt.var(bistro$r1,bistro$w)
    vars[2] <- wt.var(bistro$r2,bistro$w)
    vars[3] <- wt.var(bistro$r3,bistro$w)
    vars[4] <- wt.var(bistro$r4,bistro$w)
    vars[5] <- wt.var(bistro$r5,bistro$w)
    vars[6] <- wt.var(bistro$r6,bistro$w)
    vars[7] <- wt.var(bistro$pi1,bistro$w)
    vars[8] <- wt.var(bistro$pi2,bistro$w)
    vars[9] <- wt.var(bistro$pi3,bistro$w)
    vars[10] <- wt.var(bistro$pi4,bistro$w)

    stds <- with(bistro, sqrt(vars))

    lower <- rep(0,10)
    lower[1] <- wtd.quantile(bistro$r1,q=0.025,weight=bistro$w)
    lower[2] <- wtd.quantile(bistro$r2,q=0.025,weight=bistro$w)
    lower[3] <- wtd.quantile(bistro$r3,q=0.025,weight=bistro$w)
    lower[4] <- wtd.quantile(bistro$r4,q=0.025,weight=bistro$w)
    lower[5] <- wtd.quantile(bistro$r5,q=0.025,weight=bistro$w)
    lower[6] <- wtd.quantile(bistro$r6,q=0.025,weight=bistro$w)
    lower[7] <- wtd.quantile(bistro$pi1,q=0.025,weight=bistro$w)
    lower[8] <- wtd.quantile(bistro$pi2,q=0.025,weight=bistro$w)
    lower[9] <- wtd.quantile(bistro$pi3,q=0.025,weight=bistro$w)
    lower[10] <- wtd.quantile(bistro$pi4,q=0.025,weight=bistro$w)

    upper <- rep(0,10)
    upper[1] <- wtd.quantile(bistro$r1,q=0.975,weight=bistro$w)
    upper[2] <- wtd.quantile(bistro$r2,q=0.975,weight=bistro$w)
    upper[3] <- wtd.quantile(bistro$r3,q=0.975,weight=bistro$w)
    upper[4] <- wtd.quantile(bistro$r4,q=0.975,weight=bistro$w)
    upper[5] <- wtd.quantile(bistro$r5,q=0.975,weight=bistro$w)
    upper[6] <- wtd.quantile(bistro$r6,q=0.975,weight=bistro$w)
    upper[7] <- wtd.quantile(bistro$pi1,q=0.975,weight=bistro$w)
    upper[8] <- wtd.quantile(bistro$pi2,q=0.975,weight=bistro$w)
    upper[9] <- wtd.quantile(bistro$pi3,q=0.975,weight=bistro$w)
    upper[10] <- wtd.quantile(bistro$pi4,q=0.975,weight=bistro$w)

    med <- rep(0,10)
    med[1] <- wtd.quantile(bistro$r1,q=0.5,weight=bistro$w)
    med[2] <- wtd.quantile(bistro$r2,q=0.5,weight=bistro$w)
    med[3] <- wtd.quantile(bistro$r3,q=0.5,weight=bistro$w)
    med[4] <- wtd.quantile(bistro$r4,q=0.5,weight=bistro$w)
    med[5] <- wtd.quantile(bistro$r5,q=0.5,weight=bistro$w)
    med[6] <- wtd.quantile(bistro$r6,q=0.5,weight=bistro$w)
    med[7] <- wtd.quantile(bistro$pi1,q=0.5,weight=bistro$w)
    med[8] <- wtd.quantile(bistro$pi2,q=0.5,weight=bistro$w)
    med[9] <- wtd.quantile(bistro$pi3,q=0.5,weight=bistro$w)
    med[10] <- wtd.quantile(bistro$pi4,q=0.5,weight=bistro$w)

    samplemean <- rep(0,10)
    samplemean[1] <- mean(bistro$r1)
    samplemean[2] <- mean(bistro$r2)
    samplemean[3] <- mean(bistro$r3)
    samplemean[4] <- mean(bistro$r4)
    samplemean[5] <- mean(bistro$r5)
    samplemean[6] <- mean(bistro$r6)
    samplemean[7] <- mean(bistro$pi1)
    samplemean[8] <- mean(bistro$pi2)
    samplemean[9] <- mean(bistro$pi3)
    samplemean[10] <- mean(bistro$pi4)

    samplestd <- rep(0,10)
    samplestd[1] <- sd(bistro$r1)
    samplestd[2] <- sd(bistro$r2)
    samplestd[3] <- sd(bistro$r3)
    samplestd[4] <- sd(bistro$r4)
    samplestd[5] <- sd(bistro$r5)
    samplestd[6] <- sd(bistro$r6)
    samplestd[7] <- sd(bistro$pi1)
    samplestd[8] <- sd(bistro$pi2)
    samplestd[9] <- sd(bistro$pi3)
    samplestd[10] <- sd(bistro$pi4)

    parameter = c("r(A<->C)","r(A<->G)","r(A<->T)","r(C<->G)","r(C<->T)","r(G<->T)","pi(A)","pi(C)","pi(G)","pi(T)")
    dataOut = data.frame(parameter,means,vars,stds,lower,upper,med, samplemean, samplestd)
    write.table(dataOut, file=paste0(bistroroot,"-r.pstat"), row.names=FALSE, sep="\t")
    print(dataOut)
    return(dataOut)
}
