# r script to simulate many scenarios and determine which conditions
# yield smaller bias between true likelihood and estimated density
# Claudia February 2016

source('branch-length.r')

## simulate densities
simulateDensities = function(nsites,branch.length,nsim=10000,eta.jc=0.5,eta.tn=0.8, delta=0.001) {
    s = seq(delta,2*branch.length,delta)
    Q = randomQ(4,rescale=TRUE)
    x = simulateSequenceSummary(nsites,Q,branch.length)
    while((0 %in% rowSums(x)) || (0 %in% colSums(x))){ #discard x (counts) with rows of zero (and also columns because error in density(tn))
        Q = randomQ(4,rescale=TRUE)
        x = simulateSequenceSummary(nsites,Q,branch.length)
    }
    print(round(Q$Q,4))
    print(round(Q$p,4))
    print(min(diag(Q$Q)))
    print(x)
    ## compute the likelihood using counts x and generator Q for each s
    log.like = gtr.log.like(x,s,Q)
    ## turn into a density estimate proportional to likelihood, so truncated uniform prior in t
    log.like = log.like - max(log.like) # going to rescale anyway, may exp() accurate
    y = exp(log.like)
    ## make y sum to one
    y = y / sum(y)
    ## density at each point is y divided by delta
    y = y / delta
    ## Make a data frame
    df.true = data.frame(s,y)
    ## GTR Q optimized over data
    Qlist.gtr = optim.gtr(x,Q$r)
    #print(Qlist.gtr$branch.length)
    ## compute this density also
    log.like2 = gtr.log.like(x,s,Qlist.gtr$Q)
    log.like2 = log.like2 - max(log.like2)
    y2 = exp(log.like2)
    y2 = y2/sum(y2)
    y2 = y2 / delta
    df.gtr = data.frame(s,y2)
    ## now get a sample of JCish points to compare
    jc = simulateBranchLength.jc(nsim,x,eta.jc)
    ## density estimate
    d.jc = density(jc)
    ## another data frame
    df.jc = data.frame(x=d.jc$x,y=d.jc$y)
    ## TN estimate
    tn = simulateBranchLength.tn(nsim,x,eta.tn)
    d.tn = density(tn)
    df.tn = data.frame(x=d.tn$x,y=d.tn$y)
    ## plot from comparePlots
    ## p1 = ggplot(df.true, aes(x=s,y=y)) +
    ##     geom_line(color="blue") +
    ##         geom_line(aes(x=x,y=y),data=df.jc,color="red",linetype="dashed") +
    ##             geom_line(aes(x=x,y=y),data=df.tn,color="darkgreen",linetype="dashed") +
    ##                 geom_line(aes(x=s,y=y2),data=df.gtr,color="gold") +
    ##                     xlab('branch length') +
    ##                         ylab('densities') +
    ##                             ggtitle('Blue = TrueQ, Red = JC, Green = TN, Gold = GTR')
    ## plot(p1)
    return( list(Q=Q$Q,p=Q$p,estT=Qlist.gtr$branch.length, df.true=df.true, df.gtr=df.gtr,df.jc=df.jc, df.tn=df.tn) )
}


## Compute the bias between the JC and TN densities, and the true likelihood
## for a given value of branch.length
computeBias = function(nsites,branch.length,nsim=10000,eta.jc=0.5,eta.tn=0.8, delta=0.001) {
    D = simulateDensities(nsites,branch.length,nsim,eta.jc,eta.tn, delta)
    ## how far df.true$y from df.gtr$y
    meanSqGTR = mean((D$df.true$y-D$df.gtr$y)^2)
    #print(meanSqGTR)
    mode.true = D$df.true[which.max(D$df.true$y),1]
    mode.gtr = D$df.gtr[which.max(D$df.gtr$y),1]
    print("modes true")
    print(mode.true)
    print(mode.gtr)
    ## JC bias
    mode.jc = D$df.jc[which.max(D$df.jc$y),1]
    print("mode JC")
    print(mode.jc)
    bias.jc.true = mode.true - mode.jc
    bias.jc.gtr = mode.gtr - mode.jc
    print(bias.jc.true)
    print(bias.jc.gtr)
    ## TN bias
    mode.tn = D$df.tn[which.max(D$df.tn$y),1]
    print("mode TN")
    print(mode.tn)
    bias.tn.true = mode.true - mode.tn
    bias.tn.gtr = mode.gtr - mode.tn
    print(bias.tn.true)
    print(bias.tn.gtr)
    return( list(Q=D$Q,p=D$p,estT=D$estT, ms=meanSqGTR, bias.jc.true=bias.jc.true, bias.jc.gtr=bias.jc.gtr, bias.tn.true=bias.tn.true, bias.tn.gtr=bias.tn.gtr))
}

computeBiasSimulations = function(nsites,branch.length,nsim=10000,eta.jc=0.5,eta.tn=0.8, delta=0.001, nrep=1000) {
    Q.vec = list()
    p.vec = list()
    bias.jc.true.vec = rep(0,nrep)
    bias.jc.gtr.vec = rep(0,nrep)
    bias.tn.true.vec = rep(0,nrep)
    bias.tn.gtr.vec = rep(0,nrep)
    estT.vec = rep(0,nrep)
    for(i in 1:nrep){
        print(i)
        B = computeBias(nsites, branch.length, nsim, eta.jc, eta.tn, delta)
        Q.vec[[i]] = B$Q
        p.vec[[i]] = B$p
        bias.jc.true.vec[i] = B$bias.jc.true
        bias.jc.gtr.vec[i] = B$bias.jc.gtr
        bias.tn.true.vec[i] = B$bias.tn.true
        bias.tn.gtr.vec[i] = B$bias.tn.gtr
        estT.vec[i] = B$estT
    }
    return( list(Q.vec=Q.vec,p.vec=p.vec,estT.vec=estT.vec, bias.jc.true=bias.jc.true.vec, bias.jc.gtr=bias.jc.gtr.vec, bias.tn.true=bias.tn.true.vec, bias.tn.gtr=bias.tn.gtr.vec))
}


nrep = 100
branch.length=0.15
nsites=500
B=computeBiasSimulations(nsites,branch.length, nrep=nrep)

## need function to summarize B for a given branch.length: mean(min(Qdiag)), mean(min(p)), mean(max(p)), mean(biases), mean(estT) (and sd)
# s=vector of branch lengths
#s=seq(0.01,0.2,by=0.01)
#minDiag = function(x){
#    min(diag(x))}

# for branch length s
computeBiasSimulations.df = function(s,nsites=500,nsim=10000,eta.jc=0.5,eta.tn=0.8, delta=0.001, nrep=10) {
    df=data.frame(bl=rep(0,length(s)*nrep), minQdiag=rep(0,length(s)*nrep), biasJC.true=rep(0,length(s)*nrep), biasJC.gtr=rep(0,length(s)*nrep),  biasTN.true=rep(0,length(s)*nrep), biasTN.gtr=rep(0,length(s)*nrep), estT=rep(0,length(s)*nrep))
    B=computeBiasSimulations(nsites,s,nsim,eta.jc,eta.tn, delta, nrep)
    for(j in 1:nrep){
        df[j,1] = s
        df[j,2] = min(diag(B$Q.vec[[j]]))
        df[j,3] = B$bias.jc.true[j]
        df[j,4] = B$bias.jc.gtr[j]
        df[j,5] = B$bias.tn.true[j]
        df[j,6] = B$bias.tn.gtr[j]
        df[j,7] = B$estT.vec[j]
    }
    return ( list(df=df,B=B))
}


s=0.15
nrep=1000
nsites=500
d500=computeBiasSimulations.df(s,nrep=nrep, nsites=nsites)
save(d500,file="d500.Rda")

nsites=1000
d1000=computeBiasSimulations.df(s,nrep=nrep, nsites=nsites)
save(d1000,file="d1000.Rda")

load("df.Rda")
library('ggplot2')
d=d1000

p <- ggplot(d$df, aes(x=minQdiag, y=biasJC.true))+ggtitle("Bias True: blue JC, green TN")+
    geom_point(alpha=0.3, color='blue')+
    geom_point(aes(x=minQdiag,y=biasTN.true), data=d$df,color='green', alpha=0.2)+
    geom_hline(yintercept=mean(d$df$biasJC.true), color="blue")+
    geom_hline(yintercept=mean(d$df$biasTN.true), color="green")
pdf("../manuscript/biasTrue1000_2.pdf")
plot(p)
dev.off()

p2 <- ggplot(d$df, aes(x=minQdiag, y=biasJC.gtr))+ggtitle("Bias GTR: blue JC, green TN")+
    geom_point(alpha=0.3, color='blue')+
    geom_point(aes(x=minQdiag,y=biasTN.gtr), data=d$df,color='green', alpha=0.2)+
    geom_hline(yintercept=mean(d$df$biasJC.gtr), color="blue")+
    geom_hline(yintercept=mean(d$df$biasTN.gtr), color="green")
pdf("../manuscript/biasGTR1000_2.pdf")
plot(p2)
dev.off()

p3 <- ggplot(d$df, aes(x=biasTN.true, y=biasJC.true))+ggtitle("true blue, GTR green")+
    geom_point(alpha=0.3, color='blue')+
    geom_point(aes(x=biasTN.gtr, y=biasJC.gtr), data=d$df,color='green', alpha=0.2)+
    geom_abline(slope=1)
pdf("../manuscript/JCvsTNtrue1000_2.pdf")
plot(p3)
dev.off()

# ========================================================================================================
# coverage
source('branch-length.r')
library(ggplot2)
nsites=500
branch.length=0.15
nsim=10000
eta.jc=0.2
eta.tn=0.2
delta=0.001

computeCoverage = function(nsites,branch.length,nsim=10000,eta.jc=0.8,eta.tn=0.8, delta=0.001) {
    D = simulateDensities(nsites,branch.length,nsim,eta.jc,eta.tn, delta)
    #cutoff = c(0.01,5,10)
    cutoff = seq(0.01,5,by=0.5)
    covJC.true = 0
    covTN.true = 0
    covJC.gtr = 0
    covTN.gtr = 0
    for(c in cutoff){
        print(paste("cutoff",c,"-------"))
        subsetTrue = D$df.true[D$df.true$y>c,]
        subsetGTR = D$df.gtr[D$df.gtr$y>c,]
        subsetJC = D$df.jc[D$df.jc$y>c,]
        subsetTN = D$df.tn[D$df.tn$y>c,]
        lowerTrue = subsetTrue[1,1]
        lowerGTR = subsetGTR[1,1]
        lowerJC = subsetJC[1,1]
        lowerTN = subsetTN[1,1]
        upperTrue = subsetTrue[length(subsetTrue$s),1]
        upperGTR = subsetGTR[length(subsetGTR$s),1]
        upperJC = subsetJC[length(subsetJC$x),1]
        upperTN = subsetTN[length(subsetTN$x),1]
        print(paste("lowerTrue",lowerTrue,"upperTrue", upperTrue))
        print(paste("lowerJC",lowerJC,"upperJC", upperJC))
        print(paste("lowerTN",lowerTN,"upperTN", upperTN))
        covJC.true = covJC.true + (lowerTrue>lowerJC) +(upperTrue<lowerJC)
        covJC.gtr = covJC.gtr + (lowerGTR>lowerJC) +(upperGTR<lowerJC)
        covTN.true = covJC.true + (lowerTrue>lowerTN) +(upperTrue<lowerTN)
        covTN.gtr = covJC.gtr + (lowerGTR>lowerTN) +(upperGTR<lowerTN)
    }
    covJC.true= covJC.true/(2*length(cutoff))
    covTN.true= covTN.true/(2*length(cutoff))
    covJC.gtr= covJC.gtr/(2*length(cutoff))
    covTN.gtr= covTN.gtr/(2*length(cutoff))
    ## p1 = ggplot(D$df.true, aes(x=s,y=y)) +
    ##     geom_line(color="blue") +
    ##         geom_line(aes(x=x,y=y),data=D$df.jc,color="red",linetype="dashed") +
    ##             geom_line(aes(x=x,y=y),data=D$df.tn,color="darkgreen",linetype="dashed") +
    ##                 geom_line(aes(x=s,y=y2),data=D$df.gtr,color="gold") +
    ##                     xlab('branch length') +
    ##                         ylab('densities') +
    ##                             ggtitle('Blue = TrueQ, Red = JC, Green = TN, Gold = GTR')
    ##plot(p1)
    return( list(Q=D$Q,p=D$p,estT=D$estT,covJC.true=covJC.true, covJC.gtr=covJC.gtr,covTN.true=covTN.true, covTN.gtr=covTN.gtr))
}

computeCoverageSimulations = function(nsites,branch.length,nsim=10000,eta.jc=0.5,eta.tn=0.8, delta=0.001, nrep=1000) {
    cov.jc.true.vec = rep(0,nrep)
    cov.jc.gtr.vec = rep(0,nrep)
    cov.tn.true.vec = rep(0,nrep)
    cov.tn.gtr.vec = rep(0,nrep)
    for(i in 1:nrep){
        print(i)
        D = computeCoverage(nsites, branch.length, nsim, eta.jc, eta.tn, delta)
        cov.jc.true.vec[i] = D$covJC.true
        cov.jc.gtr.vec[i] = D$covJC.gtr
        cov.tn.true.vec[i] = D$covTN.true
        cov.tn.gtr.vec[i] = D$covTN.gtr
    }
    return( list(cov.jc.true = mean(cov.jc.true.vec),cov.jc.gtr = mean(cov.jc.gtr.vec),cov.tn.true = mean(cov.tn.true.vec),cov.tn.gtr = mean(cov.tn.gtr.vec)))
}

nrep = 10
branch.length=0.15
nsites=500
eta=0.2
B=computeCoverageSimulations(nsites,branch.length, nrep=nrep, eta.jc=eta, eta.tn=eta)
B

# eta=vector (same value for JC and TN)
computeCoverageSimulations.df = function(s,eta,nsites=500,nsim=10000, delta=0.001, nrep=1000) {
    df=data.frame(eta=eta, covJC.true=rep(0,length(eta)), covJC.gtr=rep(0,length(eta)),  covTN.true=rep(0,length(eta)), covTN.gtr=rep(0,length(eta)))
    j = 1
    for(e in eta){
        B=computeCoverageSimulations(nsites,s,nsim,eta.jc=e,eta.tn=e, delta, nrep)
        df[j,2] = B$cov.jc.true
        df[j,3] = B$cov.jc.true
        df[j,4] = B$cov.jc.gtr
        df[j,5] = B$cov.tn.true
        j = j+1
    }
    return ( df )
}

nsites=1000
nrep=100
s=0.15
eta=seq(0.1,0.9,by=0.1)
df=computeCoverageSimulations.df(s,eta, nsites=nsites, nrep=nrep)
save(df,file="covDf.Rda")

## todo: check compareCoverage, see that it works and does what we want
# do simulations and simulations.df
# change eta and see how this change
# we want to see how this changes with eta, so maybe a plot for eta, and see how coverage changes with eta





###########################################################################################################
## Errors:

## [1] 253
##         [,1]    [,2]    [,3]    [,4]
## [1,] -0.6630  0.1704  0.0054  0.4871
## [2,]  0.2807 -1.2074  0.0034  0.9232
## [3,]  1.1873  0.4577 -2.4161  0.7712
## [4,]  0.5754  0.6622  0.0042 -1.2418
## [1] 0.4068 0.2470 0.0019 0.3444
## [1] -2.416091
##      [,1] [,2] [,3] [,4]
## [1,]  176    6    1   17
## [2,]    7   97    0   10
## [3,]    1    0    0    0
## [4,]   20   14    0  151
## Error in density.default(tn) (from JCvsTN.r!829WzP#42) : 'x' contains missing values
## In addition: Warning messages:
## 1: In log(1 - p.r * prop.ag/(2 * p.a * p.g) - prop.tv/(2 * p.r)) :
##   NaNs produced
## 2: In rgamma(nsim, mu^2/v, mu/v) : NAs produced

## [1] 282
##         [,1]    [,2]    [,3]    [,4]
## [1,] -5.2268  0.0192  4.2049  1.0027
## [2,]  0.0004 -3.3838  2.3991  0.9843
## [3,]  0.0154  0.4521 -0.6174  0.1499
## [4,]  0.0087  0.4403  0.3557 -0.8048
## [1] 0.0023 0.1168 0.6198 0.2611
## [1] -5.226796
##      [,1] [,2] [,3] [,4]
## [1,]    0    0    1    0
## [2,]    0   40   13    6
## [3,]    2   20  276    5
## [4,]    0    8   12  117
## [1] "prop.ag 0.006 prop.ct 0.028 prop.tv 0.1"
## [1] "p.a 0.003 p.c 0.127 p.g 0.605 p.t 0.265 p.r 0.608 p.y 0.392"
## [1] "mu NaN v 0.00247087299017111"
## Error in density.default(tn) (from #30) : 'x' contains missing values
## In addition: Warning messages:
## 1: In log(1 - p.r * prop.ag/(2 * p.a * p.g) - prop.tv/(2 * p.r)) :
##   NaNs produced
## 2: In rgamma(nsim, mu^2/v, mu/v) : NAs produced

## log(1 - p.r*prop.ag/(2*p.a*p.g) - prop.tv/(2*p.r) )
a=0.608*0.006/(2*0.003*0.605)
b=0.1/(2*0.608)
1-a-b
#-0.08719552
log(1-a-b)
#NaN
## Tamura&Nei paper: mu formula gives good estimates unless prop.ag, prop.ct, prop.tv are very large, and n is relatively small
## In this situation, the formula may become inapplicable because of negative arguments of logarithms.
## Therefore, the application of mu formula should be confined to the case of relatively closely related sequences (mu<1) and a
## relatively large value of n




## TO DO:
# -do many many simulations with s=0.15, find bias:
# -always positive?
# -always TN<JC?
# -how close to each other?
# - we need to know how the bias compares to the spread (because saying bias=0.06 does not say anything)

## then: repeat for different values of s, to see if there is any correlation with bl

## then, do pne similar script but now for spread to choose best eta that covers the dist
# maybe add eta.tn, eta.jc vector to the function? so that for a fixed branch length, it tried many etas?
