# r script to simulate many scenarios and determine which conditions
# yield smaller bias between true likelihood and estimated density
# Claudia February 2016

source('branch-length.r')

## Compute the bias between the JC and TN densities, and the true likelihood
## for a given value of branch.length
computeBias = function(nsites,branch.length,nsim=10000,eta.jc=0.5,eta.tn=0.8, delta=0.001) {
    s = seq(delta,2*branch.length,delta)
    Q = randomQ(4,rescale=TRUE)
    x = simulateSequenceSummary(nsites,Q,branch.length)
    while(0 %in% colSums(x)){ #discard x (counts) with rows of zero
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
    ## how far df.true$y from df.gtr$y
    meanSqGTR = mean((df.true$y-df.gtr$y)^2)
    #print(meanSqGTR)
    mode.true = df.true[which.max(df.true$y),1]
    mode.gtr = df.gtr[which.max(df.gtr$y),1]
    print("modes true")
    print(mode.true)
    print(mode.gtr)
    ## JC bias
    mode.jc = df.jc[which.max(df.jc$y),1]
    print("mode JC")
    print(mode.jc)
    bias.jc.true = mode.true - mode.jc
    bias.jc.gtr = mode.gtr - mode.jc
    print(bias.jc.true)
    print(bias.jc.gtr)
    ## TN bias
    mode.tn = df.tn[which.max(df.tn$y),1]
    print("mode TN")
    print(mode.tn)
    bias.tn.true = mode.true - mode.jc
    bias.tn.gtr = mode.gtr - mode.jc
    print(bias.tn.true)
    print(bias.tn.gtr)
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
    return( list(Q=Q$Q,p=Q$p,estT=Qlist.gtr$branch.length, ms=meanSqGTR, bias.jc.true=bias.jc.true, bias.jc.gtr=bias.jc.gtr, bias.tn.true=bias.tn.true, bias.tn.gtr=bias.tn.gtr))
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

computeBiasSimulations.df = function(s,nsites=500,nsim=10000,eta.jc=0.5,eta.tn=0.8, delta=0.001, nrep=10) {
    df=data.frame(bl=rep(0,length(s)*nrep), minQdiag=rep(0,length(s)*nrep), biasJC.true=rep(0,length(s)*nrep), biasJC.gtr=rep(0,length(s)*nrep),  biasTN.true=rep(0,length(s)*nrep), biasTN.gtr=rep(0,length(s)*nrep), estT=rep(0,length(s)*nrep))
    for(i in 1:length(s)){
        print(paste0("branch length ",s[i]))
        B=computeBiasSimulations(nsites,s[i],nsim,eta.jc,eta.tn, delta, nrep)
        for(j in 1:nrep){
            df[j,1] = s[i]
            df[j,2] = min(diag(B$Q.vec[[j]]))
            df[j,3] = B$bias.jc.true[j]
            df[j,4] = B$bias.jc.gtr[j]
            df[j,5] = B$bias.tn.true[j]
            df[j,6] = B$bias.tn.gtr[j]
            df[j,7] = B$estT.vec[j]
        }
    }
    return (df)
}

s=0.15
df=computeBiasSimulations.df(s,nrep=1000)

## ERROR
## [1] 125
##         [,1]    [,2]    [,3]    [,4]
## [1,] -0.9556  0.0039  0.6559  0.2958
## [2,]  0.2686 -0.7248  0.2550  0.2012
## [3,]  1.2434  0.0070 -1.5025  0.2521
## [4,]  0.4562  0.0045  0.2051 -0.6658
## [1] 0.4565 0.0066 0.2408 0.2960
## [1] -1.502495
##      [,1] [,2] [,3] [,4]
## [1,]  194    0   20   11
## [2,]    0    0    0    0
## [3,]   20    0  113    5
## [4,]   17    1    3  116
## Error in density.default(tn) (from JCvsTN.r@14092cu#42) : 'x' contains missing values
## In addition: Warning messages:
## 1: In log(1 - p.y * prop.ct/(2 * p.c * p.t) - prop.tv/(2 * p.y)) :
##   NaNs produced
## 2: In rgamma(nsim, mu^2/v, mu/v) : NAs produced

## TO DO: figure out error (add prints), do many many simulations with s=0.15, find bias:
# -always positive?
# -always TN<JC?
# -how close to each other?

## then: repeat for different values of s, to see if there is any correlation with bl

## then, do pne similar script but now for spread to choose best eta that covers the dist
# maybe add eta.tn, eta.jc vector to the function? so that for a fixed branch length, it tried many etas?
