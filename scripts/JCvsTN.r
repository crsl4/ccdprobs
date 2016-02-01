# r script to simulate many scenarios and determine which conditions
# yield smaller bias between true likelihood and estimated density
# Claudia February 2016

source('branch-length.r')

## Compute the bias between the JC and TN densities, and the true likelihood
## for a given value of branch.length
computeBias = function(nsites,branch.length,nsim=10000,eta.jc=0.5,eta.tn=0.8, delta=0.001) {
    s = seq(delta,2*branch.length,delta)
    Q = randomQ(4,rescale=TRUE)
    print(round(Q$Q,4))
    print(round(Q$p,4))
    print(min(diag(Q$Q)))
    x = simulateSequenceSummary(nsites,Q,branch.length)
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

nrep = 1000
branch.length=0.15
nsites=500

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
    return( list(Q.vec=Q.vec,p.vec=p.vec,estT.vec=estT.vec, bias.jc.true.vec=bias.jc.true.vec, bias.jc.gtr.vec=bias.jc.gtr.vec, bias.tn.true.vec=bias.tn.true.vec, bias.tn.gtr.vec=bias.tn.gtr.vec))
}


## need function to summarize B for a given branch.length: mean(min(Qdiag)), mean(min(p)), mean(max(p)), mean(biases), mean(estT) (and sd)
## need function to repeat this for branch.length in (0,1), and create df
