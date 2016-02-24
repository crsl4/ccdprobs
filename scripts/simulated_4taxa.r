# r script to use a simulated data to compare bl summaries between mrbayes and importance sampling
# Claudia February 2016


## to do: modify below to simulate sequence at x, then at 1,2,y,3,4
## not uncertainty in topology, ccdprobs=1,0,0
## for a given tree, fix branch lengths, and simulate data
## then input data to compute branch lengths summaries, many many times
## with true sequence at x, and with estimated seq at x
## compare results
## also compare true Q with estimated Q

## this will let us know if the estimation of seqx is causing the problem


nsites = 500
branch.length = 0.15

doit = function(nsites, branch.length, eta.jc=0.5, eta.tn=0.8, nsim=10000, delta = 0.001) {
    s = seq(delta,2*branch.length,delta)
    Q = randomQ(4,rescale=TRUE)
    print(round(Q$Q,4))
    print(round(Q$p,4))
    print(min(diag(Q$Q)))
    x = simulateSequenceSummary(nsites,Q,branch.length)
    print(x)
    p1 = comparePlot(x,s,Q,nsim,eta.jc,eta.tn)
#    plot(p1)
    return (p1)
}
