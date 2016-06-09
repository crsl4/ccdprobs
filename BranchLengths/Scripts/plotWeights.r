## R script to summarize the weights.txt files after simulation
## we want plot of x=ntaxa, y=mean(ESS%) with std bars
## need to be inside Examples/Simulations/simulations1
## Claudia June 2016

library(plotrix)

## ================= simulations1 ============================================
dat1 = read.table("weights.out", header=FALSE, sep=",") ## change to txt
colnames(dat1) <- c("ntax","rep","time","seed","minW", "meanW", "maxW", "numBigW", "sample", "ess", "essp")
head(dat1)

## here we need to first do rbind for all weights.txt
dat <- dat1
nrep = nrow(dat)/10

meanDat <- aggregate(dat[, c(3,5,6,7,8,9,10,11)], list(dat$ntax), mean)
sdDat <- aggregate(dat[, c(3,5,6,7,8,9,10,11)], list(dat$ntax), sd)
fullDat <- merge(meanDat,sdDat, by="Group.1")
fullDat

## mean ESS% plot
meanplot=with(fullDat,plotCI(Group.1,essp.x*100,(100*essp.y)/sqrt(nrep),
    main ="",xlab="Number of Taxa",
    xlim=range(Group.1)*c(0.9,1.05),
    ylab="",ylim=c(0,100),axes=F,pch=16))
axis(side=1,at=fullDat$Group.1)
axis(side=2,ylim=c(0,100), at=seq(0,100,by=10), labels=seq(0,100,by=10))
mtext("ESS %", line=2.5,side=2,las=0)
lines(meanplot)
lines(x=fullDat$Group.1, y=6000*100/(20000*(2*fullDat$Group.1-3)), col="red")

## mean computer time
meanplot2=with(fullDat,plotCI(Group.1,time.x/60,time.y/(60*sqrt(nrep)),
    main ="",xlab="Number of Taxa",
    xlim=range(Group.1)*c(0.9,1.05),
    ylab="",ylim=c(0,15),axes=F,pch=16))
axis(side=1,at=fullDat$Group.1)
axis(side=2,ylim=c(0,15), at=seq(0,15,by=3), labels=seq(0,15,by=3))
mtext("Computing time (min.)", line=2.5,side=2,las=0)
lines(meanplot2)


