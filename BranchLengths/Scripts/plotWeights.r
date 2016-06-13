## R script to summarize the weights.txt files after simulation
## we want plot of x=ntaxa, y=mean(ESS%) with std bars
## need to be inside Examples/Simulations/simulations1
## Claudia June 2016

library(plotrix)

## ================= simulations1 ============================================
path = "niceBL-mvnormal/"
dat1 = read.table(paste0(path,"weights_nice_1-10.txt"), header=FALSE, sep=",")
dat2 = read.table(paste0(path,"weights_nice1-5.txt"), header=FALSE, sep=",")
dat3 = read.table(paste0(path,"weights_nice6-10.txt"), header=FALSE, sep=",")
dat4 = read.table(paste0(path,"weights_nice11-15.txt"), header=FALSE, sep=",")
dat5 = read.table(paste0(path,"weights_nice16-20.txt"), header=FALSE, sep=",")
colnames(dat1) <- c("ntax","rep","time","seed","minW", "meanW", "maxW", "numBigW", "sample", "ess", "essp")
colnames(dat2) <- c("ntax","rep","time","seed","minW", "meanW", "maxW", "numBigW", "sample", "ess", "essp")
colnames(dat3) <- c("ntax","rep","time","seed","minW", "meanW", "maxW", "numBigW", "sample", "ess", "essp")
colnames(dat4) <- c("ntax","rep","time","seed","minW", "meanW", "maxW", "numBigW", "sample", "ess", "essp")
colnames(dat5) <- c("ntax","rep","time","seed","minW", "meanW", "maxW", "numBigW", "sample", "ess", "essp")
head(dat1)

## here we need to first do rbind for all weights.txt
dat <- rbind(dat1,dat2,dat3,dat4,dat5)
head(dat)
nrep = nrow(dat)/10

## noticed problem on
dd <- subset(dat,ntax==6)
##     ntax rep time   seed         minW       meanW       maxW numBigW sample
## 4      6   1 2137  85708 2.069011e-09 0.001000000 0.04391958       8   1000
## 14     6   2  269  92520 6.578536e-10 0.001000000 0.02463126      13   1000
## 24     6   3  270  98565 2.570466e-10 0.001000000 0.02285549       4   1000
## the first replicate was when the computer fell asleep, so will substitute with the mean
dat[4,3] <- mean(dd$time[-1])

meanDat <- aggregate(dat[, c(3,5,6,7,8,9,10,11)], list(dat$ntax), mean)
sdDat <- aggregate(dat[, c(3,5,6,7,8,9,10,11)], list(dat$ntax), sd)
meanDat
sdDat
fullDat <- merge(meanDat,sdDat, by="Group.1")
fullDat

## mean ESS% plot
pdf("ess.pdf")
par(mar=c(3.1,3.6,.5,.5), mgp=c(1.5,.5,0), tck=-0.01, las=1, yaxs="r", xaxs="r")
meanplot=with(fullDat,plotCI(Group.1,essp.x*100,(100*essp.y)/sqrt(nrep),
    main ="",xlab="Number of Taxa",
    xlim=range(Group.1)*c(0.9,1.05),
    ylab="",ylim=c(0,100),axes=F,pch=16))
axis(side=1,at=fullDat$Group.1)
axis(side=2,ylim=c(0,100), at=seq(0,100,by=10), labels=seq(0,100,by=10))
mtext("ESS %", line=2.5,side=2,las=0)
lines(meanplot)
lines(x=fullDat$Group.1, y=6000*100/(20000*(2*fullDat$Group.1-3)), col="red")
dev.off()

## mean computer time
pdf("time.pdf")
par(mar=c(3.1,3.6,.5,.5), mgp=c(1.5,.5,0), tck=-0.01, las=1, yaxs="r", xaxs="r")
meanplot2=with(fullDat,plotCI(Group.1,time.x/60,time.y/(60*sqrt(nrep)),
    main ="",xlab="Number of Taxa",
    xlim=range(Group.1)*c(0.9,1.05),
    ylab="",ylim=c(0,15),axes=F,pch=16))
axis(side=1,at=fullDat$Group.1)
axis(side=2,ylim=c(0,15), at=seq(0,15,by=3), labels=seq(0,15,by=3))
mtext("Computing time (min.)", line=2.5,side=2,las=0)
lines(meanplot2)
dev.off()

