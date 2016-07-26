## R code to compare the frequency of trees
## to mrbayes

## question: do w divide by length(w) to get std error? not indep sample.

## whales --------------------------------------------------
library(ape)
dat <- read.table("run4.txt", header=TRUE)
dat <- read.table("run3.txt", header=TRUE) ## sample 100, better ess ~3%
str(dat)
logweight = dat$logWt
logw = logweight - max(logweight)
w = exp(logw) / sum(exp(logw))
dat$w <- w

ess = 1/sum(w^2)
ess
ess/length(w) ## 0.02%

means=with(dat, sapply( split( w,tree ),sum ) )
means
means[means>0.01]
which(means>0.01) ## 1902
means[1902]
## (1,2,((3,(((4,5),14),((6,7),((((8,9),11),(10,12)),13)))),(15,(((16,17),((18,19),(20,21))),((22,23),((((24,(27,28)),(25,26)),(29,30)),31))))));
##                                                                                                                                      0.986682

## in mrbayes:
##   tree tree_1 [p = 0.152, P = 0.152] = [&W 0.152173] (2,(((((31,((30,29),((26,25),((28,27),24)))),(23,22)),(((21,20),(19,18)),(17,16))),15),(((((12,10),(13,(11,(9,8)))),(7,6)),(14,(5,4))),3)),1);



dat$mean1 = means[1902]
##dat$mean2 = means[39]
head(dat)
dat$indicator1 = dat$tree == levels(dat$tree)[1902]
##dat$indicator2 = dat$tree == levels(dat$tree)[39]
v1 = sum(dat$w*(dat$indicator1 - dat$mean1)^2)
##v2 = sum(dat$w*(dat$indicator2 - dat$mean2)^2)
sqrt(v1/length(w)) ## 0.0016211
##sqrt(v2/length(w))


## cats --------------------------------------------------
library(ape)
dat <- read.table("run6.txt", header=TRUE)
str(dat)
logweight = dat$logWt
logw = logweight - max(logweight)
w = exp(logw) / sum(exp(logw))
dat$w <- w

ess = 1/sum(w^2)
ess
ess/length(w) ## 12.40%

means=with(dat, sapply( split( w,tree ),sum ) )
means

means[10]
## (1,(2,(((((7,8),9),10),12),11)),(3,((4,5),6)));
##                                        0.101369
means[39]
## (1,2,((3,((4,5),6)),(((((7,8),9),10),12),11)));
##                                        0.890654


## in mrbayes:
##   tree tree_1 [p = 0.781, P = 0.781] = [&W 0.781013] (((11,(12,(10,(9,(8,7))))),2),((6,(5,4)),3),1);
##   tree tree_2 [p = 0.144, P = 0.925] = [&W 0.143809] (2,((11,(12,(10,(9,(8,7))))),((6,(5,4)),3)),1);


dat$mean1 = means[10]
dat$mean2 = means[39]
head(dat)
dat$indicator1 = dat$tree == levels(dat$tree)[10]
dat$indicator2 = dat$tree == levels(dat$tree)[39]
v1 = sum(dat$w*(dat$indicator1 - dat$mean1)^2)
v2 = sum(dat$w*(dat$indicator2 - dat$mean2)^2)
sqrt(v1/length(w)) ## 0.004268
sqrt(v2/length(w)) ## 0.004413


## primates --------------------------------------------------
library(ape)
dat <- read.table("run2.txt", header=TRUE)
str(dat)
logweight = dat$logWt
logw = logweight - max(logweight)
w = exp(logw) / sum(exp(logw))
dat$w <- w

ess = 1/sum(w^2)
ess
ess/length(w) ## 53.37%

means=with(dat, sapply( split( w,tree ),sum ) )
means

means[4]
## (1,2,((((((3,4),5),6),7),(((8,9),10),11)),12));
##                                    9.999065e-01
means[16]
## (1,2,(((((3,(4,5)),6),7),(((8,9),10),11)),12));
##                                    9.346465e-05


## in mrbayes:
##   tree tree_1 [p = 0.995, P = 0.995] = [&W 0.995303] (2,(12,((11,(10,(9,8))),(7,(6,(5,(4,3)))))),1);
##   tree tree_2 [p = 0.005, P = 1.000] = [&W 0.004636] (2,(12,((11,(10,(9,8))),(7,(6,((5,4),3))))),1);


dat$mean1 = means[4]
dat$mean2 = means[16]
head(dat)
dat$indicator1 = dat$tree == levels(dat$tree)[4]
dat$indicator2 = dat$tree == levels(dat$tree)[16]
v1 = sum(dat$w*(dat$indicator1 - dat$mean1)^2)
v2 = sum(dat$w*(dat$indicator2 - dat$mean2)^2)
sqrt(v1/length(w)) ## 0.0001368
sqrt(v2/length(w)) ## 0.0001367


## artiodactyl----------------------------------------------------
library(ape)
dat <- read.table("run16.txt", header=TRUE)
str(dat)
logweight = dat$logWt
logw = logweight - max(logweight)
w = exp(logw) / sum(exp(logw))
dat$w <- w

ess = 1/sum(w^2)
ess
ess/length(w) ## 34.54% (artiodactyl)

means=with(dat, sapply( split( w,tree ),sum ) )
means
## (1,((2,3),(5,6)),4); (1,((2,3),5),(4,6)); (1,((2,5),(4,6)),3);
##         2.498735e-24         5.635129e-28         1.734376e-26
## (1,(2,((4,5),6)),3); (1,(2,((4,6),5)),3); (1,(2,(4,(5,6))),3);
##         2.093593e-10         2.651260e-04         5.899925e-03
## (1,(2,3),((4,5),6)); (1,(2,3),((4,6),5)); (1,(2,3),(4,(5,6)));
##         8.910928e-10         5.682623e-05         1.367697e-02
## (1,2,((3,(5,6)),4)); (1,2,(3,((4,5),6))); (1,2,(3,((4,6),5)));
##         8.127873e-23         1.933638e-06         3.471426e-01
## (1,2,(3,(4,(5,6))));
##         6.329566e-01

means[13] ## 1st tree: 0.63296 (1,2,(3,(4,(5,6))));
means[12] ## 2nd tree: 0.34714 (1,2,(3,((4,6),5)));

## in mrbayes:
## tree tree_1 [p = 0.493, P = 0.493] = [&W 0.493485] (2,(((6,5),4),3),1);
## tree tree_2 [p = 0.474, P = 0.967] = [&W 0.473759] (2,((5,(6,4)),3),1);

dat$mean1 = means[13]
dat$mean2 = means[12]
head(dat)
dat$indicator1 = dat$tree == levels(dat$tree)[13]
dat$indicator2 = dat$tree == levels(dat$tree)[12]
v1 = sum(dat$w*(dat$indicator1 - dat$mean1)^2)
v2 = sum(dat$w*(dat$indicator2 - dat$mean2)^2)
sqrt(v1/length(w)) ## 0.0068165
sqrt(v2/length(w)) ## 0.006733



