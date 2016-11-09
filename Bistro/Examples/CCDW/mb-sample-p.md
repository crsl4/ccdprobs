```shell
cp ccdw.nex.run1.p ccdw.run1.p
cp ccdw.nex.run2.p ccdw.run2.p
```
and remove first line

```R
foo = read.table("ccd.run1.p", header=TRUE)
foo = foo[2002:22001] ## remove burnin
dim(foo) ## 20000 14
require(ggplot2)
ggplot(foo,aes(x=r.A...C))+geom_density()
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
source("../../Scripts/readBistro.R")
bistro = readBistro("randomQ-1")
ggplot(foo, aes(x=pi.A.))+ geom_density(color="green")+geom_density(aes(x=p1,data=bistro,color="violet"))
ggplot(foo, aes(x=sAC))+ geom_density(color="green")+geom_density(aes(x=s1,data=bistro,color="violet"))

require(dplyr)
foo.s = select(foo, sAC,sAG,sAT,sCG,sCT,sGT)
foo.var = apply(as.matrix(foo.s),2,var)
bistro.s = select(bistro, sAC,sAG,sAT,sCG,sCT,sGT)
bistro.var = apply(as.matrix(bistro.s),2,var)
bistro.var / foo.var
```

```R
ccdw = read.dna("ccdw.fasta", format="fasta", as.character=TRUE)
table(ccdw[1,],ccd[2,])
paste(ccdw[,1],collapse="")
patterns = apply(ccdw,2,function(x) paste(x,collapse=""))
str(patterns)
unique(patterns)
table(patterns) ## very rare only cg,gt
```

maybe the problem is that the rates CG, GT are very low

partition by codon
```R
table(patterns[seq(3,1545,3)]) ## third codon position
```