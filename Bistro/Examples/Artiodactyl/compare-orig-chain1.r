orig = read.table("bistro-orig-0.par",header=FALSE)
chain1 = read.table("bistro-chain1-0.par",header=FALSE)
str(orig)
str(chain1)
library(ggplot2)
q = ggplot(orig,aes(x=V4))+geom_density(col="blue") +
    geom_density(data=chain1,aes(V2),col="green") +
    ggtitle("blue=original, green=thread")
plot(q)
