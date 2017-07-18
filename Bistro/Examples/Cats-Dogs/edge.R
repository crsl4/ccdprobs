edge = read.csv("edge.profile")
require(ggplot2)
require(dplyr)
for ( i in unique(edge$number) )
{
  p = ggplot(filter(edge,number==i), aes(x=t,y=exp(logl-max(logl)))) +
    geom_line() +
    geom_hline(yintercept=0) +
    ggtitle(paste("Edge",i)) +
    theme_bw()
  print(p)
}

