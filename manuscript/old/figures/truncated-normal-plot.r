require(ggplot2)
x0 = 0.7
x1 = seq(-4,x0,0.01)
x2 = seq(x0,4,0.01)
x3 = seq(-4,4,0.01)
df = data.frame(
  x = c(x1,x2),
  y = c( rep(0,length(x1)), dnorm(x2) )
)
df2 = data.frame(
  x=x3,
  y=dnorm(x3)
)
rm(x1,x2,x3)

pdf("truncated-normal.pdf",width=10,height=5)
ggplot(df, aes(x=x,y=y)) +
  geom_hline(yintercept=0) +
  geom_line(aes(x=x,y=y),data=df2,color="gray") +
  geom_line(color="blue") +
  geom_segment(aes(x=-1,xend=0,y=dnorm(-1),yend=dnorm(-1)),
               arrow=arrow(ends="both",type="closed",angle=20,length=unit(0.04,"npc")),
               size=0.1) +
  geom_segment(aes(x=0,xend=0,y=0,yend=dnorm(0)),color="gray",linetype="dotted") +
  annotate("text",label="sigma",parse=TRUE,x=-0.5,y=dnorm(-1)*0.9,cex=12) +
  scale_x_continuous(breaks=c(0,x0),labels=c(expression(mu),expression(x[0]))) +
  xlab('') +
  theme_bw() +
  theme(axis.text.x=element_text(size=30)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  ylab('')
dev.off()
