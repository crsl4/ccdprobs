# r script to plot df_*_mb.Rda, df_*.Rda
# with aid of bl_*.Rda

library(ggplot2)
load("cats-birds/df_birds1.Rda")
df_is = df
load("df_birds_mb.Rda")
df_mb=df
load("cats-birds/bl_birds1.Rda")

# mean for tree 1
x = as.vector(as.numeric(df_is[1,3:7][data[1:5,2]]))
y = as.vector(as.numeric(df_mb[1,3:7][data[1:5,3]]))
d=data.frame(IS=x,MB=y)
p <- ggplot(d,aes(x=IS,y=MB))+ geom_point()+geom_abline(slope=1)
plot(p)

load("cats-birds/data_birds1.Rda")
data[data$w>0.01,]
load("cats-birds/data_birds0.Rda")
data[data$w>0.01,]
