foo = read.csv("foo.csv")
require(ggplot2)
ggplot(foo, aes(x=t,y=logl)) +
  geom_line() +
  geom_smooth(method="lm",se=FALSE)

ggplot(foo, aes(x=t,y=dlogl)) +
  geom_line()

ggplot(foo, aes(x=t,y=ddlogl)) +
  geom_line()

require(dplyr)
foo = foo %>%
  mutate(y = exp(logl-max(logl))) %>%
  mutate(y = y / sum(y) / 0.001) %>%
  mutate(logy = log(y))

ggplot(foo, aes(x=t,y=logy)) +
  geom_line() +
  xlim(c(0,0.01)) +
  ylim(c(-10,10))

mu = with(foo, sum(t*y)/sum(y))
alpha = 1 + 1/mu

df = data.frame(x=seq(0,0.005,length.out=501))
df = df %>%
  mutate(y=alpha / (x+1)^(alpha+1)) %>%
  mutate(y=y/sum(y) / (0.005/500)) %>%
  mutate(z=dexp(x,1/mu))

with(df, sum(x*y)/ sum(y) )

ggplot(foo, aes(x=t,y=y)) +
  geom_point() +
  geom_point(aes(x=x,y=y), data=df, color="blue") +
  xlim(c(0,0.005))

foo2 = read.csv("foo2.csv")
ggplot(foo2, aes(x=t,y=logl)) +
  geom_line()
ggplot(foo2, aes(x=t,y=dlogl)) +
  geom_line()
ggplot(foo2, aes(x=t,y=exp(logl-max(logl)))) +
  geom_line()
