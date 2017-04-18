foo = read.csv("foo.csv")
require(ggplot2)
ggplot(foo, aes(x=t,y=logl)) +
  geom_line()
#  geom_smooth(method="lm",se=FALSE)

fit = lm(logl ~ t + I(t^2), data=foo)
fit3 = lm(logl ~ t + I(t^2) + I(t^3), data=foo)
foo = foo %>%
  mutate(fitted = fitted(fit)) %>%
  mutate(fitted3 = fitted(fit3))

ggplot(foo, aes(x=t,y=logl)) +
  geom_line() +
#  geom_line(aes(y=fitted),color="red") +
#  geom_line(aes(y=fitted3),color="gold") +
  ylim(c(max(foo$logl)-5,max(foo$logl)+0.5)) +
  theme_bw()

foo = foo %>%
  mutate(maxlogl = max(logl)) %>%
  mutate(y1 = exp(logl-maxlogl)) %>%
  mutate(y1 = y1/sum(y1)) %>%
  mutate(y2 = exp(fitted-maxlogl)) %>%
  mutate(y2 = y2/sum(y2)) %>%
  mutate(y3 = exp(fitted3-maxlogl)) %>%
  mutate(y3 = y3/sum(y3)) %>%
  select(-maxlogl)

ggplot(foo, aes(x=t)) +
  geom_line(aes(y=y1), color="black") +
#  geom_line(aes(y=y2), color="red") +
#  geom_line(aes(y=y3), color="gold") +
  geom_hline(yintercept=0)
  
ggplot(foo, aes(x=t,y=dlogl)) +
  geom_line()

ggplot(foo, aes(x=t,y=ddlogl)) +
  geom_line()

require(dplyr)
foo = foo %>%
  mutate(y = exp(logl-max(logl))) %>%
  mutate(y = y / sum(y) / 0.00005) %>%
  mutate(logy = log(y))

ggplot(foo, aes(x=t,y=logy)) +
  geom_line()
#  xlim(c(0,0.01)) +
#  ylim(c(-10,10))

mu = with(foo, sum(t*y)/sum(y))
alpha = 1 + 1/mu

df = data.frame(x=seq(0,0.005,length.out=501))
df = df %>%
  mutate(y=alpha / (x+1)^(alpha+1)) %>%
  mutate(y=y/sum(y) / (0.005/500))

with(df, sum(x*y)/ sum(y) )

ggplot(foo, aes(x=t,y=y)) +
  geom_point() +
  geom_point(aes(x=x,y=y), data=df, color="blue")

foo2 = read.csv("foo2.csv")
ggplot(foo2, aes(x=t,y=logl)) +
  geom_line()
ggplot(foo2, aes(x=t,y=dlogl)) +
  geom_line()
ggplot(foo2, aes(x=t,y=exp(logl-max(logl)))) +
  geom_line()
