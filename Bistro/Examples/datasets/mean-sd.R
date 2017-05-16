require(dplyr)
require(ggplot2)

source("../../Scripts/readBistro.R")
b024 = readBistro("bistroFixT024")
m024 = read.table("bmcmc024.par")
m024 = m024[-(1:1000),]

names(m024) = c("logl",paste0("pi",1:4),paste0("s",1:6))


meansd = function(x)
{
 df = x %>% summarize(
  n=n(),
  pi1.mean = mean(pi1),
  pi2.mean = mean(pi2),
  pi3.mean = mean(pi3),
  pi4.mean = mean(pi4),
  s1.mean = mean(s1),
  s2.mean = mean(s2),
  s3.mean = mean(s3),
  s4.mean = mean(s4),
  s5.mean = mean(s5),
  s6.mean = mean(s6),
  pi1.sd = sd(pi1),
  pi2.sd = sd(pi2),
  pi3.sd = sd(pi3),
  pi4.sd = sd(pi4),
  s1.sd = sd(s1),
  s2.sd = sd(s2),  
  s3.sd = sd(s3),
  s4.sd = sd(s4),
  s5.sd = sd(s5),
  s6.sd = sd(s6)) %>%
   mutate(pi1.scale = pi1.mean*(1-pi1.mean)/pi1.sd^2) %>%
   mutate(pi2.scale = pi2.mean*(1-pi2.mean)/pi2.sd^2) %>%
   mutate(pi3.scale = pi3.mean*(1-pi3.mean)/pi3.sd^2) %>%
   mutate(pi4.scale = pi4.mean*(1-pi4.mean)/pi4.sd^2) %>%   
   mutate(s1.scale = s1.mean*(1-s1.mean)/s1.sd^2) %>%   
   mutate(s2.scale = s2.mean*(1-s2.mean)/s2.sd^2) %>%   
   mutate(s3.scale = s3.mean*(1-s3.mean)/s3.sd^2) %>%   
   mutate(s4.scale = s4.mean*(1-s4.mean)/s4.sd^2) %>%   
   mutate(s5.scale = s5.mean*(1-s5.mean)/s5.sd^2) %>%   
   mutate(s6.scale = s6.mean*(1-s6.mean)/s6.sd^2)   
 return(df)
}

b024.msd = meansd(b024)
m024.msd = meansd(m024)

## An example tern plot
require(ggtern)
ggtern(m024,aes(s1,s2,s3)) +
  geom_point(alpha=0.2,color="blue") +
  geom_point(data=b024,color="red",alpha=0.2) +
  tern_limits(0.34,0.68,0.27)

my_tern_limits = function (T = 1, L = 1, R = 1, ...) 
{
  ret <- list()
  if (!all(sapply(list(T, L, R), function(x) {
    length(x) == 1 && is.numeric(x)
  }))) 
    stop("Arguments T, L and R must be numeric and scalar", 
         call. = FALSE)
  tryCatch({
    s <- solve(diag(-1, 3, 3) + 1, c(1 - T, 1 - L, 1 - R))
    o <- function(x) {
      x[order(x)]
    }
    T <- o(c(s[1], T))
    L = o(c(s[2], L))
    R = o(c(s[3], R))
    lims <- list(T, L, R)
    if (any(sapply(lims, function(x) {
      diff(x) == 0
    }))) 
      stop("Invalid limits, solution produces zero ranges on some scales", 
           call. = FALSE)
    if (any(sapply(lims, function(x) {
      min(x) < 0 | max(x) > 1
    }))) 
      warning("Solution to limits produces range outside of [0,1] for some scales", 
              call. = FALSE)
    ret <- list(scale_T_continuous(limits = T, ...), scale_L_continuous(limits = L, 
                                                                        ...), scale_R_continuous(limits = R, ...))
  }, error = function(e) {
    warning(e)
  })
  invisible(ret)
}

triplot = function(df1,df2,aes,color1="blue",color2="red",alpha1=0.1,alpha2=0.1)
{
  if ( nrow(df1) > 1000 )
    df1 = df1[sample(1:nrow(df1),size=1000,replace=FALSE),]
  if ( nrow(df2) > 1000 )
    df2 = df2[sample(1:nrow(df2),size=1000,replace=FALSE),]
  df = bind_rows(df1,df2)
  df = df %>%
    select_(aes$x,aes$y,aes$z)
  df = as.matrix(df)
  rsum = apply(df,1,sum)
  df = diag(1/rsum) %*% df
  m = apply(df,2,max)
  m = m + 0.05 * (1-m)
  p = ggplot(mapping=aes) +
    geom_point(data=df1,color=color1,alpha=alpha1) +
    geom_point(data=df2,color=color2,alpha=alpha2) +
    coord_tern() +
   my_tern_limits(m[2],m[1],m[3])
  return(p)
}
