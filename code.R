library(tidyverse)
library(nleqslv)
library(patchwork)

dat = read_csv("agacis.csv")

dat.long <- dat |>    
  dplyr::select(-Annual) |>                   # Remove annual column 
  pivot_longer(cols = c(Jan, Feb, Mar, Apr,   # pivot the column data into one col
                        May, Jun, Jul, Aug, 
                        Sep, Oct, Nov, Dec), 
               values_to = "Precipitation",   # store the values in Precipitation
               names_to = "Month") |>         # store the months in Month
  mutate(Precipitation = case_when(Precipitation == "M" ~ NA_character_,
                                   TRUE                 ~ Precipitation))|>
  mutate(Precipitation = as.numeric(Precipitation))

### Weibull
llweibull <- function(par, data, neg=F){
  # a <- par[1]
  # sigma <- par[2]
  a <- exp(par[1]) # go from (-inf,inf) to (0,inf)
  sigma <- exp(par[2]) # go from (-inf,inf) to (0,inf)
  
  ll <- sum(log(dweibull(x=data, shape=a, scale=sigma)), na.rm=T)
  
  return(ifelse(neg, -ll, ll))
}

weibulls = optim(fn = llweibull,
              par = c(1,1),
              data = dat.long$Precipitation,
              neg=T)

weibull.a = weibulls$par[1]
weibull.sigma = weibulls$par[2]

### Part A
gamma.MLE = function(data, par, neg=F){
  a = par[1]
  b = par[2]
  
  loglik <- sum(log(dgamma(x=data, shape = a, rate = b)), na.rm = T)
  
  return(ifelse(neg, -loglik, loglik))
}

gammas = optim(par = c(1, 1), 
              fn = gamma.MLE,
              data=dat.long$Precipitation,
              neg=T)

gamma.a = gammas$par[1]
gamma.b = gammas$par[2]

dat.gamma <- tibble(x = seq(0,15,length.out=1000)) |>
  mutate(pdf.mle = dgamma(x=x, shape=gammas$par[1], scale=gammas$par[2]))

ggplot() +
  geom_histogram(data = dat.long, aes(x = Precipitation, y = after_stat(density)),
  breaks=seq(0, 15, 1),
  color="grey")+
  geom_hline(yintercept = 0)+
  geom_line(data = dat.gamma, aes(x = x, y = pdf.mle))
  
  ###### PART B
lnorm.MLE = function(data, par, neg=F){
  mu = par[1]
  sigma = par[2]
  
  loglik <- sum(log(dlnorm(x=data, meanlog = mu, sdlog = sigma)), na.rm = T)
  
  return(ifelse(neg, -loglik, loglik))
}

lnorms = optim(par = c(1, 1), 
              fn = lnorm.MLE,
              data=dat.long$Precipitation,
              neg=T)

lnorm.mu = lnorms$par[1]
lnorm.sd = lnorms$par[2]

dat.lnorm <- tibble(x = seq(0,15,length.out=1000)) |>
  mutate(pdf.mle = dlnorm(x=x, meanlog=lnorms$par[1], sdlog=lnorms$par[2]))



ggplot() +
  geom_histogram(data = dat.long, aes(x = Precipitation, y = after_stat(density)),
                 breaks=seq(0, 15, 1),
                 color="grey")+
  geom_hline(yintercept = 0)+
  geom_line(data = dat.lnorm, aes(x = x, y = pdf.mle))

#### Part C

(weibull.gamma = llweibull(dat.long$Precipitation, par = weibulls$par)
                /gamma.MLE(dat.long$Precipitation, par = gammas$par))

#### Part D

(weibull.lnorm = llweibull(dat.long$Precipitation, par = weibulls$par)
  /lnorm.MLE(dat.long$Precipitation, par = lnorms$par))

#### Part E

(gamma.lnorm = gamma.MLE(dat.long$Precipitation, par = gammas$par)
  /lnorm.MLE(dat.long$Precipitation, par = lnorms$par))
  