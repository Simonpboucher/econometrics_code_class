
rm(list=ls())  ## clear the workspace


################ Simulate a GARCH model


TT <- 1500
yy <- matrix(0,TT,1)
hh <- matrix(0,TT,1)

e <- rnorm(TT)

mu <- 0
a0 <- 0.001
a1 <- 0.06
a2 <- 0.92

hh[1] <- a0/(1-a1-a2)

yy <- mu + sqrt(hh[1])*e[1]
for (t in 2:TT){
  hh[t] <- a0 + a1*((yy[t-1]-mu)^2) + a2*hh[t-1]
  yy[t] <- mu + sqrt(hh[t])*e[t]
}


plot(ts(yy))


like <- function(pars){

  TT <- length(yy)
  mu <- pars[1]

  a0 <- exp(pars[2])
  a1 <- exp(pars[3])
  a2 <- exp(pars[4])

  ans <- -1e10  #penalty value
  
  llike <- matrix(0,TT,1)

  if (a1 + a2 < 1){

    hh <- matrix(0,TT,1)

    hh[1] <- a0/(1-a1 - a2)

    for (t in 2:TT){
      hh[t] <- a0 + a1*((yy[t-1]-mu)^2) + a2*hh[t-1]
      llike[t] <- -0.5*log(2*pi) - 0.5*log(hh[t]) - 0.5*(yy[t]-mu)^2 /hh[t]
    }

    ans <- sum(llike[2:TT])

  }
#  print(c(mu,a0,a1,a2,ans))

  return(ans)

}




pars0 <- c(mu,log(a0),log(a1),log(a2))

pars1 <- optim(pars0,like,control=list(fnscale=-1))$par

pars1[2:4] <- exp(pars1[2:4])

print(pars1)



 

############## Comparison 

 install.packages("tseries", dependencies=TRUE)

 library("tseries")	
 
 temp <- garch(yy,trace=FALSE)$coef
 temp
 
 # see also the package "rugarch" for other GARCH models
 

############# Another way


like2 <- function(pars){

  TT <- length(yy)
  mu <- pars[1]

  a0 <- pars[2]
  a1 <- pars[3]
  a2 <- pars[4]

  ans <- -1e10  #penalty value
  
  llike <- matrix(0,TT,1)

  if (a1 + a2 < 1){

    hh <- matrix(0,TT,1)
#    hh[1] <- var(yy)
    hh[1] <- a0/(1-a1 - a2)
    for (t in 2:TT){
      hh[t] <- a0 + a1*((yy[t-1]-mu)^2) + a2*hh[t-1]
      llike[t] <- -0.5*log(2*pi) - 0.5*log(hh[t]) - 0.5*(yy[t]-mu)^2 /hh[t]
    }

    ans <- sum(llike[2:TT])

  }
#  print(c(mu,a0,a1,a2,ans))

  return(ans)

}




pars0 <- c(mu,a0,a1,a2)
temp <- optim(pars0,like2,method="L-BFGS-B",lower=c(-5,0.0001,0,0),upper=c(5,2,0.99, 0.99),control=list(fnscale=-1))

pars1 <- temp$par
pars1


 
################# Variance targeting

like3 <- function(pars){

  TT <- length(yy)
  mu <- pars[1]

  vv <- var(yy)
  a1 <- pars[2]
  a2 <- pars[3]

  ans <- -1e10  #penalty value
  
  llike <- matrix(0,TT,1)
  
  if (a1 + a2 < 1){

    hh <- matrix(0,TT,1)
    hh[1] <- vv

    for (t in 2:TT){
      hh[t] <- vv + a1*((yy[t-1]-mu)^2 - vv) + a2*( hh[t-1]-vv )
      llike[t] <- -0.5*log(2*pi) - 0.5*log(hh[t]) - 0.5*(yy[t]-mu)^2 /hh[t]
    }

    ans <- sum(llike[2:TT])

  }
#  print(c(mu,a1,a2,ans))

  return(ans)

}


pars0 <- c(mu,a1,a2)
temp <- optim(pars0,like3,method="L-BFGS-B",lower=c(-5,0,0),upper=c(5,0.99, 0.99),control=list(fnscale=-1))

pars1 <- temp$par
pars1










   











						
  
