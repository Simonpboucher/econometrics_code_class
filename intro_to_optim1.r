rm(list=ls())  ## clear the workspace

# Ordinary least squares example 

Tsize <- 500

error <- rnorm(Tsize)
x <- runif(Tsize)

y <- 0.1 + 0.6*x + 0.1*error

plot(x,y)

olsfct <- function(pars){
  b0 <- pars[1]
  b1 <- pars[2]
  ehat <- y - b0 -b1 * x
  return(sum(ehat^2))
}


pars0 <- c(0,0)
temp <- optim(pars0,olsfct,method="BFGS")

print(temp$par)


ones <- matrix(1,Tsize,1)
X <- cbind(ones,x)
betahat <- solve(t(X)%*%X)%*%t(X)%*%y

print(betahat)


# Non-linear least squares example

error <- rnorm(Tsize)
x <- runif(Tsize)

beta1 <- 0.5
beta2 <- 2.8
y <- beta1 * (x^beta2) + 0.1*error

plot(x,y)

nlsfct <- function(pars){
  b0 <- pars[1]
  b1 <- pars[2]
  ehat <- y- b0*(x^b1)
  return(sum(ehat^2))
}

pars0 <- c(1,0.5)
temp <- optim(pars0,nlsfct, method="BFGS")

print(temp$par)


################################# Maximum Likelihood Example ####################################

rm(list=ls())  ## clear the workspace

## generate artificial data

nn <- 300
X <- matrix(rnorm(nn),nn,1)
e <- rnorm(nn)
beta0 <- 0.5
beta1 <- 1.2
var0 <- 1
y <- beta0 + beta1*X + sqrt(var0)*e


## take data as given and proceed to estimation

YY <- y  

ones <- matrix(1,nn,1)
XX <- cbind(ones,X) 

betahat <- solve(t(XX)%*%XX)%*%t(XX)%*%YY

print('OLS estimate of beta')
print(c(betahat[1],betahat[2]))

yhat <- XX%*%betahat
ehat <- YY - yhat
s2 <- sum(ehat^2)/(nn-dim(XX)[2])

print('OLS estimate of error variance')
print(s2)

print('=========================================')



like <- function(pars){
  
  b0 <- pars[1]
  b1 <- pars[2]
  vv <- pars[3]

  residuals <- YY - XX %*% rbind(b0,b1)

  ssr <- sum(residuals^2)

  llike <- -(nn/2)*log(2*pi) - (nn/2)*log(vv) - (0.5*ssr/vv)  
  
  print(c(pars,llike))

  return(llike)
}



pars0 <- c(0,0,0.1)
temp <- optim(pars0,like,method="L-BFGS-B",control=list(fnscale=-1),lower=c(-100,-100,0.000001),upper=c(100,100,100))

pars1 <- temp$par

b0mle <- pars1[1]
b1mle <- pars1[2]
vmle <-  pars1[3]

print('OLS estimates')
print(c(betahat[1],betahat[2],s2))


print('ML estimates')
print(c(b0mle,b1mle,vmle))

print('check on estimate of error variance')
print((nn-2)*s2/nn)

















