rm(list=ls())  ## clear the workspace

N <- 100

error <- rnorm(N)
y <- matrix(0,N,1)

# Simulate data from an AR(1) model
phi0 <- 0.45
c0 <- 0.1
mu0 <- c0 / (1-phi0)
var0 <- 0.25

stddev0 <- sqrt(var0/(1-phi0^2)) 

y[1] <- rnorm(1, mean=mu0, sd=stddev0)
for (i in 2:N){
  y[i] <- c0 + phi0*y[i-1] + sqrt(var0)*error[i]
}	

plot(ts(y))


like_exact <- function(pars){

	  c1   <- pars[1]
	  phi1 <- pars[2]
	  var1 <- pars[3]
  
	  mu1  <- c1/(1-phi1)

	  loglike <- matrix(0,N,1)
	  loglike[1] <- -0.5*log(2*pi) - 0.5*log( var1/(1-phi1^2) ) - 0.5* (y[1] - mu1)^2 / (var1/(1-phi1^2)) 

	  for (i in 2:N){
		    loglike[i] <- -0.5*log(2*pi) - 0.5*log( var1 ) - 0.5* (y[i] - c1 - phi1*y[i-1])^2 / var1  
	  }	
	  return(sum(loglike))
}


like_cond <- function(pars){

	  c1   <- pars[1]
	  phi1 <- pars[2]
	  var1 <- pars[3]
  
	  mu1  <- c1/(1-phi1)

	  loglike <- matrix(0,N,1)

	  for (i in 2:N){
		    loglike[i] <- -0.5*log(2*pi) - 0.5*log( var1 ) - 0.5* (y[i] - c1 - phi1*y[i-1])^2 / var1  
	  }	
	  return(sum(loglike[2:N]))
}



pars0 <- c(0,0,1)
temp <- optim(pars0,like_exact,method="L-BFGS-B",control=list(fnscale=-1),lower=c(-100,-0.99,0.000001),upper=c(100,0.99,100),hessian=TRUE)

c(c0,phi0,var0)
temp$par


pars0 <- c(0,0,1)
temp <- optim(pars0,like_cond,method="L-BFGS-B",control=list(fnscale=-1),lower=c(-100,-0.99,0.000001),upper=c(100,0.99,100),hessian=TRUE)

c(c0,phi0,var0)
temp$par

# compare with OLS

ones <- matrix(1,N-1,1)
Y <- y[2:N]
X <- cbind(ones,y[1:(N-1)])

betahat <- solve( t(X)%*%X ) %*% t(X) %*% Y

betahat

Yhat <- X %*% betahat
ehat <- Y - Yhat

mean(ehat^2)
				


