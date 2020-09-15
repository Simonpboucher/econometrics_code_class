
### AR(1) example

rm(list=ls())

Tsize <- 100

y <- matrix(0, Tsize, 1)


mu <- 0.25
phi <- 0.5
sigma2 <- 0.2

e <- rnorm(Tsize)


y[1] <- mu + sqrt(sigma2)*e[1]
for (t in 2:Tsize){
	y[t] <- mu + phi*y[t-1] + sqrt(sigma2)*e[t]		
}


plot(y[1:(Tsize-1)],y[2:Tsize])
plot(ts(y))


ones <- matrix(1,Tsize-1,1)
X <- cbind(ones,y[1:(Tsize-1)])
Y <- y[2:Tsize]
betahat <- solve(t(X)%*%X)%*%t(X)%*%Y

print(betahat)

ehat <- Y - X%*%betahat


s2hat <- t(ehat) %*% ehat / ((Tsize - 1) - dim(X)[2])

V <- matrix(s2hat,2,2)*solve(t(X)%*%X)





likeAR1 <- function(pars){
  
  mu <- pars[1]
  phi <- pars[2]
  vv <- pars[3]
  N <- Tsize-1

  residuals <- Y - X %*% rbind(mu,phi)

  ssr <- sum(residuals^2)

  llike <- -(N/2)*log(2*pi) - (N/2)*log(vv) - (0.5*ssr/vv)  
  
  return(llike)
}



pars0 <- c(0,0,1)
temp <- optim(pars0,likeAR1,method="L-BFGS-B",control=list(fnscale=-1),lower=c(-100,-0.99,0.0001),upper=c(100,0.99,100), hessian=T)

pars1 <- temp$par

muhat     <- pars1[1]
phihat    <- pars1[2]
sigma2hat <-  pars1[3]

c(muhat, phihat, sigma2hat)

c(betahat, s2hat)

Hhat <- temp$hessian

Ihat <- solve(-Hhat)

Ihat
V





#### AR(2) example



Tsize <- 100

y <- matrix(0, Tsize, 1)


mu <- 0.25
phi1 <- 0.5
phi2 <- 0.25
sigma2 <- 0.2

e <- rnorm(Tsize)


y[1] <- mu + sqrt(sigma2)*e[1]
y[2] <- mu + phi1*y[1] + sqrt(sigma2)*e[1]
for (t in 3:Tsize){
	y[t] <- mu + phi1*y[t-1] + phi2*y[t-2] + sqrt(sigma2)*e[t]		
}


plot(ts(y))


check_stationarity <- function(phi1,phi2){
	flag <- 0
	if ((phi2 + phi1 < 1) & (phi2 - phi1 < 1) & (abs(phi2) < 1)) flag <- 1
	return(flag)	
}



check_stationarity(phi1,phi2)


ones <- matrix(1,Tsize-2,1)
X <- cbind(ones,y[2:(Tsize-1)],y[1:(Tsize-2)])
Y <- y[3:Tsize]
betahat <- solve(t(X)%*%X)%*%t(X)%*%Y

print(betahat)

ehat <- Y - X%*%betahat

s2hat <- t(ehat) %*% ehat / ((Tsize - 1) - dim(X)[2])

V <- matrix(s2hat,3,3)*solve(t(X)%*%X)

V



likeAR2 <- function(pars){
  
  mu <- pars[1]
  phi1 <- pars[2]
  phi2 <- pars[3]
  vv <- exp(pars[4])
  N <- Tsize-2

  residuals <- Y - X %*% rbind(mu, phi1, phi2)

  ssr <- sum(residuals^2)

  llike <- -5000
  if (check_stationarity(phi1,phi2) == 1) {

  	llike <- -(N/2)*log(2*pi) - (N/2)*log(vv) - (0.5*ssr/vv)  
  
  }
#  print(c(pars,llike))
  return(llike)
}



pars0 <- c(0,0,0,log(1))
temp <- optim(pars0, likeAR2, method="BFGS", control=list(fnscale=-1), hessian=T)


pars1 <- temp$par

muhat   <- pars1[1]
phi1hat <- pars1[2]
phi2hat <- pars1[3]
sigma2hat <-  exp(pars1[4])


c(muhat, phi1hat, phi2hat, sigma2hat)

c(betahat, s2hat)

Hhat <- temp$hessian

Ihat <- solve(-Hhat)

Ihat


tstat <- phi2hat/sqrt(Ihat[3,3])
tstat

pval <- 2*(1-pnorm(abs(tstat)))
pval




pars0 <- c(0,0.99,0,log(1))
temp <- optim(pars0, likeAR2, method="Nelder-Mead", control=list(fnscale=-1), hessian=T)






## LR test of phi2 = 0



likeAR2restricted <- function(pars){
  
  mu <- pars[1]
  phi1 <- pars[2]
  phi2 <- 0
  vv <- exp(pars[3])
  N <- Tsize-2

  residuals <- Y - X %*% rbind(mu, phi1, phi2)

  ssr <- sum(residuals^2)

  llike <- -5000
  if (check_stationarity(phi1,phi2) == 1) {

  	llike <- -(N/2)*log(2*pi) - (N/2)*log(vv) - (0.5*ssr/vv)  
  
  }
#  print(c(pars,llike))
  return(llike)
}





pars0 <- c(0,0,log(1))
temp <- optim(pars0, likeAR2restricted, method="BFGS", control=list(fnscale=-1), hessian=T)
L0 <- temp$value
L0


# same as:
pars1 <- temp$par
L0 <- likeAR2restricted(pars1)
L0



pars0 <- c(0,0,0,log(1))
temp <- optim(pars0, likeAR2, method="BFGS", control=list(fnscale=-1), hessian=T)
L1 <- temp$value
L1

LR <- 2*(L1 - L0)

pval <- 1 - pchisq(LR,1)
pval


AIC0 <- -2*L0 + 2*3
AIC1 <- -2*L1 + 2*4

AIC0
AIC1
































