rm(list=ls())  ## clear the workspace

###############

# Simulation of Gauss copula

# Begins with simulatation of multivariate normal distribution, mean zero and correlation matrix P


Tsize <- 500

P <- matrix(0,2,2)

diag(P) <- 1

P[1,2] <- 0.7
P[2,1] <- P[1,2]


Z <- cbind(rnorm(Tsize), rnorm(Tsize))

AU <- chol(P)
AL <- t(AU)

X <- matrix(0,Tsize,2)

for (t in 1:Tsize){
   X[t, ] <- AL %*% Z[t, ]		
}


cor(X[,1],X[,2])


cov(X[,1],X[,2])


tauHat <- cor.test(X[,1],X[,2],method="kendall")$estimate
show(tauHat)

rhoHat <- cor.test(X[,1],X[,2],method="spearman")$estimate
show(rhoHat)



cor(exp(X[,1]),X[,2])

tauHat <- cor.test(exp(X[,1]),X[,2],method="kendall")$estimate
show(tauHat)

rhoHat <- cor.test(exp(X[,1]),X[,2],method="spearman")$estimate
show(rhoHat)





U_Gauss <- apply(X, c(1,2), pnorm)



split.screen(c(1,2))
screen(1)
plot(X[,1],X[,2])

screen(2)
plot(U_Gauss[,1], U_Gauss[,2])




# Simulation of t copula

P <- matrix(0,2,2)

diag(P) <- 1

P[1,2] <- 0.7
P[2,1] <- P[1,2]


Z <- cbind(rnorm(Tsize), rnorm(Tsize))

nu <- 4
W <- rchisq(Tsize,df=nu)

X <- matrix(0,Tsize,2)

for (t in 1:Tsize){
   X[t, ] <-  sqrt(nu/W[t]) * AL %*% Z[t, ]		
}

U_t <- apply(X, c(1,2), pt, df=nu)



split.screen(c(1,2))
screen(1)
plot(X[,1],X[,2])

screen(2)
plot(U_t[,1], U_t[,2])




# Use of package
install.packages("copula")
library('copula')

temp <- gumbelCopula( param = 2 )
U_Gumbel = rCopula(Tsize, temp)


temp <- claytonCopula( param = 2.2 )
U_Clayton = rCopula(Tsize, temp)



split.screen(c(2,2))

screen(1)
plot(U_Gauss[,1], U_Gauss[,2])


screen(2)
plot(U_t[,1], U_t[,2])

screen(3)
plot(U_Gumbel[,1], U_Gumbel[,2])

screen(4)
plot(U_Clayton[,1], U_Clayton[,2])




# Simulation of meta-distributions

Xmeta1 <- apply(U_Gauss, c(1,2), qnorm)    # Gaussian copula, Gaussian marginals
Xmeta2 <- apply(U_t, c(1,2), qnorm)        # t copula,        Gaussian marginals
Xmeta3 <- apply(U_Gumbel, c(1,2), qnorm)   # Gumbel copula,   Gaussian marginals
Xmeta4 <- apply(U_Clayton, c(1,2), qnorm)  # Clayton copula,  Gaussian marginals



split.screen(c(2,2))

screen(1)
plot(Xmeta1[,1], Xmeta1[,2])

screen(2)
plot(Xmeta2[,1], Xmeta2[,2])

screen(3)
plot(Xmeta3[,1], Xmeta3[,2])

screen(4)
plot(Xmeta4[,1], Xmeta4[,2])



# Estimation: full maximum likelihood


like <- function(pars){
  
  mu1   <- pars[1]
  var1  <- pars[2]
  mu2   <- pars[3]
  var2  <- pars[4]
  theta <- pars[5]

  llike <- matrix(0,Tsize,1)
  clike <- matrix(0,Tsize,1)
  llike1 <- matrix(0,Tsize,1)
  llike2 <- matrix(0,Tsize,1)      
  
  temp <- normalCopula( param = theta )

  for (t in 1:Tsize){
  	llike1[t] <- -0.5*log(2*pi) -0.5*log(var1) - 0.5*(Y1[t]-mu1)^2/var1
  	llike2[t] <- -0.5*log(2*pi) -0.5*log(var2) - 0.5*(Y2[t]-mu2)^2/var2  	
  	u1 <- pnorm((Y1[t]-mu1)/sqrt(var1))
  	u2 <- pnorm((Y2[t]-mu2)/sqrt(var2))  
  	clike[t] <- log( dCopula(c(u1,u2),temp) )
  	llike[t] <- clike[t] + llike1[t] + llike2[t]
  }
  ans <- sum(llike)
  print(c(pars,ans))
  return(ans)
}



Y1 <- Xmeta1[,1]
Y2 <- Xmeta1[,2]



pars0 <- c(mean(Y1),var(Y1),mean(Y2),var(Y2),cor(Y1,Y2))

temp <- optim(pars0,like,method="L-BFGS-B",control=list(fnscale=-1),lower=c(-1,0.1,-1,0.1,-0.98),upper=c(1,5,1,5,0.98),hessian=TRUE)

pars1 <- temp$par


Ihat <- -1 * temp$hessian

SE <- sqrt( diag( solve(Ihat) ) )

print(pars1)
print(SE)


## Chech with estimation of bivariate normal


multiLike <- function(pars){
  
	mu <- pars[1:2]
	A  <- matrix(0,2,2)
 	diag(A) <- pars[3:4]
	A[1,2] <- pars[5]
 	
 	Sigma <- t(A) %*% A  
 	SigmaInv <- solve(Sigma) 	
 	SigmaDet <- det(Sigma)
 	
 	llike <- matrix(0,Tsize,1)
 	for (i in 1:Tsize){
 	
  		ssr <- t(y[i,]-mu) %*% SigmaInv %*% (y[i,]-mu)

 	    llike[i] <- -(2/2)*log(2*pi) - 0.5*log(SigmaDet) - 0.5*ssr  
  
	}  
  
    ans <- sum(llike) 
  
    #print(c(pars,ans))

    return(ans)
}


y <- cbind(Y1,Y2)
Sigma0 <- cov(y)

A0 <- chol(Sigma0)

pars0 <- c(mean(Y1),var(Y1),mean(Y2),var(Y2),A0[1,2])

temp <- optim(pars0, multiLike, method="BFGS", control=list(fnscale=-1))

muHat <- temp$par[1:2]
Ahat  <- matrix(0,2,2)
diag(Ahat) <- temp$par[3:4]
Ahat[1,2] <- temp$par[5]
 	
SigmaHat <- t(Ahat) %*% Ahat  


print(muHat)

print(SigmaHat)

print(SigmaHat[1,2]/sqrt(SigmaHat[1,1]*SigmaHat[2,2]))

                                        

#################
# Estimation: psuedo-maximum likelihood or IFM
# 1. Parametric approach


marginalLike <- function(pars){
  mu1   <- pars[1]
  var1  <- pars[2]
  mu2   <- pars[3]
  var2  <- pars[4]

  llike <- matrix(0,Tsize,1)
  llike1 <- matrix(0,Tsize,1)
  llike2 <- matrix(0,Tsize,1)      
   
  for (t in 1:Tsize){
  	llike1[t] <- -0.5*log(2*pi) -0.5*log(var1) - 0.5*(Y1[t]-mu1)^2/var1
  	llike2[t] <- -0.5*log(2*pi) -0.5*log(var2) - 0.5*(Y2[t]-mu2)^2/var2  	
  	llike[t] <- llike1[t] + llike2[t]
  }
  ans <- sum(llike)
  print(c(pars,ans))
  return(ans)
}



Y1 <- Xmeta1[,1]
Y2 <- Xmeta1[,2]




pars0 <- c(mean(Y1),var(Y1),mean(Y2),var(Y2))

temp <- optim(pars0,marginalLike,method="L-BFGS-B",control=list(fnscale=-1),lower=c(-1,0.01,-1,0.01),upper=c(1,5,1,5))

pars1 <- temp$par



mu1hat   <- pars1[1]
var1hat  <- pars1[2]
mu2hat   <- pars1[3]
var2hat  <- pars1[4]



pars0 <- c(mean(Y1),var(Y1),mean(Y2),var(Y2))


mu1hat   <- pars0[1]
var1hat  <- pars0[2]
mu2hat   <- pars0[3]
var2hat  <- pars0[4]



copulaLike <- function(pars){
  
  theta <- pars
     
  clike <- matrix(0,Tsize,1)
    
  temp <- normalCopula( param = theta )

  u1 <- pnorm((Y1-mu1hat)/sqrt(var1hat))
  u2 <- pnorm((Y2-mu2hat)/sqrt(var2hat))  
  clike <-  dCopula( matrix(cbind(u1,u2),Tsize,2), temp,log=TRUE) 

  ans <- sum(clike)
  print(c(pars,ans))
  return(ans)
}



temp <- optimize(copulaLike, interval=c(-0.98, 0.98), maximum=TRUE)

theta1hat <- temp$maximum

theta1hat


pars0 = 0.90

temp <- optim(pars0, copulaLike,method="L-BFGS-B",control=list(fnscale=-1),lower=c(-0.98),upper=c(0.98))

temp


# 1. Non-parametric approach



Fhat <- function(Y,y){
	temp <- as.numeric(Y <= y)
	sum(temp)/(length(Y) + 1)	
}




copulaLike2 <- function(pars){
  
  theta <- pars
  
  clike <- matrix(0,Tsize,1)
  
  temp <- normalCopula( param = theta )
  
  for (t in 1:Tsize){
  	u1 <- Fhat(Y1,Y1[t])
  	u2 <- Fhat(Y2,Y2[t])  	  	  	  	
  	clike[t] <- log( dCopula(c(u1,u2),temp) )
  }
  ans <- sum(clike)
  print(c(pars,ans))
  return(ans)
}



temp <- optimize(copulaLike2, interval=c(-0.98, 0.98), maximum=TRUE)

theta1hat <- temp$maximum


# Method-of-moments estimation using rank correlation, see article by Trivedi (p. 16) for mappings between rank correlations and copula parameters

tauHat <- cor.test(Y1,Y2,method="kendall")$estimate

thetahat <- sin(tauHat * pi /2)

show(thetahat)


rhoHat <- cor.test(Y1,Y2,method="spearman")$estimate

thetahat <- 2*sin(rhoHat * pi /6)

show(thetahat)



















