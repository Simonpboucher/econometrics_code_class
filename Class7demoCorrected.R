## Simulate and estimate trivariate normal 

rm(list=ls())  ## clear the workspace

library(MASS)  # needed for mvrnorm

Tsize <- 500  

muTrue <- matrix(c(0.01, -0.02, 0.03), 3, 1)

SigmaTrue <- matrix(0,3,3)

diag(SigmaTrue) <- c(6, 5, 4)

SigmaTrue[1,2] <- 1.3
SigmaTrue[1,3] <- 1.1
SigmaTrue[2,3] <- 0.9
SigmaTrue[2,1] <- SigmaTrue[1,2]
SigmaTrue[3,1] <- SigmaTrue[1,3]
SigmaTrue[3,2] <- SigmaTrue[2,3]


y <- mvrnorm(Tsize, mu = muTrue, Sigma = SigmaTrue )

pairs(y)


like <- function(pars){
  
	mu <- pars[1:3]
	A  <- matrix(0,3,3)
 	diag(A) <- pars[4:6]
	A[1,2] <- pars[7]
	A[1,3] <- pars[8]
	A[2,3] <- pars[9]
 	
 	Sigma <- t(A) %*% A  
 	SigmaInv <- solve(Sigma) 	
 	SigmaDet <- det(Sigma)
 	
 	llike <- matrix(0,Tsize,1)
 	   	
 	for (i in 1:Tsize){
 		 	
  		ssr <- t(y[i,]-mu) %*% SigmaInv %*% (y[i,]-mu)
  		
 	    llike[i] <- -(3/2)*log(2*pi) - 0.5*log(SigmaDet) - 0.5*ssr    
	}    
	
    ans <- sum(llike)   
    #print(c(pars,ans))
    return(ans)
}


mu0 <- apply(y, 2, mean)
Sigma0 <- cov(y)
A0 <- chol(Sigma0)


pars0 <- c(mu0, diag(A0), A0[1,2], A0[1,3], A0[2,3])
temp <- optim(pars0, like, method="BFGS", control=list(fnscale=-1))

muHat <- temp$par[1:3]
Ahat  <- matrix(0,3,3)
diag(Ahat) <- temp$par[4:6]
Ahat[1,2] <- temp$par[7]
Ahat[1,3] <- temp$par[8]
Ahat[2,3] <- temp$par[9]
 	
SigmaHat <- t(Ahat) %*% Ahat  


print(t(muTrue))
print(muHat)

print(t(SigmaTrue))
print(SigmaHat)

print(mu0)
print(Sigma0)




install.packages("numDeriv")   ## library needed for numerical derivatives

library(numDeriv)		 


fct <- function(pars){
  
	mu <- pars[1:3]
	Sigma  <- matrix(0,3,3)
 	diag(Sigma) <- pars[4:6]
	Sigma[1,2] <- pars[7]
	Sigma[1,3] <- pars[8]
	Sigma[2,3] <- pars[9]
	Sigma[2,1] <- Sigma[1,2]
	Sigma[3,1] <- Sigma[1,3]
	Sigma[3,2] <- Sigma[2,3]
 	
 	SigmaInv <- solve(Sigma)
 	
 	SigmaDet <- det(Sigma)
 	
 	llike <- matrix(0,Tsize,1)
 	for (i in 1:Tsize){ 	
 		
  		ssr <- t(y[i,]-mu) %*% SigmaInv %*% (y[i,]-mu)
  		
 	    llike[i] <- -(3/2)*log(2*pi) - 0.5*log(SigmaDet) - 0.5*ssr    
	}   
	 
    ans <- sum(llike)   
    return(ans)
}

pars1 <- c(muHat, diag(SigmaHat), SigmaHat[1,2], SigmaHat[1,3], SigmaHat[2,3])
Ihat <- -1 * hessian(fct, pars1)

SE <- sqrt( diag( solve(Ihat) ) )

print(SE)





