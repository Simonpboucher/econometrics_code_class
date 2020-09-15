rm(list=ls())  ## clear the workspace


totrep <- 1000

stats <- matrix(0, totrep, 1)

rej <- 0


Tsize <- 30

for (irep in 1:totrep){
	
	y <- matrix(0, Tsize, 1)
	
	e <- rnorm(Tsize)
	
	mu <- 0.1
	phi <- 0.95
	
	y[1] <- mu + sqrt(1/(1-phi^2))*e[1]
	
	
	x <- runif(Tsize)
	
	for (t in 2:Tsize){

		y[t] <- mu + phi*x[t-1] + e[t]					

#		y[t] <- mu + phi*y[t-1] + e[t]					

	}
	


	xx   <- cbind(1,x[ 1 : (Tsize-1)]  )
	
#	xx   <- cbind(1,y[ 1 : (Tsize-1)]  )	
	
	
	
	yy   <- y[2:Tsize]
			
	betahat <- solve( t(xx) %*% xx) %*% t(xx) %*% yy	
		
	ehat <- yy - xx %*% betahat								
				
	df <- dim(xx)[1] - dim(xx)[2]
	sigmahat<- t(ehat)%*%ehat/df
	
	
	varbeta <- matrix(sigmahat,2,2) *  solve( t(xx) %*% xx)
	
	tstat <- (betahat[2]-phi)/sqrt(varbeta[2,2])
	

	pval <- 2*(1- pt(abs(tstat),df))
	
	if (pval < 0.05) rej <- rej +1
	

	stats[irep] <- betahat[2]
		
}


print(mean(stats)-phi)

print(100*rej/totrep)