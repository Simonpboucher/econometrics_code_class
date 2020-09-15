#install.packages('tseries', dependencies=TRUE)

library('tseries')

rm(list=ls())

# S&P 500
x <- get.hist.quote(instrument = "^gspc", start = "2004-01-02", quote = "AdjClose") 


plot(x)

y <- 100*diff(log(x))

plot(ts(y))

ret <- as.numeric(y$Adjusted)


acf(ret)

acf(ret^2)




###################################################
###################################################
# Example of Monte Carlo simulation


rm(list=ls())



totrep <- 1000

stats <- matrix(0, totrep, 1)

rej <- 0


Tsize <- 30

for (irep in 1:totrep){
	
	y <- matrix(0, Tsize, 1)
	
	e <- rnorm(Tsize)
	
	mu <- 0.1
	phi <- 0.99
	
	y[1] <- mu + sqrt(1/(1-phi^2))*e[1]
	
	
	for (t in 2:Tsize){

		y[t] <- mu + phi*y[t-1] + e[t]					

	}
	
	
	X   <- cbind(1, y[ 1 : (Tsize-1)] )	
	
	Y   <- y[2:Tsize]

			
	betahat <- solve( t(X) %*% X) %*% t(X) %*% Y	
		
	ehat <- Y - X %*% betahat								
				
	df <- dim(X)[1] - dim(X)[2]

	sigmahat<- t(ehat)%*%ehat/df
	
	
	varbeta <- matrix(sigmahat,2,2) *  solve( t(X) %*% X)
	
	tstat <- (betahat[2]-phi)/sqrt(varbeta[2,2])
	

	pval <- 2*(1- pt(abs(tstat),df))
	
	if (pval < 0.05) rej <- rej +1
	

	stats[irep] <- betahat[2]
		
}


print(mean(stats)-phi)

print(100*rej/totrep)









