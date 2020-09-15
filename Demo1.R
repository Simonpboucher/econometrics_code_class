rm(list=ls())
#install.packages("tseries")
library('tseries')



x1 <- get.hist.quote(instrument = "MSFT", start = "1990-01-02", quote = "AdjClose",compression = "m") 


x <-   as.numeric(x1$Adjusted)


Tsize <- length(x)



ret <- (x[2:Tsize] - x[1:(Tsize-1)])/x[1:(Tsize-1)]


split.screen(c(2,1))
screen(1)
plot(ts(x[2:Tsize], frequency = 12, start=c(1990,02)), ylab="Prices", main="Microsoft")

screen(2)
plot(ts(ret, frequency = 12, start=c(1990,02)), ylab="Returns")




h <- hist(ret, breaks=50, main="Histogram with Normal Curve", xlab="Returns on Microsoft")

xfit <-seq(min(ret),max(ret),length=40) 
yfit<-dnorm(xfit,mean=mean(ret),sd=sd(ret)) 
yfit <- yfit*diff(h$mids[1:2])*length(ret) 
lines(xfit, yfit, col="blue", lwd=2)



mean(ret)
sqrt(var(ret))


skewness <- function(y){
  TT <- length(y)
  mm <- mean(y)
  vv <- sum((y-mm)^2)/TT
  z <- (y-mm)/sqrt(vv)
  Sk <- mean(z^3)
  return(Sk)
}


kurtosis <- function(y){
  TT <- length(y)
  mm <- mean(y)
  vv <- sum((y-mm)^2)/TT
  z <- (y-mm)/sqrt(vv)
  Ku <- mean(z^4)
  return(Ku)
}


myfun <- function(y){
  TT <- length(y)
  mm <- mean(y)
  vv <- sum((y-mm)^2)/TT
  z <- (y-mm)/sqrt(vv)
  Sk <- mean(z^3)  
  Ku <- mean(z^4)
  return(list(mean=mm,variance=vv,skewness=Sk,kurtosis=Ku))
}



myfun(ret)


Tsize <- length(ret)

temp <- myfun(ret)
shat <- temp$skewness
khat <- temp$kurtosis
JB <- Tsize*(shat^2/6 + (khat-3)^2/24)
JB

pvalue <- 1- pchisq(JB,2)
print(pvalue)

