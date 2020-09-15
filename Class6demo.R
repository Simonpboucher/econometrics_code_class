rm(list=ls())
install.packages("tseries")
library('tseries')


x1 <- get.hist.quote(instrument = "ge", start = "1970-01-01", quote = "Close") 

ret1 <- 100*diff(log(as.numeric(x1$Close)))

split.screen(c(2,1))
screen(1)
plot(x1)

screen(2)
plot(ts(ret1))





x2 <- get.hist.quote(instrument = "ibm", start = "1970-01-01", quote = "Close") 

ret2 <- 100*diff(log(as.numeric(x2$Close)))

split.screen(c(2,1))
screen(1)
plot(x2)

screen(2)
plot(ts(ret2))




x3 <- get.hist.quote(instrument = "xom", start = "1970-01-01", quote = "Close") 

ret3 <- 100*diff(log(as.numeric(x3$Close)))


split.screen(c(2,1))
screen(1)
plot(x3)

screen(2)
plot(ts(ret3))


x4 <- get.hist.quote(instrument = "^gspc", start = "1970-01-01", quote = "Close") 

ret4 <- 100*diff(log(as.numeric(x4$Close)))


split.screen(c(2,1))
screen(1)
plot(x4)

screen(2)
plot(ts(ret4))



######

split.screen(c(2,2))

screen(1)
plot(ts(ret1), ylab="GE")

screen(2)
plot(ts(ret2), ylab="IBM")

screen(3)
plot(ts(ret3), ylab="Mobil")


screen(4)
plot(ts(ret4), ylab="SP500")

######

Y <-data.frame(GE=ret1,IBM=ret2,Mobil=ret3,SP500=ret4)

cov(Y)

cor(Y)

# Scatterplot matrix
pairs(Y)

###########################

rm(list=ls())

library(mvtnorm)

x1 <-seq(-2, 2, by=0.1)
x2 <-seq(-2, 2, by=0.1)

density <-function(x){
  sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  z <- dmvnorm(x, mean=c(0,0), sigma=sigma)
  return(z)
}

fgrid <- function(x, y, f){
    z <- matrix(nrow=length(x), ncol=length(y))
    for(i in 1:length(x)){
        for(j in 1:length(y)){
            z[i,j] <- f(c(x[i], y[j]))
        }
    }
    return(z)
}

d <- fgrid(x1, x2, density)

contour(x1, x2, d, nlevels=5, main="Bivariate Normal Contours", xlab=expression(x[1]), ylab=expression(x[2]))


###########################

Tsize <- 1000
sigma <- diag(2)
Y1sim <- rmvt(Tsize,sigma=sigma,df=3)

Y2sim <- cbind(rt(Tsize,df=3)*sqrt((3-2)/3), rt(Tsize,df=3)*sqrt((3-2)/3))

split.screen(c(1,2))

screen(1)
plot(Y1sim, xlab=expression(x[1]), ylab=expression(x[2]),main="Bivariate-t")
abline(v=0)
abline(h=0)

screen(2)
plot(Y2sim, xlab=expression(x[1]), ylab=expression(x[2]),main="Independent t")
abline(v=0)
abline(h=0)





