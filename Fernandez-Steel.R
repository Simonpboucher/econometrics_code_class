x <- seq(-10,10,0.01)
x <- matrix(x,length(x),1)


fstar <- function(x){
  if (x < 0){
    ans <- dt(xi*x,df)*2/(xi + 1/xi)
  } 
  else
  {
    ans <- dt(x/xi,df)*2/(xi + 1/xi)
  }
  return(ans)
}


df <- 5
xi <- 2

y <- apply(x,1,fstar)

plot(x,dt(x,df),type="l",col="blue")
lines(x,y,col="red")



