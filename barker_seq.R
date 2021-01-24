N <- 1e4 #sample size for estimating the expectation
n <- 1e4 #number of values for theta to be evaluated between 0 and 10.
r <- 1:10 #for different values of r
set.seed(16)

vals <- list()
theta <- seq(0.1, 10, length.out = n)

g <- function(z, r){   # Calculates g_{r}(e^z)
  s <- exp(z)
  temp <- (s^(r+1) - s)/(s^(r+1) - 1)
  
  return(temp)
}

for(i in 1:length(r)){
  
  vals[[i]] <- list()   #  Store values for each r -- h(theta), 
  z <- rnorm(N, mean = 0, sd = 1)
  h <- numeric(length = n) 
  
  for(j in 1:n){
    foo <- sqrt(theta[j])*z - theta[j]/2
    est <- theta[j]*g(foo, r[i])
    h[j] <- mean(est)
  }
  
  k <- which.max(h)
  theta.hat <- theta[k]
  aoar <- h[k]/theta[k]
  
  vals[[i]][[1]] <- h
  vals[[i]][[2]] <- theta.hat
  vals[[i]][[3]] <- aoar
}

####################################
### Collecting AOAR

acc <- numeric(length = length(r))
for(i in 1:length(r)){
  acc[i] <- vals[[i]][[3]]
}

### Collecting theta.hat and taking its square root
lhat <- numeric(length = length(r))
for(i in 1:length(r)){
  lhat[i] <- sqrt(vals[[i]][[2]])
}

### Values
acc
lhat

### Plot
plot(r, acc, type = "b", pch = 16, ylim = c(0.155, 0.24), ylab = "acceptance rate", xlab = "r")
abline(h = 0.234, lty = 2)
text(x = 1.5, y = 0.23, labels = "aoar = 0.234", cex = 0.8)




#####################################
## Miscellaneous
###################################
x <- 3
z <- numeric(length = N)
gz <- numeric(length = N)
est <- numeric(length = N)

for(i in 1:N){
  z[i] <- rnorm(1, 0, 1)
  foo <- sqrt(x)*z[i] - x/2
  gz[i] <- x*g(foo, 5)
  est[i] <- mean(gz[1:i])
}



