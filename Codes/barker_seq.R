#################################################################
## Finding AOAR and optimal variance for Generalized Barker's  
## acceptance probability corresponding to different values of r. 
#################################################################

# This computation is based on Monte Calro aproximation
# since numerical integration is somewhat unstable.

set.seed(16)

N <- 1e5 #sample size for estimating the expectation
n <- 1e4 #number of values for theta to be evaluated between 0 and 10.
r <- 1:10 #for different values of r

vals <- list() 
theta <- seq(5.5, 6.5, length.out = n) #reparameterize l as theta where theta = (l^2)*I

g <- function(s, r){
  # Calculates the r-Barker's balancing function for a given "s".
  # s - argument for the balancing function g_r
  
  temp <- (s^(r+1) - s)/(s^(r+1) - 1)
  
  return(temp)
}

## Finding optimal values for each h.
z <- rnorm(N, mean = 0, sd = 1)
for(i in 1:length(r)){
  print(i)
  vals[[i]] <- list()   #  Store values for each r
  
  h <- numeric(length = n) 
  
  for(j in 1:n){
    foo <- exp(sqrt(theta[j])*z - theta[j]/2)
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


########### Exporting values
barker <- list(r, lhat, acc)
# save(barker, file = "barker.result")


########### Plotting optimal acceptance
pdf("Plots/Gen_barker.pdf", width = 6, height = 6)
plot(r, acc, type = "b", pch = 16, ylim = c(0.155, 0.24), ylab = "optimal acceptance", xlab = "r")
abline(h = 0.234, lty = 2)
text(x = 1.1, y = 0.23, labels = "0.234", cex = 0.8)
dev.off()




