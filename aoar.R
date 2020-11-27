N = 1e4    ## Sample size
n = 1e2    ## No. of samples
x = 0

i.grand <- function(z, theta){              ## Integrand in h(theta)
  temp <- -z*sqrt(theta) + theta/2
  ig <- (1 + exp(temp))^(-1)
  
  return(ig)
}

mc <- function(seed, N, theta){             ## Estimates integral for different theta 
  set.seed(seed + x)
  est <- numeric(length = length(theta))
  samp <- rnorm(N, mean = 0, sd = 1)
  for(i in 1:length(theta)){
    vals = theta[i]*i.grand(samp, theta[i])
    est[i] = mean(vals)
  }
  return(est)
}

par(mfrow = c(2, 2))

##########################################################################
######  GLOBAL 
##########################################################################

theta = seq(0, 10, by = 0.1)
est <- matrix(0, nrow = n, ncol = length(theta))

for(i in 1:n){
  est[i, ] <- mc(i, N = N, theta = theta)
  cat("\r", i)
}

final <- colMeans(est)
plot(theta, final, type = "l", col = "blue", xlab = "theta", ylab = "h(theta)")
abline(v = 4)
abline(v = 7)

##########################################################################
######  ZOOMING IN
##########################################################################

theta1 = seq(4, 7, by = 0.01)
est1 <- matrix(0, nrow = n, ncol = length(theta1))
for(i in 1:n){
  est1[i, ] <- mc(i, N = N, theta = theta1)
  cat("\r", i)
}

final1 <- colMeans(est1)
plot(theta1, final1, type = "l", col = "blue", xlab = "theta", ylab = "h(theta)")
abline(v = 5.75)
abline(v = 6.25)

##########################################################################
######  ZOOMING IN
##########################################################################

theta2 = seq(5.751, 6.249, by = 0.001)
est2 <- matrix(0, nrow = n, ncol = length(theta2))
for(i in 1:n){
  est2[i, ] <- mc(i, N = N, theta = theta2)
  cat("\r", i)
}

final2 <- colMeans(est2)
plot(theta2, final2, type = "l", col = "blue", xlab = "theta", ylab = "h(theta)")
abline(v = 5.975)
abline(v = 6.075)

##########################################################################
######  Final
##########################################################################

theta3 = seq(5.9751, 6.0749, by = 0.0005)
est3 <- matrix(0, nrow = n, ncol = length(theta3))
for(i in 1:n){
  est3[i, ] <- mc(i, N = N, theta = theta3)
  cat("\r", i)
}

final3 <- colMeans(est3)
plot(theta3, final3, type = "l", col = "blue", xlab = "theta", ylab = "h(theta)")
t.max <- theta3[which.max(final3)]                        ## Optimum value of theta. l is then sqrt(theta/I) for any I
abline(v = t.max)
print(t.max)

l <- sqrt(t.max)
print(l)

##########################################################################
######  Average Optimal Acceptance rate (AOAR)
##########################################################################

# Theoretical
foo <- rnorm(1e6, mean = 0, sd = 1)
vals <- i.grand(foo, t.max)
aoar <- mean(vals)
print(aoar)


