N = 1e6   ## Sample size
n = 1e2    ## No. of samples
x = 0

i.grand <- function(z, theta){              ## Integrand in h(theta)
  temp <- -z*sqrt(theta) + theta/2
  ig <- (1 + exp(temp))^(-1)
  
  return(mean(ig))
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

pdf("finding_l.pdf", height = 6, width = 10)
par(mfrow = c(1, 2))

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
abline(v = 6.02)

# ##########################################################################
# ######  ZOOMING IN
# ##########################################################################

# theta1 = seq(4, 7, by = 0.01)
# est1 <- matrix(0, nrow = n, ncol = length(theta1))
# for(i in 1:n){
#   est1[i, ] <- mc(i, N = N, theta = theta1)
#   cat("\r", i)
# }

# final1 <- colMeans(est1)
# plot(theta1, final1, type = "l", col = "blue", xlab = "theta", ylab = "h(theta)")
# abline(v = 5.75)
# abline(v = 6.25)

# ##########################################################################
# ######  ZOOMING IN
# ##########################################################################

# theta2 = seq(5.751, 6.249, by = 0.001)
# est2 <- matrix(0, nrow = n, ncol = length(theta2))
# for(i in 1:n){
#   est2[i, ] <- mc(i, N = N, theta = theta2)
#   cat("\r", i)
# }

# final2 <- colMeans(est2)
# plot(theta2, final2, type = "l", col = "blue", xlab = "theta", ylab = "h(theta)")
# abline(v = 5.975)
# abline(v = 6.075)

##########################################################################
######  Final
##########################################################################

theta3 = seq(6.0315, 6.033, length = 5e2)
est3 <- matrix(0, nrow = n, ncol = length(theta3))
for(i in 1:n){
  est3[i, ] <- mc(i, N = N, theta = theta3)
  cat("\r", i)
}

final3 <- colMeans(est3)
plot(theta3, final3, type = "l", col = "blue", xlab = expression(theta), ylab = expression(h(theta)))
t.max <- theta3[which.max(final3)]                        ## Optimum value of theta. l is then sqrt(theta/I) for any I
abline(v = t.max)
print(t.max)

l <- sqrt(t.max)
print(l) # 2.456109

##########################################################################
######  Average Optimal Acceptance rate (AOAR)
##########################################################################

# Theoretical
aoar <- numeric(length = 1e2)
for(t in 1:1e2)
{
  cat("\r", t)
  foo <- rnorm(1e8, mean = 0, sd = 1)
  aoar[t] <- i.grand(foo, t.max)
}

# .1589796
print(mean(aoar))


