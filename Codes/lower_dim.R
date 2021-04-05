#############################################
###  Optimal acceptance rate versus dimension
#############################################
library(Rfast)
library(mcmcse)

AR.step <- function(current, sigma, alpha = "mh"){         
  # Accept - reject step
  # current = current state of the chain
  # sigma   = standard deviation of the gaussian proposal
  # alpha   = If "barker" specified then Barker's acceptance function, else M-H acceptance function
  
  d <- length(current)
  y <- rnorm(d, mean = current, sd = sigma)
  s <- prod(dnorm(y, mean = 0, sd = 1)/dnorm(current, mean = 0, sd = 1))
  
  if(alpha == "barker"){
    a = s/(1+s)
  }else{
    a = min(1, s)
  }
 
  u <- runif(1)
  if(u <= a){
    return(c(y, 1))
  }
  else{
    return(c(current, 0))
  }
}

AR.sample <- function(init, N, sigma, K, alpha = "mh"){
  # Accept reject MCMC sample: Returns a list of batch means and acceptance function
  # init   = initial value X0
  # N      = size of the sample
  # sigma  = standard deviation of the gaussian proposal
  # K      = batch-size
  # alpha  = If "barker" specified then Barker's acceptance function, else M-H acceptance function
  
  d <- length(init)
  step <- matrix(0, nrow = N, ncol = d)   # batch means
  a <- 0
  
  step[1, ] <- init
  for(i in 2:N){
     foo <-  AR.step(step[i-1, ], sigma, alpha)
     step[i,] <- foo[1:d]
     a <- a + foo[d+1]

  }
  return(list(step, a/N))
}


N = 1e6
d = 1:10
K = 100      #batch size
###########################################################
#### Minimizing variance between Batch means (M-H) Algorithm
###########################################################

sigma_mh <- numeric(length = length(d))   # optimal sigma
a_mh <- numeric(length = length(d))
eff_mh <- numeric(length = length(d))       # optimal acceptance 

sigmaC_mh <- numeric(length = length(d))   # optimal sigma
aC_mh <- numeric(length = length(d))
effC_mh <- numeric(length = length(d)) 

init <- numeric()

for(j in 1:length(d)){
  print(paste0("Doing for d = ", d[j]))
  sigma <- seq(2/sqrt(d[j]), 3/sqrt(d[j]), length.out = 21)
  eff_j <- numeric(length = length(sigma))
  a_j <- numeric(length = length(sigma))

  first_corr <- numeric(length = length(sigma))

  init <- c(init, rnorm(1, 0, 1))
  
  for(i in 1:length(sigma)){
    samp <- AR.sample(init, N, sigma[i], K)
    e <- diag(mcse.multi(samp[[1]], r = 1, size = K)$cov)
    a_j[i] <- samp[[2]]
    eff_j[i] <- mean(e)

    # minimum autocorrelation
    first_corr[i] <- mean(diag(matrix(cor(samp[[1]][-1,], samp[[1]][-N,]), ncol = d[j], nrow = d[j]) ))
    cat("\r", i)
  }
  eff_mh[j] <- which.min(eff_j)
  sigma_mh[j] <- sigma[eff_mh[j]]
  a_mh[j] <- a_j[eff_mh[j]]

  effC_mh[j] <- which.min(first_corr)
  sigmaC_mh[j] <- sigma[effC_mh[j]]
  aC_mh[j] <- a_j[effC_mh[j]]
  print(paste0("Done for d = ", j))
}

pdf("acc_mh.pdf", height = 6, width = 6)
plot(d, a_mh, xlab = "Dimensions", ylab = "Optimal acceptance (M-H)", type = "b", pch = 16)
dev.off()

pdf("acc_mh_Corr.pdf", height = 6, width = 6)
plot(d, aC_mh, xlab = "Dimensions", ylab = "Optimal acceptance (M-H)", type = "b", pch = 16)
dev.off()
save(a_mh,sigma_mh, eff_mh, aC_mh, sigmaC_mh, effC_mh, file = "opt_MH")


###########################################################
#### Minimizing variance between Batch means (Barker Algorithm)
###########################################################

sigma_b <- numeric(length = length(d))     # Optimal Sigma
a_b <- numeric(length = length(d))        # Optimal acceptance
eff_b <- numeric(length = length(d)) 

sigmaC_b <- numeric(length = length(d))     # Optimal Sigma
aC_b <- numeric(length = length(d))        # Optimal acceptance
effC_b <- numeric(length = length(d)) 
init <- numeric()

for(j in 1:length(d)){
  print(paste0("Doing for d = ", d[j]))
  sigma <- seq(2/sqrt(d[j]), 3/sqrt(d[j]), length.out = 21)
  eff_jb <- numeric(length = length(sigma))
  a_j <- numeric(length = length(sigma))
  
  first_corr <- numeric(length = length(sigma))
  init <- c(init, rnorm(1, 0, 1))
  
  for(i in 1:length(sigma)){
    samp <- AR.sample(init, N, sigma[i], K, alpha = "barker")
    e <- diag(mcse.multi(samp[[1]], r = 1, size = K)$cov)
    a_j[i] <- samp[[2]]
    eff_jb[i] <- mean(e)

    # minimum autocorrelation
    first_corr[i] <- mean(diag(matrix(cor(samp[[1]][-1,], samp[[1]][-N,]), ncol = d[j], nrow = d[j]) ))
    cat("\r", i)
  }
  eff_b[j] <- which.min(eff_jb)  
  sigma_b[j] <- sigma[ eff_b[j] ]
  a_b[j] <- a_j[ eff_b[j] ]

  effC_b[j] <- which.min(first_corr)
  sigmaC_b[j] <- sigma[effC_b[j]]
  aC_b[j] <- a_j[effC_b[j]]  
  print(paste0("Done for d = ", j))
}

# pdf("acc_bark.pdf", height = 6, width = 6)
# plot(d, a_b, xlab = "Dimensions", ylab = "Optimal acceptance (Barker's)", type = "b", pch = 16)
# dev.off()

# pdf("acc_bark_Corr.pdf", height = 6, width = 6)
# plot(d, aC_b, xlab = "Dimensions", ylab = "Optimal acceptance (Barker's)", type = "b", pch = 16)
# dev.off()

# save(a_b, sigma_b, eff_b, effC_b, sigmaC_b, aC_b, file = "opt_bark")
# a_b

# pdf("acc_both.pdf", height = 6, width = 6)
# plot(d, a_b, ylim = c(.14,.5), xlab = "Dimensions", ylab = "Optimal acceptance (Barker's)", type = "b")
# lines(d, a_mh, type = "b", pch = 16, lty = 2)
# legend("topright", legend = c("Barker's", "MH"), lty = 1:2)
# dev.off()

pdf("acc_both_Corr.pdf", height = 6, width = 6)
plot(d, aC_b, ylim = c(.14,.5), xlab = "Dimensions", ylab = "Optimal acceptance (Barker's)", type = "b")
lines(d, aC_mh, type = "b", pch = 16, lty = 2)
legend("topright", legend = c("Barker's", "MH"), lty = 1:2)
dev.off()
