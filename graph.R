#############################################
###  Optimal acceptance rate versus dimension
#############################################
library(Rfast)

AR.step <- function(current, sigma, alpha = "mh"){         
  # Accept - reject step
  # current = current state of the chain
  # sigma   = standard deviation of the gaussian proposal
  # alpha   = If "barker" specified then Barker's acceptance function, else M-H acceptance function
  
  d <- length(current)
  y <- rmvnorm(1, mu = current, sigma = diag(sigma^2, d))
  s <- dmvnorm(y, mu = rep(0, d), sigma = diag(d))/dmvnorm(current, mu = rep(0, d), sigma = diag(d))
  
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
  b <- numeric(d)
  bm <- matrix(0, nrow = N/K, ncol = d)   # batch means
  a <- 0
  
  for(i in 1:N){
    step = AR.step(init, sigma, alpha)
    b = b + step[1:d]
    a = a + step[d+1]
    init = step[1:d]
    
    if(i%%K == 0){
      bm[i/K, ] <- b/K
      b <- numeric(d)
    }
  }
  return(list(bm, a/N))
}


N = 1e6
d = 1:10
K = 1000      #batch size
###########################################################
#### Minimizing variance between Batch means (M-H) Algorithm
###########################################################

sigma_mh <- numeric(length = length(d))   # optimal sigma
a_mh <- numeric(length = length(d))       # optimal acceptance 
init <- numeric()

for(j in 1:length(d)){
  print(paste0("Doing for d = ", d[j]))
  sigma <- seq(2/sqrt(d[j]), 3/sqrt(d[j]), length.out = 21)
  eff_j <- numeric(length = length(sigma))
  a_j <- numeric(length = length(sigma))
  
  init <- c(init, rnorm(1, 0, 1))
  
  for(i in 1:length(sigma)){
    samp <- AR.sample(init, N, sigma[i], K)
    e <- colVars(samp[[1]])
    a_j[i] <- samp[[2]]
    eff_j[i] <- mean(e)
    cat("\r", i)
  }
  sigma_mh[j] <- sigma[which.min(eff_j)]
  a_mh[j] <- a_j[which.min(eff_j)]
  print(paste0("Done for d = ", j))
}

plot(d, a_mh, xlab = "Dimensions", ylab = "Optimal acceptance (M-H)", type = "b")
a_mh


###########################################################
#### Minimizing variance between Batch means (Barker Algorithm)
###########################################################

sigm_b <- numeric(length = length(d))     # Optimal Sigma
a_b <- numeric(length = length(d))        # Optimal acceptance
init <- numeric()

for(j in 1:length(d)){
  print(paste0("Doing for d = ", d[j]))
  sigma <- seq(2/sqrt(d[j]), 3/sqrt(d[j]), length.out = 21)
  eff_j <- numeric(length = length(sigma))
  a_j <- numeric(length = length(sigma))
  
  init <- c(init, rnorm(1, 0, 1))
  
  for(i in 1:length(sigma)){
    samp <- AR.sample(init, N, sigma[i], K, alpha = "barker")
    e <- colVars(samp[[1]])
    a_j[i] <- samp[[2]]
    eff_j[i] <- mean(e)
    cat("\r", i)
  }
  sigma_b[j] <- sigma[which.min(eff_j)]
  a_b[j] <- a_j[which.min(eff_j)]
  print(paste0("Done for d = ", j))
}

plot(d, a_b, xlab = "Dimensions", ylab = "Optimal acceptance (Barker's)", type = "b")
a_b
