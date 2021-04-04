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
d = 100
K = 100      #batch size


###########################################################
#### Minimizing variance between Batch means (Barker Algorithm)
###########################################################

sigma <- c(2.38/sqrt(d), 2.46/sqrt(d))
eff_jb <- numeric(length = length(sigma))
a_j <- numeric(length = length(sigma))
  
first_corr <- numeric(length = length(sigma))
init <- rnorm(d, mean = 0, sd = 1) 
for(i in 1:length(sigma)){
  samp <- AR.sample(init, N, sigma[i], K, alpha = "barker")
  e <- diag(mcse.multi(samp[[1]], r = 1, size = K)$cov)
  a_j[i] <- samp[[2]]
  eff_jb[i] <- mean(e)
    
  first_corr[i] <- mean(diag(matrix(cor(samp[[1]][-1,], samp[[1]][-N,]), ncol = d, nrow = d) ))
  cat("\r", i)
}

