##########################################################################################
#     Running Barker's alogrithm to find the optimal proposal variance. The target is 
#     a 50 dimensional Gaussian distribution with mean 0 and identity covariance matrix.
#     Proposals are also Gaussian with iid components with variance l^2/d. The aim will 
#     be to find the optimal value of l by minimizing the first order auto-correlations.
##########################################################################################


library(mcmcse)
library(tictoc)

sampler <- function(samp, ar_step, sigma){
 
  # samp : the matrix containing draws from N(0, 1). Used to generate proposals and store the MC.
  # ar_step : a random uniform draw used at the accpet-reject step
  # sigma : proposal standard deviation

  acc_prob <- 0        # keeps track of the no. of acceptances

  for(i in 2:N){

    curr <- samp[i-1, ]                       # current state
    prop <- samp[i-1, ] + sigma*samp[i, ]     # proposed state
    temp <- sum(dnorm(prop, log = TRUE) - dnorm(samp[i-1, ], log = TRUE))
    one_by_a <- 1 + exp(-temp)

    if(1/ar_step >= one_by_a){
      samp[i, ] <- prop
      acc_prob <- acc_prob + 1
    }
    else{
      samp[i, ] <- curr
    }
  }

  return(list(samp, acc_prob/N))
}


#############################################
# Parameters

M <- 1e3  # no. of iterations
N <- 1e6  # length of the chain
d <- 50   # dimensions
K <- 100  # batch size
sigma <- seq(2/sqrt(d), 3/sqrt(d), length.out = 51)


##############################################
# Variables to store data

eff_bm <- matrix(0, nrow = M, ncol = length(sigma))      # Store the batch means se calculated using mcmcse package 
eff_fc <- matrix(0, nrow = M, ncol = length(sigma))      # Stores first order autocorrelation
acc_rate <- matrix(0, nrow = M, ncol = length(sigma))    # Stores acceptance probabilities

tic()

for(j in 1:M){

  print(paste0("Doing for m = ", j))

  set.seed(j)
  xi <- matrix(rnorm(d*N, mean = 0, sd = 1), ncol = d)
  prob <- runif(1)

  bm_j <- numeric(length = length(sigma))
  fc_j <- numeric(length = length(sigma))
  a_j <- numeric(length = length(sigma))

  for(i in 1:length(sigma)){

    samp <- sampler(samp = xi, ar_step = prob, sigma = sigma[i])
    e <- diag(mcse.multi(samp[[1]], r = 1, size = K)$cov)
    a_j[i] <- samp[[2]]
    bm_j[i] <- mean(e)

    # minimum autocorrelation
    fc_j[i] <- mean(diag(matrix(cor(samp[[1]][-1,], samp[[1]][-N,]), ncol = d, nrow = d) ))
    cat("\r", i)
  }
  
  eff_bm[j, ] <- bm_j  
  eff_fc[j, ] <- fc_j
  acc_rate[j, ] <- a_j
 
  print(paste0("Done for m = ", j))
}

toc()


#########################################
# Save the results
res <- list(sigma, eff_bm, eff_fc, acc_rate)
save <- (res, file = "multi_gaussian")


############################################################################################################
# END