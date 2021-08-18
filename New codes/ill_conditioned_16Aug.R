##############################################################################################
#     Running Barker's alogrithm to find the optimal proposal variance. The target is 
#     an ill-conditioned 50 dimensional Gaussian distribution with mean 0. The eigenvalues
#     of the covariance matrix are evenly spaced over the interval [0.1, 100]. Proposals are
#     also Gaussian with iid components with variance l^2/d. The aim will be to find the
#     optimal value of l by minimizing the first order auto-correlations.
##############################################################################################


########### Defining Target ##################

#install.packages("randcorr")
library(doParallel)

library(randcorr)   # Generating a random Correlation matrix

d = 50
v <- seq(0.1, 1, length.out = 30)
v <- c(v , seq(2, 100, length.out = 20))

set.seed(0)
P <- randcorr(d)
E <- diag(sqrt(v))%*%P%*%diag(sqrt(v))   # E is the covariance matrix

S <- eigen(E)         

# Suppose E_inv = t(L)%*%L is the inverse of the civariance matrix of our target distribution
# then L is given in the next line. We only store L since that is all we would need.

L <- diag(sqrt(1/S$values))%*%t(S$vectors)


library(mcmcse)
library(mvtnorm)
library(tictoc)

sampler <- function(samp, ar_step, init, sigma){
 
  # samp : the matrix containing draws from N(0, 1). Used to generate proposals and store the MC.
  # ar_step : a random uniform draw used at the accpet-reject step
  # sigma : proposal standard deviation

  acc_prob <- 0        # keeps track of the no. of acceptances
  samp[1, ] <- init

  for(i in 2:N){

    curr <- samp[i-1, ]                       # current state
    prop <- samp[i-1, ] + sigma*samp[i, ]     # proposed state
    temp <- sum(dnorm(L%*%prop, log = TRUE) - dnorm(L%*%curr, log = TRUE))
    one_by_a <- exp(temp)/(1 + exp(temp))

    if(ar_step[i] <= one_by_a){
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

M <- 1e2   # no. of iterations
N <- 1e6  # length of the chain
K <- 2000  # batch size
sigma <- seq(0.3/sqrt(d), 0.5/sqrt(d), length.out = 51)


##############################################
# Variables to store data

eff_bm <- matrix(0, nrow = M, ncol = length(sigma))      # Store the batch means se calculated using mcmcse package 
eff_fc <- matrix(0, nrow = M, ncol = length(sigma))      # Stores first order autocorrelation
eff_ct <- matrix(0, nrow = M, ncol = length(sigma))      # Stores estimate of convergence time (R&R 2001)
eff_ess <- matrix(0, nrow = M, ncol = length(sigma))     # Stores multiESS()
acc_rate <- matrix(0, nrow = M, ncol = length(sigma))    # Stores acceptance probabilities

# Number of cores
detectCores()
registerDoParallel(cores = detectCores()-2)


doingReps <- function(j)
{
  print(paste0("Doing for m = ", j))

  set.seed(j)
  xi <- matrix(rnorm(d*N, mean = 0, sd = 1), ncol = d)
  prob <- runif(N)
  init <- rmvnorm(1, mean = rep(0, d), sigma = E)

  bm_j <- numeric(length = length(sigma))
  fc_j <- numeric(length = length(sigma))
  ct_j <- numeric(length = length(sigma))
  ess_j <- numeric(length = length(sigma))
  a_j <- numeric(length = length(sigma))

  for(i in 1:length(sigma)){

    samp <- sampler(samp = xi, ar_step = prob, init = init, sigma = sigma[i])
    a_j[i] <- samp[[2]]

    e <- mcse.multi(samp[[1]], r = 1, size = K)$cov
    bm_j[i] <- mean(diag(e))
    ess_j[i] <- multiESS(samp[[1]], covmat = e)

    # minimum autocorrelation
    C <- matrix(cor(samp[[1]][-1,], samp[[1]][-N,]), ncol = d, nrow = d)
    fc_j[i] <- mean(diag(C))
    ct_j[i] <- -1/log(C[1,1])

    cat("\r", i)
  }
  
  print(paste0("Done for m = ", j))
  return(cbind(bm_j, fc_j, ct_j, ess_j, a_j))
 
}

tic()
foo <- foreach(j = 1:M) %dopar% 
{
  doingReps(j)
}

toc()

final.out <- array(unlist(foo), dim = c(length(sigma), 5, M))


eff_bm <- t(final.out[,1, ])
eff_fc <- t(final.out[,2, ])
eff_ct <- t(final.out[,3, ])
eff_ess <- t(final.out[,4, ])
acc_rate <- t(final.out[,5, ])

#########################################
# Save the results

res <- list(sigma, eff_bm, eff_fc, eff_ess, eff_ct, acc_rate)
save(res, file = "ill_gaussian_16Aug")

# Plots
pdf(file = "ill_gaussian_plots_16Aug.pdf")
plot(sigma, colMeans(acc_rate), type = "l")
plot(colMeans(acc_rate), colMeans(eff_bm), type = "l")
plot(colMeans(acc_rate), colMeans(eff_fc), type = "l")
plot(colMeans(acc_rate), colMeans(eff_ess), type = "l")
plot(colMeans(acc_rate), colMeans(eff_ct), type = "l")
dev.off()


############################################################################################################
# END