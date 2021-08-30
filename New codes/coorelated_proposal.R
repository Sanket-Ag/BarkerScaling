##############################################################################################
#     A simple covariance matrix
#     We will do two target distributions. Both are zero mean Gaussian with each component 
#     having unit variance and the correlation between any two components is rho. 
#     The two cases we do are rho = 0.4 and rho = 0.85
##############################################################################################


library(doParallel)
library(tictoc)
library(mvtnorm)

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





############################################################################################################
############################################################################################################
############################################################################################################
# First we will do for rho = 0.85
# Same code will then be used for rho = 0.4 case.


########### Defining Target (rho = 0.85) ##################

d = 50
E <- matrix(0.85, nrow = d, ncol = d)
diag(E) <- 1
S <- eigen(E)

# Suppose E_inv = t(L)%*%L is the inverse of the civariance matrix of our target distribution
# then L is given in the next line. We only store L since that is all we would need.

L <- diag(sqrt(1/S$values))%*%t(S$vectors)

#############################################
# Parameters

M <- 1e3   # no. of iterations
N <- 1e5   # length of the chain
sigma <- seq(2/sqrt(d), 3/sqrt(d), length.out = 51)


##############################################
# Variables to store data

eff_fc <- matrix(0, nrow = M, ncol = length(sigma))      # Stores estimated convergence time in bar{x}
eff_ct <- matrix(0, nrow = M, ncol = length(sigma))      # Stores estimated convergence time in x_1 - bar{x}
acc_rate <- matrix(0, nrow = M, ncol = length(sigma))    # Stores acceptance probabilities


# Number of cores
detectCores()
registerDoParallel(cores = detectCores()-2)


doingReps <- function(j)
{
  print(paste0("Doing for m = ", j))

  set.seed(j)
  xi <- rmvnorm(N, mean = rep(0, d), sigma = E)
  prob <- runif(N)
  init <- rmvnorm(1, mean = rep(0, d), sigma = E)

  fc_j <- numeric(length = length(sigma))
  ct_j <- numeric(length = length(sigma))
  a_j <- numeric(length = length(sigma))

  for(i in 1:length(sigma)){

    samp <- sampler(samp = xi, ar_step = prob, init = init, sigma = sigma[i])
    a_j[i] <- samp[[2]]

    xbar <- rowMeans(samp[[1]])
    fc_j[i] <- cor(xbar[-1], xbar[-N])

    x <- samp[[1]][, 1] - xbar
    ct_j[i] <- cor(x[-1], x[-N])
   
  }
  
  return(cbind(fc_j, ct_j, a_j))
 
}

tic()
foo <- foreach(j = 1:M) %dopar% 
{
  doingReps(j)
}
toc()

final.out <- array(unlist(foo), dim = c(length(sigma), 3, M))

eff_fc <- t(final.out[ ,1, ])
eff_ct <- t(final.out[ ,2, ])
acc_rate <- t(final.out[ ,3, ])

#########################################
# Save the results

res <- list(sigma, eff_fc, eff_ct, acc_rate)
save(res, file = "CoorProp_085")

# Plots
pdf(file = "CoorProp_085.pdf")
plot(sigma, colMeans(acc_rate), type = "l")
abline(h = 0.158)
plot(colMeans(acc_rate), -1/log(colMeans(eff_ct)), type = "l", main = "x_i - bar{x}", ylab = "convergence time")
plot(colMeans(acc_rate), -1/log(colMeans(eff_fc)), type = "l", main = "bar{x}", ylab = "convergence time")
dev.off()




############################################################################################################
############################################################################################################
############################################################################################################
# Now we will do for rho = 0.4
# The code remains same as above. Only the description of the target changes.
# Also the values of sigma.


########### Defining Target (rho = 0.85) ##################

d = 50
E <- matrix(0.4, nrow = d, ncol = d)
diag(E) <- 1
S <- eigen(E)

# Suppose E_inv = t(L)%*%L is the inverse of the civariance matrix of our target distribution
# then L is given in the next line. We only store L since that is all we would need.

L <- diag(sqrt(1/S$values))%*%t(S$vectors)

#############################################
# Parameters

M <- 1e3    # no. of iterations
N <- 1e5    # length of the chain
sigma <- seq(2/sqrt(d), 3/sqrt(d), length.out = 51)


##############################################
# Variables to store data

eff_fc <- matrix(0, nrow = M, ncol = length(sigma))      # Stores estimated convergence time in bar{x}
eff_ct <- matrix(0, nrow = M, ncol = length(sigma))      # Stores estimated convergence time in x_1 - bar{x}
acc_rate <- matrix(0, nrow = M, ncol = length(sigma))    # Stores acceptance probabilities


# Number of cores
detectCores()
registerDoParallel(cores = detectCores()-2)


doingReps <- function(j){

  print(paste0("Doing for m = ", j))

  set.seed(j)
  xi <- rmvnorm(N, mean = rep(0,d), sigma = E)
  prob <- runif(N)
  init <- rmvnorm(1, mean = rep(0, d), sigma = E)

  fc_j <- numeric(length = length(sigma))
  ct_j <- numeric(length = length(sigma))
  a_j <- numeric(length = length(sigma))

  for(i in 1:length(sigma)){

    samp <- sampler(samp = xi, ar_step = prob, init = init, sigma = sigma[i])
    a_j[i] <- samp[[2]]

    xbar <- rowMeans(samp[[1]])
    fc_j[i] <- cor(xbar[-1], xbar[-N])

    x <- samp[[1]][, 1] - xbar
    ct_j[i] <- cor(x[-1], x[-N])
   
  }

  return(cbind(fc_j, ct_j, a_j))

}

tic()
foo <- foreach(j = 1:M) %dopar% 
{
  doingReps(j)
}
toc()

final.out <- array(unlist(foo), dim = c(length(sigma), 3, M))

eff_fc <- t(final.out[ ,1, ])
eff_ct <- t(final.out[ ,2, ])
acc_rate <- t(final.out[ ,3, ])


#########################################
# Save the results

res <- list(sigma, eff_fc, eff_ct, acc_rate)
save(res, file = "CoorProp_04")

# Plots
pdf(file = "CoorProp_04.pdf")
plot(sigma, colMeans(acc_rate), type = "l")
abline(h = 0.158)
plot(colMeans(acc_rate), -1/log(colMeans(eff_ct)), type = "l", main = "x_i - bar{x}", ylab = "convergence time")
plot(colMeans(acc_rate), -1/log(colMeans(eff_fc)), type = "l", main = "bar{x}", ylab = "convergence time")
dev.off()

############################################################################################################
# END