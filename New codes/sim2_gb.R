
##############################################################################################
#     A simple covariance matrix
#     Target distribution is zero mean Gaussian with each component having unit variance and 
#     the correlation between any two components is rho = 0.85. 
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
    a <- (exp(temp) + exp(2*temp))/(1 + exp(temp) + exp(2*temp))

    if(ar_step[i] <= a){
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

M <- 1e3    # no. of iterations
N <- 1e6    # length of the chain
sigma <- seq(0.8/sqrt(d), 1.2/sqrt(d), length.out = 51)


##############################################
# Variables to store data

eff_fc <- matrix(0, nrow = M, ncol = length(sigma))      # Stores estimated convergence time in bar{x}
eff_ct <- matrix(0, nrow = M, ncol = length(sigma))      # Stores estimated convergence time in x_1 - bar{x}
acc_rate <- matrix(0, nrow = M, ncol = length(sigma))    # Stores acceptance probabilities


# Number of cores
detectCores()
registerDoParallel(cores = detectCores()-6)


doingReps <- function(j)
{
  print(paste0("Doing for m = ", j))

  set.seed(j)
  xi <- matrix(rnorm(d*N, mean = 0, sd = 1), ncol = d)
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
foo <- foreach(j = 1:M, .packages = "mvtnorm") %dopar% 
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
save(res, file = "sim2_gb")

# Plots
pdf(file = "sim2_gb.pdf")
par(mar = c(5.1, 5, 2.1, 2.1))
plot(sigma, colMeans(acc_rate), type = "l", ylab = "acceptance rate", xlab = expression(sigma), cex.main = 2.25, cex.lab = 1.75, cex.axis = 1.75)
abline(h = 0.197)
plot(colMeans(acc_rate), -1/log(colMeans(eff_ct)), type = "l", main = expression('x'[1] - bar(x)), ylab = "convergence time", xlab = "acceptance rate", cex.main = 2.25, cex.lab = 1.75, cex.axis = 1.75)
plot(colMeans(acc_rate), -1/log(colMeans(eff_fc)), type = "l", main = expression(bar(x)), ylab = "convergence time", xlab = "acceptance rate", cex.main = 2.25, cex.lab = 1.75, cex.axis = 1.75)
dev.off()
