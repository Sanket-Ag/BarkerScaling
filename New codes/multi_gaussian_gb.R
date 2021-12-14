library(doParallel)
library(tictoc)
##########################################################################################
#     Running Generalized Barker's alogrithm (with r = 2) to find the optimal proposal  
#     variance. The target is a 50 dimensional Gaussian distribution with mean 0 and 
#     identity covariance matrix. Proposals are also Gaussian with iid components with  
#     variance l^2/d. The aim will be to find the optimal value of l by minimizing the 
#     first order auto-correlations.
##########################################################################################

set.seed(678)   

sampler <- function(samp, ar_step, init, sigma){

  # samp : the matrix containing draws from N(0, 1). Used to generate proposals and store the MC.
  # ar_step : a random uniform draw used at the accpet-reject step
  # sigma : proposal standard deviation

  acc_prob <- 0        # keeps track of the no. of acceptances
  samp[1, ] <- init

  for(i in 2:N){

    curr <- samp[i-1, ]                       # current state
    prop <- samp[i-1, ] + sigma*samp[i, ]     # proposed state
    temp <- sum(dnorm(prop, log = TRUE) - dnorm(curr, log = TRUE))
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


#############################################
# Parameters

M <- 1e3  # no. of iterations
N <- 1e6  # length of the chain
d <- 50   # dimensions
sigma <- seq(2/sqrt(d), 3/sqrt(d), length.out = 51)


##############################################
# Variables to store data


eff_fc <- matrix(0, nrow = M, ncol = length(sigma))      # Stores estimated convergence time in bar{x}
eff_ct <- matrix(0, nrow = M, ncol = length(sigma))      # Stores estimated convergence time in x_1 - bar{x}
eff_ff <- matrix(0, nrow = M, ncol = length(sigma))      # Stores estimated convergence time in x_1
acc_rate <- matrix(0, nrow = M, ncol = length(sigma))    # Stores acceptance probabilities


# Number of cores
detectCores()
registerDoParallel(cores = detectCores() - 2)


doingReps <- function(j)
{
  print(paste0("Doing for m = ", j))

  xi <- matrix(rnorm(N*d), ncol = d)
  prob <- runif(N)
  init <- rnorm(d)

  fc_j <- numeric(length = length(sigma))
  ct_j <- numeric(length = length(sigma))
  ff_j <- numeric(length = length(sigma))
  a_j <- numeric(length = length(sigma))

  for(i in 1:length(sigma)){

    samp <- sampler(samp = xi, ar_step = prob, init = init, sigma = sigma[i])
    x <- samp[[1]][, 1]
    a_j[i] <- samp[[2]]

    ff_j[i] <- cor(x[-1], x[-N])

    xbar <- rowMeans(samp[[1]])
    fc_j[i] <- cor(xbar[-1], xbar[-N])

    y <- x - xbar
    ct_j[i] <- cor(y[-1], y[-N])
   
  }
  
  return(cbind(fc_j, ct_j, a_j, ff_j))
 
}

tic()
foo <- foreach(j = 1:M) %dopar% 
{
  doingReps(j)
}
toc()

final.out <- array(unlist(foo), dim = c(length(sigma), 4, M))

eff_fc <- t(final.out[ ,1, ])
eff_ct <- t(final.out[ ,2, ])
eff_ff <- t(final.out[ ,4, ])
acc_rate <- t(final.out[ ,3, ])

##################################3######
# Save the results

res <- list(sigma, eff_fc, eff_ct, eff_ff, acc_rate)
save(res, file = "multi_gaussian_2")

# Plots
pdf(file = "muti_gaussian_2.pdf")
par(mar = c(5.1, 5, 2.1, 2.1))
plot(sigma, colMeans(acc_rate), type = "l", ylab = "acceptance rate", xlab = expression(sigma), cex.main = 2.25, cex.lab = 1.75, cex.axis = 1.75)
abline(h = 0.197)
plot(colMeans(acc_rate), -1/log(colMeans(eff_ct)), type = "l", main = expression('x'[1] - bar(x)), ylab = "convergence time", xlab = "acceptance rate", cex.main = 2.25, cex.lab = 1.75, cex.axis = 1.75)
plot(colMeans(acc_rate), -1/log(colMeans(eff_fc)), type = "l", main = expression(bar(x)), ylab = "convergence time", xlab = "acceptance rate", cex.main = 2.25, cex.lab = 1.75, cex.axis = 1.75)
plot(colMeans(acc_rate), -1/log(colMeans(eff_ff)), type = "l", main = expression('x'[1]), ylab = "convergence time", xlab = "acceptance rate", cex.main = 2.25, cex.lab = 1.75, cex.axis = 1.75)
dev.off()
