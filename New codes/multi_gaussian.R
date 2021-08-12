##########################################################################################
#     Running Barker's alogrithm to find the optimal proposal variance. The target is 
#     a 50 dimensional Gaussian distribution with mean 0 and identity covariance matrix.
#     Proposals are also Gaussian with iid components with variance l^2/d. The aim will 
#     be to find the optimal value of l by minimizing the first order auto-correlations.
##########################################################################################


library(mcmcse)


AR.sample <- function(init, sigma){
  

  acc_prob <- 0
  #xi[1, ] <- init

  for(i in 2:N){
    temp <- (2*sigma*sum(xi[i-1, ]*xi[i, ]) + sigma*sigma*sum(xi[i, ]*xi[i, ]))/2
    one_by_a <- 1 + exp(temp)

    if(1/runif(1) >= one_by_a){
      xi[i, ] <- xi[i-1, ] + sigma*xi[i, ]
      acc_prob <- acc_prob + 1
    }
    else{
      xi[i, ] <- xi[i-1, ]
    }
  }

  return(list(xi, acc_prob/N))
}


N = 1e6
d = 50
K = 100      #batch size

set.seed(567)
xi <- matrix(rnorm(d*N, mean = 0, sd = 1), ncol = d, byrow = T)


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

  #print(paste0("Doing for d = ", d[j]))

  sigma <- seq(2/sqrt(d[j]), 3/sqrt(d[j]), length.out = 51)
  eff_jb <- numeric(length = length(sigma))
  a_j <- numeric(length = length(sigma))
  
  first_corr <- numeric(length = length(sigma))
  init <- c(init, rnorm(1, 0, 1))
  
  for(i in 1:length(sigma)){
  
    print(paste0("Doing for sigma = ", sigma[i]))

    samp <- AR.sample(init, sigma[i])
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

pdf("acc_bark.pdf", height = 6, width = 6)
plot(sigma, eff_jb, xlab = "sigma", ylab = "1/Efficiency", type = "l", col = "blue")
plot(sigma, first_corr, xlab = "sigma", ylab = "1/Efficiency", type = "l", col = "blue")
dev.off()


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