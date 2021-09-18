####################################################
####################################################


library(doParallel)
library(tictoc)
#library(mvtnorm)

set.seed(345)   

sampler <- function(y, X, samp, ar_step, init, sigma)
{

	# samp : the matrix containing draws from N(0, 1). Used to generate proposals and store the MC.
	# ar_step : a random uniform draw used at the accpet-reject step
	# sigma : proposal standard deviation

	one.minus.yx <- (1 - y)*X

	# log likelihood
	logf <- function(beta)
	{	
		-sum(beta^2)/2 - sum(log(1 + exp(-X%*%beta))) - sum(one.minus.yx%*%beta)
	}

	acc_prob <- 0        # keeps track of the no. of acceptances
	samp[1, ] <- init
	curr <- samp[1, ]

	for(i in 2:N)
	{
		prop <- samp[i-1, ] + sigma*samp[i, ]     # proposed state

		log.rat <- logf(prop) - logf(curr)
		a <- exp(log.rat)/(1 + exp(log.rat))

		if(ar_step[i] < a)
		{
			curr <- prop
			acc_prob <- acc_prob + 1
		}
		samp[i, ] <- curr
	}

	return(list("beta" = samp, "accept.prob" = acc_prob/N))
}


############################################################
############################################################
# Bayesian Logistic Regression with Titanic Data
############################################################
############################################################


############################################################
# Getting started.
# install.packages(titanic)
library(titanic)
library(mcmcse)


# Taking the relevant parameters in the data. 
# Leaving out Embarked since not important
titanic_sub <- titanic_train[ ,c(2,3,5,6,7,8,10,12)]
titanic_sub$Pclass <- as.factor(titanic_sub$Pclass)
titanic_sub$Sex <- as.factor(titanic_sub$Sex)
titanic_sub$Embarked <- as.factor(titanic_sub$Embarked)

# Original data has missing age
# colSums(is.na(titanic_sub))

# remove missing data
titanic_comp <- titanic_sub[complete.cases(titanic_sub), ]

# Embark should be only 3 levels, but there is a 4th, which is due to data entry error
# which(titanic_comp$Embarked == "")
titanic_comp <- titanic_comp[-which(titanic_comp$Embarked == ""),  ]

glmFit <- glm(Survived ~., data = titanic_comp, family = "binomial")

X <- model.matrix(glmFit)
y <- titanic_comp$Survived


load("bayes_logit_learn")

d <- dim(X)[2]
E <- res[[2]]

S <- eigen(E)
E_half <- S$vectors%*%diag(sqrt(S$values))


# Suppose E_inv = t(L)%*%L is the inverse of the civariance matrix of our target distribution
# then L is given in the next line. We only store L since that is all we would need.

L <- diag(sqrt(1/S$values))%*%t(S$vectors)

#############################################
# Parameters

M <- 1e2   # no. of iterations
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

  xi <- t(E_half%*%matrix(rnorm(N*d), nrow = d))
  prob <- runif(N)
  init <- E_half%*%rnorm(d)

  fc_j <- numeric(length = length(sigma))
  ct_j <- numeric(length = length(sigma))
  a_j <- numeric(length = length(sigma))

  for(i in 1:length(sigma)){

    samp <- sampler(y = y, X = X, samp = xi, ar_step = prob, init = glmFit$coeff, sigma = sigma[i])
    a_j[i] <- samp$accept.prob

    xbar <- rowMeans(samp$beta)
    fc_j[i] <- cor(xbar[-1], xbar[-N])

    x <- samp$beta[, 1] - xbar
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
save(res, file = "bayes_logit")

# Plots
pdf(file = "bayes_logit.pdf")
plot(sigma, colMeans(acc_rate), type = "l")
abline(h = 0.183)
plot(colMeans(acc_rate), -1/log(colMeans(eff_ct)), type = "l", main = "x_i - bar{x}", ylab = "convergence time")
plot(colMeans(acc_rate), -1/log(colMeans(eff_fc)), type = "l", main = "bar{x}", ylab = "convergence time")
dev.off()
