

# Barker's algorithm with multivariate normal proposal
bayes_logit_mcmc <- function(y, X, N = 1e7, prop.sd = .3, start = NULL)
{
	p <- dim(X)[2]
	one.minus.yx <- (1 - y)*X	


	# log likelihood
	logf <- function(beta)
	{
		-sum(beta^2)/2 - sum(log(1 + exp(-X%*%beta))) - sum(one.minus.yx%*%beta)
	}

	# random initialization
	beta.mat <- matrix(0, nrow = N, ncol = p)
	beta.mat[1, ] <- start
	beta <- beta.mat[1, ]
	accept <- 0

	for(i in 2:N)
	{
		U <- runif(1)
		prop <- rnorm(p, mean = beta, sd = 0.4*prop.sd)

		log.rat <- logf(prop) - logf(beta)
		a <- exp(log.rat)/(1 + exp(log.rat))

		if(U < a)
		{
			beta <- prop
			accept <- accept + 1
		}
		beta.mat[i, ] <- beta
	}
	return(list("beta" = beta.mat, "accept.prob" = accept/N))
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
std.dev <- sqrt(diag(summary(glmFit)$cov.scaled))

X <- model.matrix(glmFit)
y <- titanic_comp$Survived

output <- bayes_logit_mcmc(y = y, X = X, start = glmFit$coeff, prop.sd = std.dev)
print(output$accept.prob)

plot.ts(output$beta)
multiESS(output$beta, r = 1)


cov.mat <- cov(output$beta)
res <- list(output$accept.prob, cov.mat)
save(res, file = "bayes_logit_learn")

