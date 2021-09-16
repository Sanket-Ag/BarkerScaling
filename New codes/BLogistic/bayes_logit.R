
# Random walk M-H algorithm with multivariate normal proposal
bayes_logit_mh <- function(y, X, N = 1e4, prop.sd = .35, start = NULL)
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
		prop <- rnorm(p, mean = beta, sd = .1*dum2)

		log.rat <- logf(prop) - logf(beta)
		if(log(U) < log.rat)
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

output <- bayes_logit_mh(y =y, X = X, start = glmFit$coeff, prop.sd = .01)
plot.ts(output$beta)
multiESS(output$beta, r = 1)
