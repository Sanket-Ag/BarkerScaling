#####################################################################
#           Script 2: Normal Linear Mixed Models (NLMM)
# 
# Want to compare the approximate expression for likelihood obtained
# using second order Taylor series exapnasion with the exact
# expression in the simple case of NLMM. Here (p = 3)
#####################################################################

set.seed(11)   #For reproducibility

n <- 100                               #number of observations
s0 <- 1                                 #error variance
sb0 <- 1                               #random effects variance

x1 <- rnorm(n, 0, 1)                   #fixed effect predictors 
x2 <- rnorm(n, 0, 1)
X <- matrix(c(rep(1, n), x1, x2), nrow = n)

a0 <- c(3, 2, 5)                       #fixed effect coeff.
z <- rnorm(n, 0, 1)                    #random effect predictors

# Calculating Likelihoods
Z <- diag(z)
I_n <- diag(rep(1,n))

exact_L <- function(a2 = a0[2], s = s0, sb = sb0){
  a <- c(a0[1], a2, a0[3])
  
  Sigma <- sb*(t(Z)%*%Z) + s*I_n
  c <- ((2*pi)^(-n/2))/sqrt(det(Sigma))
  d <- crossprod(y - X%*%a, solve(Sigma))%*%(y - X%*%a)
  L <- c*exp(-d/2)
}

approx_L <- function(a2 = a0[2], s = s0, sb = sb0){       
  a <- c(a0[1], a2, a0[3])
  
  L0 <- as.numeric(((2*pi*s)^(-n/2))*exp(-crossprod(y - X%*%a)/(2*s)))
  d <- 1 + as.numeric((sb/(2*s*s))*(crossprod(crossprod(Z,y - X%*%a)) - s*crossprod(z)))
  L1 <- L0*d
}


#########################################################################
# First we generate various line plots for different values of error variance
# and random effects variance. Each time the response variable is generated 
# separately.
############################### Varying a2 #############################

a2 <- seq(-5, 5, 0.1)
el <- numeric(length(a2))
al <- numeric(length(a2))
s1 <- c(0.05, 0.5, 1, 5)
sb1 <- c(0.05, 0.5, 1, 5)

# Using model equation
par(mfrow = c(2, 4))
for(s in s1){
  b <- rnorm(n, 0, sd = sqrt(sb0))
  e <- rnorm(n, 0, sd = sqrt(s))
  y <- X%*%a0 + z*b + e
  
  for(i in 1:length(a2)){
    el[i] <- exact_L(a2 = a2[i], s = s)
    al[i] <- approx_L(a2 = a2[i], s = s)
  }
  
  plot(a2, el, type="l", ylab = "L", main = paste0(s))
  plot(a2, al, type="l", ylab = "L", main = "Approx")
}

# Using marginal likelihood
par(mfrow = c(2, 4))
for(s in s1){
  y <- rnorm(100, mean = X%*%a0, sd = sqrt(s + sb0*(z^2)))
  
  for(i in 1:length(a2)){
    el[i] <- exact_L(a2 = a2[i], s = s)
    al[i] <- approx_L(a2 = a2[i], s = s)
  }
  
  plot(a2, el, type="l", ylab = "L", main = paste0(s))
  plot(a2, al, type="l", ylab = "L", main = "Approx")
}

# For different values of sb
par(mfrow = c(2, 4))
for(sb in sb1){
  b <- rnorm(n, 0, sd = sqrt(sb))
  e <- rnorm(n, 0, sd = sqrt(s0))
  y <- X%*%a0 + z*b + e
  
  for(i in 1:length(a2)){
    el[i] <- exact_L(a2 = a2[i], sb = sb)
    al[i] <- approx_L(a2 = a2[i], sb = sb)
  }
  
  plot(a2, el, type="l", ylab = "L", main = paste0(sb))
  plot(a2, al, type="l", ylab = "L", main = "Approx")
}

################################ Cross-combination of s and sb #####################
par(mfrow = c(4,4))
for(s in s1){
  for(sb in sb1){
    b <- rnorm(n, 0, sd = sqrt(sb))
    e <- rnorm(n, 0, sd = sqrt(s))
    y <- X%*%a0 + z*b + e
    
    for(i in 1:length(a2)){
      el[i] <- exact_L(a2 = a2[i], s = s, sb = sb)
      al[i] <- approx_L(a2 = a2[i], s = s, sb = sb)
    }
    
    plot(a2, el, col = "blue", type="l", main = paste0("s = ", s, ", sb = ", sb))
    lines(a2, al, col = "red")
  }
}

#choose sb = 0.05 and s = 5 because we got the best performance there.
s <- 5
sb <- 0.05
b <- rnorm(n, 0, sd = sqrt(sb))
e <- rnorm(n, 0, sd = sqrt(s))
y <- X%*%a0 + z*b + e

for(i in 1:length(a2)){
  el[i] <- exact_L(a2 = a2[i], sb = sb, s = s)
  al[i] <- approx_L(a2 = a2[i], sb = sb, s = s)
}

print(max(el/al))

par(mfrow = c(1, 1))
plot(a2, el, col = "blue", type="l", ylab = "L", lwd=2)
lines(a2, al, col = "red", lty=2, lwd=2)
legend("topleft", c("Exact", "Approx"), col = c("blue", "red"), lty=c(1,2), lwd=c(2,2), cex = 0.7, bty = "n")

####################################################################################################
## Visualizing the slicings for different values of s and sb
# Now, we generate y using model equation by fixing s = 5 and sb = 1.
# Further we make line plots similar to before for the fixed response 
# variables but by changing s and sb. The difference between what is done above
# and this one is that we fix the response in this case.

s <- 5
sb <- 0.05
b <- rnorm(n, 0, sd = sqrt(sb))
e <- rnorm(n, 0, sd = sqrt(s))
y <- X%*%a0 + z*b + e

par(mfrow = c(4,4))
for(s in s1){
  for(sb in sb1){
    for(i in 1:length(a2)){
      el[i] <- exact_L(a2 = a2[i], s = s, sb = sb)
      al[i] <- approx_L(a2 = a2[i], s = s, sb = sb)
    }
    
    plot(a2, el, col = "blue", type="l", main = paste0("s = ", s, ", sb = ", sb))
    lines(a2, al, col = "red")
  }
}