#####################################################################
#####################################################################
##                                                                 ##
##  M.Sc. Project (MTH598A):  Sanket Agrawal                       ##
##                Guided By:  Dootika Vats                         ##
##                                                                 ##
#####################################################################
#####################################################################

#####################################################################
#           Script 1: Normal Linear Mixed Models (NLMM)
# 
# Want to compare the approximate expression for likelihood obtained
# using second order Taylor series exapnasion with the exact
# expression in the simple case of NLMM.
#####################################################################

#####################################################################
# MODEL:  y[i] = x[i]*a + z[i]*b[i] + e[i],  i = 1, ..., n
# Assmp.  e[i] ~ N(0, s) (iid); b[i] ~ N(0, sb) (iid) 
#####################################################################

# For our demonstration, we fix s = s_b = 0, a = 5 and generate n = 100 obs.

set.seed(1)   #For reproducibility

n <- 100                               #number of observations
s0 <- 1                                 #error variance
sb0 <- 1                                #random effects variance
x <- rnorm(n, 0, 1)                    #fixed effect predictors
a0 <- 1                                 #fixed effect coeff.
z <- rnorm(n, 0, 1)                    #random effect predictors
y <- numeric(length = n)               #response

for(i in 1:n){
  y[i] <- rnorm(1, mean = x[i]*a0, sd = sqrt(s0 + sb0*(z[i]^2)))
}

# Calculating Likelihoods
Z <- diag(z)
I_n <- diag(rep(1,n))

exact_L <- function(a = a0, s = s0){                      
  Sigma <- sb0*(t(Z)%*%Z) + s*I_n
  c <- ((2*pi)^(-n/2))/sqrt(det(Sigma))
  d <- crossprod(y - a*x, solve(Sigma))%*%(y - a*x)
  L <- c*exp(-d/2)
}
  
approx_L <- function(a = a0, s = s0){                    
  L0 <- as.numeric(((2*pi*s)^(-n/2))*exp(-crossprod(y - a*x)/(2*s)))
  d <- 1 + (sb0/(2*s*s))*(crossprod(crossprod(Z,y - a*x)) - s*crossprod(z))[1]
  L1 <- L0*d
  L2 <- L0*exp(d-1)
  L <- c(L1, L2)
}

######################## When varying both a and s ########################################
a <- seq(-5, 5, 0.1)
s <- seq(0.1, 5, 0.1)

el <- matrix(data = NA, nrow = length(a), ncol = length(s))
for(i in 1:length(a)){
  for(j in 1:length(s)){
    el[i, j] = exact_L(a = a[i], s = s[j])
  }
}

al <- matrix(data = NA, nrow = length(a), ncol = length(s))
for(i in 1:length(a)){
  for(j in 1:length(s)){
    al[i, j] = approx_L(a[i], s[j])[1]
  }
}

#library(plotly)

fig <- plot_ly(x = s, y = a, z = el)
fig <- fig %>% add_surface()
fig

fig <- plot_ly(x = s, y = a, z = al)
fig <- fig %>% add_surface()
fig
##############################################################################################



################################### When varying only a ######################################
a <- seq(-5, 5, 0.1)
el <- numeric(length(a))
al <- numeric(length(a))
al1 <- numeric(length(a))

for(i in 1:length(a)){
  el[i] <- exact_L(a = a[i])
  al[i] <- approx_L(a = a[i])[1]
}

par(mfrow = c(1, 1))
plot(a, el, type="l", ylab = "L")
plot(a, al, type="l", ylab = "L")


plot(a, 1e5*al, col = "red", type="l", ylab = "L")
lines(a, el, col = "blue")
legend("topleft", c("Exact", "1e5*Approx"), col = c("blue", "red"), lty=c(1,1), cex = 0.7, bty = "n")
##############################################################################################