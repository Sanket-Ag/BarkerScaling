theta <- seq(0, 10, length.out = 1e3)
h <- seq(0, 10, length.out = 1e3)
lhat <- numeric(length(h))
aoar <- numeric(length(h))

for(i in 1:length(h)){
  q <- -sqrt(h[i] + theta)/2  
  y <- theta*2*pnorm(q)
  t <- theta[which.max(y)]
  lhat[i] <- sqrt(t)
  aoar[i] <- 2*pnorm(-sqrt(h[i] + t)/2)    
}

plot(h, lhat, type = "l")
abline(h = 2.46)
plot(h, aoar, type = "l")
abline(h = 0.159)
