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
abline(v = 1.9)
plot(h, aoar, type = "l")
abline(h = 0.159)
abline(v = 1.9)

s <- seq(0.1, 100, length.out = 1e4)
h <- 1.9
barker <- s/(1+s)
bedard <- pnorm((log(s) - h/2)/sqrt(h)) + s*pnorm((-log(s) - h/2)/sqrt(h))

plot(s, barker, type = "l", col = "blue", ylab="acceptance prob.")
lines(s, bedard, col = "red")
legend("bottomright", legend = c("Barker", "Bedard (h = 1.9)"), col = c("blue", "red"), lty = 1)
