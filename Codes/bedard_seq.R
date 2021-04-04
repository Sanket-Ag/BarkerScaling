#################################################################
## Finding AOAR and optimal variance for Bedard's acceptance 
## probability corresponding to different values of h. 
#################################################################

set.seed(1)

theta <- seq(0, 10, length.out = 1e4) #reparameterize l as theta where theta = (l^2)*I
h <- seq(0, 10, length.out = 1e4)
lhat <- numeric(length(h))
aoar <- numeric(length(h))

## Finding optimal values for each h.
for(i in 1:length(h)){
  q <- -sqrt(h[i] + theta)/2  
  y <- theta*2*pnorm(q)
  t <- theta[which.max(y)]
  lhat[i] <- sqrt(t)
  aoar[i] <- 2*pnorm(-sqrt(h[i] + t)/2)    
}


########### Exporting values
bedard <- list(h, lhat, aoar)
save(bedard, file = "bedard.result")

########### Plotting optimal acceptance
#pdf("bedard.pdf", width = 6, height = 6)
plot(h, aoar, type = "l", xlab = expression(h), ylab = "optimal acceptance")
abline(h = 0.158, lty = 2)
abline(v = 1.913, lty = 2)
text(x = 2.1, y = 0.21, labels = expression(h ~ "= 1.913"), srt = 90, cex = 0.8)
text(x = 9.6, y = 0.165, labels = expression(alpha ~ "= 0.158"), cex = 0.8)
#dev.off()
