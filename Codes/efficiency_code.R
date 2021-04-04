##################################################
### Plotting efficiency versus acceptance rate
### Plotting efficiency versus proposal scaling
##################################################


func_barker <- function(x, l)
{
  # Calculates integrand for finding M_B(l)
  foo <- x - log(1 + exp(x)) - log(2*pi*l^2)/2 - (x + l^2/2)^2/(2*l^2)
  return(exp(foo))
}

l <- seq(0.01, 10, length = 1e3)
bark_seq <- numeric(length(l))
mh_seq <- numeric(length(l))

for(i in 1:length(l))
{
  bark_seq[i] <- integrate(func_barker, lower = -Inf, upper = Inf, l = l[i])$value
  mh_seq[i] <- 2*pnorm(-l[i]/2)
}


# Plotting h(l) versus acceptance rate

#pdf("eff_vs_acc.pdf", height = 4, width = 6)
plot(mh_seq, l^2*mh_seq, type = "l", lty = 2, xlab = "Acceptance rate", ylab = "Efficiency")
lines(bark_seq, l^2*bark_seq,  lty = 1)
abline(v = c(bark_seq[which.max(l^2*bark_seq)], mh_seq[which.max(l^2*mh_seq)]), lty = c(1,2))
legend("topright", legend = c("Barker's", "MH"), lty = c(1, 2), cex = 0.8)
text(x = 0.1277295, y = 0.2, labels = expression(alpha ~ "= 0.158"), srt = 90, cex = 0.7)
text(x = 0.266722, y = 0.2, labels = expression(alpha ~ "= 0.234"), srt = 90, cex = 0.7)
#dev.off()

# Plotting h(l) versus l
# Comparison for barker's and MH.
#pdf("ratio.pdf", height = 4, width = 6)
plot(l, bark_seq/mh_seq, type = "l", xlab = expression(l), ylab = "Relative efficiency")
#dev.off()
