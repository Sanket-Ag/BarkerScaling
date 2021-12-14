load("multi_gaussian")
mg_barker <- res

load("multi_gaussian_2")
mg_genbarker <- res 

acc_b <- colMeans(mg_barker[[5]])
acc_gb <- colMeans(mg_genbarker[[5]])
ff_b <- colMeans(mg_barker[[4]])
ff_gb <- colMeans(mg_genbarker[[4]])
fc_b <- colMeans(mg_barker[[2]])
fc_gb <- colMeans(mg_genbarker[[2]])


pdf(file = "sim1_barker.pdf")
par(mar = c(5.1, 5, 2.1, 2.1))
plot(acc_b, -1/log(ff_b), type = "l", main = expression('x'[1]), ylab = "convergence time", xlab = "acceptance rate", cex.main = 2.25, cex.lab = 1.75, cex.axis = 1.75, xlim = c(0.1, 0.28), ylim = c(86, 114))
lines(acc_gb, -1/log(ff_gb), type = "l", lty = 2)
legend("topright", legend = c("r = 1", "r = 2"), lty = 1:2)
plot(acc_b, -1/log(fc_b), type = "l", main = expression(bar(x)), ylab = "convergence time", xlab = "acceptance rate", cex.main = 2.25, cex.lab = 1.75, cex.axis = 1.75, xlim = c(0.1, 0.28), ylim = c(86, 114))
lines(acc_gb, -1/log(fc_gb), type = "l", lty = 2)
legend("topright", legend = c("r = 1", "r = 2"), lty = 1:2)
dev.off()