# Plots for the isotropic simulations

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


pdf(file = "sim1_x1.pdf", height = 4, width = 5)
par(mar = c(5.1, 5, 2.1, 2.1))
plot(acc_b, -1/log(ff_b), type = "l", main = expression('x'[1]), ylab = "convergence time", xlab = "acceptance rate", cex.main = 1.25, cex.lab = 1.25, cex.axis = 1.25, xlim = c(0.1, 0.28), ylim = c(86, 114))
lines(acc_gb, -1/log(ff_gb), type = "l", lty = 2)
legend("topright", legend = c("r = 1", "r = 2"), lty = 1:2)
dev.off()

pdf(file = "sim1_mean.pdf", height = 4, width = 5)
par(mar = c(5.1, 5, 2.1, 2.1))
plot(acc_b, -1/log(fc_b), type = "l", main = expression(bar(x)), ylab = "convergence time", xlab = "acceptance rate", cex.main = 1.25, cex.lab = 1.25, cex.axis = 1.25, xlim = c(0.1, 0.28), ylim = c(86, 114))
lines(acc_gb, -1/log(fc_gb), type = "l", lty = 2)
legend("topright", legend = c("r = 1", "r = 2"), lty = 1:2)
dev.off()



# Plots for the non-isotropic simulations

load("25Aug_rho085")
mg_barker <- res

load("sim2_gb")
mg_genbarker <- res


acc_b <- colMeans(mg_barker[[4]])
acc_gb <- colMeans(mg_genbarker[[4]])
ff_b <- colMeans(mg_barker[[3]])
ff_gb <- colMeans(mg_genbarker[[3]])
fc_b <- colMeans(mg_barker[[2]])
fc_gb <- colMeans(mg_genbarker[[2]])

pdf(file = "sim2_x1.pdf", height = 4, width = 5)
par(mar = c(5.1, 5, 2.1, 2.1))
plot(acc_b, -1/log(ff_b), type = "l", main = expression('x'[1] - bar(x)), ylab = "convergence time", xlab = "acceptance rate", cex.main = 1.25, cex.lab = 1.25, cex.axis = 1.25, xlim = c(0.1, 0.28), ylim = range(-1/log(ff_b), -1/log(ff_gb)))
lines(acc_gb, -1/log(ff_gb), type = "l", lty = 2)
legend("topright", legend = c("r = 1", "r = 2"), lty = 1:2)
dev.off()

pdf(file = "sim2_mean.pdf", height = 4, width = 5)
par(mar = c(5.1, 5, 2.1, 2.1))
plot(acc_b, -1/log(fc_b), type = "l", main = expression(bar(x)), ylab = "convergence time", xlab = "acceptance rate", cex.main = 1.25, cex.lab = 1.25, cex.axis = 1.25, xlim = c(0.1, 0.28), ylim = range(-1/log(fc_b), -1/log(fc_gb)))
lines(acc_gb, -1/log(fc_gb), type = "l", lty = 2)
legend("topright", legend = c("r = 1", "r = 2"), lty = 1:2)
dev.off()