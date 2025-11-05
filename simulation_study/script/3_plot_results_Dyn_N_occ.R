load("simulation_study/script/sim_output_v3.RData")

# true values
 gam = 1
 omega = 0.8
 lambda = 5
 p = 0.3

par(mfrow=c(2,2))
boxplot(exp(mod[, "gam"]), ylim=c(0,3),  main="gam")
abline(h=gam, col=2, lwd=2)
#abline(h=gam/2, lty=2, col=2, lwd=1)

boxplot(mod[, "omega"], ylim=c(0.65,0.9), main="omega")
abline(h=omega, col=2, lwd=2)

boxplot(mod[, "lam"], ylim=c(3,14), main="lambda")
abline(h=lambda, col=2, lwd=2)

boxplot(mod[, "rho"], ylim=c(0.1,0.40), main="rho")
abline(h=p, col=2, lwd=2)

