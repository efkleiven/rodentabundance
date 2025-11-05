load("bestandsindekser/script/sim_output_base.RData")

# true values
 gam = 1
 omega = 0.8
 lambda = 5
 p = 0.3

 mod_base[1,1:20]
par(mfrow=c(2,2))
boxplot(mod_base[,"gamma"], ylim=c(0.4,1.5),  main="gamma")
abline(h=gam, col=2, lwd=2)
#abline(h=gam/2, lty=2, col=2, lwd=1)

boxplot(mod_base[, "omega"], ylim=c(0.65,0.9), main="omega")
abline(h=omega, col=2, lwd=2)

boxplot(mod_base[, "lam"], ylim=c(3,14), main="lambda")
abline(h=lambda, col=2, lwd=2)

boxplot(mod_base[, "rho"], ylim=c(0.1,0.40), main="rho")
abline(h=p, col=2, lwd=2)

