load("simulation_study/output/sim_output_seasonalgains.RData")

# true values
 gam = 1
 omega = 0.8
 lambda = 5
 p = 0.3

par(mfrow=c(1,2))
boxplot(exp(mod[, "gam[1]"]), ylim=c(-1,1),  main="gam")
abline(h=0, col=2, lwd=2)

boxplot(exp(mod[, "gam[2]"]), ylim=c(0,3),  main="gam")
abline(h=1, col=2, lwd=2)

boxplot(mod[, "omega"], ylim=c(0.65,0.9), main="omega")
abline(h=omega, col=2, lwd=2)

boxplot(mod[, "lam"], ylim=c(3,14), main="lambda")
abline(h=lambda, col=2, lwd=2)

boxplot(mod[, "rho"], ylim=c(0.1,0.40), main="rho")
abline(h=p, col=2, lwd=2)

# violine plot
library(ggplot2)
library(reshape2)

# Example: assuming 'mod' is your posterior sample matrix or data frame
# with columns like: "gam[1]", "gam[2]", "omega", "lam", "rho"

# Prepare parameters
params <- data.frame(
  gam1   = exp(mod[, "gam[1]"]),
  gam2   = exp(mod[, "gam[2]"]),
  omega  = mod[, "omega"],
  lambda = mod[, "lam"],
  rho    = mod[, "rho"]
)

# Reshape to long format
params_long <- melt(params, variable.name = "Parameter", value.name = "Value")

# Ensure consistent order for plotting
params_long$Parameter <- factor(params_long$Parameter,
                                levels = c("gam1", "gam2", "omega", "lambda", "rho"))

# Reference values
refs <- data.frame(
  Parameter = factor(c("gam1", "gam2", "omega", "lambda", "rho"),
                     levels = levels(params_long$Parameter)),
  ref_value = c(0, 1, 0.8, 5, 0.3)
)

# Create the plot

p <- ggplot(params_long, aes(x = Parameter, y = Value)) +
  geom_violin(fill = "lightblue", color = "grey30", alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.8) +
  geom_segment(data = refs,
               aes(x = as.numeric(Parameter) - 0.2,
                   xend = as.numeric(Parameter) + 0.2,
                   y = ref_value, yend = ref_value),
               color = "red", linewidth = 1.2) +
  facet_wrap(~Parameter, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior distributions of parameters",
       y = "Value", x = "")

# Display the plot
print(p)

# Save as PNG
ggsave("simulation_study/plot/posterior_violin_seasonal.png", plot = p, width = 8, height = 6, dpi = 300)
