# Set parameters
set.seed(1)  # For reproducibility
n_families <- 10000  # Number of families
n_parents <- 2       # twin parents
beta <- 0.3          # Coefficient for education transmission
alpha <- 0.2         # Coefficient for social factors
gamma <- 0.8
phi3_sq <- 0.3          # Contribution of Unique experiences
phi2_sq <- 0.2          # Contribution of Common family influence
phi1_sq <- 1-phi3_sq-phi2_sq # Contribution of Heritability

ntime <- 100 # time of reproducibility

# Store results for 100 runs
results <- data.frame(beta0 = numeric(ntime), beta0_lower = numeric(ntime), beta0_upper = numeric(ntime),
                      beta1 = numeric(ntime), beta1_lower = numeric(ntime), beta1_upper = numeric(ntime),
                      beta2 = numeric(ntime), beta2_lower = numeric(ntime), beta2_upper = numeric(ntime))
                      

# Run the simulation
for (i in 1:ntime) {
  
  # Simulate data
  n <- n_families * n_parents
  family_ids <- rep(1:n_families, each = n_parents)
  Z1 <- rnorm(n)  # Social factors
  Z2 <- rnorm(n)  # Social factors
  H1 <- rnorm(n_families)  # Heritability
  H2 <- rnorm(n)          # Heritability
  C <- rnorm(n_families)  # Common family influence
  E <- rnorm(n)           # Unique experiences
  e1 <- rnorm(n)                 # Individual-level noise
  e2 <- rnorm(n)                 # Individual-level noise
  
  X <- sqrt(phi1_sq)*H1[family_ids] + sqrt(phi2_sq)*C[family_ids] + sqrt(phi3_sq)*E  # Parent education outcome
  
  Y1 <- beta * (sqrt(phi1_sq)*H1[family_ids] +  sqrt(phi3_sq)*E) +
                alpha * Z1 + gamma*sqrt(phi2_sq)*C[family_ids] + e1  # Biological child outcome
  Y2 <- beta * (sqrt(phi1_sq)*H2 + sqrt(phi3_sq)*E) +
                alpha * Z2 + gamma*sqrt(phi2_sq)*C[family_ids] + e2 # Adopted child outcome
  
  # OLS model
  biog <- lm(Y1 ~ X + Z1)
  beta0 <- coef(biog)["X"]
  se0 <- summary(biog)$coefficients["X", "Std. Error"]  # Standard error for beta2
  ci0_lower <- beta0 - 1.96 * se0  # Lower bound of CI
  ci0_upper <- beta0 + 1.96 * se0  # Upper bound of CI
 
  # Twin-parent model
  X_centered <- X - ave(X, family_ids)
  Z1_centered <- Z1 - ave(Z1, family_ids)
  Y1_centered <- Y1 - ave(Y1, family_ids)
  
  twin <- lm(Y1_centered ~ X_centered + Z1_centered)
  beta1 <- coef(twin)["X_centered"]
  se1 <- summary(twin)$coefficients["X_centered", "Std. Error"]  # Standard error for beta1
  ci1_lower <- beta1 - 1.96 * se1  # Lower bound of CI
  ci1_upper <- beta1 + 1.96 * se1  # Upper bound of CI
  
  # Adoption model
  adopt <- lm(Y2 ~ X + Z2)
  beta2 <- coef(adopt)["X"]
  se2 <- summary(adopt)$coefficients["X", "Std. Error"]  # Standard error for beta2
  ci2_lower <- beta2 - 1.96 * se2  # Lower bound of CI
  ci2_upper <- beta2 + 1.96 * se2  # Upper bound of CI
  
  # Store results

  results$beta0[i] <- beta0
  results$beta0_lower[i] <- ci0_lower
  results$beta0_upper[i] <- ci0_upper
  
  results$beta1[i] <- beta1
  results$beta1_lower[i] <- ci1_lower
  results$beta1_upper[i] <- ci1_upper
  
  results$beta2[i] <- beta2
  results$beta2_lower[i] <- ci2_lower
  results$beta2_upper[i] <- ci2_upper
}

# Load required libraries
library(ggplot2)
library(patchwork)

# Beta estimates plot
beta_plot <- ggplot(results, aes(x = 1:nrow(results))) +
  geom_line(aes(y = beta0, color = "Beta0"), linewidth = 1) + 
  geom_ribbon(aes(ymin = beta0_lower, ymax = beta0_upper, fill = "Beta0 CI"), alpha = 0.2) +
  geom_line(aes(y = beta1, color = "Beta1"), linewidth = 1) + 
  geom_ribbon(aes(ymin = beta1_lower, ymax = beta1_upper, fill = "Beta1 CI"), alpha = 0.2) +
  geom_line(aes(y = beta2, color = "Beta2"), linewidth = 1) + 
  geom_ribbon(aes(ymin = beta2_lower, ymax = beta2_upper, fill = "Beta2 CI"), alpha = 0.2) +
  labs(title = "A. Estimates of Beta Over Time", x = "", y = "") + 
  geom_hline(yintercept = beta, color = "black", linetype = "dashed") + 
  geom_hline(yintercept = beta*phi3_sq + gamma*phi2_sq, color = "black", linetype = "dashed") +
  geom_hline(yintercept = beta*(phi1_sq + phi3_sq) + gamma*phi2_sq, color = "black", linetype = "dashed") +
  theme_minimal()

# Calculate phi1_sq_hat
results$phi1_sq_hat <- (results$beta0 - results$beta2) / results$beta1

# Density plot for phi1_sq_hat
density_plot <- ggplot(results, aes(x = phi1_sq_hat)) +
  geom_density(fill = "blue", alpha = 0.5) +
  geom_vline(xintercept = phi1_sq, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "B. Smooth Density of Heritability Contribution", x = "", y = "") +
  theme_minimal()

# Combine the plots
simulation <- beta_plot + density_plot

ggsave("figs/simulation.svg", plot = simulation, width=7, height=4)
