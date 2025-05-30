# Load required libraries
library(dplyr)
library(coda)
library(MASS)
library(gridExtra)
library(ggplot2)
# Load simulated data
matches <- read.csv("simulated_data.csv")

# Set up constants
num_teams <- 20
n_matches <- nrow(matches)
goal_diff <- matches$GoalDiff
home_team <- matches$HomeTeam
away_team <- matches$AwayTeam

# MCMC settings
n_iter <- 5000
burn_in <- 1000

# Initial values
mu <- 0
H <- 0.1
A <- rep(0, num_teams)
D <- rep(0, num_teams)
p <- 0.3

# Store results
samples <- matrix(NA, nrow = n_iter, ncol = 2 + 2*(num_teams-1) + 1)
colnames(samples) <- c("mu", "H", paste0("A_", 2:num_teams), paste0("D_", 2:num_teams), "p")

# Prior variance
sigma2 <- 1e4

# Set seed for reproducibility
set.seed(123)

# MCMC sampling loop
for (iter in 1:n_iter) {
  
  # Compute λ1 and λ2
  lambda1 <- exp(pmin(pmax(mu + H + A[home_team] - D[away_team], -20), 20))
  lambda2 <- exp(pmin(pmax(mu + A[away_team] - D[home_team], -20), 20))
  
  # Skellam zero probability
  s <- 2 * sqrt(lambda1 * lambda2)
  f_PD0 <- exp(-(lambda1 + lambda2)) * besselI(s, nu = 0, expon.scaled = TRUE) * exp(s)
  f_PD0[is.na(f_PD0) | is.infinite(f_PD0)] <- 1e-10
  
  # Zero-inflation probabilities for draws
  draw_probs <- rep(NA, n_matches)
  draw_ids <- which(goal_diff == 0)
  draw_probs[draw_ids] <- p / (p + (1-p)*f_PD0[draw_ids])
  draw_probs[draw_ids] <- pmin(pmax(draw_probs[draw_ids], 1e-8), 1-1e-8)
  
  # Sample δ_i
  delta <- rep(0, n_matches)
  delta[draw_ids] <- rbinom(length(draw_ids), 1, draw_probs[draw_ids])
  
  # Generate latent w1, w2
  w1 <- numeric(n_matches)
  w2 <- numeric(n_matches)
  
  for (i in 1:n_matches) {
    z <- goal_diff[i]
    if (delta[i] == 1) {
      w1[i] <- 0
      w2[i] <- 0
    } else {
      if (z > 0) {
        w2[i] <- rpois(1, lambda2[i])
        w1[i] <- w2[i] + z
      } else if (z < 0) {
        w1[i] <- rpois(1, lambda1[i])
        w2[i] <- w1[i] - z
      } else {
        w1[i] <- rpois(1, lambda1[i])
        w2[i] <- w1[i]
      }
    }
  }
  
  # Update p (Beta posterior)
  a_post <- 1 + sum(delta)
  b_post <- 1 + n_matches - sum(delta)
  p <- rbeta(1, a_post, b_post)
  
  # Design matrix X
  X <- matrix(0, nrow = n_matches, ncol = 2 + 2*(num_teams-1))
  X[,1] <- 1
  X[,2] <- 1
  
  for (i in 1:n_matches) {
    if (home_team[i] != 1) X[i, 2 + home_team[i]-1] <- 1
    if (away_team[i] != 1) X[i, 2 + away_team[i]-1] <- 1
    if (away_team[i] != 1) X[i, 2+(num_teams-1) + away_team[i]-1] <- -1
    if (home_team[i] != 1) X[i, 2+(num_teams-1) + home_team[i]-1] <- -1
  }
  
  # Outcome vector Y
  Y <- c(w1, w2)
  X_full <- rbind(X, X)
  
  # Normal approximation for parameter update
  V_beta <- solve(t(X_full) %*% X_full / sigma2 + diag(1/sigma2, ncol(X_full)))
  logY <- log(Y + 1e-3)  # avoid log(0)
  m_beta <- V_beta %*% (t(X_full) %*% logY / sigma2)
  
  beta <- MASS::mvrnorm(1, m_beta, V_beta)
  
  mu <- beta[1]
  H  <- beta[2]
  A[-1] <- beta[3:(2+num_teams-1)]
  D[-1] <- beta[(2+num_teams):(2+2*(num_teams-1))]
  A[1] <- -sum(A[-1])
  D[1] <- -sum(D[-1])
  
  # Save sample
  samples[iter,] <- c(mu, H, A[-1], D[-1], p)
  
  if (iter %% 500 == 0) cat("Iteration:", iter, "\n")
}

# Post-processing
samples_mcmc <- mcmc(samples[burn_in:n_iter,])

# Traceplots
plot(samples_mcmc)

# Posterior summary
summary(samples_mcmc)
# Extract key parameters for plotting
params_to_plot <- c("mu", "H", "A_2", "D_2", "p")

# Convert to a data frame
samples_df <- as.data.frame(samples_mcmc)

# Traceplots using base R
par(mfrow=c(3,2))  # arrange 3 rows, 2 cols
for (param in params_to_plot) {
  plot(samples_df[[param]], type="l", main=paste("Traceplot of", param),
       xlab="Iteration", ylab=param, col="blue")
}
par(mfrow=c(1,1))
par(mfrow=c(3,2))
for (param in params_to_plot) {
  plot(density(samples_df[[param]]), main=paste("Density of", param),
       xlab=param, col="darkgreen", lwd=2)
}
par(mfrow=c(1,1))
par(mfrow=c(3,2))
for (param in params_to_plot) {
  acf(samples_df[[param]], main=paste("ACF of", param), col="purple")
}
par(mfrow=c(1,1))
# Calculate ESS for all parameters
ess_values <- effectiveSize(samples_mcmc)


#########################
# Calculate posterior means
posterior_means <- colMeans(samples_mcmc)

# Calculate 95% credible intervals (2.5% and 97.5% quantiles)
credible_intervals <- apply(samples_mcmc, 2, quantile, probs = c(0.025, 0.975))

# Combine means and credible intervals into a clean data frame
results_df <- data.frame(
  Parameter = colnames(samples_mcmc),
  Mean = posterior_means,
  `2.5%` = credible_intervals[1,],
  `97.5%` = credible_intervals[2,]
)

# Print results for key parameters first
key_params <- c("mu", "H", "A_2", "D_2", "p")
print(results_df[results_df$Parameter %in% key_params,])

# If you want to see all parameters
print(results_df)
