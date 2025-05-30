library(dplyr)
library(Bessel)

set.seed(123)

# Simulation parameters
num_teams <- 20
num_matches <- num_teams * (num_teams - 1)

# True parameters with completely uninformative priors (as per paper)
p <- runif(1, 0, 1) # Uniform(0,1) for zero-inflation
mu <- rnorm(1, 0, sqrt(1e4)) # N(0,10^4)
H <- rnorm(1, 0, sqrt(1e4)) # N(0,10^4)
A <- rnorm(num_teams, 0, sqrt(1e4)) # N(0,10^4) for all teams
D <- rnorm(num_teams, 0, sqrt(1e4)) # N(0,10^4) for all teams

# Generate match combinations with numerical stability safeguards
matches <- expand.grid(
  HomeTeam = 1:num_teams,
  AwayTeam = 1:num_teams
) %>% 
  filter(HomeTeam != AwayTeam) %>%
  mutate(
    # Calculate log-lambdas first (with stability bounds)
    log_lambda1 = mu + H + A[HomeTeam] - D[AwayTeam],
    log_lambda2 = mu + A[AwayTeam] - D[HomeTeam],
    
    # Apply reasonable bounds (-20,20) to prevent numerical overflow
    log_lambda1 = pmax(pmin(log_lambda1, 20), -20),
    log_lambda2 = pmax(pmin(log_lambda2, 20), -20),
    
    # Convert to expected goals
    lambda1 = exp(log_lambda1),
    lambda2 = exp(log_lambda2),
    
    # Apply realistic football constraints (optional but recommended)
    lambda1 = pmin(pmax(lambda1, 0.05), 10), # Min 0.05, max 10 goals expected
    lambda2 = pmin(pmax(lambda2, 0.05), 10),
    
    # Vectorized Skellam PMF at zero with numerical safeguards
    fPD0 = {
      s <- 2 * sqrt(lambda1 * lambda2)
      case_when(
        s > 700 ~ 0,  # Handle extremely large values
        is.na(s) ~ 0,  # Handle NA cases
        TRUE ~ exp(-(lambda1 + lambda2)) * besselI(s, nu = 0)
      )
    },
    
    # Zero-inflated draw probability
    is_draw_prob = p + (1 - p) * ifelse(is.na(fPD0), 0, fPD0),
    is_draw_prob = pmin(pmax(is_draw_prob, 0), 1), # Ensure valid probability
    
    # Simulate match outcome
    is_draw = rbinom(n(), 1, is_draw_prob),
    HomeGoals = ifelse(is_draw == 1, 0, rpois(n(), lambda1)),
    AwayGoals = ifelse(is_draw == 1, 0, rpois(n(), lambda2)),
    GoalDiff = HomeGoals - AwayGoals
  )

# Verify realistic outcomes despite uninformative priors
outcome_stats <- matches %>%
  summarise(
    AvgHomeGoals = mean(HomeGoals),
    AvgAwayGoals = mean(AwayGoals),
    DrawRate = mean(GoalDiff == 0),
    HomeWinRate = mean(GoalDiff > 0),
    AwayWinRate = mean(GoalDiff < 0),
    MaxGoalDiff = max(abs(GoalDiff)),
    PropLargeWins = mean(abs(GoalDiff) > 3)
  )

print(outcome_stats)

# View score distribution
cat("\nMost common scorelines:\n")
print(sort(table(paste(matches$HomeGoals, "-", matches$AwayGoals)), decreasing = TRUE)[1:10])

# Save results
write.csv(matches, "simulated_data.csv", row.names = FALSE)
