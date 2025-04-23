# Load libraries
library(ggplot2)
library(dplyr)

# Set seed for reproducibility
set.seed(42)

# Generate samples
n_samples <- 100000

# 1. True Beta(1,1) ~ Uniform(0,1)
beta_samples <- runif(n_samples)

# 2. "Soft" Beta(1,1): inv_logit(normal(0, 1))
logit <- function(x) log(x / (1 - x))
inv_logit <- function(x) 1 / (1 + exp(-x))
soft_samples <- inv_logit(rnorm(n_samples, mean = 0, sd = 1))

# Create a dataframe for ggplot
df <- data.frame(
  value = c(beta_samples, soft_samples),
  distribution = rep(c("Beta(1,1) - Uniform", "Soft Beta(1,1)"), each = n_samples)
)

# Plot densities
ggplot(df, aes(x = value, fill = distribution)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = c("skyblue", "tomato")) +
  labs(title = "Beta(1,1) vs. Soft Beta(1,1) (inv_logit(N(0,1)))",
       x = "rho",
       y = "Density") +
  theme_minimal()