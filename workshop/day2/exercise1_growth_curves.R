# Day 2 exercises
library(tidyverse)
library(here)
library(stats4)

data1 <- read.csv(here("day2", "ex1.dat"), sep = " ")

data1 <- data1 %>%
  select(Age, X)

Age <- data1[,1]
Length <- data1[,2]

minusLL1 <- function(logLinf, loga50, logdelta, logSigma){
  # Extract params
  Linf <- exp(logLinf); a50 <- exp(loga50); delta <- exp(logdelta); Sigma <- exp(logSigma)
  
  # Make model predictions
  Pred <- Linf*(1 + exp(-log(19)*(Age - a50)/delta))^(-1)
  
  # Compute the negative log-likelihood
  NegLogL <- (-1)*sum(dnorm(Length, Pred, Sigma, log = TRUE))
  
  return(NegLogL)
}

# Plot data to find starting values
plot(Age, Length)

# Fit model 

starts <- list(logLinf = log(100), loga50 = log(12), logdelta = log(10), logSigma = 1)
mleOutput <- mle(minusLL1, start = starts)
print(summary(mleOutput))

# Extract parameter estimates to check fit
Linf <- exp(coef(mleOutput)[1])
a50 <- exp(coef(mleOutput)[2])
delta <- exp(coef(mleOutput)[3])
Sigma <- exp(coef(mleOutput)[4])

best_LL <- -as.numeric(logLik(mleOutput))

# Calculate AICc
AIC1 <- 2*best_LL + (2*length(Age)*4)/(length(Age) - 4 - 1)

# Plot fit
Pred <- Linf*(1 + exp(-log(19)*(Age - a50)/delta))^(-1)

results <- tibble(Age, Length, Pred)
ggplot(data = results) +
  geom_point(aes(Age, Length)) +
  geom_line(aes(Age, Pred), colour = "blue")

# Second function ------------------------------------------------------------
minusLL2 <- function(logLinf, a0, logkappa, logSigma){
  # Extract params
  Linf <- exp(logLinf); kappa <- exp(logkappa); Sigma <- exp(logSigma)
  
  # Make model predictions
  Pred <- Linf*(1 - exp(-kappa*(Age - a0)))
  
  # Compute the negative log-likelihood
  NegLogL <- (-1)*sum(dnorm(Length, Pred, Sigma, log = TRUE))
  
  return(NegLogL)
}

# Plot data to find starting values
kappa <- 0.1
a0 <- 1
Linf <- 110

Pred_test <- Linf*(1 - exp(-kappa*(Age - a0)))
results <- tibble(Age, Length, Pred_test)
ggplot(data = results) +
  geom_point(aes(Age, Length)) +
  geom_line(aes(Age, Pred_test), colour = "blue")

# Fit model 
starts <- list(logLinf = log(100), a0 = 0, logkappa = log(0.2), logSigma = log(50))
mleOutput2 <- mle(minusLL2, start = starts, method = "Nelder-Mead")
print(summary(mleOutput2))

# Extract parameter estimates to check fit
Linf_2 <- exp(coef(mleOutput2)[1])
a0 <- coef(mleOutput2)[2]
kappa <- exp(coef(mleOutput2)[3])
Sigma_2 <- exp(coef(mleOutput2)[4])

best_LL_2 <- -as.numeric(logLik(mleOutput2))

# Calculate AICc
AIC2 <- 2*best_LL_2 + (2*length(Age)*4)/(length(Age) - 4 - 1)

# Plot fit
Pred2 <- Linf_2*(1 - exp(-kappa*(Age - a0)))

results2 <- tibble(Age, Length, Pred2)
ggplot(data = results2) +
  geom_point(aes(Age, Length)) +
  geom_line(aes(Age, Pred2), colour = "blue")

# Calculate Model Weights -----------------------------------------------------



# Plot both lines ----------------------------------------------------------
total_results <- tibble(Age, Length, Pred, Pred2)
ggplot(data = total_results) +
  geom_point(aes(Age, Length)) +
  geom_line(aes(Age, Pred), colour = "blue") +
  geom_line(aes(Age, Pred2), colour = "red")