# Homework 1
library(tidyverse)
library(here)
library(stats4)

# Load Data ---------------------------------------------------------------
data <- read_delim(here("homework1", "HWK1_data.txt"))

catch <- data$Catch
obs <- c(1096, 178, 980)
cv <- c(0.35, 0.5, 0.48)

# Log Likelihood ----------------------------------------------------------

negloglik <- function(logK, logr, logtau){
  # extract parameters
  K <- exp(logK); r <- exp(logr); tau <- exp(logtau); N <- rep(0, length(catch))
  
  # generate biomass series
  N[1] <- K
  for(i in 2:length(N)){
    N[i] = max((N[i-1] + r*N[i-1]*(1 - (N[i-1]/K)^2.39) - catch[i-1]), 0.01)
  }
  
  # pull out specific abundance predictions for 1988, 1993, 2003
  preds <- N[c(67, 72, 82)]
  
  # calculate sigma for each year with abundance data
  sig <- sqrt(cv^2 + tau^2)
  
  # calculate likelihood
  
  nLL <- sum(-1*dlnorm(obs, log(preds), sig, log = T))
  
  return(nLL)
}


# Fit Model ---------------------------------------------------------------

# given the abundance estimates increase from 1993 to 2003, I'll use
# a starting value for r > 1

starts <- list(logK = log(2000), logr = log(1.3), logtau = log(2))
mle_output <- mle(negloglik, start = starts)

summary(mle_output)



# Plot Results ------------------------------------------------------------

# extract parameters from mle output
K <- exp(mle_output@coef[1])
r <- exp(mle_output@coef[2])
tau <- exp(mle_output@coef[3])

# create abundance time series
pop <- rep(0, length(catch))
pop[1] <- K

for(i in 2:length(pop)){
  pop[i] = pop[i-1] + r*pop[i-1]*(1 - (pop[i-1]/K)^2.39) - catch[i-1]
}

data <- cbind(data, pop)

# plot data
ggplot(data = data) +
  geom_line(aes(x = Year, y = pop)) +
  geom_point(aes(x = Year, y = Abundance)) +
  ylab("Population Size")

# It seems unlikely that the Fin whale population would bounce around this much,
# given how long lived the species is. This makes me think that mle found a local
# minimum negative log likelihood, not the true global minimum. This is likely because
# the abundance observations are very different. If I change the starting 
# r value to below 1, I see if that result is better than my first model fit.

starts <- list(logK = log(1000), logr = log(0.9), logtau = log(1))
mle_output2 <- mle(negloglik, start = starts)
summary(mle_output2)

if(logLik(mle_output2) < logLik(mle_output)){
  print("yes")
}

# Here is a plot of the population trajectory using the better parameter estimates

# extract parameters
K <- exp(mle_output2@coef[1])
r <- exp(mle_output2@coef[2])
tau <- exp(mle_output2@coef[3])

# create abundance time series
pop2 <- rep(0, length(catch))
pop2[1] <- K

for(i in 2:length(pop2)){
  pop2[i] = pop2[i-1] + r*pop2[i-1]*(1 - (pop2[i-1]/K)^2.39) - catch[i-1]
}

data <- cbind(data, pop2)


# plot data
ggplot(data = data) +
  geom_line(aes(x = Year, y = pop2)) +
  geom_point(aes(x = Year, y = Abundance)) +
  ylab("Population Size")

# The negative log-likelihood is lower and the population trajectory appears more
# consistent with Fin whale biology and exploitation history of whales.
# (Question: Are the differences in abundance estimates related to misidentification
# by observers?)


# Likelihood Profile ------------------------------------------------------
best_LL <- -1*as.numeric(logLik(mle_output2))

# since MSYR is approximately 0.705r, first will find the likelihood for different
# r values, then transform 

r_vals <- seq(from = 0.1, to = 0.7, by = 0.005)
negloglikes <- rep(0, length(r_vals))

for(i in 1:length(r_vals)){
  starts <- list(logK = log(2000), logtau = log(2))
  fixed <- list(logr = log(r_vals[i]))
  
  output <- mle(negloglik, start = starts, fixed = fixed)
  negloglikes[i] <- -as.numeric(logLik(output)) - best_LL
  
}

MSYR_vals <- 0.705*r_vals

likelihood <- as.data.frame(cbind(MSYR_vals, negloglikes))

ggplot() + geom_line(data = likelihood, aes(x = MSYR_vals, y = negloglikes)) +
  xlab("MSYR") + ylab("Negative Log Likelihood")

# The MSYR corresponding to the minimum negative log-likelihood is ~ 0.11
# This ratio is higher than the ratio of catch to abundance estimates for 1988, 1993, and 2003,
# as well as the ratio of catch to model population size. If the population trajectory is
# reasonable, then potentially more fin whales could be harvested sustainably.
# However, the discrepancy between the population model estimates and the abundance data seems
# to suggest that better survey methods are needed to confirm this.