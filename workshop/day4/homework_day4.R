# Homework Day 4

data <- read.csv(here("day4", "Homework Day4.csv"), header = FALSE)
colnames(data) <- c("year", "catch", "catch_rate")


# Population Model Function -----------------------------------------------------------------
PopMod <- function(catch, r, K, phi){
  # pull out parameters
  C <- catch; r <- r; k <- K; phi <- phi
  
  # calculate derived parameters
  n_0 <- phi*K # initial population
  n <- length(catch)
  
  # vector for predictions
  pop <- rep(NA, length(C))
  
  for(i in 1:length(C)){
    if(i == 1){
      pop[i] = n_0
      
    }else if(i > 1){
      pop[i] = pop[i-1] + r*pop[i-1]*(1 - (pop[i-1]/k)) - C[i-1]
      
    }
  }
  
  return(pop)
}

# Negative log-likelihood function -----------------------------------------------------------
NegLogLike <- function(pop, catch_rate){
  
  # extract observations and predictions
  c_r <- catch_rate; pop <- pop
  
  # calculate catchability coefficient
  ratio <- c_r/pop
  for(i in ratio){
    if(i < 0.0001){
      ratio[i] <- 0.0002 
    }
  }
  
  q <- exp(mean(ratio)) 
  
  # calculate sigma
  sig <- sqrt(mean((log(c_r) - log(q*pop))^2)) 
  
  catch_rate_preds <- log(q*pop)
  

  likelihood <- -sum(dlnorm(x = c_r, meanlog = catch_rate_preds, sdlog = sig, log = TRUE))
  return(likelihood)
  
}

# Threshold creation function --------------------------------------------------------
# need to have the population prediction and the likelihood in one function call for optim
# optim is needed to set up the threshold
threshold_func <- function(log_r, log_K, log_phi){
  r <- exp(log_r); K <- exp(log_K); phi <- exp(log_phi)
  
  pop <- PopMod(catch, r, K, phi)
  
  loglike <- NegLogLike(pop, catch_rate)
  
  return(loglike)
}

starts <- list(log_r = log(0.3), log_K = log(6000), log_phi = log(0.8))
output <- mle(threshold_func, start = starts)
# SIR Algorithm -----------------------------------------------------------------------
DoSir <- function(Nout=1000,data)
{
  
  # Read in the basic data
  data <- data
  catch <- data[,2]
  catch_rate <- data[,3]
  
  
  # Storage for the parameters and the final depletion
  Vals <- matrix(0,ncol=4,nrow=Nout)
  colnames(Vals) <- c("K", "r", "phi", "End_Deplet")
  
  # Storage for the total likelihood encountered in the first stage sampling
  # and the number of draws
  AveLike <- 0
  Ntest <- 0
  
  # Reset parameters for SIR
  Threshold <- exp(0) # Markus' estimate, still need to get threshold to work   
  Cumu <- 0
  Ndone <- 0
  while (Ndone < Nout)
  {
    # Generate from priors
    K <- runif(1, min = 5000, max = 8000)
    r <- runif(1, min = 0.15, max = 0.45)
    phi <- runif(1, min = 0.8, max = 1.2)
    
    # Predict the population biomass
    pop <- PopMod(catch, r, K, phi)
    
    
    # Compute the negative log-likelihood
    neg_log <- NegLogLike(pop, catch_rate)
    TheLike <- exp(neg_log-77.48) #this correction is Markus', haven't figured out my own
    AveLike <- AveLike + TheLike
    Ntest <- Ntest + 1
    
    # Determine if a parameter vector is to be saved
    Cumu <- Cumu + TheLike
    while (Cumu > Threshold & Ndone < Nout)
    {
      Ndone <- Ndone + 1
      Cumu <- Cumu - Threshold
      # save posterior info
      Vals[Ndone, 1] <- K
      Vals[Ndone, 2] <- r
      Vals[Ndone, 3] <- phi
      Vals[Ndone, 4] <- pop[27]/K
    }
  }
  
  return(Vals)
}