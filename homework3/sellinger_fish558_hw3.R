# Homework 3
library(tidyverse)
library(here)
library(gridExtra)

##################### PART A ##############################################
# load data ---------------------------------------------------------------
data <- read_csv(here("homework3", "bowhead_whale_catch.csv"))

obs_years <- c(1978, 1980, 1981, 1982, 1983, 1985, 1986, 1987, 1988, 1993)
obs <- c(4820, 3900, 4389, 6572, 6268, 5132, 7251, 5151, 6609, 7778)
cv <- c(0.273, 0.314, 0.253, 0.311, 0.321, 0.269, 0.186, 0.298, 0.113, 0.071)

obs_data <- data.frame(obs_years, obs, cv)

# initiate parameters -----------------------------------------------------
s0 <- 0.9
s1 <- 0.95
K1plus <- 15000
fmax <- 0.29

Nyears <- nrow(data) # number of catch observations

params <- c(s0, s1, K1plus, fmax) # put parameters in vector to pass to pop model

# define population model -------------------------------------------------
pop_model <- function(params, data, Nyears){
  
  # extract catch data 
  Catch <- data$Catch
  Nages <- 14 # number of age classes
  
  # extract parameters
  s0 <- params[1] 
  s1 <- params[2] 
  K1plus <- params[3]
  fmax <- params[4]
  
  # calculate derived parameters
  f0 <- (1-s1) / (s0*s1^(12))
  N_initial <- K1plus
  z <- 2.39
  
  # create transition matrix for population projection
  transition <- matrix(0, nrow = 14, ncol = 14) 
  transition[1,14] <- f0  
  transition[2,1] <- s0
  
  for(i in 2:13){
    transition[i+1,i] <- s1
  }
  
  transition[14,14] <- s1
  
  # before projecting the population need to determine the # of animals in each age class
  ee <- eigen(transition) #eigen 
  age_struct <- as.numeric(ee$vectors[,1]/sum(ee$vectors[,1]))#population by age class
  age_struct <- as.numeric(ee$vectors[,1])
  
  # create population matrix and catch matrix
  N <- matrix(NA, nrow = Nages, ncol = Nyears+1) # rows are ages, columns years (+1 because we need 2003 estimates for part B)
  C <- matrix(NA, nrow = 13, ncol = (Nyears+1)) # nrow = 13 because age 0 whales aren't caught
  
  # fill in age and catch data for 1848
  N[,1] <- age_struct*K1plus/sum(age_struct[-1]) # need to take out age 0 whales from the total, since they don't contribute to carrying capacity
  C[,1]  <- Catch[1]*N[2:Nages,1]/sum(N[2:Nages,1]) 
  
  for(t in 2:(Nyears+1)){
    Ntmp <- c(N[1, t-1],N[2:Nages, t-1] - C[,t-1]) # make a temp vector with population # by age class after yearly catch
    Ntmp <- mapply(max, Ntmp, 0.0001) #to make sure not negative
    
    N[,t] <- transition %*% Ntmp # calculate the number of whales that survive
    N[,t] <- mapply(max, N[, t], 0.0001) # remove any negatives
    Ntot <- sum(N[2:Nages,t]) # find N1+ (# of whales greater than age 0)
    
    if(!is.na(Catch[t])){
      C[,t] <- (N[2:Nages,t]/Ntot) * Catch[t] #Calculate catches per group for years where there's catch
    } else {C[,t] <- 0} #For years with no catches
    
    N[1,t] <- N[Nages,t]*(f0 + (fmax - f0)*(1- (Ntot/K1plus)^z)) # calculate # of age 0 whales
  }
  
  return(N)
}
  
# define log likelihood  ---------------------------------------------------------
minusLL <- function(params, preds, observed){
  
  #likelihood estimates
  sigma <- observed$cv
  N_observed <- observed$obs
  
  L <- -1*sum(-(log(preds)-log(N_observed))^2 / (2*sigma^2))
  
  return(L)
}

# Question A2 -------------------------------------------------------------
predictions <- pop_model(params, data, Nyears)

colnames(predictions) <- 1848:2003

# calculate population in 2002
pop_2002 <- sum(predictions[,"2002"])
print(pop_2002)

# calculate the negative log likelihood
N_predicted <- predictions[-1, c("1978", "1980", "1981", "1982", "1983", "1985", "1986", "1987", "1988", "1993")]
N_predicted <- colSums(N_predicted)

negloglike <- minusLL(params, N_predicted, obs_data)
print(negloglike)


# Question A3 -------------------------------------------------------------
DoSIR1 <- function(Nout, observed){
  
  # Create matrix for parameter estimates
  Vals <- matrix(0, ncol=4, nrow=Nout)
  
  # Create matrix for the total likelihood encountered in the first stage sampling
  # and the number of draws
  total_like <- 0
  Ntest <- 0
  
  # Reset parameters for SIR
  Threshold <- 1 
  Cumu <- 0
  Ndone <- 0
  while (Ndone < Nout){
    
    # Counter for Ntest
    Ntest <- Ntest + 1
    
    # Priors
    s0 <- runif(1, 0.8, 1)
    s1 <- runif(1, 0.9, 1)
    K1plus <- runif(1, 10000, 20000)
    fmax <- runif(1, 0.25, 0.33)
    pars <- c(s0, s1, K1plus, fmax)
    
    # Call the population model
    pop <- pop_model(pars, data, Nyears)
    pop_total <- sum(pop[,156])
    
    #Determine if extinct or not
    TheLike  <- ifelse(pop_total>1, 1,0)
    
    # Determine if a parameter vector is to be saved
    Cumu <- Cumu + TheLike
    while (Cumu >= Threshold & Ndone < Nout)
    {
      Ndone <- Ndone + 1
      Cumu <- Cumu - Threshold
      # save values
      Vals[Ndone,] <- c(s0, s1,K1plus,fmax)
    }
  }
  
  return(Vals)
  
}

SIR1_output <- DoSIR1(Nout=1000)

colnames(SIR1_output) <- c("s0", "s1", "K1plus", "fmax")
SIR1_output <- as.data.frame(SIR1_output)

# Question A4 -------------------------------------------------------------
DoSIR2 <- function(Nout,minusLL, popmodel, data, observed){
  
  # Data stuff
  data <- data
  Catch <- data$Catch
  Nyears<- nrow(data)
  
  # Create matrix of parameters
  Vals <- matrix(0,ncol= 4,nrow=Nout)
  
  # Create matrix for the total likelihood encountered in the first stage sampling
  # and the number of draws
  total_like <- 0
  Ntest <- 0
  
  # Reset parameters for SIR
  Threshold <- exp(0)
  Cumu <- 0
  Ndone <- 0
  while (Ndone < Nout){
    # Counter for Ntest
    Ntest <- Ntest + 1
    # Generate from priors
    S0 <- runif(1, 0.8, 1)
    S1plus <- runif(1, 0.9, 1)
    K1plus <- runif(1, 10000, 20000)
    fmax <- runif(1, 0.25, 0.33)
    pars <- c(S0, S1plus, K1plus, fmax)
    
    # Call the population model
    N_projected <- popmodel(pars, data, Nyears)
    colnames(N_projected) <- 1847:2003
    
    
    Years_with_est <- c(1978, 1980, 1981, 1982, 1983, 1985, 1986, 1987, 1988, 1993)
    N_predicted <- N_projected[, c("1978", "1980", "1981", "1982", "1983", "1985", "1986", "1987", "1988", "1993")]
    N_predicted <- colSums(N_predicted)
    
    # Compute the negative log-likelihood and hence the likelihood
    NegLogLike <- minusLL(params = pars, preds = N_predicted, obs = observed)
    TheLike <- exp(-1*NegLogLike+2.22)
    #print(TheLike)
    # Determine if a parameter vector is to be saved
    Cumu <- Cumu + TheLike
    #print(Cumu)
    while (Cumu >= Threshold & Ndone < Nout)
    {
      Ndone <- Ndone + 1
      Cumu <- Cumu - Threshold
      # need to save posterior informative here - parameters
      Vals[Ndone,] <- c(S0, S1plus,K1plus,fmax)
      print(Ndone)
    }
  }
  
  return(Vals)
  
}

A4 <- DoSIR2(Nout=200, minusLL = minusLL, popmodel = pop_model, data = data, observed = obs_data)

# saving results from the SIR2 algorithm into a csv file so I don't have to rerun it 
SIR2_results <- as.data.frame(A4)
colnames(SIR2_results) <- c("s0", "s1", "K1plus", "fmax")
write.csv(SIR2_results, "homework3/sir_output.csv")


# Question A5 -------------------------------------------------------------
vecs <- nrow(A4)

pop <- matrix(data = NA, nrow = vecs, ncol = Nyears) # population estimates for each of the parameter vectors

for(i in 1:vecs){
  for(t in 1:Nyears){
    whales <- pop_model(A4[i,], data, Nyears)[2:14,t] # remove age 0 whales
    pop[i,t]<- sum(whales) #getting 1+ year data for each year
  }
}

# saving population estimates to a csv as well
colnames(pop) <- paste("y_", 1848:2002, sep = '')
pop <- as.data.frame(pop)
write.csv(pop, "homework3/whale_sim_estimates.csv")

# calculate the median and 95% CI for each of the years
quant_data <- matrix(data = NA, nrow = Nyears, ncol = 3)


for(i in 1:Nyears){
  quant_data[i,] <- quantile(pop[,i], c(0.05, 0.5, 0.95))
}

quant_data <- as.data.frame(quant_data)
colnames(quant_data) <- c("5%", "median", "95%")
quant_data$year <- seq(1848,2002)

# create plots for SIR1 results, the post-model, pre-data distributions for our four parameters
p1 <- ggplot(data = SIR1_output, aes(x = s0)) + geom_histogram(bins = 30)
p2 <- ggplot(data = SIR1_output, aes(x = s1)) + geom_histogram(bins = 30)
p3 <- ggplot(data = SIR1_output, aes(x = K1plus)) + geom_histogram(bins = 30)
p4 <- ggplot(data = SIR1_output, aes(x = fmax)) + geom_histogram(bins = 30)

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
# compare our SIR1 results with our SIR2 results
p5 <- ggplot(data = SIR2_results, aes(x = s0)) + geom_histogram(bins = 30)
p6 <- ggplot(data = SIR2_results, aes(x = s1)) + geom_histogram(bins = 30)
p7 <- ggplot(data = SIR2_results, aes(x = K1plus)) + geom_histogram(bins = 30)
p8 <- ggplot(data = SIR2_results, aes(x = fmax)) + geom_histogram(bins = 30)

grid.arrange(p5, p6, p7, p8, ncol = 2, nrow = 2)
# adding in the priors for s1 and K1plus seemed to improve our estimates for those parameters
# to get even better estimates, we could increase the number of parameter vectors we sample

# create plot for our whales' time trajectory from 1848 to 2002

ggplot(data = quant_data) +
  geom_line(aes(x = year, y = median)) +
  geom_line(aes(x = year, y = `5%`), linetype = 2) +
  geom_line(aes(x = year, y = `95%`), linetype = 2) + 
  ylab("N+1 Bowheads") + xlab("Year")

##################### PART B ##############################################

# re-load data if necessary------------------------------------------------------------
SIR2_results <- read_csv(here("homework3", "sir_output.csv"))
SIR2_results <- SIR2_results[,-1] # loaded with extra column, deleted it

pop <- read_csv(here("homework3", "whale_sim_estimates.csv"))

# calculate states of nature ----------------------------------------------
# want to know the percentage of simulations that fall within each of our 
# states of nature
pop2002 <- pop[,155]

pop2002 %>% 
  filter(y_2002 < 7000) # 54 sims 
pop2002 %>% 
  filter(y_2002 >= 7000 & y_2002 < 8000) # 102 sims
pop2002 %>% 
  filter(y_2002 >= 8000) # 44 sims

# create harvest rules ----------------------------------------------------

catch1 <- rep(67, 20)
catch2 <- rep(134, 20)
catch3 <- rep(201, 20)

Nyears <- 20

# update pop model --------------------------------------------------------
pop_model2 <- function(params, data, Nyears){
  
  # extract catch data 
  Catch <- data
  Nages <- 14 # number of age classes
  
  # extract parameters
  s0 <- params[1] 
  s1 <- params[2] 
  K1plus <- params[3]
  fmax <- params[4]
  N_initial <- params[5]
  
  # calculate derived parameters
  f0 <- (1-s1) / (s0*s1^(12))
  z <- 2.39
  
  # create transition matrix for population projection
  transition <- matrix(0, nrow = 14, ncol = 14) 
  transition[1,14] <- f0  
  transition[2,1] <- s0
  
  for(i in 2:13){
    transition[i+1,i] <- s1
  }
  
  transition[14,14] <- s1
  
  # before projecting the population need to determine the # of animals in each age class
  ee <- eigen(transition) #eigen 
  age_struct <- as.numeric(ee$vectors[,1]/sum(ee$vectors[,1]))#population by age class
  age_struct <- as.numeric(ee$vectors[,1])
  
  # create population matrix and catch matrix
  N <- matrix(NA, nrow = Nages, ncol = Nyears+1) # rows are ages, columns years (+1 because we need 2023 estimates for part B)
  C <- matrix(NA, nrow = 13, ncol = (Nyears+1)) # nrow = 13 because age 0 whales aren't caught
  
  # fill in age and catch data for 2002
  N[,1] <- age_struct*N_initial/sum(age_struct[-1]) # need to take out age 0 whales from the total, since they don't contribute to carrying capacity
  C[,1]  <- Catch[1]*N[2:Nages,1]/sum(N[2:Nages,1]) 
  
  for(t in 2:(Nyears+1)){
    Ntmp <- c(N[1, t-1],N[2:Nages, t-1] - C[,t-1]) # make a temp vector with population # by age class after yearly catch
    Ntmp <- mapply(max, Ntmp, 0.0001) #to make sure not negative
    
    N[,t] <- transition %*% Ntmp # calculate the number of whales that survive
    N[,t] <- mapply(max, N[, t], 0.0001) # remove any negatives
    Ntot <- sum(N[2:Nages,t]) # find N1+ (# of whales greater than age 0)
    
    if(!is.na(Catch[t])){
      C[,t] <- (N[2:Nages,t]/Ntot) * Catch[t] #Calculate catches per group for years where there's catch
    } else {C[,t] <- 0} #For years with no catches
    
    N[1,t] <- N[Nages,t]*(f0 + (fmax - f0)*(1- (Ntot/K1plus)^z)) # calculate # of age 0 whales
  }
  
  year_ests <- colSums(N[2:14,])
  targ_years <- c(year_ests[1], year_ests[2], year_ests[21])
  return(targ_years)
}


# run simulations ---------------------------------------------------------

par_data <- as.matrix(cbind(SIR2_results, pop2002))

# first harvest rule
results1 <- matrix(data = NA, nrow = 200, ncol = 3)

for(i in 1:nrow(par_data)){
  pars <- par_data[i,]
  results1[i,] <- pop_model2(params = pars, data = catch1, Nyears = Nyears)
}
colnames(results1) <- c("y_2002", "y_2003", "y_2023")
results1 <- as.data.frame(results1)

# second harvest rule
results2 <- matrix(data = NA, nrow = 200, ncol = 3)

for(i in 1:nrow(par_data)){
  pars <- par_data[i,]
  results2[i,] <- pop_model2(params = pars, data = catch2, Nyears = Nyears)
}
colnames(results2) <- c("y_2002", "y_2003", "y_2023")
results2 <- as.data.frame(results2)

# third harvest rule
results3 <- matrix(data = NA, nrow = 200, ncol = 3)

for(i in 1:nrow(par_data)){
  pars <- par_data[i,]
  results3[i,] <- pop_model2(params = pars, data = catch3, Nyears = Nyears)
}
colnames(results3) <- c("y_2002", "y_2003", "y_2023")
results3 <- as.data.frame(results3)
# summarize results -------------------------------------------------------

# first harvest rule
results1 %>% 
  filter(y_2002 < 7000) %>% 
  filter(y_2023 > y_2003) # 52
results1 %>% 
  filter(y_2002 >= 7000 & y_2002 < 8000) %>% 
  filter(y_2023 > y_2003) # 102
results1 %>% 
  filter(y_2002 > 8000) %>% 
  filter(y_2023 > y_2003) # 44

# second harvest rule
results2 %>% 
  filter(y_2002 < 7000) %>% 
  filter(y_2023 > y_2003) # 21
results2 %>% 
  filter(y_2002 >= 7000 & y_2002 < 8000) %>% 
  filter(y_2023 > y_2003) # 93
results2 %>% 
  filter(y_2002 > 8000) %>% 
  filter(y_2023 > y_2003) # 44

# third harvest control rule
results3 %>% 
  filter(y_2002 < 7000) %>% 
  filter(y_2023 > y_2003) # 2
results3 %>% 
  filter(y_2002 >= 7000 & y_2002 < 8000) %>% 
  filter(y_2023 > y_2003) # 46
results3 %>% 
  filter(y_2002 > 8000) %>% 
  filter(y_2023 > y_2003) # 44
