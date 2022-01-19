#install necessary libraries
library(tidyverse)
library(stats4)
library(here)



#Download data
data <- read_csv(here("homework3", "bowhead_whale_catch.csv"))


###Part A2
#Parameters
S0 <- 0.9 #.7?
S1plus <-0.95 #
K1plus <- 15000
N_initial <- K1plus #Assume initial population  size is at pre-exploitation carrying capacity
fmax <- 0.29
f0 <- (1-S1plus) / (S0*S1plus^(13)) 
pars <- c(S0, S1plus, K1plus, fmax)
Nyears <- as.numeric(nrow(data))

#Function of population model
popmodel <- function(pars, data, Nyears){
  #data 
  Catch <- data$Catch
  Nages <- 14
  
  #parameters
  S0 <- pars[1] #0.7
  S1plus <- pars[2] #0.9
  K1plus <- pars[3] #15000
  fmax <- pars[4]
  f0 <- (1-S1plus) / (S0*S1plus^(12)) 
  N_initial <- K1plus
  z <- 2.39
  
  #parameter matrix:
  parmat <- matrix(0, nrow = 14, ncol = 14) #Set up age class x age class matrix
  parmat[1,14] <- f0 + (fmax - f0)*(1- (N_initial/K1plus)^z) #Just for this cell, use initial  population size
  parmat[2,1] <- S0
  for(i in 2:13){
    parmat[i+1,i] <- S1plus
  }
  
  parmat[14,14] <- S1plus
  ee <- eigen(parmat) #eigen 
  age_struct <- as.numeric(ee$vectors[,1]/sum(ee$vectors[,1]))#population by age class
  
  #Abundance
  N <- matrix(NA, nrow = Nages, ncol = Nyears+1) #Abundance matrix w/ ages down rows and years across columns
  C <- matrix(NA, nrow = 13, ncol = (Nyears+1)) #Catch matrix
  N[,1] <- age_struct*K1plus
  C[,1] <- 0
  
  for(t in 2:(Nyears+1)){
    N[2:Nages,t-1] <- N[2:Nages, t-1] - C[,t-1] 
    N[2:Nages,t-1] <- mapply(max, N[2:Nages, t-1], 0.0001) #to  make sure not negative
    N[,t] <- parmat %*% N[,t-1]
    N[,t] <- mapply(max, N[, t], 0.0001)
    Ntot <- sum(N[2:Nages,t])
    
    if(!is.na(Catch[i-1])){
      C[,t] <- (N[2:Nages,t]/Ntot) * Catch[t-1] #Calculate catches per group for years where there's catch
    } else {C[,t]} <- 0 #For years with no catches
    
    z <- 2.39
    parmat[1,14] <- f0 + (fmax - f0)*(1- (Ntot/K1plus)^z) 
  }
  
  return(N)
  print(N)
}

#Get predicted abundance values from population model
N_projected <- popmodel(pars=pars, data=data, Nyears=Nyears)

#Label columns with year names
colnames(N_projected) <- 1847:2002

#Abundance in 2002
Abundance_2002 <- N_projected[,156]
print(Abundance_2002)
Abundance_2002_total <- sum(Abundance_2002)
print(Abundance_2002_total)

#Likelihood part

#List of data w/ abundance estimates & CV available
Years_with_est <- c(1978, 1980, 1981, 1982, 1983, 1985, 1986, 1987, 1988, 1993)
N_observed <- c(4820, 3900, 4389, 6572, 6268, 5132, 7251, 5151, 6609, 7778)
CV <- c(0.273, 0.314, 0.253, 0.311, 0.321, 0.269, 0.186, 0.298, 0.113, 0.071)

observed <- data.frame(Years_with_est, N_observed, CV)
N_predicted <- N_projected[, c("1978", "1980", "1981", "1982", "1983", "1985", "1986", "1987", "1988", "1993")]
N_predicted <- colSums(N_predicted)

minusLL <- function(pars, N_predicted, observed){
  
  #likelihood estimates
  sigma <- observed$CV
  N_observed <- observed$N_observed
  L <- -1*sum(-(log(N_predicted)-log(N_observed))^2 / (2*sigma^2))
  
  return(L)
}

#Get negative log likelihood
negative_log_likelihood <- minusLL(pars=pars, N_predicted=N_predicted, observed=observed)
print(negative_log_likelihood)



#Part A3

DoSir <- function(Nout, observed){
  
  # Read in the basic data
  
  # Create matrix of parameters
  Vals <- matrix(0,ncol=4,nrow=Nout)
  
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
    # Generate from priors
    S0 <- runif(1, 0.8, 1)
    S1plus <- runif(1, 0.9, 1)
    K1plus <- runif(1, 10000, 20000)
    fmax <- runif(1, 0.25, 0.33)
    pars <- c(S0, S1plus, K1plus, fmax)
    
    # Call the population model
    Pop <- popmodel(pars, data, Nyears)
    Pop_total <- sum(Pop[,156])
    
    #Determine if extinct or not
    TheLike  <- ifelse(Pop_total>1, 1,0)
    
    # Determine if a parameter vector is to be saved
    Cumu <- Cumu + TheLike
    while (Cumu >= Threshold & Ndone < Nout)
    {
      Ndone <- Ndone + 1
      Cumu <- Cumu - Threshold
      # need to save posterior informative here - params plus ratio of depletion
      Vals[Ndone,] <- c(S0, S1plus,K1plus,fmax)
    }
  }
  
  return(Vals)
  
}

test <- DoSir(Nout=1000, observed=observed)


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#A4: Use  the  SIR  algorithm  to  sample  200  parameter  vectors  from  the  posterior distribution.

DoSir2 <- function(Nout,minusLL, popmodel, data, observed){
  
  # Data stuff
  data <- data
  Catch <- data$Catch
  Nyears<- nrow(data)
  
  # Create matrix of parameters
  Vals <- matrix(0,ncol=4,nrow=Nout)
  
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
    colnames(N_projected) <- 1847:2002
    Years_with_est <- c(1978, 1980, 1981, 1982, 1983, 1985, 1986, 1987, 1988, 1993)
    N_predicted <- N_projected[, c("1978", "1980", "1981", "1982", "1983", "1985", "1986", "1987", "1988", "1993")]
    N_predicted <- colSums(N_predicted)
    
    # Compute the negative log-likelihood and hence the likelihood
    NegLogLike <- minusLL(pars=pars, N_predicted=N_predicted, observed=observed)
    TheLike <- exp(-1*NegLogLike-2.22)
    
    # Determine if a parameter vector is to be saved
    Cumu <- Cumu + TheLike
    while (Cumu >= Threshold & Ndone < Nout)
    {
      Ndone <- Ndone + 1
      Cumu <- Cumu - Threshold
      # need to save posterior informative here - params plus ratio of depletion
      Vals[Ndone,] <- c(S0, S1plus,K1plus,fmax)
      print(Ndone)
    }
  }
  
  return(Vals)
  
}

test2 <- DoSir2(Nout=10, minusLL=minusLL, popmodel=popmodel, data=data, observed=observed)

#Summarize  the  results  of  1)  and  2) :

#by  plots  of  distributions  for S0, S1plus, K, fmax

#posterior  for  the  time-trajectory  of  1+  population  size from 1848 to 2002 (you can summarize the distribution of 1+ population size for each year by its 5th, median and 95th percentiles).



# Task B ------------------------------------------------------------------
# Create a decision table

# goal will be to use the median parameter values from the A4 question, in the meantime
# I'll just pick pre-set values for now to get the code working

# define parameters
S0 <- 0.9 
S1plus <-0.95 
K1plus <- 15000
N_initial <- K1plus #Assume initial population  size is at pre-exploitation carrying capacity
fmax <- 0.29
f0 <- (1-S1plus) / (S0*S1plus^(13)) 
pars <- c(S0, S1plus, K1plus, fmax)

Nyears <- 20

# will simulate 100 different starting values for 2003 (not sure if this is the correct method, but will
# use it for now)

sim1 <- runif(100, min = 100, max = 7000)
sim2 <- runif(100, min = 7000, max = 8000)
sim3 <- runif(100, min = 8000, max = 20000)

# simulate the catch data for each of our harvest options
catch1 <- rep(67, 20)
catch2 <- rep(134, 20)
catch1 <- rep(67, 20)
pop_start <- 