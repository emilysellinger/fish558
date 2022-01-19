# Gray Whale exercise

library(tidyverse)
library(here)

# 
Ex3 <- function()
{
  # call for the logistic model
  log_mod_vals <- DoSir(Nout=200,Model=F)
  
  # call for the exponential model
  exp_mod_vals <- DoSir(Nout=200,Model=T) 
  
  values <- list(log_mod_vals, exp_mod_vals)
  return(values)
}



# Read Data -------------------------------------------------------------------------------------
ReadData <- function()
{
  FileName <- "day3/Ex3a.csv"
  TheData1 <- matrix(scan(FileName,skip=1,n=22*3,sep=','),ncol=3,byrow=T)
  FileName <- "day3/Ex3b.csv"
  TheData2 <- matrix(scan(FileName,skip=1,n=38*2,sep=','),ncol=2,byrow=T)
  
  Outs <- NULL
  Outs$SurveyYr <- TheData1[,1]
  Outs$SurveyEst <- TheData1[,2]
  Outs$SurveyCV <- TheData1[,3]
  Outs$CatchYr <- TheData2[,1]
  Outs$CatchVal <- TheData2[,2]
  return(Outs)
}

# Population Function ----------------------------------------------------------------------------
PopModel <- function(Catch,r,K,InitPop,ExponModel)
{
  # Conduct projections and return a vector called "pop" that includes population numbers
  # Choose between the exponential and logistic models depends on ExponModel
  C <- Catch; r <- r; k <- K; n_0 <- InitPop;
  pop <- rep(NA, length(C))
  
  if(ExponModel == T){
    
    # calculate projections
    for(i in 1:length(C)){
      if(i == 1){
        pop[i] = (1 + r)*n_0
      }else if(i > 1){
        pop[i] = (1 + r)*pop[i-1] - C[i-1]
      }
    }
  }else if(ExponModel == F){
    
    # calculate projections
    for(i in 1:length(C)){
      if(i == 1){
        pop[i] = n_0 + r*n_0*(1 - (n_0/k))
        
      }else if(i > 1){
        pop[i] = pop[i-1] + r*pop[i-1]*(1 - (pop[i-1]/k)) - C[i-1]
        
      }
    }
    
  }
  
  return(pop)
}

# Likelihood Function --------------------------------------------------------------------------
Likelihood <- function(Pop,SurveyYr,SurveyEst,SurveyCV,AddCV)
{
  # Account for the additional CV  
  UseCV <- sqrt(SurveyCV^2+AddCV^2)
  
  # Extract the predictions corresponding to the observations and compute the negatuve log-likelihood
  Preds <- Pop[SurveyYr]
  Residuals <- log(UseCV)+0.5*(log(Preds)-log(SurveyEst))^2/UseCV^2
  LogLike <- sum(Residuals)
  
  return(LogLike)
}

# SIR Algorithm ----------------------------------------------------------------------------------
DoSir <- function(Nout=1000,Model)
{
  
  # Read in the basic data
  TheData <- ReadData()
  Yr1 <- TheData$CatchYr[1]
  Catch <- TheData$CatchVal
  SurveyEst <- TheData$SurveyEst
  SurveyCV <- TheData$SurveyCV
  SurveyYr <- TheData$SurveyYr
  Nyears <- length(Catch)
  
  # Storage for the parameters and the final depletion
  Vals <- matrix(0,ncol=5,nrow=Nout)
  colnames(Vals) <- c("K", "r", "Pop1965", "AddCV", "End_Deplet")
  
  # Storage for the total likelihood encountered in the first stage sampling
  # and the number of draws
  AveLike <- 0
  Ntest <- 0
  
  # Reset parameters for SIR
  Threshold <- exp(0)   
  Cumu <- 0
  Ndone <- 0
  while (Ndone < Nout)
  {
    # Generate from priors
    K <- runif(1, min = 20000, max = 50000)
    r <- runif(1, min = 0, max = 0.15)
    Pop1965 <- runif(1, min = 10000, max = 15000)
    AddCV <- runif(1, min = 0.1, max = 0.2)
    
    # Call the population model
    Pop <- PopModel(Catch,r,K,Pop1965,Model)
    
    # Compute the negative log-likelihood and hence the likelihood
    NegLogLike <- Likelihood(Pop,SurveyYr-Yr1+1,SurveyEst,SurveyCV,AddCV)
    TheLike <- exp(-1*NegLogLike-32.19)
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
      Vals[Ndone, 3] <- Pop1965
      Vals[Ndone, 4] <- AddCV
      Vals[Ndone, 5] <- Pop[38]/K
    }
  }
  
  BayesFact <- AveLike/Ntest
  return(list(Vals, BayesFact))
}

values <- Ex3()

# Calculate Bayes Factor ---------------------------------------------------------------------------
exp_mod_results <- DoSir(200, T)
log_mod_results <- DoSir(200, F)

bayes_fac_exp <- exp_mod_results[[2]]/(log_mod_results[[2]] + exp_mod_results[[2]])
bayes_fac_log <- log_mod_results[[2]]/(log_mod_results[[2]] + exp_mod_results[[2]])
# Task B --------------------------------------------------------------------------------------------
ReadData2 <- function()
{
  FileName <- "day3/Ex3a.csv"
  TheData1 <- matrix(scan(FileName,skip=1,n=22*3,sep=','),ncol=3,byrow=T)
  FileName <- "day3/Ex3b.csv"
  TheData2 <- matrix(scan(FileName,skip=1,n=38*2,sep=','),ncol=2,byrow=T)
  
  Outs <- NULL
  Outs$SurveyYr <- TheData1[,1]
  Outs$SurveyYr <- Outs$SurveyYr[-c(21,22)]
  Outs$SurveyEst <- TheData1[,2]
  Outs$SurveyEst <- Outs$SurveyEst[-c(21,22)]
  Outs$SurveyCV <- TheData1[,3]
  Outs$SurveyCV <- Outs$SurveyCV[-c(21,22)]
  Outs$CatchYr <- TheData2[,1]
  Outs$CatchVal <- TheData2[,2]
  return(Outs)
}

DoSir2 <- function(Nout=1000,Model)
{
  
  # Read in the basic data
  TheData <- ReadData2()
  Yr1 <- TheData$CatchYr[1]
  Catch <- TheData$CatchVal
  SurveyEst <- TheData$SurveyEst
  SurveyCV <- TheData$SurveyCV
  SurveyYr <- TheData$SurveyYr
  Nyears <- length(Catch)
  
  # Storage for the parameters and the final depletion
  Vals <- matrix(0,ncol=5,nrow=Nout)
  colnames(Vals) <- c("K", "r", "Pop1965", "AddCV", "End_Deplet")
  
  # Storage for the total likelihood encountered in the first stage sampling
  # and the number of draws
  AveLike <- 0
  Ntest <- 0
  
  # Reset parameters for SIR
  Threshold <- exp(0)   
  Cumu <- 0
  Ndone <- 0
  while (Ndone < Nout)
  {
    # Generate from priors
    K <- runif(1, min = 20000, max = 50000)
    r <- runif(1, min = 0, max = 0.15)
    Pop1965 <- runif(1, min = 10000, max = 15000)
    AddCV <- runif(1, min = 0.1, max = 0.2)
    
    # Call the population model
    Pop <- PopModel(Catch,r,K,Pop1965,Model)
    
    # Compute the negative log-likelihood and hence the likelihood
    NegLogLike <- Likelihood(Pop,SurveyYr-Yr1+1,SurveyEst,SurveyCV,AddCV)
    TheLike <- exp(-1*NegLogLike-32.19)
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
      Vals[Ndone, 3] <- Pop1965
      Vals[Ndone, 4] <- AddCV
      Vals[Ndone, 5] <- Pop[38]/K
    }
  }
  
  BayesFact <- AveLike/Ntest
  return(list(Vals, BayesFact))
}


# Calculate Bayes Factor 
exp_mod_results <- DoSir2(200, T)
log_mod_results <- DoSir2(200, F)

bayes_fac_exp <- exp_mod_results[[2]]/(log_mod_results[[2]] + exp_mod_results[[2]])
bayes_fac_log <- log_mod_results[[2]]/(log_mod_results[[2]] + exp_mod_results[[2]])
# much higher support for exponential model when last two years are removed

# Task C -----------------------------------------------------------------------------------------
# going to write a new SIR function where last column is the ending population
# then will calculate the number of the 200 runs that are above 90% of K

DoSir3 <- function(Nout=1000,Model)
{
  
  # Read in the basic data
  TheData <- ReadData()
  Yr1 <- TheData$CatchYr[1]
  Catch <- TheData$CatchVal
  SurveyEst <- TheData$SurveyEst
  SurveyCV <- TheData$SurveyCV
  SurveyYr <- TheData$SurveyYr
  Nyears <- length(Catch)
  
  # Storage for the parameters and the final depletion
  Vals <- matrix(0,ncol=5,nrow=Nout)
  colnames(Vals) <- c("K", "r", "Pop1965", "AddCV", "End_Pop")
  
  # Storage for the total likelihood encountered in the first stage sampling
  # and the number of draws
  AveLike <- 0
  Ntest <- 0
  
  # Reset parameters for SIR
  Threshold <- exp(0)   
  Cumu <- 0
  Ndone <- 0
  while (Ndone < Nout)
  {
    # Generate from priors
    K <- runif(1, min = 20000, max = 50000)
    r <- runif(1, min = 0, max = 0.15)
    Pop1965 <- runif(1, min = 10000, max = 15000)
    AddCV <- runif(1, min = 0.1, max = 0.2)
    
    # Call the population model
    Pop <- PopModel(Catch,r,K,Pop1965,Model)
    
    # Compute the negative log-likelihood and hence the likelihood
    NegLogLike <- Likelihood(Pop,SurveyYr-Yr1+1,SurveyEst,SurveyCV,AddCV)
    TheLike <- exp(-1*NegLogLike-32.19)
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
      Vals[Ndone, 3] <- Pop1965
      Vals[Ndone, 4] <- AddCV
      Vals[Ndone, 5] <- Pop[38]
    }
  }
  
  BayesFact <- AveLike/Ntest
  return(Vals)
}

values <- DoSir3(Nout = 200, F) # logistic population
Ks <- values[,1]
Ends <- values[,5]

count <-  0
for(i in 1:200){
  if(Ends[i] > 0.9*Ks[i]){
    count <- count + 1
  }
}

prob <- count/200 # in this run, the probability the population is close to the carrying capacity is 0.17


# check solution file for Andre's code plotting the results