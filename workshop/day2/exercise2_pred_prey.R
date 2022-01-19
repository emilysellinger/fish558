# Exercise 2

data <- read.csv(here("day2", "EX2R.DAT"), sep = "")

# Notes - going to need to use optim
# predators - Nc

predators <- data[,1]
consumptions <- data[,c(5,6,7)]
prey <- data[,c(2, 3, 4)]

SSQ <- function(observed, predicted){
  sum((observed - predicted)^2)
}

minSSQ <- function(param_vect, Type = "Linear"){
  if(Type == "Linear"){
    alphas <- exp(param_vect) # input must be logged
    # matrix for predictions
    Pred <- matrix(NA, nrow = length(predators), ncol = length(alphas))
    
    # calculate predictions
    for(i in 1:length(alphas)){
      Pred[,i] <- alphas[i]*predators
    }
  }else if(Type == "Holling"){
    alphas <- exp(param_vect[1:3]) # vector of vector notation
    betas <- exp(param_vect[4:6])
    
    # matrix for predictions
    Pred <- matrix(NA, nrow = length(predators), ncol = length(alphas))
    
    # calculate predictions
    for(i in 1:length(alphas)){
      Pred[,i] <- alphas[i]*predators/(1 + betas[i]*prey[i])
    }
  }else if(Type == "Sigmoid"){
    alphas <- exp(param_vect[1:3])
    betas <- exp(param_vect[4:6])
    gamma <- length(alphas) + length(betas)
    
    # matrix for predictions
    Pred <- matrix(NA, nrow = length(predators), ncol = length(alphas))
    
    # calculate predictions
    for(i in 1:length(alphas)){
      Pred[,i] <- (alphas[i]*predators*(prey[,i])^(gamma - 1))/(1 + betas[i]*(prey[,i])^gamma)
    }
    
  }else if(Type == "Pre-emption"){
    alphas <- exp(param_vect[1:3])
    betas <- exp(param_vect[4:6])
    lambdas <- exp(param_vect[7:9])
    
    # matrix for predictions
    Pred <- matrix(NA, nrow = length(predators), ncol = length(alphas))
    
    # calculate predictions
    for(i in 1:length(alphas)){
      Pred[,i] <- (alphas[i]*predators)/(1 + betas[i]*prey[,i] + lambdas[i]*predators)
    }
  }
  
  sum_of_sqs <- SSQ(log(consumptions), log(Pred))
  
  return(sum_of_sqs)
}

# parameter vectors
initial_alphas <- c(1,1,1)
initial_betas <- c(1,1,1)
initial_lambdas <- c(1,1,1)


