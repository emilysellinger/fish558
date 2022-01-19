# Homework 5
library(tidyverse)
library(here)

##################### TASK B ####################################################

# beta rv generator -------------------------------------------------------

beta_rv <- function(mu, cv){
  
  if(mu == 1){
    num <- 1
  }else{
    alpha <- (1 - mu*(1 + cv^2))/(cv^2)
    beta <- alpha*(1 - mu)/mu
    num <- rbeta(1, alpha, beta)
  }
  
  return(num)
}

########################### TASK C ###############################################


# population model --------------------------------------------------------

pop_mod <- function(sj, area){
  pop_inits <- c(98, 68, 48, 38, 27, 24, 16, 11, 12, 13, 3, 11, 1, 5, 2, 113)
  
  preds <- matrix(0, nrow = 16, ncol = 101)
  preds[,1] <- pop_inits
  b <- 0.75
  sa <- 0.95
  
  for(i in 2:ncol(preds)){
    # calculate the phi's for each management scenario
    if(area == 0){
      phi1 <- 1; phi2 <- 1; phi3 <- 1
    }
    if(area == 1){
      phi1 <- beta_rv(0.7, 0.1); phi2 <- 1; phi3 <- 1
    }
    if(area == 2){
      phi1 <- 1; phi2 <- beta_rv(0.95,0.1); phi3 <- 1
    }
    if(area == 3){
      phi1 <- 1; phi2 <- 1; phi3 <- beta_rv(0.97, 0.1)
    }
    
    # calculate the total number of adults
    ntilda <- sum(preds[10:16, i-1])
    
    preds[1, i] <- rbinom(1, rbinom(1, ntilda, b), phi1)
    preds[2, i] <- rbinom(1, rbinom(1, preds[1, i-1], sj), phi2)
    preds[3, i] <- rbinom(1, rbinom(1, preds[2, i-1], sj), phi2)
    preds[4, i] <- rbinom(1, rbinom(1, preds[3, i-1], sj), phi2)
    preds[5, i] <- rbinom(1, rbinom(1, preds[4, i-1], sj), phi2)
    preds[6, i] <- rbinom(1, rbinom(1, preds[5, i-1], sj), phi2)
    preds[7, i] <- rbinom(1, rbinom(1, preds[6, i-1], sj), phi2)
    preds[8, i] <- rbinom(1, rbinom(1, preds[7, i-1], sj), phi2)
    preds[9, i] <- rbinom(1, rbinom(1, preds[8, i-1], sj), phi2)
    preds[10, i] <- rbinom(1, rbinom(1, preds[9, i-1], sa), phi3)
    preds[11, i] <- rbinom(1, rbinom(1, preds[10, i-1], sa), phi3)
    preds[12, i] <- rbinom(1, rbinom(1, preds[11, i-1], sa), phi3)
    preds[13, i] <- rbinom(1, rbinom(1, preds[12, i-1], sa), phi3)
    preds[14, i] <- rbinom(1, rbinom(1, preds[13, i-1], sa), phi3)
    preds[15, i] <- rbinom(1, rbinom(1, preds[14, i-1], sa), phi3)
    preds[16, i] <- rbinom(1, rbinom(1, preds[15, i-1], sa), phi3) + rbinom(1, rbinom(1, preds[16, i-1], sa), phi3)
  }
  
  # calculate the total number of adults after 100 years
  total <- sum(preds[10:16, 101])
  return(total)
}


# simulation function -----------------------------------------------------

run_sims <- function(nsims, sj, area){
  results <- rep(0, nsims)
  for(i in 1:nsims){
    results[i] <- pop_mod(sj, area)
  }
  return(results)
}

# results for harvest strategy 1 -------------------------------------------
sj <- c(0.75, 0.8, 0.82)

harvest1 <- matrix(data = 0, nrow = 1000, ncol = 3)

for(j in 1:3){
  harvest1[,j] <- run_sims(1000, sj[j], 1)
}
harvest1 <- as.data.frame(harvest1)
write_csv(harvest1, file = "homework5/harvest_strat1.csv")

# results for harvest strategy 2 ------------------------------------------
harvest2 <- matrix(data = 0, nrow = 1000, ncol = 3)

for(j in 1:3){
  harvest2[,j] <- run_sims(1000, sj[j], 2)
}

harvest2 <- as.data.frame(harvest2)
write_csv(harvest2, file = "homework5/harvest_strat2.csv")

# results for harvest strategy 3 ------------------------------------------
harvest3 <- matrix(data = 0, nrow = 1000, ncol = 3)

for(j in 1:3){
  harvest3[,j] <- run_sims(1000, sj[j], 3)
}

harvest3 <- as.data.frame(harvest3)
write_csv(harvest3, file = "homework5/harvest_strat3.csv")
# results for no harvest strategy -----------------------------------------
no_harvest <- matrix(data = 0, nrow = 1000, ncol = 3)

for(j in 1:3){
  no_harvest[,j] <- run_sims(1000, sj[j], 0)
}

no_harvest <- as.data.frame(no_harvest)
write_csv(no_harvest, file = "homework5/no_harvest_strat.csv")


# calculate expected performance metrics ----------------------------------

# read in data if necessary
harvest1 <- read_csv(here("homework5", "harvest_strat1.csv"))
harvest2 <- read_csv(here("homework5", "harvest_strat2.csv"))
harvest3 <- read_csv(here("homework5", "harvest_strat3.csv"))
no_harvest <- read_csv(here("homework5", "no_harvest_strat.csv"))

# expected number of adults after 100 years
apply(harvest1, 2, mean)
apply(harvest2, 2, mean)
apply(harvest3, 2, mean)
apply(no_harvest, 2, mean)

# probability the number of adults is less than 1000
prob <- function(vec){
  n <- 0
  for(i in 1:length(vec)){
    if(vec[i] < 1000){
      n <- n + 1
    }
  }
  
  prop <- n/length(vec)
  return(prop)
}

apply(harvest1, 2, prob)
apply(harvest2, 2, prob)
apply(harvest3, 2, prob)
apply(no_harvest, 2, prob)

################################### PART D #####################################################################

# harvesting model --------------------------------------------------------

# There are definitely faster ways to code this model
hunt_mod <- function(sj, area){
  pop_inits <- c(98, 68, 48, 38, 27, 24, 16, 11, 12, 13, 3, 11, 1, 5, 2, 113)
  
  preds <- matrix(0, nrow = 16, ncol = 101)
  preds[,1] <- pop_inits
  
  # create matrix for number of hunted individuals
  harvests <- matrix(0, nrow = 16, ncol = 101)
  
  b <- 0.75
  sa <- 0.95
  
  for(i in 2:ncol(preds)){
    # calculate the phi's for each management scenario
    if(area == 1){
      phi1 <- beta_rv(0.7, 0.1); phi2 <- 1; phi3 <- 1
    }
    if(area == 2){
      phi1 <- 1; phi2 <- beta_rv(0.95,0.1); phi3 <- 1
    }
    if(area == 3){
      phi1 <- 1; phi2 <- 1; phi3 <- beta_rv(0.97, 0.1)
    }
    
    # calculate the total number of adults
    ntilda <- sum(preds[10:16, i-1])
    
    x1 <- rbinom(1, ntilda, b)
    preds[1, i] <- rbinom(1, x1 , phi1)
    harvests[1, i] <- x1 - preds[1, i]
    #
    x2 <- rbinom(1, preds[1, i-1], sj)
    preds[2, i] <- rbinom(1, x2, phi2)
    harvests[2, i] <- x2 - preds[2, i]
    #
    x3 <- rbinom(1, preds[2, i-1], sj)
    preds[3, i] <- rbinom(1, x3, phi2)
    harvests[3, i] <- x3 - preds[3, i]
    #
    x4 <- rbinom(1, preds[3, i-1], sj)
    preds[4, i] <- rbinom(1, x4, phi2)
    harvests[4, i] <- x4 - preds[4, i]
    #
    x5 <- rbinom(1, preds[4, i-1], sj)
    preds[5, i] <- rbinom(1, x5, phi2)
    harvests[5, i] <- x5 - preds[5, i]
    #
    x6 <- rbinom(1, preds[5, i-1], sj)
    preds[6, i] <- rbinom(1, x6, phi2)
    harvests[6, i] <- x6 - preds[6, i]
    #
    x7 <- rbinom(1, preds[6, i-1], sj)
    preds[7, i] <- rbinom(1, x7, phi2)
    harvests[7, i] <- x7 - preds[7, i]
    #
    x8 <- rbinom(1, preds[7, i-1], sj)
    preds[8, i] <- rbinom(1, x8, phi2)
    harvests[8, i] <- x8 - preds[8, i]
    #
    x9 <- rbinom(1, preds[8, i-1], sj)
    preds[9, i] <- rbinom(1, x9, phi2)
    harvests[9, i] <- x9 - preds[9, i]
    #
    x10 <- rbinom(1, preds[9, i-1], sa)
    preds[10, i] <- rbinom(1, x10, phi3)
    harvests[10, i] <- x10 - preds[10, i]
    #
    x_11 <- rbinom(1, preds[10, i-1], sa)
    preds[11, i] <- rbinom(1, x_11, phi3)
    harvests[11, i] <- x_11 - preds[11, i]
    #
    x12 <- rbinom(1, preds[11, i-1], sa)
    preds[12, i] <- rbinom(1, x12, phi3)
    harvests[12, i] <- x12 - preds[12, i]
    #
    x13 <- rbinom(1, preds[12, i-1], sa)
    preds[13, i] <- rbinom(1, x13, phi3)
    harvests[13, i] <- x13 - preds[13, i]
    #
    x14 <- rbinom(1, preds[13, i-1], sa)
    preds[14, i] <- rbinom(1, x14, phi3)
    harvests[14, i] <- x14 - preds[14, i]
    #
    x15 <- rbinom(1, preds[14, i-1], sa)
    preds[15, i] <- rbinom(1, x15, phi3)
    harvests[15,i] <- x15 - preds[15, i]
    #
    x16a <- rbinom(1, preds[15, i-1], sa); x16b <- rbinom(1, preds[16, i-1], sa); x16 <- x16a + x16b
    preds[16, i] <- rbinom(1, x16a, phi3) + rbinom(1, x16b, phi3)
    harvests[16, i] <- x16 - preds[16, i]
  }
  
  # calculate the total number of animals harvested over 100 years
  total <- sum(harvests)
  return(total)
}
# simulation function 2 -----------------------------------------------------

run_sims2 <- function(nsims, sj, area){
  results <- rep(0, nsims)
  for(i in 1:nsims){
    results[i] <- hunt_mod(sj, area)
  }
  return(results)
}
# run simulations for harvest strategy 1 ----------------------------------
sj <- c(0.75, 0.8, 0.82)

hunts1 <- matrix(data = 0, nrow = 1000, ncol = 3)

for(j in 1:3){
  hunts1[,j] <- run_sims2(1000, sj[j], 1)
}

hunts1 <- as.data.frame(hunts1)
write_csv(hunts1, file = "homework5/hunts_strat1.csv")


# run simulations for harvest strategy 2 ----------------------------------
hunts2 <- matrix(data = 0, nrow = 1000, ncol = 3)

for(j in 1:3){
  hunts2[,j] <- (run_sims2(1000, sj[j], 2))*5
}

hunts2 <- as.data.frame(hunts2)
write_csv(hunts2, file = "homework5/hunts_strat2.csv")

# run simulations for harvest strategy 3 ----------------------------------
hunts3 <- matrix(data = 0, nrow = 1000, ncol = 3)

for(j in 1:3){
  hunts3[,j] <- (run_sims2(1000, sj[j], 3))*30
}

hunts3 <- as.data.frame(hunts3)
write_csv(hunts3, file = "homework5/hunts_strat3.csv")

# summarize results -------------------------------------------------------
apply(hunts1, 2, mean)
apply(hunts2, 2, mean)
apply(hunts3, 2, mean)

# load in data if needed --------------------------------------------------
hunts1 <- read_csv(here("homework5", "hunts_strat1.csv"))
hunts2 <- read_csv(here("homework5", "hunts_strat2.csv"))
hunts3 <- read_csv(here("homework5", "hunts_strat3.csv"))