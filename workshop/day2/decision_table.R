
rvals <- c(0.1, 0.2, 0.3)
kvals <- c(1000, 2000, 3000)
harvest_rates <- c(0.05, 0.1, 0.15)

# want to compare B/K in 20 years

# try something else




projection <- function(carrying_cap = 1000, rate = 0.1, harvest = 0.05){
  K <- carrying_cap
  r <- rate
  h <- harvest
  s_0 <- 0.2*K
  
  biomass <- rep(NA, 21)
  
  for(i in 1:21){
    if(i == 1){
      biomass[i] <- s_0
    }else if(i > 1){
      biomass[i] <- biomass[i-1] + r*biomass[i-1]*(1 - (biomass[i-1]/K)) - h*biomass[i-1]
    }
  }
  
  #print(biomass)
  b_over_k <- biomass[21]/K
  
  return(b_over_k)
}


projection(1000, 0.3, 0.15)

# carrying capacity does not affect the depletion levels, only "state of nature"
# is the growth rate of the stock

depletion_level <- tibble{
  harvest_lev = c("h_0.05", "h_0.1", "h_0.15")
}