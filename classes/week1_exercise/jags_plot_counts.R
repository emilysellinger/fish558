library(here)
library(tidyverse)
library(runjags)
Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.0")

# load data ---------------------------------------------------------------

data <- read_delim(here("classes/week1_exercise", "Lecture3.dat"))


# write model -------------------------------------------------------------

model <- "model {
  D[1] <- 0
  for(i in 2 : N){ 
 
 D[i] ~ dlnorm(mu, precision)
 
 C[i] ~ dpois(D[i] * E[i])
 
 }
 
 # Priors
 mu ~ dunif(-1000, 1000)
 precision ~ dunif(0, 1000)
  
}"


# data and inits ----------------------------------------------------------

C <- data$Count
E <- data$Effort
N <- length(C)
model_data <- list(C = C, E = E, N = N)

# starting values for chains
inits1 <- list(mu = 0, precision = 500)
inits2 <- list(mu = 0, precision = 200)
inits3 <- list(mu = 1, precision = 400,
               .RNG.name="base::Wichmann-Hill", .RNG.seed=9)
inits4 <- list(mu = 0, precision = 800,
               .RNG.name="base::Wichmann-Hill", .RNG.seed=4)
inits5 <- list(mu = 0, precision = 600,
               .RNG.name="base::Wichmann-Hill", .RNG.seed=3)

# combine chains
inits <- list(inits1,inits2,inits3,inits4,inits5)


# run model ---------------------------------------------------------------
results <- run.jags(model = model, monitor = c("mu", "precision", "D"), 
                    data = model_data, n.chains = 5, method="rjags", inits = inits,
                    plots=T,silent.jag=T)

print(results)
plot(results, layout = runjags.getOption("plot.layout"),
     new.windows = runjags.getOption("new.windows"), file = "testjags.pdf")


# plot posteriors ---------------------------------------------------------

mcmc <- rbind(results$mcmc[[1]],results$mcmc[[2]], results$mcmc[[3]])


hist(mcmc[,1],xlab="mu",main="")
hist(mcmc[,2],xlab="precision",main="")

# get distributions for each plot
data %>%
  group_by(Plot) %>%
  count()


ResSum <- matrix(0,ncol=23,nrow=5)
for (Iyr in 1:23)
  ResSum[,Iyr] <- quantile(mcmc[,Iyr+3],probs=c(0.05,0.25,0.5,0.75,0.95))  
xx <- seq(1:23)
plot(xx,ResSum[3,],xlab="Year",ylab="Depletion",type='l',lwd=3,ylim=c(0,1.3))
xx2 <- c(xx,rev(xx))
polygon(xx2,c(ResSum[1,],rev(ResSum[5,])),col="gray50")
polygon(xx2,c(ResSum[2,],rev(ResSum[4,])),col="gray95")
lines(xx,ResSum[3,],lwd=3,lty=1)

# Plot the posterior predictive distribution
par(mfrow=c(2,2))
plot(xx,I,pch=16,ylim=c(0,100))
ResSum <- matrix(0,ncol=23,nrow=5)
for (Iyr in 1:23)
  ResSum[,Iyr] <- quantile(mcmc[,26+Iyr],probs=c(0.05,0.25,0.5,0.75,0.95))  
lines(xx,ResSum[3,],lwd=3,lty=1) 
lines(xx,ResSum[1,],lwd=1,lty=2) 
lines(xx,ResSum[5,],lwd=1,lty=2) 

