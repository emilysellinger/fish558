library(nimble)
library(igraph)

# note: ~ and <- have meanings same as JAGS/BUGS
pumpCode <- nimbleCode({
  
  for (i in 1:N) {
    theta[i] ~ dgamma(alpha,beta)
    lambda[i] <- theta[i]*t[i]
    x[i] ~ dpois(lambda[i])
  }
  
  alpha ~ dexp(1.0)
  beta ~ dgamma(0.1,1.0)
})

# we need to divide the data (something that appears in the likelihood function) and
# constants - two separate lists (warnings are issued if constants are included in data)
pumpConsts <- list(N = 10, t = c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5))

pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4,  22))


pumpInits <- list(alpha = 1, beta = 1, theta = rep(0.1, pumpConsts$N))

# call to set up model - compiles and checks for correct syntax
pump <- nimbleModel(code = pumpCode, name = "pump", constants = pumpConsts, data = pumpData, inits = pumpInits)

pump$getNodeNames()

pump$logProb_x
pump$lifted_d1_over_beta


pump$modelDef

pump$modelDef$BUGScode

# Plot the graphs
pump$plotGraph()

pump$getDependencies(c("alpha"))
pump$getDependencies(c("beta"))

# Check the way the samplers are applied.
pump$checkConjugacy()

# generate from the distribution for theta (can use nimble to simulate data)
set.seed(0)
simulate(pump,"theta") # takes current alpha and betas to create new thetas
print(pump$theta)

# calculate the log probabilities (log posterior distribution)
pump$calculate(pump$getDependencies(c("theta")))
# can use function like you would use a function in R


mcmc.out <- nimbleMCMC(code = pumpCode, constants = pumpConsts,
                       data = pumpData, inits = pumpInits, 
                       monitors = c("alpha","beta","theta"),
                       nchains = 2, niter = 10000,thin=1,nburnin=2000,
                       samplesAsCodaMCMC = TRUE,
                       summary = TRUE, WAIC = TRUE)

# Compile the model
Cpump <- compileNimble(pump,showCompilerOutput = TRUE)

# smart to compile model first if you plan to run the model multiple times
# as opposed to the code above, where it is compiled & run in the same function call
mcmc.out <- nimbleMCMC(model=Cpump,
                       monitors = c("alpha","beta","theta"),
                       nchains = 2, niter = 10000, thin=1,nburnin=2000,
                       samplesAsCodaMCMC = TRUE,
                       summary = TRUE, WAIC = TRUE)

print(str(mcmc.out))
mcmc.out$summary


pumpConf <- configureMCMC(pump, print = TRUE)

