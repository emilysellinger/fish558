library(stats4)

# the negative log-likelihood for example 1 in production model where no info on biology and all fleets have
# same selectivity area
minusLL <- function(logK,logr, Type=1)
{
  # Extract the parameters
  K <- exp(logK); r <- exp(logr); Nyear <- length(C); B <- rep(0,Nyear+1)
  
  # Calculate the biomass series
  B[1] <- K
  for (Iyear in 1:Nyear)
   {
    B[Iyear+1] <- B[Iyear] + r*B[Iyear]*(1.0-B[Iyear]/K) - C[Iyear]
    if (B[Iyear] < 0.01) B[Iyear] <- 0.01
   }  

  # ML estimate for q
  qhat <- 0; bot <- 0;
  for (Iyear in 1:Nyear)
   { qhat <- qhat + log(I[Iyear]/B[Iyear]); bot <- bot + 1 }  
  qhat <- exp(qhat/bot)
  
  # Compute the sum of squared errors and hence the ML estimate for Sigma
  SS <-  0;
  for (Iyear in 1:Nyear)
    SS <- SS + (log(I[Iyear]) - log(qhat*B[Iyear]))^2
  Sigma <- sqrt(SS/Nyear)
  
  # Concentrated likelihood
  NegLogLike <- Nyear*log(Sigma);
  #cat(r,K,NegLogLike,"\n")

  # material to return
  Outs <- NULL
  Outs$B <- B
  Outs$CpueHat <-B*qhat
  Outs$r <- r
  Outs$K <- K
  
  if (Type==1) return(NegLogLike)
  if (Type==2) return(Outs)
 }  

# the negative log-likelihood for age structured model
minusLL2 <- function(logMSY,logr, Type=1)
{
  # Extract the parameters
  Nyear <- length(C); B <- rep(0,Nyear+1); MSY <- exp(logMSY); r <- exp(logr);
  K <- 4*MSY/r # this is from the logistic growth equation (b/c biomass at MSY - max sustainable yield - is K/2)

  # biomass series
  B[1] <- K
  for (Iyear in 1:Nyear)
  {
    B[Iyear+1] <- B[Iyear] + r*B[Iyear]*(1.0-B[Iyear]/K) - C[Iyear]
    if (B[Iyear] < 0.01) B[Iyear] <- 0.01
  }  
  
  # ML estimate for q
  qhat <- 0; bot <- 0;
  for (Iyear in 1:Nyear)
  { qhat <- qhat + log(I[Iyear]/B[Iyear]); bot <- bot + 1 }  
  qhat <- exp(qhat/bot)
  
  # Concentrated likelihood
  SS <-  0;
  for (Iyear in 1:Nyear)
    SS <- SS + (log(I[Iyear]) - log(qhat*B[Iyear]))^2
  Sigma <- sqrt(SS/Nyear)
  
  NegLogLike <- Nyear*log(Sigma);
  #cat(r,K,NegLogLike,"\n")
  
  # material to return
  Outs <- NULL
  Outs$B <- B
  Outs$CpueHat <-B*qhat
  Outs$r <- r
  Outs$K <- K
  
  if (Type==1) return(NegLogLike)
  if (Type==2) return(Outs)
}  

par(mfrow=c(2,2))


# Catch and effort data
C <- c(15.9,25.7,28.5,23.7,25.0,33.3,28.2,19.7,17.5,19.3,21.6,23.1,22.5,22.5,23.6,29.1,14.4,13.2,28.4,34.6,37.5,25.9,25.3) 
I <- c(61.89,78.98,55.59,44.61,56.89,38.27,33.84,36.13,41.95,36.63,36.33,38.82,34.32,37.64,34.01,32.16,26.88,36.61,30.07,30.75,23.36,22.36,21.91)
Yr <- seq(from=1967,to=1989)


# Part 1 (fit the model and plot diagnostics)
# ===========================================

# fit the model
start <- list(logK=log(400),logr=log(0.3))
fixed <- list(Type=1)
mleOutput <- mle(minusLL,start=start,fixed=fixed)
print(summary(mleOutput))

# Extract the model parameters
logK <- coef(mleOutput)[1]
logr <- coef(mleOutput)[2]
Output <- minusLL(logK,logr,Type=2)

# Plot the fit
plot(c(Yr,1990),Output$B,ylim=c(0,max(Output$B)*1.05),type="l",xlab="Year",ylab="Biomass ('000t)")
plot(c(Yr,1990),Output$CpueHat,ylim=c(0,max(I)*1.05),type="l",xlab="Year",ylab="Biomass ('000t)")
points(Yr,I,pch=16)

# Part 2 (profile likelihood for MSY)
# ===================================

# fit the model
start <- list(logMSY=log(400*0.3/2),logr=log(0.01))
fixed <- list(Type=1)
mleOutput <- mle(minusLL2,start=start,fixed=fixed)
print(summary(mleOutput))

# create likelihood profile for MSY (means MSY needs to vary)
BestLL <- -as.numeric(logLik(mleOutput))

MSYs <- seq(from=10,to=30,by=0.1)
Nprof <- length(MSYs)
NegLogLiks <- rep(0,Nprof)
for (IMSY in 1:Nprof)
{
  # fit the model but FIX logSigma
  start <- list(logr=log(0.01))
  fixed <- list(logMSY=log(MSYs[IMSY]), Type = 1)
  mleOutput <- mle(minusLL2,start=start,fixed=fixed)
  NegLogLiks[IMSY]<- -as.numeric(logLik(mleOutput))-BestLL
} 
print(NegLogLiks)

plot(MSYs,NegLogLiks,xlab="MSY",ylab="Negative log-likelihood",type="l",lty=1,ylim=c(0,10))
# Indicate 95% CI
abline(h=3.84/2)
points(MSYs,0,pch=16)
