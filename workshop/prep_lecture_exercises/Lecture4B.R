

library(stats4)

C <- NULL; I <- NULL

# the negative log-likelihood 
minusLL <- function(logK,logr, Type=1)
{
  # Extract the parameters
  K <- exp(logK); r <- exp(logr); Nyear <- length(C); B <- rep(0,Nyear+1)
  
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
  
  SS <-  0;
  for (Iyear in 1:Nyear)
    SS <- SS + (log(I[Iyear]) - log(qhat*B[Iyear]))^2
  Sigma <- sqrt(SS/Nyear)
  
  NegLogLike <- Nyear*log(Sigma);
  #cat(r,K,NegLogLike,"\n")

  Outs <- NULL
  Outs$B <- B
  Outs$CpueHat <-B*qhat
  Outs$r <- r
  Outs$K <- K
  
  if (Type==1) return(NegLogLike)
  if (Type==2) return(Outs)
 }  

# the negative log-likelihood 
minusLL2 <- function(logMSY,logr, Type=1)
{
  # Extract the parameters
  MSY <- exp(logMSY); r <- exp(logr); Nyear <- length(C); B <- rep(0,Nyear+1)
  K <- 4*MSY/r
  
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
  
  SS <-  0;
  for (Iyear in 1:Nyear)
    SS <- SS + (log(I[Iyear]) - log(qhat*B[Iyear]))^2
  Sigma <- sqrt(SS/Nyear)
  
  NegLogLike <- Nyear*log(Sigma);
  #cat(r,K,NegLogLike,"\n")
  
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
Yr <- seq(from=1965,to=1987)


# Part 1 (fit the model and plot diagnostics)
# ===========================================

# fit the model
start <- list(logK=log(400),logr=log(0.3))
fixed<- list(Type=1)
mleOutput <- mle(minusLL,start=start,fixed=fixed)
print(summary(mleOutput))

# Extract the model parameters
logK <- coef(mleOutput)[1]
logr <- coef(mleOutput)[2]
Output <- minusLL(logK,logr,Type=2)

plot(c(Yr,1988),Output$B,ylim=c(0,max(Output$B)*1.05),type="l",xlab="Year",ylab="Biomass ('000t)")
plot(c(Yr,1988),Output$CpueHat,ylim=c(0,max(I)*1.05),type="l",xlab="Year",ylab="Biomass ('000t)")
points(Yr,I,pch=16)

# Part 2 (profile likelihood for MSY)
# ===================================

# fit the model
start <- list(logMSY=log(400*0.3/2),logr=log(0.3))
mleOutput <- mle(minusLL2,start=start,fixed=fixed)
print(summary(mleOutput))

# Extract the model parameters
logMSY <- coef(mleOutput)[1]
logr <- coef(mleOutput)[2]
Output <- minusLL2(logMSY,logr,Type=2)
print(Output)
print(length(Output$B))

plot(c(Yr,1988),Output$B,ylim=c(0,max(Output$B)*1.05),type="l",xlab="Year",ylab="Biomass ('000t)")
plot(c(Yr,1988),Output$CpueHat,ylim=c(0,max(I)*1.05),type="l",xlab="Year",ylab="Biomass ('000t)")
points(Yr,I,pch=16)





