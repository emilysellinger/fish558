library(mvtnorm)
library(stats4)
library(MASS)
library(coda)


# Define priors ---------------------------------------------------------------


# Likelihood function ================================================================================

NegLogLike2 <- function(pars,DataUsed,Print=F)
{
  b <- 15.5; S1 <- 0.7;SJ <- 0.8; SA <- 0.95; gamma <- 0.2
  Nint <- c(100,10,5,1) # initial numbers at each stage
  
  Nyear <- DataUsed$Nyear

  S0 <- pars[1]
  Lambda <- pars[2]
  q <- pars[3]
  EggVal <- exp(pars[4:(4+Nyear-1-1)])
  
  Effort = DataUsed$Effort
  EggObs <- DataUsed$EggObs
  SDEgg <- DataUsed$SDEgg
  AdultObs <- DataUsed$AdultObs
  SDAdult <- DataUsed$SDAdult
  CatchObs <- DataUsed$CatchObs
  SDCatch <- DataUsed$SDCatch
  Ayears <- DataUsed$Ayears
  AgeDataUse <- DataUsed$AgeDataUse
    
  N <- array(0,dim=c(Nyear+1,4))
  
  # PopModel
  for(i in 1:(Nyear+1)){
    if(i == 1){
      N[i,1] <- Nint[1]
      N[i,2] <- Nint[2]
      N[i,3] <- Nint[3]
      N[i,4] <- Nint[4]
      
    }else if(i > 1){
      N[i,1] <-  EggVal[i-1]
      N[i,2] <-  S0*N[i-1,1]
      N[i,3] <-  S1*N[i-1,2] + (1-gamma)*N[i-1,3]*SJ
      N[i,4] <-  gamma*N[i-1,3]*SJ + SA*N[i-1,4] - (1-Lambda)*q*Effort[i-1]*N[i-1,4]
    }
  }
  
  
  # Egg negative log likelihood
  EggPreds <- N[,1]
  
  Egg_NegLL <- -sum(dlnorm(EggObs, log(EggPreds[1:(length(EggPreds)-1)]), sd = SDEgg, log = TRUE))
  
  
  # Age Structure log likelihood (multinomial)
  age_structure_pred <- N[c(2,5,10),]
  
  age_struct_ll <- vector(length = 3)
  for(i in 1:3){
    age_struct_ll[i] <- -sum(dmultinom(x = AgeDataUse[i,], prob = age_structure_pred[i,], log = TRUE))
  }
  
  AgeStruc_NegLL <- sum(age_struct_ll)
  
  # Number of Adults interacting with fishery likelihood
  CatchPred <- log(q*N[-c(1),4]*Effort)
  
  Catch_NegLL <- -sum(dlnorm(CatchObs, CatchPred, sd = SDCatch, log = TRUE))
  
  # Counts of adults
  AdultPreds <- log(N[,4])
  Adult_NegLL <- -sum(dlnorm(AdultObs, AdultPreds[1:(length(EggPreds)-1)], sd = SDAdult, log = TRUE))
  
  total_ll <- sum(Egg_NegLL, AgeStruc_NegLL, Catch_NegLL, Adult_NegLL)
  
  return(total_ll)

}

# MCMC =================================================================================================================

DoMCMC<-function(Xinit,DataUsed,Ndim,covar,Nsim=1000,Nburn=0,Nthin=1)
{
  Xcurr <- Xinit
  Fcurr <- -1*NegLogLike2(Xcurr,DataUsed)
  Outs <- matrix(0,nrow=(Nsim-Nburn),ncol=(Ndim+1))
  Ipnt <- 0; Icnt <- 0
  for (Isim in 1:Nsim)
  {
    Xnext <- rmvnorm(1, mean=Xcurr, sigma=covar)
    Fnext <- -1*NegLogLike2(Xnext,DataUsed)
    Rand1 <- log(runif(1,0,1))
    if (Fnext > Fcurr+Rand1)
    {Fcurr <- Fnext; Xcurr <- Xnext }   
    if (Isim > Nburn & Isim %% Nthin == 0)
    {
      Icnt <- Icnt + 1; Outs[Icnt,] <- c(Xcurr,Fcurr); cat("saving",Icnt,"\n")     
    }
  } 
  xx <- seq(1:Icnt)
  par(mfrow=c(4,4),mar=c(3,4,2,1))
  for (II in 1:(Ndim+1))
  {
    yy <- Outs[,II][1:Icnt]
    if (II <= Ndim)
      lab1 <- paste("Parameter ",II)
    else
      lab1 <- "Posterior density"
    if (II <= Ndim) yy <- exp(yy)
    plot(xx,yy,xlab="Cycle number",ylab=lab1,type='b',pch=16,cex=0.02)
  }
  par(mfrow=c(4,4),mar=c(3,4,2,1))
  for (II in 1:(Ndim+1))
  {
    yy <- Outs[,II][1:Icnt]
    if (II <= Ndim)
      lab1 <- paste("Parameter ",II)
    else
      lab1 <- "Posterior density"
    if (II <= Ndim) yy <- exp(yy)
    hist(yy,ylab=lab1,main="")
  }
  return(Outs[1:Icnt,])
}

# ================================================================================


OFile <- "day4/Ex4a.dat"
Nyear <- scan(OFile,skip=1,n=1,quiet=T)
Effort <- scan(OFile,skip=5,n=Nyear,quiet=T)
EggsObs <- scan(OFile,skip=7,n=Nyear,quiet=T)
SDEgg <- scan(OFile,skip=9,n=1,quiet=T)
AdultObs  <-scan(OFile,skip=11,n=Nyear,quiet=T)
SDAdult <- scan(OFile,skip=13,n=1,quiet=T)
CatchObs  <-scan(OFile,skip=15,n=Nyear,quiet=T)
SDCatch <- scan(OFile,skip=17,n=1,quiet=T)
print(SDCatch)
Nayears <- scan(OFile,skip=20,n=1,quiet=T)
Ayears <- scan(OFile,skip=22,n=Nayears,quiet=T)
AgeComp <- matrix(scan(OFile,skip=24,n=Nayears*4,quiet=T),ncol=4,byrow=T)

OFile <- "day4/Ex4b.dat"
Npars <- 14
Vectors <- scan(OFile,skip=2,n=Npars,quiet=T)
SDVec <- scan(OFile,skip=4,n=Npars,quiet=T)
Varco <- matrix(scan(OFile,skip=6,Npars*Npars,quiet=T),ncol=Npars)

DataUsed <- list(EggObs=EggsObs,SDEgg=SDEgg,AdultObs=AdultObs,SDAdult=SDAdult,CatchObs=CatchObs,SDCatch=SDCatch,Ayears=Ayears,AgeDataUse=AgeComp,Effort=Effort,Nyear=Nyear)
print(DataUsed)
print(Vectors)
NegLogLike2(Vectors,DataUsed=DataUsed,Print=T)

# use
set.seed(123456)
Outs <- DoMCMC(Vectors, DataUsed, 14, Varco, Nsim = 300000, Nburn = 0, Nthin = 10)
write(t(Outs), "ParVec.out", ncolumn = 15) 

# Coda

Outs <- matrix(scan("ParVec.Out"),ncol=15,byrow=T)  
MCMC1 <- mcmc(Outs)
ggmcmc(ggs(MCMC1),file="E:\\output.pdf")

# Tabular Diagnostics
print(geweke.diag(MCMC1))
print(heidel.diag(MCMC1))
print(autocorr.diag(MCMC1))

# Graphical diagnostics
par(mfrow=c(4,4),mar=c(3,4,2,1))
traceplot(MCMC1)
autocorr.plot(MCMC1)

# 
