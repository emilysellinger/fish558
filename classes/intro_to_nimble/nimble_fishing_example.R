library(nimble)
library(igraph)
library(coda)
library(lattice)

set.seed(1991)

TrueSpatial <- c(0,1,2,3,4,5)*0.4
SigmaProc <- 0.25
TempPar <- 0.25

# Generate data ============================================================

GenData <- function(Ndata=100,Nind=10)
 {
  Data <- NULL
  Nspatial <- length(TrueSpatial)
  Individual <- matrix(NA,nrow=Nind,ncol=Nspatial)
  
  # Loop over individuals
  for (Ind in 1:Nind)
   {  
    # Inidvidual effects
    Individual[Ind,] <- rnorm(Nspatial,0,SigmaProc)
    for (Idata in 1:(Ndata/Nind))
     {
      Temp <- rnorm(Nspatial,0,3)
      LinearPred <- exp(TrueSpatial + Temp*TempPar + Individual[Ind,])
      LinearPred <- LinearPred/sum(LinearPred)
      Actual <- rmultinom(1,1:Nspatial,LinearPred)
      Actual <- which(Actual==1)
      Out <- c(Ind,Idata,Actual,Temp)
      Data <- rbind(Data,Out)
     }
  }  
  colnames(Data) <- c("Ind","Repl","Selected",rep(paste("Temp",1:Nspatial,sep="")))
  print(head(Data))
  
  Outs <- NULL
  Outs$Data <- Data
  Outs$Individial <- Individual
  Outs$TempPar <- TempPar
  Outs$TrueSpatial <- TrueSpatial
  print(str(Outs))
  return(Outs)
}

# ============================================================

# generate the data
Ndata <- 1000; Nind <- 10
GenStuff <- GenData(Ndata=Ndata,Nind=Nind)  
Data <- GenStuff$Data
#print(Data) 

 # ============================================================
 # Estimation method
 # ============================================================


 TheCode <- nimbleCode({ 
   
   for(k in 1:Ncat){
     for(i in 1:Nind){
       for(j in 1:Ndata)
     }
   }
   Ptilda = exp(alpha + beta*TempPar + gamma)
   
   
   
   alpha ~ dnorm(0,1000)
   alpha[1] <- 0
   beta ~ dnorm(0,1000)
   gamma ~ dnorm(0, sigma)
   sigma ~ dunif(0,5)
 })

 TheConsts <- list("Ncat" = 6, "Ndata" = Ndata, "Nind" = Nind, "IndexI" = Data[,1])
 #print(TheConsts)
                
 TheData <- list("Actual" = Data[,3],"Temp"=Data[,4:9])
 #print(TheData)
 
 #initial values
 inits1 <- list(Individual=matrix(rnorm(Nind*6,0,1),nrow=Nind,ncol=6),SDInd=0.1,TrueSpatial=rnorm(5,0,1),TempPar=0,Pred=matrix(rnorm(Ndata*6,0,1),nrow=Ndata,ncol=6))
 inits2 <- list(Individual=matrix(rnorm(Nind*6,0,1),nrow=Nind,ncol=6),SDInd=2.0,TrueSpatial=rnorm(5,0,1),TempPar=0,Pred=matrix(rnorm(Ndata*6,0,1),nrow=Ndata,ncol=6))
 TheInits <- list(inits1,inits2)
 print(str(inits1))

 # Create the model
 test <- nimbleModel(code = TheCode, name = "code", constants = TheConsts, data = TheData, inits = inits1)

 # Apply MCMC
 mcmc.out <- nimbleMCMC(code = TheCode, constants = TheConsts,
                        data = TheData, inits = TheInits, 
                        monitors = c("Individual","SDInd","TrueSpatial","TempPar"),
                        nchains = 2, niter = 100,thin=8,nburnin=20,
                        samplesAsCodaMCMC = TRUE,
                        summary = TRUE, WAIC = TRUE)
 
 print(mcmc.out$summary$all.chains)

 # Gelman-Rubin
 #print(gelman.diag(mcmc.out$samples))

 # graphs
 par(mfrow=c(3,3))
 Index <- which(colnames(mcmc.out$samples$chain1)=="SDInd")
 Ncycle <- length(mcmc.out$samples$chain1[,1])
 plot(1:Ncycle,mcmc.out$samples$chain1[,Index],xlab="cycle",ylab="SD",col="red",lty=1,type='l')
 lines(1:Ncycle,mcmc.out$samples$chain2[,Index],col="blue")
 
 Index <- which(colnames(mcmc.out$samples$chain1)=="TempPar")
 Ncycle <- length(mcmc.out$samples$chain1[,1])
 plot(1:Ncycle,mcmc.out$samples$chain1[,Index],xlab="cycle",ylab="Tempar",col="red",lty=1,type='l')
 lines(1:Ncycle,mcmc.out$samples$chain2[,Index],col="blue")
 
 Index <- which(colnames(mcmc.out$samples$chain1)=="TrueSpatial[1]")
 Ncycle <- length(mcmc.out$samples$chain1[,1])
 plot(1:Ncycle,mcmc.out$samples$chain1[,Index],xlab="cycle",ylab="True Spatial 1",col="red",lty=1,type='l')
 lines(1:Ncycle,mcmc.out$samples$chain2[,Index],col="blue")
 Index <- which(colnames(mcmc.out$samples$chain1)=="TrueSpatial[2]")
 Ncycle <- length(mcmc.out$samples$chain1[,1])
 plot(1:Ncycle,mcmc.out$samples$chain1[,Index],xlab="cycle",ylab="True Spatial 2",col="red",lty=1,type='l')
 lines(1:Ncycle,mcmc.out$samples$chain2[,Index],col="blue")
 Index <- which(colnames(mcmc.out$samples$chain1)=="TrueSpatial[3]")
 Ncycle <- length(mcmc.out$samples$chain1[,1])
 plot(1:Ncycle,mcmc.out$samples$chain1[,Index],xlab="cycle",ylab="True Spatial 3",col="red",lty=1,type='l')
 lines(1:Ncycle,mcmc.out$samples$chain2[,Index],col="blue")
 Index <- which(colnames(mcmc.out$samples$chain1)=="TrueSpatial[4]")
 Ncycle <- length(mcmc.out$samples$chain1[,1])
 plot(1:Ncycle,mcmc.out$samples$chain1[,Index],xlab="cycle",ylab="True Spatial 4",col="red",lty=1,type='l')
 lines(1:Ncycle,mcmc.out$samples$chain2[,Index],col="blue")
 Index <- which(colnames(mcmc.out$samples$chain1)=="TrueSpatial[5]")
 Ncycle <- length(mcmc.out$samples$chain1[,1])
 plot(1:Ncycle,mcmc.out$samples$chain1[,Index],xlab="cycle",ylab="True Spatial 5",col="red",lty=1,type='l')
 lines(1:Ncycle,mcmc.out$samples$chain2[,Index],col="blue")
 
