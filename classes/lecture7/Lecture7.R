library(runjags)
Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.0")
library("rstan")
library(devtools)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(nimble)

# ==========================================================================================
# ==========================================================================================


DoFitJAGS <- function(C,I,nsample=5000,thin=5,Plot=F,silent=T)
{
  model <- "model 
  {
  
   #time step [1] conditions (note: T for truncation)
   Pmed[1] <-0
   P[1]~dlnorm(Pmed[1], isigma2)T(0.05,1.6)
  
   #time steps of model  
   for( t in 2 : N )
   {
    Pmed[t] <- log(max(P[t - 1] + (r * P[t - 1]) * (1 - P[t - 1]) - C[t - 1] / K, 0.001) )	
    FF[t-1] <- C[t-1]/P[t-1]/K
    P[t] ~ dlnorm(Pmed[t],isigma2)T(0.05,1.5)
   }
  
   # Likelihood
   for( t in 1 : N )
    {
     Imed[t] <- log((q * K) * P[t])
     I[t] ~ dlnorm(Imed[t],itau2)
    
     # posterior predictions of the index (mean)
     index[t]<-(q*K*P[t])
     # posterior predictive distrubution (hint, the parameterization of dlnorm is not the same as in R)
     I.new[t]~dlnorm(log(index[t]), itau2)
    }
  
   #priors 
   r ~ dlnorm( -1.38, 3.845)T(0.01,1.2)
   isigma2 ~ dgamma(3.785,0.0102)
   itau2 ~ dgamma(1.709,0.00861)
   iq ~ dgamma(0.001,0.001)T( 0.5,100)
   K ~ dlnorm(5.0429,3.7603)T(10,1000)
 
   sigma2 <- 1/isigma2
   tau2 <- 1/itau2
   q <- 1/iq
   sigma <- sqrt(sigma2)
   tau <- sqrt(tau2)
  
   #additional parameters and preditions
   MSP <-  r*K/4
   EMSP <-  r/(2*q)
   P1990 <-  P[N] + r*P[N]*(1-P[N]) - C[N]/K
   B1990 <-  P
 }"
  
 cat("calling JAGS\n")
 N <- length(C)
 data <- list(C=C,I=I,N=N)
 
 set.seed(1901)
 
 #initial values
 inits1 <- list(r=0.8, K=200, iq=0.5, isigma2=100, itau2=100, P=c(0.99,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60,0.58,0.56))
 inits2 <- list(r=0.5, K=300, iq=0.8, isigma2=200, itau2=200, P=c(0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99))
 inits <- list(inits1,inits2)
 
 # Run Jags
 results <- run.jags(model=model, monitor=c("r","K","q","sigma","tau","P","FF","I.new"), 
                     data=data, n.chains=2, method="rjags", inits=inits,
                     plots=T,sample=nsample,thin=thin,silent.jags=silent)
 
 if (Plot==T) print(results)
 
 # Combine results from multiple chains
 mcmc <- rbind(results$mcmc[[1]],results$mcmc[[2]])
 
 # Plots
 if (Plot==T)
  {
   ResSum <- matrix(0,ncol=23,nrow=5)
   for (Iyr in 1:23)
    ResSum[,Iyr] <- quantile(mcmc[,Iyr+5],probs=c(0.05,0.25,0.5,0.75,0.95))  
   xx <- seq(1:23)
   plot(xx,ResSum[3,],xlab="Year",ylab="Depletion",type='l',lwd=3,ylim=c(0,1.3))
   xx2 <- c(xx,rev(xx))
   polygon(xx2,c(ResSum[1,],rev(ResSum[5,])),col="gray50")
   polygon(xx2,c(ResSum[2,],rev(ResSum[4,])),col="gray95")
   lines(xx,ResSum[3,],lwd=3,lty=1)
 
   # Plot the posterior predictive distribution
   plot(xx,I,pch=16,ylim=c(0,100))
   ResSum <- matrix(0,ncol=23,nrow=5)
   for (Iyr in 1:23)
    ResSum[,Iyr] <- quantile(mcmc[,50+Iyr],probs=c(0.05,0.25,0.5,0.75,0.95))  
   lines(xx,ResSum[3,],lwd=3,lty=1) 
   lines(xx,ResSum[1,],lwd=1,lty=2) 
   lines(xx,ResSum[5,],lwd=1,lty=2) 
 }
 
 # Results
 res <- apply(mcmc,2,median)
 #print(res)
 return(res)
}

# ====================================================================================================

DoFitStan <- function(C,I,nsample=5000,thin=5,Plot=F,silent=T)
{
  cat("calling Stan\n")
  N <- length(C)
  tuna_data <- list(C=C,I=I,N=N)
  
  #set.seed(1901)
  
  #initial values
  inits1 <- list(r=0.8, K=200, iq=0.5, isigma2=100, itau2=100, P=c(0.99,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60,0.58,0.56))
  inits2 <- list(r=0.8, K=200, iq=0.5, isigma2=100, itau2=100, P=c(0.99,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60,0.58,0.56))
  inits <- list(inits1,inits2)
  
  fit <- stan(file = 'Schaefer.stan', data = tuna_data, 
              iter = nsample, chains = 2,init=inits, warmup=0.1*nsample,verbose=F)
  #print(fit)
  la <- extract(fit, permuted = TRUE) # return a list of arrays 
 
  if (Plot==T)
   {
    quants <- matrix(0,nrow=3,ncol=N)
    for (II in 1:N)
     quants[,II] <- quantile(la$P[,II],probs=,c(0.05,0.5,0.95))
    ymax <- max(quants)*1.05
    plot(1:N,quants[2,],xlab="Year",ylab="Depletion",type="b",lty=1,pch=16,ylim=c(0,ymax))
    lines(1:N,quants[1,],lwd=1,lty=2)
    lines(1:N,quants[3,],lwd=1,lty=2)
  
    quants <- matrix(0,nrow=3,ncol=N)
    for (II in 1:N)
     quants[,II] <- quantile(la$Inew[,II],probs=,c(0.05,0.5,0.95))
    ymax <- max(quants)*1.05
    plot(1:N,I,xlab="Year",ylab="Index",type="b",lty=1,pch=16,ylim=c(0,ymax))
    lines(1:N,quants[2,],lwd=2,lty=1)
    lines(1:N,quants[1,],lwd=1,lty=2)
    lines(1:N,quants[3,],lwd=1,lty=2)
    #print(str(la))
   }  
  
  res <- NULL
  res <- c(median(la$r),median(la$K),median(la$q),median(la$sigma),median(la$tau))
  res <- c(res,apply(la$P,2,median))
  res <- c(res,apply(la$FF[,1:(N-1)],2,median))
  res <- c(res,apply(la$Inew,2,median))
  print(res)
  return(res)
}

# ====================================================================================================
# ====================================================================================================

DoFitNimble <- function(C,I,nsample=5000,thin=5,Plot=F,silent=T)
{
  
 TheCode <- nimbleCode({ 
  #time step [1] conditions (note: T for truncation)
  Pmed[1] <- 0
  #P[1] ~ T(dlnorm(Pmed[1], isigma2),0.05,1.6)
  P[1] ~ dlnorm(Pmed[1], isigma2)
  
  #time steps of model  
  for( t in 2 : N )
  {
    Pmed[t] <- log(max(P[t - 1] + (r * P[t - 1]) * (1 - P[t - 1]) - C[t - 1] / K, 0.001) )	
    #P[t] ~ T(dlnorm(Pmed[t],isigma2),0.05,1.5) 
    FF[t-1] <- C[t-1]/(P[t-1]*K)
    P[t] ~ dlnorm(Pmed[t],isigma2)
  }
  
  # Likelihood
  for( t in 1 : N )
  {
    Imed[t] <- log((q * K) * P[t])
    I[t] ~ dlnorm(Imed[t],itau2)
    
    #posterior predictions (hint, the parameterization of dlnorm is not the same as in R)
    index[t]<-(q*K*P[t])
    I.new[t]~dlnorm(log(index[t]), itau2)
  }
  
  #priors 
  r ~ T(dlnorm( -1.38, 3.845),0.01,1.2)
  isigma2 ~ dgamma(3.785,0.0102)
  itau2 ~ dgamma(1.709,0.00861)
  iq ~ T(dgamma(0.001,0.001),0.5,100)
  K ~ T(dlnorm(5.0429,3.7603),10,1000)
  
  sigma2 <- 1/isigma2
  tau2 <- 1/itau2
  q <- 1/iq
  sigma <- sqrt(sigma2)
  tau <- sqrt(tau2)
  
  #additional parameters and preditions
  MSP <-  r*K/4
  EMSP <-  r/(2*q)
  P1990 <-  P[N] + r*P[N]*(1-P[N]) - C[N]/K
  B1990 <-  P[N]
 })

 cat("calling Nimble\n")
 N <- length(C)
 data <- list(C=C,I=I,N=N)
 
 TheConsts <- list("N" = N, "C" = C)
 TheData <- list("I" = I)
 
 #initial values
 inits1 <- list(r=0.8, K=200, iq=0.5, isigma2=100, itau2=100, P=c(0.99,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60,0.58,0.56))
 inits2 <- list(r=0.5, K=300, iq=0.8, isigma2=200, itau2=200, P=c(0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99,0.99))
 TheInits <- list(inits1,inits2)

 # Create the model
 #test <- nimbleModel(code = TheCode, name = "code", constants = TheConsts,
#                     data = TheData, inits = inits1)
 
 # Conduct MCMC sampling
 mcmc.out <- nimbleMCMC(code = TheCode, constants = TheConsts,
                        data = TheData, inits = TheInits,
                        nchains = 2, niter = nsample,thin=thin,nburnin=0.1*nsample,
                        #                       nchains = 2, niter = 2000,
                        summary = TRUE, WAIC = TRUE,
                        samplesAsCodaMCMC=TRUE,progressBar=FALSE,
                        monitors = c('r','K','q','sigma','tau','P','FF','I.new'))
 
 # Output statistics
 mcmc <- rbind(mcmc.out$samples$chain1,mcmc.out$samples$chain2)
 
 if (Plot==T)
  {
   Index <-which(colnames(mcmc)=="P[1]")
   ResSum <- matrix(0,ncol=23,nrow=5)
   for (Iyr in 1:23)
    ResSum[,Iyr] <- quantile(mcmc[,Iyr+Index-1],probs=c(0.05,0.25,0.5,0.75,0.95))  
   xx <- seq(1:23)
   plot(xx,ResSum[3,],xlab="Year",ylab="Depletion",type='l',lwd=3,ylim=c(0,1.3))
   xx2 <- c(xx,rev(xx))
   polygon(xx2,c(ResSum[1,],rev(ResSum[5,])),col="gray50")
   polygon(xx2,c(ResSum[2,],rev(ResSum[4,])),col="gray95")
   lines(xx,ResSum[3,],lwd=3,lty=1)
 
   # Posterior predictive
   Index <-which(colnames(mcmc)=="I.new[1]")
   plot(xx,I,pch=16,ylim=c(0,100))
   ResSum <- matrix(0,ncol=23,nrow=5)
   for (Iyr in 1:23)
    ResSum[,Iyr] <- quantile(mcmc[,Iyr+Index-1],probs=c(0.05,0.25,0.5,0.75,0.95))  
   lines(xx,ResSum[3,],lwd=3,lty=1) 
   lines(xx,ResSum[1,],lwd=1,lty=2) 
   lines(xx,ResSum[5,],lwd=1,lty=2) 
 }
 
 res <- NULL
 res <- median(mcmc[,which(colnames(mcmc)=="r")])
 res <- c(res,median(mcmc[,which(colnames(mcmc)=="K")]))
 res <- c(res,median(mcmc[,which(colnames(mcmc)=="q")]))
 res <- c(res,median(mcmc[,which(colnames(mcmc)=="sigma")]))
 res <- c(res,median(mcmc[,which(colnames(mcmc)=="tau")]))
 Index <-which(colnames(mcmc)=="P[1]")
 for (Iyr in 1:23)
  res <- c(res,median(mcmc[,Iyr+Index-1]))
 Index <-which(colnames(mcmc)=="FF[1]")
 for (Iyr in 1:22)
  res <- c(res,median(mcmc[,Iyr+Index-1]))
 Index <-which(colnames(mcmc)=="I.new[1]")
 for (Iyr in 1:23)
   res <- c(res,median(mcmc[,Iyr+Index-1]))
 
 print(res)
 return(res)
}

# ====================================================================================================
# ==========================================================================================


DoSims <- function(Nsim=100,Type="JAGS",r,q,K,sigma,tau,FF,Plot=F,FileOut)
{
 # Space for results
 ErrMat <- matrix(0,nrow=Nsim,ncol=Noutput)

 set.seed(7810)
 Nyear <- 23
 for (Isim in 1:Nsim)
  {
   cat("Running simulation",Isim,"\n")
   B <- rep(0,Nyear+1)  
   # Project forward with process error
   B[1] <- rlnorm(1,0,sigma)*K
   for (Iyear in 1:Nyear)
    {
     Bdet <- as.numeric(B[Iyear] + r*B[Iyear]*(1.0-B[Iyear]/K) - B[Iyear]*FF[Iyear])
     B[Iyear+1] <- rlnorm(1,log(Bdet),sigma)
     C[Iyear] <-  B[Iyear]*FF[Iyear]
     I[Iyear] <- rlnorm(1,log(q*B[Iyear]),tau)
   }
   
   if (Type=="JAGS") Estim <- DoFitJAGS(C,I,nsample=1000,thin=1,Plot=F)
   if (Type=="Stan") Estim <- DoFitStan(C,I,nsample=1000,thin=1,Plot=F)
   if (Type=="Nimble") Estim <- DoFitNimble(C,I,nsample=1000,thin=1,Plot=F)
   
    ErrMat[Isim,] <- (Estim-Opmodel)/Opmodel 
  }

 # Plot histograms of results
 if (Plot==T)
  {pdf(file = "Lecture7_plots.pdf")
   par(mfrow=c(5,5))
   UseStats <- c(1:28)
   for (II in UseStats)  
    {
     hist(ErrMat[,II],main=Stats[II])
     abline(v=0,lty=2,lwd=4,col="red")  
   }
   dev.off()
  } 

 # Summary stat
 colnames(ErrMat) <- Stats
 Noutputs <- length(ErrMat[1,])
 Bias <- rep(0,Noutputs)
 names(Bias) <- Stats
 MARE <- rep(0,Noutputs)
 names(MARE) <- Stats
 for (Istat in 1:Noutputs)
 {
  Bias[Istat] <- mean(ErrMat[,Istat])
  MARE[Istat] <- median(abs(ErrMat[,Istat]))
 }
 #print(summary(ErrMat))
 ResFile <- cbind(Stats,Bias,MARE)
 #print(ResFile)
 write(t(ResFile),FileOut,ncol=3)
 

 ReturnObj <- NULL
 ReturnObj$Bias <- Bias
 ReturnObj$MARE <- MARE
 return(ReturnObj)
}

# ==========================================================================================
# ==========================================================================================

# Fit the operating model and extract the parameters (which can be changed)
C<-c(15.9,25.7,28.5,23.7,25.0,33.3,28.2,19.7,17.5,19.3,21.6,23.1,22.5,22.5,23.6,29.1,14.4,13.2,28.4,34.6,37.5,25.9,25.3) 
I<-c(61.89,78.98,55.59,44.61,56.89,38.27,33.84,36.13,41.95,36.63,36.33,38.82,34.32,37.64,34.01,32.16,26.88,36.61,30.07,30.75,23.36,22.36,21.91)

par(mfrow=c(3,2))
Opmodel <- DoFitJAGS(C,I,nsample=50000,thin=1,Plot=F,silent=F)
K <- Opmodel[2]; r <- Opmodel[1]; q <- Opmodel[3]; sigma <- Opmodel[4]; tau <- Opmodel[5]
FF <- c(Opmodel[29:50],Opmodel[50])
Noutput <- length(Opmodel)
print(Opmodel)
Stats <- names(Opmodel)

Opmodel <- DoFitStan(C,I,nsample=1000,thin=1,Plot=F,silent=F)
print(Opmodel)
#Opmodel <- DoFitNimble(C,I,nsample=1000,thin=1,Plot=F,silent=F)
#print(Opmodel)
#AA


# note, make sure to do plot = T 
# Do the actual simulations
Resu1 <- DoSims(Nsim=20,Type="JAGS",r=r,K=K,q=q,sigma=sigma,tau=tau,FF=FF,Plot=T,FileOut="Save1.out")


# ideally you want to run for more than 10 simulations and more the nsamples
# if you increase the number of sims, it will take longer to run it, but not to recompile it

