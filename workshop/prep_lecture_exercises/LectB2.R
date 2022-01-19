library(stats4)

# Global variables (not good practice to use them)
Nages <- 11
Nyear <- 20
Ages <- seq(from=0,to=Nages-1)

# Set the variables 
SigmaC <- 0.1 # increases relative error CIs
SigmaI <- 0.2 # increases relative error CIs
SigmaR <- 0.8 # worse predictions, increases relative error CIs
M <- 0.2

# the negative log-likelihood 
minusLL <- function(pars, Type=1)
 {
  # Define the matrices and vectors (N goes beyond Nyear)
  N <- matrix(0,nrow=Nyear+1,ncol=Nages)
  CAAPred <- matrix(0,nrow=Nyear,ncol=Nages)
  Nvul <- rep(0,Nyear)
  Expl <- rep(0,Nyear)
  
  # Fill in the first row and first column of the N matrix
  for (Iages in 1:Nages) N[1,Iages]  <- exp(pars[Iages])*exp(-M*(Iages-1))
  for (Iyear in 2:Nyear) N[Iyear,1] <- exp(pars[Nages+Iyear-1])

  # set an offset so it easy to find the rest of the parameters
  Offset <- Nages+Nyear-1
  
  # Compute selectivity as a function of age
  Sel50 <- pars[Offset+1]
  Sel95 <- pars[Offset+2]
  Sel <- 1.0/(1.0+exp(-log(19)*(Ages-Sel50)/(Sel95-Sel50)))

  # Fill in the N-matrix
  for (Iyear in 1:Nyear)
   {
    # Find vulnerable numbers and the exploitation rate (never exceeds 0.99). Vulnerable numbers are halfway through the year
    Nvul[Iyear] <- sum(N[Iyear,]*Sel)*exp(-M/2)
    Expl[Iyear] <- min(C[Iyear]/Nvul[Iyear],0.99)
    
    # Now update the numbers (no plus-group)
    for (Iage in 1:Nages) CAAPred[Iyear,Iage] <- N[Iyear,Iage]*Sel[Iage]*Expl[Iyear]*exp(-M/2)
    for (Iage in 2:Nages) N[Iyear+1,Iage] <- (N[Iyear,Iage-1]*exp(-M/2)-CAAPred[Iyear,Iage-1])*exp(-M/2)

   }  
  
  # material to return
  Outs <- NULL
  
  # How to proceed depends on "Type" (1 computes the Negloglikehood; 2=computes the predicted values)
  if (Type==1 || Type==2)
   {
    
    # Negative Log-likelihood for the index data
    NegLogLike <- 0
    for (Iyear in 1:Nyear)
      NegLogLike <- NegLogLike - dnorm(I[Iyear],Nvul[Iyear],SigmaI,TRUE) 
    
    # Negative Log-likelihood for the CAA data
    for (Iyear in 1:Nyear)
     for (Iage in 1:Nages)
      NegLogLike <- NegLogLike - dnorm(CAA[Iyear,Iage],CAAPred[Iyear,Iage],SigmaC,TRUE) 

    # Save the results
    Outs$N <- N
    Outs$Nvul <- Nvul
    Outs$CAA <- CAAPred
    #if (Type==1) cat(NegLogLike,"\n")

   }  

  # How to proceed depends on "Type" (3 generates a new data set)
  if (Type==3)
   {
    #plot(1:20,Nvul,type="l",ylim=c(0,max(Nvul)*1.2))
    Nvul <- Nvul*exp(rnorm(Nyear,0,SigmaI)-SigmaI^2/2)
    #points(1:20,Nvul,pch=16)
    for (Iyear in 1:Nyear) CAAPred[Iyear,] <- CAAPred[Iyear,]*exp(rnorm(Nages,0,SigmaI)-SigmaI^2/2)
    Outs$N <- N
    Outs$Nvul <- Nvul
    Outs$CAA <- CAAPred
   }  

  if (Type==1) return(NegLogLike)
  if (Type==2 || Type == 3) return(Outs)
 }  

# Generate process error
set.seed(78671)

# Create a psuedo set
truepars <- c(rnorm(Nyear+Nages-1,0,SigmaR),5,9)
C <- c(0.1,0.15,0.2,0.25,0.3,0.35,0.5,0.55,0.7,0.7,0.7,0.7,0.6,0.2,0.1,0.1,0.1,0.2,0.2,0.2)
I <- rep(0,Nyear); CAA <- matrix(0,nrow=Nyear,ncol=Nages)
True <- minusLL(truepars,Type=2)

# I want to test the estimator (based on 200 data sets)
Nsim <- 200 # increasing sims reduces error
RelErrors <- matrix(0,nrow=Nsim,ncol=Nyear)
for (Isim in 1:Nsim)
 {  
  print(Isim)
  # Create a psuedo set (note that regenerates the process error)
  truepars <- c(rnorm(Nyear+Nages-1,0,SigmaR),5,9)
  C <- c(0.1,0.15,0.2,0.25,0.3,0.35,0.5,0.55,0.7,0.7,0.7,0.7,0.6,0.2,0.1,0.1,0.1,0.2,0.2,0.2)
  I <- rep(0,Nyear); CAA <- matrix(0,nrow=Nyear,ncol=Nages)
  True <- minusLL(truepars,Type=2)
  
  # Generate a dataset and save the results
  Generated <- minusLL(truepars,Type=3)
  I <- Generated$Nvul
  CAA <- Generated$CAA
 
  # Use optim to fit the model
  pars <- c(rep(0,Nyear+Nages-1),5,9)
  fit <- optim(pars,minusLL,Type=1)
  Estimated <- minusLL(fit$par,Type=2)

  # Store relative errors (percentages)
  RelErrors[Isim,] <- (Estimated$Nvul-True$Nvul)/True$Nvul*100
 }  

# Now summarize the results as relative error distributions
Quants <- matrix(0,nrow=5,ncol=Nyear)
for (Iyear in 1:Nyear)
 Quants[,Iyear] <- quantile(RelErrors[,Iyear],prob=c(0.05,0.25,0.5,0.75,0.95))  
ymax <- max(abs(Quants),10)

# plot the relative error distributions
par(mfrow=c(1,1),oma=c(1,1,1,1))
plot(1:20,Quants[3,],xlab="Year",ylab="Relative error (%)",ylim=c(-ymax,ymax))
xx <- c(1:20,20:1)
yy <- c(Quants[1,],rev(Quants[5,]))
polygon(xx,yy,col="gray10")
yy <- c(Quants[2,],rev(Quants[4,]))
polygon(xx,yy,col="gray90")
lines(1:20,Quants[3,],lwd=4,col="red")
abline(h=0,lty=2,lwd=2,col="blue")



