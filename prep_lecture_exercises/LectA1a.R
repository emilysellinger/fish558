library(stats4)

# note mle is a wrapper for optim - mle has better syntax for fixed parameters than optim
# mle is nice for a likelihood profile where you are planning on fixing one of the parameters
# hessian  = TRUE in optim and invert it to get variance covariance matrix
Ages=NULL
Lengths=NULL

# the negative log-likelihood 
minusLL <- function(logLinf,logK,t0,logSigma)
 {
  # Extract the parameters
  Linf <- exp(logLinf); Kappa <- exp(logK); Sigma <- exp(logSigma)
  
  # make the model predictions
  Pred <- Linf*(1.0-exp(-Kappa*(Ages-t0))) # Von Bonferrony equation
  
  # Compute the negative log-likelihood
  NegLogL <- -1*sum(dnorm(Lengths,Pred,Sigma,TRUE))
  
  return(NegLogL)

 }  

 # Psuedo data (always good for testing)
 set.seed(781010)
 Nsample <- 5000
 Ages <<- sample(c(1:20),Nsample,replace=TRUE) # double arrow assigns global variable
 Lengths <<- 100*(1.0-exp(-0.2*(Ages-0.1)))+rnorm(Nsample,0,5)
 plot(Ages,Lengths,xlab="Ages",ylab="Lengths")
 
 # fit the model
 start <- list(logLinf=log(80),logK=log(0.15),t0=0,logSigma=1)
 fixed <- NULL
 mleOutput <- mle(minusLL,start=start)
 print(summary(mleOutput))
 
 # Extract the parameters
 Linf <- exp(coef(mleOutput)[1])
 Kappa <- exp(coef(mleOutput)[2])
 t0 <- coef(mleOutput)[3]
 Sigma <- exp(coef(mleOutput)[4])
 cat("parameter estimates\n")
 cat("Linf = ",Linf,"\n")
 cat("Kappa = ",Kappa,"\n")
 cat("t0 = ",t0,"\n")
 cat("Sigma = ",Sigma,"\n")
 
