# Homework 4
library(mvtnorm)
library(stats4)
library(MASS)
library(coda)
library(ggplot2)
library(ggmcmc)
library(here)
library(tidyverse)

# load data ---------------------------------------------------------------
counts <- read_csv(here("homework4", "stream_count_data.csv"))
counts <- counts[,-1]
inits <- read_csv(here("homework4", "stream_initial_params.csv"))
var_covar <- read_csv(here("homework4", "stream_variance_covariance.csv"))

params <- unname(as_vector(inits[,2]))
var_covar <- unname(as.matrix(var_covar))


################################## PART A ##########################################################

# population model ----------------------------------------------------------
pop_mod <- function(mus, alphas){
  distance <- seq(0, 100, 10)
  preds <- matrix(0, nrow = 20, ncol = 11)
  sites <- nrow(preds)
  dists <- ncol(preds)
  
  for(i in 1:sites){
    for(j in 1:dists){
      preds[i, j] <- mus[i]*exp(alphas[i]*distance[j])
    }
  }
  
  return(preds)
  
}

mu_inits <- unname(as_vector(inits[6:25, 2]))
alpha_inits <- unname(as_vector(inits[26:45, 2]))

# test population model
test <- pop_mod(mu_inits, alpha_inits)

# negative log likelihood -------------------------------------------------

likelihood <- function(predicted, data){
  
  data <- as.vector(as.matrix(data))
  predicted <- as.vector(as.matrix(predicted))
  sumLL<-sum(dpois(data, predicted, log = T))
  
  return(sumLL)
}

#test likelihood
likelihood(test, counts)


# total density ---------------------------------------------------------------

PostDens<-function(pars, data){
  
  MuF <- pars[1]
  MuL <- pars[2]
  SigmaMu <- pars[3]
  alphabar <- pars[4]
  SigmaAlpha <- pars[5]
  muis<-pars[6:25]
  alphais<-pars[26:45]
  
  # get predictions
  preds <- pop_mod(mus = muis, alphas = alphais)
  
  # calculate prediction likelihood
  like <- likelihood(preds, data)
  mu_means <- rep(NA, length(muis))
  
  # calculate mu tilda i for each stream
  for (i in 1:length(muis)){
    mu_means[i]<- MuF + (((MuL-MuF)*(i-1))/19)
  }
  
  
  mu_prior<-sum(dnorm(muis, mu_means, sd = SigmaMu, log = T))
  alpha_prior<-sum(dnorm(alphais, alphabar, sd = SigmaAlpha, log = T))
  SigmaMu_prior <- dnorm(SigmaMu, 1, sd = 0.1, log = T)
  SigmaAlpha_prior<-dnorm(SigmaAlpha, 1, sd = 0.1, log = T)
  
  # calculate total density
  postdens<-like + mu_prior + alpha_prior + SigmaMu_prior + SigmaAlpha_prior
  return(postdens)
}

PostDens(params, counts)

# MCMC --------------------------------------------------------------------
DoMCMC<-function(Xinit,DataUsed,Ndim,covar,Nsim=1000,Nburn=0,Nthin=1)
{
  Xcurr <- Xinit
  Fcurr <- PostDens(Xinit, DataUsed)
  Outs <- matrix(0,nrow=(Nsim-Nburn),ncol=(Ndim+1))
  Ipnt <- 0; Icnt <- 0
  for (Isim in 1:Nsim)
  {
    Xnext <- rmvnorm(1, mean = Xcurr, sigma = covar)
    Fnext <- PostDens(Xnext, DataUsed)
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


outs <- DoMCMC(params, counts, Ndim = 45, covar = 0.05*var_covar, Nsim = 500000, Nburn = 100000, Nthin = 1000)

# label columns to make interpretation easier
colnames(outs) <- c("muF", "muL", "sig_mu", "alpha_bar", "sig_alpha", "mu1", "mu2", "mu3", "mu4", "mu5", "mu6",
                    "mu7", "mu8", "mu9", "mu10", "mu11", "mu12", "mu13", "mu14", "mu15", "mu16", "mu17", "mu18",
                    "mu19", "mu20", "alpha1", "alph2", "alpha3", "alpha4", "alpha5", "alpha6", "alpha7", "alpha8",
                    "alpha9", "alpha10", "alpha11", "alpha12", "alpha13", "alpha14", "alpha15", "alpha16","alpha17",
                    "alpha18", "alpha19", "alpha20", "posterior_density")

# save MCMC output for later analysis
write_csv(as.data.frame(outs),"homework4/mcmc_output.csv")

# read in MCMC results ----------------------------------------------------
outs <- read_csv(here("homework4", "mcmc_output.csv"))


outs2 <- mcmc(outs) # to better look at trace plots/distributions
summary(outs2)
plot(outs2)

##################################### PART B #################################################################
# extract parameter vectors -----------------------------------------------

mu12 <- outs[, "mu12"]
alpha12 <- outs[, "alpha12"]


# predict counts ----------------------------------------------------------

predict <- function(mus, alphas){
  preds <- mus*exp(alphas*15)
  
  return(preds)
}

preds <- mapply(predict, mu12, alpha12)


# plot results ------------------------------------------------------------

df <- as.data.frame(preds)
ggplot(data = df, aes(x = preds)) +
  geom_histogram(bins = 20, colour = "black", fill = "gray") +
  labs(title = "Posterior for stream 12, distance 15km") + xlab("expected count") + ylab("frequency")

############################################## PART C ########################################################

# predict mu[7,j] for each distance ------------------------------

mu7 <- outs[,"mu7"]
alpha7 <- outs[, "alpha7"]

# prediction function for mu[7,j]
predict2 <- function(mus, alphas){
  distance <- seq(0, 100, 10)
  preds <- matrix(0, nrow = length(mus), ncol = 11)
  vars <- nrow(preds)
  dists <- ncol(preds)
  
  for(i in 1:vars){
    for(j in 1:dists){
      preds[i, j] <- mus[i]*exp(alphas[i]*distance[j])
    }
  }
  
  return(preds)
}

# store results

preds2 <- predict2(mu7, alpha7)


# predict counts at each distance -----------------------------------------

ndists <- ncol(preds2)
nsims <- nrow(preds2)

pred_counts <- matrix(0, nrow = nsims, ncol = ndists)

for(i in 1:nsims){
  for(j in 1:ndists){
    pred_counts[i,j] <- rpois(1, preds2[i,j])
  }
}

# plot results ------------------------------------------------------------
df <- matrix(0, nrow = 11, ncol = 5)
df[,1] <- seq(0, 100, 10)
df[,2] <- unname(as_vector(counts[7,]))

for(i in 1:ncol(pred_counts)){
  df[i,3] <- unname(quantile(pred_counts[,i], probs = c(0.05)))
}

for(i in 1:ncol(pred_counts)){
  df[i,4] <- unname(quantile(pred_counts[,i], probs = c(0.5)))
}

for(i in 1:ncol(pred_counts)){
  df[i,5] <- unname(quantile(pred_counts[,i], probs = c(0.95)))
}

colnames(df) <- c("distance", "observed", "five_perc", "median", "ninetyfive_perc")
df <- as.data.frame(df)

ggplot(data = df, aes(x = distance, y = median)) +
  geom_line() + geom_line(data = df, aes(x = distance, y = five_perc), linetype = 2) +
  geom_line(data = df, aes(x = distance, y = ninetyfive_perc), linetype = 2) +
  geom_point(data = df, aes(x = distance, y = observed), colour = "red") + ylab("mosquito counts") +
  labs(title = "Posterior predictive distribution for stream 7")
