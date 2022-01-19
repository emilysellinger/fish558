library(devtools)
library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) # can run STAN in parallel, don't really need for this prob

set.seed(666)
x <- rnorm(100,1,2)
y <- 5 + 2*x + rnorm(100,0,2)
plot(x,y)

# note: parameter/variable names NEED to match 

Lecture1a <- function()
{
 
 # two chains 
 inits1 <- list(alpha=2,beta=1,sigma=0.5) # initial values of the parameters for chain 1
 inits2 <- list(alpha=1,beta=2,sigma=1) # initial values of paramas for chain 2
 inits <- list(inits1,inits2) # note: this is a list of lists
 reg_data <- list(x=x,y=y,N=100)
  
 fit <- stan(file = 'classes/intro_to_stan/Regression.stan', data = reg_data, 
              iter = 2000, chains = 2,init=inits,warmup=500,verbose=F)
 print(fit)
 la <- extract(fit, permuted = FALSE) # return a list of arrays 
 pairs(la)  
 
  
 # one chain - note the list of lists for inits
 inits1 <- list(alpha=2,beta=1,sigma=0.5)
 inits <- list(inits1)
 reg_data <- list(x=x,y=y,N=100)
 
 #fit <- stan(file = 'Regression.stan', data = reg_data, 
 #           iter = 2000, chains = 1,init=inits,warmup=500,verbose=F,
 #            diagnostic_file = "lect1.diag",sample_file="lect1.samp")
 #print(fit)
 
 fit <- stan(file = 'classes/intro_to_stan/Regression.stan', data = reg_data, 
             iter = 100, chains = 2,verbose=F)
 print(fit)
 # Extract just alpha and beta
 la2 <-extract(fit, pars=c("alpha","beta"),permuted = FALSE)
 print(la2)
 plot(fit)
 #print(summary(fit))

}  

Lecture1a()
