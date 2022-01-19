
# Sophia's code ----------------------------------------------------------------- 
gray <- read.table(here("day3", "GRAY.TXT"), header = F)
colnames(gray) <- c("year", "abundance")

sigma <- 0.1
year <- gray[, 1]
grid <- 100

# Population model ------------------------------------------------------------
pop_model <- function(p1968, alpha, year, sigma) {
  epsilon <- rnorm(n=length(year), mean=0, sd=sigma^2)
  pop = p1968 * exp(alpha * (year - 1968)) * exp(epsilon)
  return(pop)
}

# Likelihood ------------------------------------------------------------------
likelihood <- function(pred, abundance) {
  logLike <- sum(dpois(abundance, pred, log=TRUE))
  return(logLike)
}

# Grid search -----------------------------------------------------------------
xp <- seq(from=10000, to=15000, length.out=grid)
p1968 <- dnorm(xp, mean=12000, sd=(1000))
  
xa <- seq(from=0, to=0.05, length.out=grid)
alpha <- dnorm(xa, mean=0.04, sd=(0.01))
  
priors <- like <- matrix(data=NA, nrow=grid, ncol=grid)
pred <- c()
abundance <- gray[, 2]
  
for(i in 1:grid) {
  for(j in 1:grid) {
    pred <- pop_model(p1968=xp[i], alpha=xa[j], year=year, sigma=sigma)
    like[i, j] <- likelihood(pred=pred, abundance=abundance)
    priors[i, j] <- p1968[i] * alpha[j]
  }
}


post <- like*priors

PofD1 <- sum(post)

jointpos <- post / PofD1
mar_p1968 <- p1968 / PofD1
mar_alpha <- alpha / PofD1

# Plot ------------------------------------------------------------------------
persp(xp, xa, jointpos, 
      xlab=expression(xp), ylab=expression(xa), zlab="Density", 
      r=8, theta=45)
plot(xp, mar_p1968, type="l")
plot(xa, mar_alpha, type="l")
