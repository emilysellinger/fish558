Nsize <- 3; Narea <- 2;Ndim <- Nsize*Narea

# Read in the marices 
DataFile <- "LectB3.txt"
S <- matrix(scan(DataFile,skip=1,n=Ndim*Ndim,quiet=T),nrow=Ndim,ncol=Ndim,byrow=T)
A <- matrix(scan(DataFile,skip=Ndim+2,n=Ndim*Ndim,quiet=T),nrow=Ndim,ncol=Ndim,byrow=T)
X <- matrix(scan(DataFile,skip=2*Ndim+3,n=Ndim*Ndim,quiet=T),nrow=Ndim,ncol=Ndim,byrow=T)

# Identity matrix - will be used often
I <- matrix(0,ncol=Ndim,nrow=Ndim); diag(I) = 1

FindEqn <- function(FF=0)
 {
  
  # Specify the H matrix
  H <- matrix(0,ncol=Ndim,nrow=Ndim); diag(H) = 1
  H[Nsize,Nsize] <- 1.0 - 0.2/0.8*FF
  H[Nsize+Nsize,Nsize+Nsize] <- 1.0 - FF
  
  # Multiply the matrices
  Mat1 <- S
  Mat2 <- (X %*% (A %*% (S %*% (H %*% S))))

  # Matrix inversion multiplied by a recruitment vector with a 1 in the first row
  Neqn <- solve(I-Mat2)[,1]
  #print(Neqn)
  
  # Code to chck the equilibrium is correct
  #Test <- S %*% Neqn
  #Test <- H %*% Test
  #Test <- S %*% Test
  #Test <- A %*% Test
  #Test <- X %*% Test
  #Test[1] <- Test[1] + 1
  #print(Test-Neqn)
  
  # Compute the yield
  Catch <- sum((I-H) %*% (Mat1 %*% Neqn))
  
  # Return
  Outs <- NULL
  Outs$FF <- FF
  Outs$Neqn <- Neqn
  Outs$Catch <- Catch
  return(Outs)
 }  


# =================================================================================

# Part 1 (yield vs F)
FFs <- seq(from=0,to=1,by=0.1)
Yields <- rep(0,length(FFs))
Spawn <- rep(0,length(FFs))
for (II in 1:length(FFs))
 {
  ModelPrj <- FindEqn(FFs[II])
  Yields[II] <- ModelPrj$Catch
  Spawn[II] <- sum(ModelPrj$Neqn[c(3,6)])
 }
par(mfrow=c(2,1),oma=c(1,1,1,1))
plot(FFs,Yields,xlab="Fishing effort",ylab="Yield",type="l",lty=1)
plot(FFs,Spawn,xlab="Fishing effort",ylab="Spawners",type="l",lty=1)



# Bonus Q -------------------------------------------------------------------------------
# goal is to find what level of effort that results in the number of mature individuals 
# (age class 3) at 35% of their unfished state
func <- function(FF,Spawn0)
{
  func <-  sum(FindEqn(FF)$Neqn[c(3,6)]) -0.35*Spawn0
  return(func)
}  


# Part 2 find F such that mature numbers are 35% of that in equilirbium
Spawn0 <- sum(FindEqn(0)$Neqn[c(3,6)])
soln <- uniroot(func,lower=0,upper=1,Spawn0=Spawn0)$root
cat("Optimal effort = ",soln,"; check =",sum(FindEqn(soln)$Neqn[c(3,6)])/Spawn0,"\n")



