#############################################################
## Estimating a confidence interval for the PAR where the  ##
## Delta method has being used to estimate the variance.   ##
## A Fisher-z-transformation is also provided by this      ## 
## function. Note that before this function will run the   ##
## user must first source Useful_Functions.R and install   ##
## the package Matrix.                                     ##  
##                                                         ##
## Function inputs: a, b, c and d which represent observed ##
## values in the following 2 x 2 table                     ##
##       |D+ | D-                                          ##
##    -------------                                        ##
##    E+ | a | b                                           ##
##    E- | c | d                                           ##
##                                                         ##  
## Function output: PAR, Lower & Upper CI (untransformed), ##
## LowerF & UpperF CI (Fisher transformed)                 ##
#############################################################

# DeltaCI(a=22, b=25, c=82, d=251) 

DeltaCI <- function(a,b,c,d, conf =0.95){
  require(Matrix)
  n= a + b + c + d
  PAR = ((a + c)/n) - (c/(c + d))
  
  # Partial derivatives of the PAR function #
  dPARda <- (1/n)-((a + c)/n^2)
  dPARdb <- -((a + c)/n^2)
  dPARdc <- (c/(c + d)^2) - ((a + c)/n^2) + (1/n) - (1/(c + d))
  dPARdd <-  (c/(c + d)^2) - ((a + c)/n^2)
  
  # Calculating the SE for untransformed Delta Method #
  f <- matrix(c(dPARda, dPARdb, dPARdc, dPARdd), nrow=1, ncol=4)
  ftrans <- t(f)
  V <- Diagonal(n=4, x=c(a,b,c,d))
  
  VAR = as.vector((f)%*%(V)%*%(ftrans))
  
  SE = sqrt(VAR)
  Res = cbind(PAR, SE)
  # CIs for Delta Method
  multiple <- qnorm(0.5 + conf/2)
  Lower <- Res[1] - multiple*Res[2]
  Upper <- Res[1] + multiple*Res[2]
  
  # Calculating the SE for the Fisher-Z-transformed Delta Method #
  K <- 0.5*((1/(1+PAR))+(1/(1-PAR)))      # K = Constant
  dgPARda <- K*dPARda
  dgPARdb <- K*dPARdb 
  dgPARdc <- K*dPARdc
  dgPARdd <- K*dPARdd   
  
  fK <- matrix(c(dgPARda, dgPARdb, dgPARdc, dgPARdd), nrow=1, ncol=4)
  fKtrans <- t(fK)
  VK <- Diagonal(n=4, x=c(a,b,c,d))
  Var.transA = as.vector((fK)%*%(V)%*%(fKtrans))
  
  SE.trans = sqrt(Var.transA)
  Restrans= cbind(Partrans(PAR), SE.trans)
  
  # CIs for Fisher-transformed Delta Method
  Lowert <- Restrans[1] - multiple*Restrans[2]
  Uppert <- Restrans[1] + multiple*Restrans[2]
  # Back-transformed CIs
  LowerF <- Untrans(Lowert)
  UpperF <- Untrans(Uppert)
  
  CI <- data.frame(cbind(PAR,Lower, Upper, LowerF, UpperF), row.names=""); CI
}
