##################################################################
## Useful functions needed to run the developed methods in this ## 
## collection of code                                           ##
##################################################################

#############################################
# PAR Function Bootstrap & Jackknife Method #
#############################################

PAR2 <- function(a, b, c, d){
  n= a + b + c + d
  ((a + c)/n) - (c/(c + d))
}

#############################################################################
# Transforming by Fisher Z PAR function for Bootstrap and Jackknife Methods #
#############################################################################

Partrans <- function(PAR){
  TPar <- (1/2)*log((1+PAR)/(1-PAR))
}

######################################
# Undoing Fisher's Z transformation  #
######################################

Untrans <- function(z){
  un.trans <- (exp(2*z)-1)/(exp(2*z)+1)
}

##############################################################
# Function for creating multinomial contingency tables for   #
# cross-Sectional studies. Function required for simulations #
##############################################################

ctable <- function(p,q,e,n){
  p1 <- p*e          # a 
  p2 <- q*(1-e)      # c
  p3 <- (1-p)*e      # b
  p4 <- (1-q)*(1-e)  # d  
  probs <- c(p1,p2,p3,p4)
  
  # This repeat function has been added to avoid situations where both c and d are #
  # 0 resulting in an NA when calculating the PAR. This type of situation in       #
  # reality is unlikely.                                                           #
  
  repeat{
    c.table <- rmultinom(n = 1, size = n, prob = probs)
    if(c.table[3,] + c.table[4,] > 0){
      break
    }
  }
  c.table
}

##################################################################
# Estimation of the PAR and its variance using the delta method. #
# This function is required to perform simulations               #
##################################################################

DeltaPAR <- function(a,b,c,d){
  library(Matrix)
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
  colnames(Restrans) <- c("Par.trans", "SE.trans")
  
  AllRes <- data.frame(cbind(Res, Restrans)); AllRes
}

